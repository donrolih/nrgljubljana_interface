#include "./output.hpp"

#include <charconv>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <system_error>

namespace nrgljubljana_interface::nrg {

  namespace {

    constexpr std::string_view whitespace = " \t\r\n\f\v";

    std::string_view trim(std::string_view s) {
      auto const first = s.find_first_not_of(whitespace);
      if (first == std::string_view::npos) return {};
      auto const last = s.find_last_not_of(whitespace);
      return s.substr(first, last - first + 1);
    }

    std::vector<std::string_view> fields(std::string_view s) {
      std::vector<std::string_view> result;
      std::size_t pos = 0;
      while ((pos = s.find_first_not_of(whitespace, pos)) != std::string_view::npos) {
        auto const end = s.find_first_of(whitespace, pos);
        result.push_back(s.substr(pos, end == std::string_view::npos ? std::string_view::npos : end - pos));
        pos = end;
      }
      return result;
    }

    // Lines of a text with their numbers
    class lines {
      public:
      lines(std::string_view text, std::string_view source) : text_(text), source_(source) {}

      bool next(std::string_view &line) {
        if (text_.empty()) return false;
        auto const eol = text_.find('\n');
        line           = text_.substr(0, eol);
        text_          = eol == std::string_view::npos ? std::string_view{} : text_.substr(eol + 1);
        ++number_;
        return true;
      }

      [[noreturn]] void fail(std::string const &what) const {
        std::string msg;
        if (not source_.empty()) msg += std::string(source_) + ": ";
        if (number_ > 0) msg += "line " + std::to_string(number_) + ": ";
        throw std::invalid_argument(msg + what);
      }

      private:
      std::string_view text_, source_;
      std::size_t number_ = 0;
    };

    template <typename T> bool parse_number(std::string_view s, T &value) {
      auto const [end, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
      return ec == std::errc{} and end == s.data() + s.size();
    }

  } // namespace

  step_table_t parse_step_table(std::string_view text, std::string_view source) {
    lines in(text, source);
    std::vector<std::string> columns;
    std::vector<double> values;
    std::size_t rows = 0;
    std::string_view raw;
    while (in.next(raw)) {
      auto const line = trim(raw);
      if (line.empty()) continue;
      if (line.front() == '#') {
        if (rows > 0) in.fail("header line after data");
        columns.clear();
        for (auto f : fields(line.substr(1))) columns.emplace_back(f);
        continue;
      }
      if (columns.empty()) in.fail("data before a header line with column names");
      auto const f = fields(line);
      if (f.size() != columns.size()) in.fail("expected " + std::to_string(columns.size()) + " values, got " + std::to_string(f.size()));
      for (auto v : f) {
        double x = 0;
        if (not parse_number(v, x)) in.fail("not a number: \"" + std::string(v) + "\"");
        values.push_back(x);
      }
      ++rows;
    }
    if (columns.empty()) in.fail("no header line with column names");

    step_table_t t{columns, nda::array<double, 2>(static_cast<long>(rows), static_cast<long>(columns.size()))};
    for (std::size_t r = 0; r < rows; ++r)
      for (std::size_t c = 0; c < columns.size(); ++c) t.values(static_cast<long>(r), static_cast<long>(c)) = values[r * columns.size() + c];
    return t;
  }

  std::map<std::string, double> parse_final_values(std::string_view text, std::string_view source) {
    auto const t = parse_step_table(text, source);
    if (t.values.extent(0) != 1) {
      std::string msg = source.empty() ? std::string() : std::string(source) + ": ";
      throw std::invalid_argument(msg + "expected exactly one row, got " + std::to_string(t.values.extent(0)));
    }
    std::map<std::string, double> result;
    for (std::size_t c = 0; c < t.columns.size(); ++c) result[t.columns[c]] = t.values(0, static_cast<long>(c));
    return result;
  }

  energy_flow_t parse_energy_flow(std::string_view text, std::string_view source) {
    lines in(text, source);
    energy_flow_t flow;
    std::vector<energy_level_t> step;
    std::string_view raw;
    while (in.next(raw)) {
      auto line = trim(raw);
      if (line.empty()) { // end of a step; a step without levels is a lone empty line
        flow.push_back(std::move(step));
        step.clear();
        continue;
      }
      energy_level_t level;
      auto const first_end = line.find_first_of(whitespace);
      if (first_end == std::string_view::npos or not parse_number(line.substr(0, first_end), level.energy))
        in.fail("expected an energy, then (quantum numbers) and [degeneracy]");
      line = trim(line.substr(first_end));
      while (not line.empty() and line.front() == '(') {
        auto const close = line.find(')');
        if (close == std::string_view::npos) in.fail("unclosed '('");
        std::vector<int> qn;
        for (auto v : fields(line.substr(1, close - 1))) {
          int x = 0;
          if (not parse_number(v, x)) in.fail("not an integer quantum number: \"" + std::string(v) + "\"");
          qn.push_back(x);
        }
        level.subspaces.push_back(std::move(qn));
        line = trim(line.substr(close + 1));
      }
      if (level.subspaces.empty()) in.fail("no (quantum numbers)");
      if (line.size() < 3 or line.front() != '[' or line.back() != ']' or not parse_number(line.substr(1, line.size() - 2), level.degeneracy)
          or level.degeneracy < 1)
        in.fail("expected [degeneracy] at the end");
      step.push_back(std::move(level));
    }
    if (not step.empty()) flow.push_back(std::move(step));
    return flow;
  }

  std::string spectrum_file(std::string const &nrg_correlator, char type) {
    std::string prefix;
    switch (type) {
      case 's': prefix = "corr"; break;
      case 't': prefix = "spin"; break;
      case 'o': prefix = "orbspin"; break;
      case 'd': prefix = "spec"; break;
      case 'q': prefix = "specq"; break;
      default: throw std::invalid_argument(std::string("no spectrum list of type '") + type + "'");
    }
    return prefix + "_FDM_dens_" + nrg_correlator + ".bin";
  }

  std::string read_text_file(std::string const &path) {
    std::ifstream f(path, std::ios::binary);
    if (not f) throw std::runtime_error(path + ": cannot be read");
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
  }

} // namespace nrgljubljana_interface::nrg
