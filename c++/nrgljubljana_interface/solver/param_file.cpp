#include "./param_file.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace nrgljubljana_interface {

  namespace {

    std::string join(std::vector<std::string> const &items) {
      std::string result;
      for (auto const &item : items) result += (result.empty() ? "" : " ") + item;
      return result;
    }

  } // namespace

  std::string render_param_file(constr_params_t const &cp, solve_params_t const &sp, nrg_params_t const &np) {
    std::ostringstream F;
    F << std::boolalpha; // important: we want to output true/false strings
    F << std::setprecision(std::numeric_limits<double>::max_digits10); // number of decimal digits necessary to differentiate all values of this type
    F << "[extra]" << std::endl;
    for (const auto &i : sp.model_parameters) F << i.first << "=" << i.second << std::endl;
    F << "[param]" << std::endl;
    F << "bandrescale=" << sp.bandrescale << std::endl; // !
    F << "model=" << cp.model << std::endl;
    F << "symtype=" << cp.symtype << std::endl;
    F << "Lambda=" << sp.Lambda << std::endl;
    F << "xmax=" << np.xmax << std::endl;
    F << "Nmax=" << np.Nmax << std::endl;
    F << "mMAX=" << np.mMAX << std::endl;
    F << "keep=" << sp.keep << std::endl;
    F << "keepenergy=" << sp.keepenergy << std::endl;
    F << "keepmin=" << sp.keepmin << std::endl;
    F << "T=" << sp.T << std::endl;
    F << "ops=" << cp.ops << std::endl;
    F << "specs=" << cp.specs << std::endl;
    F << "specd=" << cp.specd << std::endl;
    F << "spect=" << cp.spect << std::endl;
    F << "specq=" << cp.specq << std::endl;
    F << "specot=" << cp.specot << std::endl;
    F << "specgt=" << np.specgt << std::endl;
    F << "speci1t=" << np.speci1t << std::endl;
    F << "speci2t=" << np.speci2t << std::endl;
    F << "specchit=" << cp.specchit << std::endl;
    F << "specv3=" << cp.specv3 << std::endl;
    F << "v3mm=" << np.v3mm << std::endl;
    F << "dmnrg=" << np.dmnrg << std::endl;
    F << "cfs=" << np.cfs << std::endl;
    F << "fdm=" << np.fdm << std::endl;
    F << "fdmexpv=" << np.fdmexpv << std::endl;
    F << "dmnrgmats=" << np.dmnrgmats << std::endl;
    F << "fdmmats=" << np.fdmmats << std::endl;
    F << "mats=" << np.mats << std::endl;
    F << "alpha=" << sp.alpha << std::endl;
    F << "gamma=" << sp.gamma << std::endl; // ?
    F << "discretization=" << np.discretization << std::endl;
    F << "z=" << np.z << std::endl;
    F << "Nz=" << sp.Nz << std::endl; // !
    F << "polarized=" << cp.polarized << std::endl;
    F << "pol2x2=" << cp.pol2x2 << std::endl;
    F << "rungs=" << cp.rungs << std::endl;
    F << "tri=" << np.tri << std::endl;
    F << "preccpp=" << np.preccpp << std::endl;
    F << "diag=" << np.diag << std::endl;
    F << "diagratio=" << np.diagratio << std::endl;
    F << "dsyevrlimit=" << np.dsyevrlimit << std::endl;
    F << "zheevrlimit=" << np.zheevrlimit << std::endl;
    F << "restart=" << np.restart << std::endl;
    F << "restartfactor=" << np.restartfactor << std::endl;
    F << "safeguard=" << np.safeguard << std::endl;
    F << "safeguardmax=" << np.safeguardmax << std::endl;
    F << "fixeps=" << np.fixeps << std::endl;
    F << "betabar=" << np.betabar << std::endl;
    F << "gtp=" << np.gtp << std::endl;
    F << "chitp=" << np.chitp << std::endl;
    F << "finite=" << np.finite << std::endl;
    F << "cfsgt=" << np.cfsgt << std::endl;
    F << "cfsls=" << np.cfsls << std::endl;
    F << "fdmgt=" << np.fdmgt << std::endl;
    F << "fdmls=" << np.fdmls << std::endl;
    F << "fdmexpvn=" << np.fdmexpvn << std::endl;
    F << "finitemats=" << np.finitemats << std::endl;
    F << "dm=" << np.dm << std::endl;
    F << "broaden_max=" << cp.mesh_max << std::endl;                // !
    F << "broaden_min=" << cp.mesh_min << std::endl;                // !
    F << "broaden_ratio=" << cp.mesh_ratio << std::endl;            // !
    F << "broaden_min_ratio=" << np.broaden_min_ratio << std::endl; // keep this one?
    F << "omega0=" << np.omega0 << std::endl;
    F << "omega0_ratio=" << np.omega0_ratio << std::endl;
    F << "diagth=" << np.diagth << std::endl;
    F << "substeps=" << np.substeps << std::endl;
    F << "strategy=" << np.strategy << std::endl;
    F << "Ninit=" << np.Ninit << std::endl;
    F << "reim=" << np.reim << std::endl;
    F << "dumpannotated=" << np.dumpannotated << std::endl;
    F << "dumpabs=" << np.dumpabs << std::endl;
    F << "dumpscaled=" << np.dumpscaled << std::endl;
    F << "dumpprecision=" << np.dumpprecision << std::endl;
    F << "dumpgroups=" << np.dumpgroups << std::endl;
    F << "grouptol=" << np.grouptol << std::endl;
    F << "dumpdiagonal=" << np.dumpdiagonal << std::endl;
    F << "savebins=" << np.savebins << std::endl;
    F << "broaden=" << np.broaden << std::endl;
    F << "emin=" << np.emin << std::endl;
    F << "emax=" << np.emax << std::endl;
    F << "bins=" << np.bins << std::endl;
    F << "accumulation=" << np.accumulation << std::endl;
    F << "linstep=" << np.linstep << std::endl;
    F << "discard_trim=" << np.discard_trim << std::endl;
    F << "discard_immediately=" << np.discard_immediately << std::endl;
    F << "goodE=" << np.goodE << std::endl;
    F << "NN1=" << np.NN1 << std::endl;
    F << "NN2even=" << np.NN2even << std::endl;
    F << "NN2avg=" << np.NN2avg << std::endl;
    F << "NNtanh=" << np.NNtanh << std::endl;
    F << "width_td=" << np.width_td << std::endl;
    F << "width_custom=" << np.width_custom << std::endl;
    F << "prec_td=" << np.prec_td << std::endl;
    F << "prec_custom=" << np.prec_custom << std::endl;
    F << "prec_xy=" << np.prec_xy << std::endl;
    F << "resume=" << np.resume << std::endl;
    F << "log=" << np.log << std::endl;
    F << "logall=" << np.logall << std::endl;
    F << "done=" << np.done << std::endl;
    F << "calc0=" << np.calc0 << std::endl;
    F << "lastall=" << np.lastall << std::endl;
    F << "lastalloverride=" << np.lastalloverride << std::endl;
    F << "dumpsubspaces=" << np.dumpsubspaces << std::endl;
    F << "dump_f=" << np.dump_f << std::endl;
    F << "dumpenergies=" << np.dumpenergies << std::endl;
    F << "logenumber=" << np.logenumber << std::endl;
    F << "stopafter=" << np.stopafter << std::endl;
    F << "forcestop=" << np.forcestop << std::endl;
    F << "removefiles=" << np.removefiles << std::endl;
    F << "noimag=" << np.noimag << std::endl;
    F << "checksumrules=" << np.checksumrules << std::endl;
    F << "checkdiag=" << np.checkdiag << std::endl;
    F << "checkrho=" << np.checkrho << std::endl;
    F << "data_has_rescaled_energies=false" << std::endl; // NRGLj-TRIQS interface uses unscaled energies
    return F.str();
  }

  void apply_request(solve_request_t const &r, templates::manifest_t const &m, std::map<std::string, char> const &correlator_type,
                     constr_params_t &cp, nrg_params_t &np) {
    std::map<char, std::vector<std::string>> lists;
    auto const names = requested_correlators(r, m);
    for (auto const &name : names)
      for (auto const &entry : m.correlators.at(name)) {
        auto &list = lists[correlator_type.at(entry.nrg_name)];
        if (std::find(list.begin(), list.end(), entry.nrg_name) == list.end()) list.push_back(entry.nrg_name);
      }
    cp.specs  = join(lists['s']);
    cp.specd  = join(lists['d']);
    cp.spect  = join(lists['t']);
    cp.specq  = join(lists['q']);
    cp.specot = join(lists['o']);

    np.savebins = not names.empty();
    if (r.fdm) {
      np.fdm     = true;
      np.fdmexpv = true;
    } else {
      np.fdm     = false;
      np.fdmexpv = false;
      np.dmnrg   = false;
      np.cfs     = false;
    }
    if (r.energies > 0) np.dumpannotated = static_cast<size_t>(r.energies);
  }

  std::vector<double> twist_values(solve_request_t const &r, int Nz) {
    if (not r.z.empty()) return r.z;
    std::vector<double> z;
    for (int i = 1; i <= Nz; ++i) z.push_back(double(i) / Nz);
    return z;
  }

  bool complete_z_average(std::vector<double> const &z, int Nz) {
    if (Nz < 1 or z.size() != static_cast<std::size_t>(Nz)) return false;
    auto sorted = z;
    std::sort(sorted.begin(), sorted.end());
    for (int i = 1; i <= Nz; ++i)
      if (std::abs(sorted[i - 1] - double(i) / Nz) > 1e-12) return false;
    return true;
  }

  bool check_method(solve_params_t const &sp, solve_request_t const &r) {
    if (sp.method.empty()) return false;
    if (sp.method != "fdm") throw std::invalid_argument("solve_params.method=\"" + sp.method + "\" is not supported; only \"fdm\"");
    if (not r.fdm) throw std::invalid_argument("solve_params.method=\"fdm\" contradicts a request with fdm=false");
    return true;
  }

} // namespace nrgljubljana_interface
