#include <c2py/c2py.hpp>

#include <nda/c2py/converters.hpp>
#include <triqs/c2py_converters/gf.hpp>
#include <triqs/c2py_converters/mesh.hpp>

#include <nrgljubljana_interface/dmft/causality.hpp>
#include <nrgljubljana_interface/dmft/self_energy.hpp>
#include <nrgljubljana_interface/dmft/measurements.hpp>
#include <nrgljubljana_interface/dmft/lattice/lattice.hpp>
#include <nrgljubljana_interface/dmft/lattice/bethe_lattice.hpp>
#include <nrgljubljana_interface/dmft/lattice/tabulated_lattice.hpp>
#include <nrgljubljana_interface/dmft/lattice/callable_lattice.hpp>
#include <nrgljubljana_interface/dmft/loop/mixing.hpp>
#include <nrgljubljana_interface/dmft/loop/mu_control.hpp>
#include <nrgljubljana_interface/dmft/loop/convergence.hpp>
#include <nrgljubljana_interface/dmft/loop/dmft_loop.hpp>
#include <nrgljubljana_interface/solver_core.hpp>

#include "../solver_core.wrap.hxx"

#include <memory>
#include <optional>
#include <type_traits>
#include <variant>

namespace nrgljubljana_interface::dmft {

  /// The lattices and the mixers a loop can take from Python
  using lattice_variant_t = std::variant<bethe_lattice_t, tabulated_lattice_t, callable_lattice_t>;
  using mixer_variant_t   = std::variant<linear_mixer, broyden_mixer>;

  namespace detail {

    template <typename Base, typename Variant> std::unique_ptr<Base> copy_of(Variant const &v) {
      return std::visit([](auto const &x) -> std::unique_ptr<Base> { return std::make_unique<std::decay_t<decltype(x)>>(x); }, v);
    }

    // Copies of the lattice and the mixer, constructed before the loop that refers to them
    struct loop_parts_t {
      std::unique_ptr<lattice_t const> lattice_copy;
      std::unique_ptr<mixer_t> mixer_copy;
      loop_parts_t(lattice_variant_t const &lattice, mixer_variant_t const &mixer)
         : lattice_copy(copy_of<lattice_t const>(lattice)), mixer_copy(copy_of<mixer_t>(mixer)) {}
    };

  } // namespace detail

  /// The DMFT loop with the solver of the package; the loop holds copies of the lattice and the mixer
  class py_dmft_loop : private detail::loop_parts_t, public dmft_loop {
   public:
    py_dmft_loop(solver_core &solver, lattice_variant_t const &lattice, mixer_variant_t const &mixer, convergence_t &convergence, loop_params_t params,
                 std::optional<mu_control_t> const &mu_control = {})
       : detail::loop_parts_t(lattice, mixer), dmft_loop(solver, *lattice_copy, *mixer_copy, convergence, std::move(params), mu_control) {}
  };

  /// μ for the next iteration, with Σ held fixed (see MuControlT)
  C2PY_WRAP_AS_METHOD mu_result_t adjust(mu_control_t const &self, mu_input_t const &input, lattice_variant_t const &lattice) {
    return std::visit([&](auto const &l) { return self.adjust(input, l); }, lattice);
  }

  /// Only the filling at μ, one lattice evaluation (see MuControlT)
  C2PY_WRAP_AS_METHOD mu_result_t filling(mu_control_t const &self, mu_input_t const &input, lattice_variant_t const &lattice) {
    return std::visit([&](auto const &l) { return self.filling(input, l); }, lattice);
  }

  /// Call the function with the StepResult of each step, after it is committed; what it returns is ignored
  C2PY_WRAP_AS_METHOD void set_callback(py_dmft_loop &self, std::function<c2py::pyref(step_result_t const &)> callback) {
    self.set_callback([callback](step_result_t const &r) { (void)callback(r); });
  }

  /// Δ of the next step; throws before initialize
  C2PY_WRAP_AS_METHOD g_w_t Delta(py_dmft_loop const &self) { return self.Delta(); }

  /// Σ of the last step; None before the first
  C2PY_WRAP_AS_METHOD std::optional<g_w_t> Sigma(py_dmft_loop const &self) { return self.Sigma(); }

  /// G_loc of the last step; None before the first
  C2PY_WRAP_AS_METHOD std::optional<g_w_t> G_loc(py_dmft_loop const &self) { return self.G_loc(); }

  /// The result of the last step; None before the first
  C2PY_WRAP_AS_METHOD std::optional<step_result_t> last_step(py_dmft_loop const &self) { return self.last_step(); }

} // namespace nrgljubljana_interface::dmft

#include "_dmft.wrap.cxx"
