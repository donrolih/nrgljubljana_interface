
// C.f. https://numpy.org/doc/1.21/reference/c-api/array.html#importing-the-api
#define PY_ARRAY_UNIQUE_SYMBOL _cpp2py_ARRAY_API
#ifndef CLAIR_C2PY_WRAP_GEN
#ifdef __clang__
// #pragma clang diagnostic ignored "-W#warnings"
#endif
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wcast-function-type"
#pragma GCC diagnostic ignored "-Wcpp"
#endif

#define C2PY_VERSION_MAJOR 0
#define C2PY_VERSION_MINOR 1

#include <c2py/c2py.hpp>
#include <c2py/serialization/h5.hpp>

using c2py::operator""_a;

// ==================== enums =====================

// ==================== module classes =====================

// --------- class _c2py_cls_0 -----------
using _c2py_cls_0                                            = nrgljubljana_interface::dmft::causality_report_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_0>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_0> = "nrgljubljana_interface.dmft._dmft.CausalityReportT";

static int synth_constructor_0(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(
       PyExc_RuntimeError,
       ("Error in constructing nrgljubljana_interface::dmft::causality_report_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_0> *)self)->_c = new _c2py_cls_0{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::causality_report_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_0> *)self)->_c);
  de("projected_points", self_c.projected_points, true);
  de("largest_violation", self_c.largest_violation, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_0> = synth_constructor_0;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_0> = c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
projected_points : {par_0}, default=0

largest_violation : {par_1}, default=0

)DOC",
                                                                      "par", {c2py::python_typename<long>(), c2py::python_typename<double>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_0>[] = {
   {"__write_hdf5__", c2py::tpxx_write_h5<_c2py_cls_0>, METH_VARARGS, "  "},
   {"__getstate__", c2py::getstate_h5<_c2py_cls_0>, METH_NOARGS, ""},
   {"__setstate__", c2py::setstate_h5<_c2py_cls_0>, METH_O, ""},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_0 = R"DOC(Number of (block, mesh point) pairs at which an eigenvalue of S was raised)DOC";
constexpr auto _c2py_doc_member_1 = R"DOC(The largest amount by which an eigenvalue was raised (floor − λ); 0 if none)DOC";
static PyObject *prop_get_dict_0(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_0> *)self)->_c);
  c2py::pydict dic;
  dic["projected_points"]  = self_c.projected_points;
  dic["largest_violation"] = self_c.largest_violation;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_0>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_0::projected_points, _c2py_cls_0>("projected_points", _c2py_doc_member_0),
   c2py::getsetdef_from_member<&_c2py_cls_0::largest_violation, _c2py_cls_0>("largest_violation", _c2py_doc_member_1),
   {"__dict__", (getter)prop_get_dict_0, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_0> =
   R"DOC(What causal_projection changed)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_0>;
// --------- class _c2py_cls_1 -----------
using _c2py_cls_1                                            = nrgljubljana_interface::dmft::occupancy_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_1>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_1> = "nrgljubljana_interface.dmft._dmft.OccupancyT";

static int synth_constructor_1(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::occupancy_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_1> *)self)->_c = new _c2py_cls_1{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::occupancy_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_1> *)self)->_c);
  de("matrix", self_c.matrix, true);
  de("total", self_c.total, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_1> = synth_constructor_1;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_1> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
matrix : {par_0}, default={}

total : {par_1}, default=0

)DOC",
                      "par", {c2py::python_typename<nrgljubljana_interface::dmft::block_matrices_t>(), c2py::python_typename<double>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_1>[] = {

   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_2 = R"DOC(N_ij = ∫ f(ω) A_ij(ω) dω per block, Hermitian)DOC";
constexpr auto _c2py_doc_member_3 = R"DOC(The sum over the blocks of degeneracy × Re tr N, e.g. the spin-summed filling)DOC";
static PyObject *prop_get_dict_1(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_1> *)self)->_c);
  c2py::pydict dic;
  dic["matrix"] = self_c.matrix;
  dic["total"]  = self_c.total;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_1>[] = {c2py::getsetdef_from_member<&_c2py_cls_1::matrix, _c2py_cls_1>("matrix", _c2py_doc_member_2),
                                                        c2py::getsetdef_from_member<&_c2py_cls_1::total, _c2py_cls_1>("total", _c2py_doc_member_3),
                                                        {"__dict__", (getter)prop_get_dict_1, nullptr, "", nullptr},
                                                        {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_1> =
   R"DOC(The occupancy of a Green's function)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_1>;
// --------- class _c2py_cls_2 -----------
using _c2py_cls_2                                            = nrgljubljana_interface::dmft::distance_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_2>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_2> = "nrgljubljana_interface.dmft._dmft.DistanceT";

static int synth_constructor_2(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::distance_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_2> *)self)->_c = new _c2py_cls_2{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::distance_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_2> *)self)->_c);
  de("elements", self_c.elements, true);
  de("max", self_c.max, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_2> = synth_constructor_2;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_2> = c2py::replace_tags(
   R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
elements : {par_0}, default={}

max : {par_1}, default=0

)DOC",
   "par",
   {c2py::python_typename<
       std::vector<nda::basic_array<double, 2, nda::C_layout, 'M', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>>>(),
    c2py::python_typename<double>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_2>[] = {

   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_4 = R"DOC(Per block and element: ∫ |Re ΔA_ij| dω + ∫ |Im ΔA_ij| dω, which is ∫ |ΔA_ii| dω on the diagonal)DOC";
constexpr auto _c2py_doc_member_5 = R"DOC(The largest element over all blocks)DOC";
static PyObject *prop_get_dict_2(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_2> *)self)->_c);
  c2py::pydict dic;
  dic["elements"] = self_c.elements;
  dic["max"]      = self_c.max;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_2>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_2::elements, _c2py_cls_2>("elements", _c2py_doc_member_4),
   c2py::getsetdef_from_member<&_c2py_cls_2::max, _c2py_cls_2>("max", _c2py_doc_member_5),
   {"__dict__", (getter)prop_get_dict_2, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_2> =
   R"DOC(The distance between two spectral functions)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_2>;
// --------- class _c2py_cls_3 -----------
using _c2py_cls_3                                            = nrgljubljana_interface::dmft::bethe_lattice_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_3>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_3> = "nrgljubljana_interface.dmft._dmft.BetheLatticeT";
static const auto _c2py_init_0 = c2py::dispatcher_c_kw_t{c2py::c_constructor<_c2py_cls_3, double, double>("half_bandwidth", "sigma_floor"_a = 1e-12)};
template <> constexpr initproc c2py::tp_init<_c2py_cls_3> = c2py::pyfkw_constructor<_c2py_init_0>;
template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_3> = _c2py_init_0.doc(R"DOC(
Throws std::invalid_argument for a half-bandwidth that is not positive and finite, and for an invalid floor
)DOC");
// description
static auto const _c2py_fun_0 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_3 const &self) -> decltype(auto) { return self.description(); }, "self")};

// half_bandwidth
static auto const _c2py_fun_1 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_3 const &self) -> decltype(auto) { return self.half_bandwidth(); }, "self")};

// local_gf
static auto const _c2py_fun_2 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_3 const &self, const nrgljubljana_interface::g_w_t &Sigma, double mu) -> decltype(auto) { return self.local_gf(Sigma, mu); }, "self",
   "Sigma", "mu")};

// moments
static auto const _c2py_fun_3 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_3 const &self, const nrgljubljana_interface::g_w_t &Sigma, double mu) -> decltype(auto) { return self.moments(Sigma, mu); }, "self",
   "Sigma", "mu")};

// onsite_energy
static auto const _c2py_fun_4 = c2py::dispatcher_f_kw_t{
   c2py::cmethod([](_c2py_cls_3 const &self, const triqs::gfs::gf_struct_t &structure) -> decltype(auto) { return self.onsite_energy(structure); },
                 "self", "structure")};

// self_consistency
static auto const _c2py_fun_5 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_3 const &self, const nrgljubljana_interface::g_w_t &Sigma, double mu) -> decltype(auto) { return self.self_consistency(Sigma, mu); },
   "self", "Sigma", "mu")};

// sigma_floor
static auto const _c2py_fun_6 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_3 const &self) -> decltype(auto) { return self.sigma_floor(); }, "self")};

static const auto _c2py_doc_0 = _c2py_fun_0.doc(R"DOC(
The type and the parameters, e.g. to check on resume that the lattice is the same
)DOC");
static const auto _c2py_doc_1 = _c2py_fun_1.doc(R"DOC()DOC");
static const auto _c2py_doc_2 = _c2py_fun_2.doc(R"DOC(
G0 alone; by default the first of moments
)DOC");
static const auto _c2py_doc_3 = _c2py_fun_3.doc(R"DOC(
(G0, G1)
)DOC");
static const auto _c2py_doc_4 = _c2py_fun_4.doc(R"DOC(
⟨h⟩ per block of the given structure
)DOC");
static const auto _c2py_doc_5 = _c2py_fun_5.doc(R"DOC(
(G_loc, Δ) with G_loc = G0 and Δ = G1 G0⁻¹, which avoids subtracting two large numbers where G0 is small.
Throws std::invalid_argument for a singular G0, naming the block and the frequency.
)DOC");
static const auto _c2py_doc_6 = _c2py_fun_6.doc(R"DOC()DOC");

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_3>[] = {
   {"description", (PyCFunction)c2py::pyfkw<_c2py_fun_0>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_0.c_str()},
   {"half_bandwidth", (PyCFunction)c2py::pyfkw<_c2py_fun_1>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_1.c_str()},
   {"local_gf", (PyCFunction)c2py::pyfkw<_c2py_fun_2>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_2.c_str()},
   {"moments", (PyCFunction)c2py::pyfkw<_c2py_fun_3>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_3.c_str()},
   {"onsite_energy", (PyCFunction)c2py::pyfkw<_c2py_fun_4>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_4.c_str()},
   {"self_consistency", (PyCFunction)c2py::pyfkw<_c2py_fun_5>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_5.c_str()},
   {"sigma_floor", (PyCFunction)c2py::pyfkw<_c2py_fun_6>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_6.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

template <>
const std::string c2py::tp_doc<_c2py_cls_3> = R"DOC(The Bethe lattice of half-bandwidth D, ρ(ε) = 2 √(D² − ε²) / (π D²), the same for every orbital.

G0(ζ) = 2 / (ζ + √(ζ − D) √(ζ + D)) with principal roots, the physical branch for Im ζ > 0; G1 = (D²/4) G0²; the
on-site energy is 0. ζ = ω + μ − Re σ + i max(−Im σ, floor), as in the tabulated lattice; a matrix Σ goes through its
eigenvectors.)DOC"
   + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_3>;
// --------- class _c2py_cls_4 -----------
using _c2py_cls_4                                            = nrgljubljana_interface::dmft::dos_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_4>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_4> = "nrgljubljana_interface.dmft._dmft.DosT";

static int synth_constructor_3(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::dos_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_4> *)self)->_c = new _c2py_cls_4{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError, ("Error in constructing nrgljubljana_interface::dmft::dos_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_4> *)self)->_c);
  de("energy", self_c.energy, false);
  de("density", self_c.density, false);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_4> = synth_constructor_3;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_4> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
energy : {par_0}

density : {par_1}

)DOC",
                      "par", {c2py::python_typename<std::vector<double>>(), c2py::python_typename<std::vector<double>>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_4>[] = {

   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_6 = R"DOC()DOC";
constexpr auto _c2py_doc_member_7 = R"DOC()DOC";
static PyObject *prop_get_dict_3(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_4> *)self)->_c);
  c2py::pydict dic;
  dic["energy"]  = self_c.energy;
  dic["density"] = self_c.density;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_4>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_4::energy, _c2py_cls_4>("energy", _c2py_doc_member_6),
   c2py::getsetdef_from_member<&_c2py_cls_4::density, _c2py_cls_4>("density", _c2py_doc_member_7),
   {"__dict__", (getter)prop_get_dict_3, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_4> =
   R"DOC(A tabulated density of states: strictly increasing energies and non-negative densities, normalized to 1)DOC"
   + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_4>;
// --------- class _c2py_cls_5 -----------
using _c2py_cls_5                                            = nrgljubljana_interface::dmft::tabulated_lattice_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_5>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_5> = "nrgljubljana_interface.dmft._dmft.TabulatedLatticeT";
static const auto _c2py_init_1 =
   c2py::dispatcher_c_kw_t{c2py::c_constructor<_c2py_cls_5, nrgljubljana_interface::dmft::dos_t, double>("dos", "sigma_floor"_a = 1e-12),
                           c2py::c_constructor<_c2py_cls_5, std::map<std::string, std::vector<nrgljubljana_interface::dmft::dos_t>>, double>(
                              "per_element", "sigma_floor"_a = 1e-12)};
template <> constexpr initproc c2py::tp_init<_c2py_cls_5>    = c2py::pyfkw_constructor<_c2py_init_1>;
template <> const std::string c2py::tp_ctor_doc<_c2py_cls_5> = _c2py_init_1.doc(R"DOC()DOC");
// description
static auto const _c2py_fun_7 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_5 const &self) -> decltype(auto) { return self.description(); }, "self")};

// local_gf
static auto const _c2py_fun_8 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_5 const &self, const nrgljubljana_interface::g_w_t &Sigma, double mu) -> decltype(auto) { return self.local_gf(Sigma, mu); }, "self",
   "Sigma", "mu")};

// moments
static auto const _c2py_fun_9 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_5 const &self, const nrgljubljana_interface::g_w_t &Sigma, double mu) -> decltype(auto) { return self.moments(Sigma, mu); }, "self",
   "Sigma", "mu")};

// onsite_energy
static auto const _c2py_fun_10 = c2py::dispatcher_f_kw_t{
   c2py::cmethod([](_c2py_cls_5 const &self, const triqs::gfs::gf_struct_t &structure) -> decltype(auto) { return self.onsite_energy(structure); },
                 "self", "structure")};

// read_dos
static auto const _c2py_fun_11 = c2py::dispatcher_f_kw_t{c2py::cfun([](const std::string &path) { return _c2py_cls_5::read_dos(path); }, "path")};

// self_consistency
static auto const _c2py_fun_12 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_5 const &self, const nrgljubljana_interface::g_w_t &Sigma, double mu) -> decltype(auto) { return self.self_consistency(Sigma, mu); },
   "self", "Sigma", "mu")};

// sigma_floor
static auto const _c2py_fun_13 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_5 const &self) -> decltype(auto) { return self.sigma_floor(); }, "self")};

static const auto _c2py_doc_7  = _c2py_fun_7.doc(R"DOC(
The kind (shared or per element), the number of points and a checksum of each table, and the floor
)DOC");
static const auto _c2py_doc_8  = _c2py_fun_8.doc(R"DOC(
G0 alone; by default the first of moments
)DOC");
static const auto _c2py_doc_9  = _c2py_fun_9.doc(R"DOC(
(G0, G1)
)DOC");
static const auto _c2py_doc_10 = _c2py_fun_10.doc(R"DOC(
⟨h⟩ per block of the given structure
)DOC");
static const auto _c2py_doc_11 = _c2py_fun_11.doc(R"DOC(
A DOS from a file of two columns, energy and density; blank lines and lines starting with # are skipped
)DOC");
static const auto _c2py_doc_12 = _c2py_fun_12.doc(R"DOC(
(G_loc, Δ) with G_loc = G0 and Δ = G1 G0⁻¹, which avoids subtracting two large numbers where G0 is small.
Throws std::invalid_argument for a singular G0, naming the block and the frequency.
)DOC");
static const auto _c2py_doc_13 = _c2py_fun_13.doc(R"DOC()DOC");

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_5>[] = {
   {"description", (PyCFunction)c2py::pyfkw<_c2py_fun_7>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_7.c_str()},
   {"local_gf", (PyCFunction)c2py::pyfkw<_c2py_fun_8>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_8.c_str()},
   {"moments", (PyCFunction)c2py::pyfkw<_c2py_fun_9>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_9.c_str()},
   {"onsite_energy", (PyCFunction)c2py::pyfkw<_c2py_fun_10>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_10.c_str()},
   {"read_dos", (PyCFunction)c2py::pyfkw<_c2py_fun_11>, METH_VARARGS | METH_KEYWORDS | METH_STATIC, _c2py_doc_11.c_str()},
   {"self_consistency", (PyCFunction)c2py::pyfkw<_c2py_fun_12>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_12.c_str()},
   {"sigma_floor", (PyCFunction)c2py::pyfkw<_c2py_fun_13>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_13.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

template <>
const std::string c2py::tp_doc<_c2py_cls_5> = R"DOC(A lattice given by one tabulated DOS for all orbitals, or one per diagonal element of each block.

The moments are Hilbert transforms of the Steffen interpolant of the DOS (NRG Ljubljana numerics), at
z = ω + μ − Re Σ + i max(−Im Σ, floor). With one DOS for all orbitals a matrix Σ goes through its eigenvectors;
with a DOS per element an off-diagonal Σ is an error. The on-site energy is ∫ ε ρ(ε) dε on the diagonal.

Each DOS is checked at construction: at least 3 finite points, strictly increasing energies, ρ ≥ 0, and
∫ ρ dε = 1 within 1e-8. It is not rescaled. Errors are std::invalid_argument.)DOC"
   + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_5>;
// --------- class _c2py_cls_6 -----------
using _c2py_cls_6                                            = nrgljubljana_interface::dmft::callable_lattice_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_6>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_6> = "nrgljubljana_interface.dmft._dmft.CallableLatticeT";
static const auto _c2py_init_2                               = c2py::dispatcher_c_kw_t{
   c2py::c_constructor<_c2py_cls_6, nrgljubljana_interface::dmft::callable_lattice_t::moments_t, nrgljubljana_interface::dmft::block_matrices_t,
                                                     std::optional<nrgljubljana_interface::dmft::callable_lattice_t::local_gf_t>, std::string>(
      "moments", "onsite_energy", "local_gf"_a = std::optional<nrgljubljana_interface::dmft::callable_lattice_t::local_gf_t>{},
      "description"_a = "callable")};
template <> constexpr initproc c2py::tp_init<_c2py_cls_6>    = c2py::pyfkw_constructor<_c2py_init_2>;
template <> const std::string c2py::tp_ctor_doc<_c2py_cls_6> = _c2py_init_2.doc(R"DOC()DOC");
// description
static auto const _c2py_fun_14 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_6 const &self) -> decltype(auto) { return self.description(); }, "self")};

// local_gf
static auto const _c2py_fun_15 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_6 const &self, const nrgljubljana_interface::g_w_t &Sigma, double mu) -> decltype(auto) { return self.local_gf(Sigma, mu); }, "self",
   "Sigma", "mu")};

// moments
static auto const _c2py_fun_16 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_6 const &self, const nrgljubljana_interface::g_w_t &Sigma, double mu) -> decltype(auto) { return self.moments(Sigma, mu); }, "self",
   "Sigma", "mu")};

// onsite_energy
static auto const _c2py_fun_17 = c2py::dispatcher_f_kw_t{
   c2py::cmethod([](_c2py_cls_6 const &self, const triqs::gfs::gf_struct_t &structure) -> decltype(auto) { return self.onsite_energy(structure); },
                 "self", "structure")};

// self_consistency
static auto const _c2py_fun_18 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_6 const &self, const nrgljubljana_interface::g_w_t &Sigma, double mu) -> decltype(auto) { return self.self_consistency(Sigma, mu); },
   "self", "Sigma", "mu")};

static const auto _c2py_doc_14 = _c2py_fun_14.doc(R"DOC(
The type and the parameters, e.g. to check on resume that the lattice is the same
)DOC");
static const auto _c2py_doc_15 = _c2py_fun_15.doc(R"DOC(
G0 alone; by default the first of moments
)DOC");
static const auto _c2py_doc_16 = _c2py_fun_16.doc(R"DOC(
(G0, G1)
)DOC");
static const auto _c2py_doc_17 = _c2py_fun_17.doc(R"DOC(
⟨h⟩ per block of the given structure
)DOC");
static const auto _c2py_doc_18 = _c2py_fun_18.doc(R"DOC(
(G_loc, Δ) with G_loc = G0 and Δ = G1 G0⁻¹, which avoids subtracting two large numbers where G0 is small.
Throws std::invalid_argument for a singular G0, naming the block and the frequency.
)DOC");

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_6>[] = {
   {"description", (PyCFunction)c2py::pyfkw<_c2py_fun_14>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_14.c_str()},
   {"local_gf", (PyCFunction)c2py::pyfkw<_c2py_fun_15>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_15.c_str()},
   {"moments", (PyCFunction)c2py::pyfkw<_c2py_fun_16>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_16.c_str()},
   {"onsite_energy", (PyCFunction)c2py::pyfkw<_c2py_fun_17>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_17.c_str()},
   {"self_consistency", (PyCFunction)c2py::pyfkw<_c2py_fun_18>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_18.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

template <>
const std::string c2py::tp_doc<_c2py_cls_6> =
   R"DOC(A lattice given by functions: the moments (G0, G1) of Σ and μ, optionally G0 alone, and a fixed on-site energy
⟨h⟩ per block. The results must have the block structure and mesh of Σ and be finite; ⟨h⟩ must match the
structure it is asked for. The description is what resume compares, so it should identify the lattice. Errors
are std::invalid_argument.)DOC"
   + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_6>;
// --------- class _c2py_cls_7 -----------
using _c2py_cls_7                                            = nrgljubljana_interface::dmft::mix_input_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_7>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_7> = "nrgljubljana_interface.dmft._dmft.MixInputT";

static int synth_constructor_4(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::mix_input_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_7> *)self)->_c = new _c2py_cls_7{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::mix_input_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_7> *)self)->_c);
  de("iteration", self_c.iteration, true);
  de("Gamma", self_c.Gamma, false);
  de("Gamma_out", self_c.Gamma_out, false);
  de("mu", self_c.mu, true);
  de("filling", self_c.filling, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_7> = synth_constructor_4;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_7> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
Gamma : {par_0}

Gamma_out : {par_1}

iteration : {par_2}, default=0

mu : {par_3}, default=0

filling : {par_4}, default=0

)DOC",
                      "par",
                      {c2py::python_typename<nrgljubljana_interface::g_w_t>(), c2py::python_typename<nrgljubljana_interface::g_w_t>(),
                       c2py::python_typename<long>(), c2py::python_typename<double>(), c2py::python_typename<double>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_7>[] = {

   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_8  = R"DOC(The iteration, from 1)DOC";
constexpr auto _c2py_doc_member_9  = R"DOC(The spectral parts of Δ and of Δ_out from the lattice; same structure and mesh)DOC";
constexpr auto _c2py_doc_member_10 = R"DOC(The spectral parts of Δ and of Δ_out from the lattice; same structure and mesh)DOC";
constexpr auto _c2py_doc_member_11 = R"DOC(μ of the solve, and the filling n(μ), for a mixer that also moves μ)DOC";
constexpr auto _c2py_doc_member_12 = R"DOC(μ of the solve, and the filling n(μ), for a mixer that also moves μ)DOC";
static PyObject *prop_get_dict_4(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_7> *)self)->_c);
  c2py::pydict dic;
  dic["iteration"] = self_c.iteration;
  dic["Gamma"]     = self_c.Gamma;
  dic["Gamma_out"] = self_c.Gamma_out;
  dic["mu"]        = self_c.mu;
  dic["filling"]   = self_c.filling;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_7>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_7::iteration, _c2py_cls_7>("iteration", _c2py_doc_member_8),
   c2py::getsetdef_from_member<&_c2py_cls_7::Gamma, _c2py_cls_7>("Gamma", _c2py_doc_member_9),
   c2py::getsetdef_from_member<&_c2py_cls_7::Gamma_out, _c2py_cls_7>("Gamma_out", _c2py_doc_member_10),
   c2py::getsetdef_from_member<&_c2py_cls_7::mu, _c2py_cls_7>("mu", _c2py_doc_member_11),
   c2py::getsetdef_from_member<&_c2py_cls_7::filling, _c2py_cls_7>("filling", _c2py_doc_member_12),
   {"__dict__", (getter)prop_get_dict_4, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_7> =
   R"DOC(What a mixer is given in one iteration)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_7>;
// --------- class _c2py_cls_8 -----------
using _c2py_cls_8                                            = nrgljubljana_interface::dmft::mix_report_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_8>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_8> = "nrgljubljana_interface.dmft._dmft.MixReportT";

static int synth_constructor_5(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::mix_report_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_8> *)self)->_c = new _c2py_cls_8{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::mix_report_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_8> *)self)->_c);
  de("method", self_c.method, true);
  de("residual_norm", self_c.residual_norm, true);
  de("history_size", self_c.history_size, true);
  de("condition", self_c.condition, true);
  de("fallback", self_c.fallback, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_8> = synth_constructor_5;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_8> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
method : {par_0}, default={}

residual_norm : {par_1}, default=0

history_size : {par_2}, default=0

condition : {par_3}, default=0

fallback : {par_4}, default=false

)DOC",
                      "par",
                      {c2py::python_typename<std::string>(), c2py::python_typename<double>(), c2py::python_typename<long>(),
                       c2py::python_typename<double>(), c2py::python_typename<bool>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_8>[] = {
   {"__write_hdf5__", c2py::tpxx_write_h5<_c2py_cls_8>, METH_VARARGS, "  "},
   {"__getstate__", c2py::getstate_h5<_c2py_cls_8>, METH_NOARGS, ""},
   {"__setstate__", c2py::setstate_h5<_c2py_cls_8>, METH_O, ""},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_13 = R"DOC("linear", "broyden" or "linear (fallback)")DOC";
constexpr auto _c2py_doc_member_14 = R"DOC(‖Γ_out − Γ‖ = √(Σ_k w_k Σ_ij |ΔΓ_ij(ω_k)|²) over all blocks, w_k the integration weights of the mesh)DOC";
constexpr auto _c2py_doc_member_15 =
   R"DOC(Broyden: the number of stored iterations, the condition number of its step, whether it fell back to linear)DOC";
constexpr auto _c2py_doc_member_16 = R"DOC()DOC";
constexpr auto _c2py_doc_member_17 = R"DOC()DOC";
static PyObject *prop_get_dict_5(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_8> *)self)->_c);
  c2py::pydict dic;
  dic["method"]        = self_c.method;
  dic["residual_norm"] = self_c.residual_norm;
  dic["history_size"]  = self_c.history_size;
  dic["condition"]     = self_c.condition;
  dic["fallback"]      = self_c.fallback;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_8>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_8::method, _c2py_cls_8>("method", _c2py_doc_member_13),
   c2py::getsetdef_from_member<&_c2py_cls_8::residual_norm, _c2py_cls_8>("residual_norm", _c2py_doc_member_14),
   c2py::getsetdef_from_member<&_c2py_cls_8::history_size, _c2py_cls_8>("history_size", _c2py_doc_member_15),
   c2py::getsetdef_from_member<&_c2py_cls_8::condition, _c2py_cls_8>("condition", _c2py_doc_member_16),
   c2py::getsetdef_from_member<&_c2py_cls_8::fallback, _c2py_cls_8>("fallback", _c2py_doc_member_17),
   {"__dict__", (getter)prop_get_dict_5, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_8> =
   R"DOC(What a mixer did; the same for all mixers)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_8>;
// --------- class _c2py_cls_9 -----------
using _c2py_cls_9                                            = nrgljubljana_interface::dmft::mix_result_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_9>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_9> = "nrgljubljana_interface.dmft._dmft.MixResultT";

static int synth_constructor_6(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::mix_result_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_9> *)self)->_c = new _c2py_cls_9{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::mix_result_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_9> *)self)->_c);
  de("Gamma_mix", self_c.Gamma_mix, false);
  de("mu_next", self_c.mu_next, true);
  de("report", self_c.report, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_9> = synth_constructor_6;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_9> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
Gamma_mix : {par_0}

mu_next : {par_1}, default={}

report : {par_2}, default={}

)DOC",
                      "par",
                      {c2py::python_typename<nrgljubljana_interface::g_w_t>(), c2py::python_typename<std::optional<double>>(),
                       c2py::python_typename<nrgljubljana_interface::dmft::mix_report_t>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_9>[] = {

   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_18 = R"DOC(The mixed Γ, before the causal projection and the rebuild of Δ)DOC";
constexpr auto _c2py_doc_member_19 = R"DOC(μ for the next iteration, if the mixer moves μ)DOC";
constexpr auto _c2py_doc_member_20 = R"DOC()DOC";
static PyObject *prop_get_dict_6(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_9> *)self)->_c);
  c2py::pydict dic;
  dic["Gamma_mix"] = self_c.Gamma_mix;
  dic["mu_next"]   = self_c.mu_next;
  dic["report"]    = self_c.report;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_9>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_9::Gamma_mix, _c2py_cls_9>("Gamma_mix", _c2py_doc_member_18),
   c2py::getsetdef_from_member<&_c2py_cls_9::mu_next, _c2py_cls_9>("mu_next", _c2py_doc_member_19),
   c2py::getsetdef_from_member<&_c2py_cls_9::report, _c2py_cls_9>("report", _c2py_doc_member_20),
   {"__dict__", (getter)prop_get_dict_6, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_9> = R"DOC(The result of mixing)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_9>;
// --------- class _c2py_cls_10 -----------
using _c2py_cls_10                                            = nrgljubljana_interface::dmft::linear_mixer;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_10>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_10> = "nrgljubljana_interface.dmft._dmft.LinearMixer";
static const auto _c2py_init_3                                = c2py::dispatcher_c_kw_t{c2py::c_constructor<_c2py_cls_10, double>("alpha")};
template <> constexpr initproc c2py::tp_init<_c2py_cls_10>    = c2py::pyfkw_constructor<_c2py_init_3>;
template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_10> = _c2py_init_3.doc(R"DOC(
Throws std::invalid_argument unless 0 < α ≤ 1
)DOC");
// alpha
static auto const _c2py_fun_19 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_10 const &self) -> decltype(auto) { return self.alpha(); }, "self")};

// description
static auto const _c2py_fun_20 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_10 const &self) -> decltype(auto) { return self.description(); }, "self")};

// mix
static auto const _c2py_fun_21 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_10 &self, const nrgljubljana_interface::dmft::mix_input_t &input) -> decltype(auto) { return self.mix(input); }, "self", "input")};

// moves_mu
static auto const _c2py_fun_22 = c2py::dispatcher_f_kw_t{
   c2py::cmethod([](_c2py_cls_10 const &self, long iteration) -> decltype(auto) { return self.moves_mu(iteration); }, "self", "iteration")};

// reset
static auto const _c2py_fun_23 = c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_10 &self) -> decltype(auto) { return self.reset(); }, "self")};

static const auto _c2py_doc_19 = _c2py_fun_19.doc(R"DOC()DOC");
static const auto _c2py_doc_20 = _c2py_fun_20.doc(R"DOC(
The type and the parameters, e.g. to check on resume that the mixer is the same
)DOC");
static const auto _c2py_doc_21 = _c2py_fun_21.doc(R"DOC(
Throws std::invalid_argument if Γ and Γ_out differ in structure or mesh, or are not finite
)DOC");
static const auto _c2py_doc_22 = _c2py_fun_22.doc(R"DOC(
Whether the mixer sets μ in this iteration (mix_result_t::mu_next); then the loop takes μ from it
)DOC");
static const auto _c2py_doc_23 = _c2py_fun_23.doc(R"DOC(
Forget the history, for a new run
)DOC");

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_10>[] = {
   {"alpha", (PyCFunction)c2py::pyfkw<_c2py_fun_19>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_19.c_str()},
   {"description", (PyCFunction)c2py::pyfkw<_c2py_fun_20>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_20.c_str()},
   {"mix", (PyCFunction)c2py::pyfkw<_c2py_fun_21>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_21.c_str()},
   {"moves_mu", (PyCFunction)c2py::pyfkw<_c2py_fun_22>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_22.c_str()},
   {"reset", (PyCFunction)c2py::pyfkw<_c2py_fun_23>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_23.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

template <>
const std::string c2py::tp_doc<_c2py_cls_10> = R"DOC(Linear mixing: Γ_mix = α Γ_out + (1 − α) Γ, element by element over all blocks)DOC"
   + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_10>;
// --------- class _c2py_cls_11 -----------
using _c2py_cls_11                                            = nrgljubljana_interface::dmft::mu_control_params_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_11>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_11> = "nrgljubljana_interface.dmft._dmft.MuControlParamsT";

static int synth_constructor_7(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(
       PyExc_RuntimeError,
       ("Error in constructing nrgljubljana_interface::dmft::mu_control_params_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_11> *)self)->_c = new _c2py_cls_11{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::mu_control_params_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_11> *)self)->_c);
  de("goal", self_c.goal, true);
  de("maxdx", self_c.maxdx, true);
  de("under", self_c.under, true);
  de("root_tol", self_c.root_tol, true);
  de("mu_tol", self_c.mu_tol, true);
  de("max_evaluations", self_c.max_evaluations, true);
  de("weight_tol", self_c.weight_tol, true);
  de("epsabs", self_c.epsabs, true);
  de("epsrel", self_c.epsrel, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_11> = synth_constructor_7;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_11> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
goal : {par_0}, default={}

maxdx : {par_1}, default=0.05

under : {par_2}, default=1.0

root_tol : {par_3}, default=1e-6

mu_tol : {par_4}, default=1e-7

max_evaluations : {par_5}, default=12

weight_tol : {par_6}, default=1e-4

epsabs : {par_7}, default=1e-10

epsrel : {par_8}, default=1e-9

)DOC",
                      "par",
                      {c2py::python_typename<std::optional<double>>(), c2py::python_typename<double>(), c2py::python_typename<double>(),
                       c2py::python_typename<double>(), c2py::python_typename<double>(), c2py::python_typename<long>(),
                       c2py::python_typename<double>(), c2py::python_typename<double>(), c2py::python_typename<double>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_11>[] = {
   {"__write_hdf5__", c2py::tpxx_write_h5<_c2py_cls_11>, METH_VARARGS, "  "},
   {"__getstate__", c2py::getstate_h5<_c2py_cls_11>, METH_NOARGS, ""},
   {"__setstate__", c2py::setstate_h5<_c2py_cls_11>, METH_O, ""},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_21 = R"DOC(The target filling (the total of occupancy); unset: μ is kept fixed)DOC";
constexpr auto _c2py_doc_member_22 = R"DOC(The largest step of μ in one iteration)DOC";
constexpr auto _c2py_doc_member_23 = R"DOC(The damping of the step to the root, 0 < under ≤ 1; the bracket is searched up to maxdx / under)DOC";
constexpr auto _c2py_doc_member_24 = R"DOC(The allowed |n − goal| of the root, and the accuracy of the root in μ)DOC";
constexpr auto _c2py_doc_member_25 = R"DOC(The allowed |n − goal| of the root, and the accuracy of the root in μ)DOC";
constexpr auto _c2py_doc_member_26 = R"DOC(Lattice evaluations per iteration, at least 4)DOC";
constexpr auto _c2py_doc_member_27 = R"DOC(The allowed |∫ A_ii dω − 1| of the lattice spectral function)DOC";
constexpr auto _c2py_doc_member_28 = R"DOC(Tolerances of the Fermi integral of the filling)DOC";
constexpr auto _c2py_doc_member_29 = R"DOC(Tolerances of the Fermi integral of the filling)DOC";
static PyObject *prop_get_dict_7(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_11> *)self)->_c);
  c2py::pydict dic;
  dic["goal"]            = self_c.goal;
  dic["maxdx"]           = self_c.maxdx;
  dic["under"]           = self_c.under;
  dic["root_tol"]        = self_c.root_tol;
  dic["mu_tol"]          = self_c.mu_tol;
  dic["max_evaluations"] = self_c.max_evaluations;
  dic["weight_tol"]      = self_c.weight_tol;
  dic["epsabs"]          = self_c.epsabs;
  dic["epsrel"]          = self_c.epsrel;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_11>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_11::goal, _c2py_cls_11>("goal", _c2py_doc_member_21),
   c2py::getsetdef_from_member<&_c2py_cls_11::maxdx, _c2py_cls_11>("maxdx", _c2py_doc_member_22),
   c2py::getsetdef_from_member<&_c2py_cls_11::under, _c2py_cls_11>("under", _c2py_doc_member_23),
   c2py::getsetdef_from_member<&_c2py_cls_11::root_tol, _c2py_cls_11>("root_tol", _c2py_doc_member_24),
   c2py::getsetdef_from_member<&_c2py_cls_11::mu_tol, _c2py_cls_11>("mu_tol", _c2py_doc_member_25),
   c2py::getsetdef_from_member<&_c2py_cls_11::max_evaluations, _c2py_cls_11>("max_evaluations", _c2py_doc_member_26),
   c2py::getsetdef_from_member<&_c2py_cls_11::weight_tol, _c2py_cls_11>("weight_tol", _c2py_doc_member_27),
   c2py::getsetdef_from_member<&_c2py_cls_11::epsabs, _c2py_cls_11>("epsabs", _c2py_doc_member_28),
   c2py::getsetdef_from_member<&_c2py_cls_11::epsrel, _c2py_cls_11>("epsrel", _c2py_doc_member_29),
   {"__dict__", (getter)prop_get_dict_7, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_11> =
   R"DOC(Parameters of mu_control_t)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_11>;
// --------- class _c2py_cls_12 -----------
using _c2py_cls_12                                            = nrgljubljana_interface::dmft::mu_input_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_12>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_12> = "nrgljubljana_interface.dmft._dmft.MuInputT";

static int synth_constructor_8(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::mu_input_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_12> *)self)->_c = new _c2py_cls_12{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::mu_input_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_12> *)self)->_c);
  de("Sigma", self_c.Sigma, false);
  de("mu", self_c.mu, true);
  de("T", self_c.T, true);
  de("degeneracy", self_c.degeneracy, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_12> = synth_constructor_8;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_12> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
Sigma : {par_0}

mu : {par_1}, default=0

T : {par_2}, default=0

degeneracy : {par_3}, default={}

)DOC",
                      "par",
                      {c2py::python_typename<nrgljubljana_interface::g_w_t>(), c2py::python_typename<double>(), c2py::python_typename<double>(),
                       c2py::python_typename<std::map<std::string, long>>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_12>[] = {

   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_30 = R"DOC(The self-energy, held fixed)DOC";
constexpr auto _c2py_doc_member_31 = R"DOC(μ of the solve)DOC";
constexpr auto _c2py_doc_member_32 = R"DOC(The temperature, and the degeneracy of each block (the manifest's))DOC";
constexpr auto _c2py_doc_member_33 = R"DOC()DOC";
static PyObject *prop_get_dict_8(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_12> *)self)->_c);
  c2py::pydict dic;
  dic["Sigma"]      = self_c.Sigma;
  dic["mu"]         = self_c.mu;
  dic["T"]          = self_c.T;
  dic["degeneracy"] = self_c.degeneracy;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_12>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_12::Sigma, _c2py_cls_12>("Sigma", _c2py_doc_member_30),
   c2py::getsetdef_from_member<&_c2py_cls_12::mu, _c2py_cls_12>("mu", _c2py_doc_member_31),
   c2py::getsetdef_from_member<&_c2py_cls_12::T, _c2py_cls_12>("T", _c2py_doc_member_32),
   c2py::getsetdef_from_member<&_c2py_cls_12::degeneracy, _c2py_cls_12>("degeneracy", _c2py_doc_member_33),
   {"__dict__", (getter)prop_get_dict_8, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_12> =
   R"DOC(What the μ control is given in one iteration)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_12>;
// --------- class _c2py_cls_13 -----------
using _c2py_cls_13                                            = nrgljubljana_interface::dmft::mu_report_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_13>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_13> = "nrgljubljana_interface.dmft._dmft.MuReportT";

static int synth_constructor_9(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::mu_report_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_13> *)self)->_c = new _c2py_cls_13{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::mu_report_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_13> *)self)->_c);
  de("status", self_c.status, true);
  de("filling", self_c.filling, true);
  de("filling_next", self_c.filling_next, true);
  de("weight_error", self_c.weight_error, true);
  de("weight_error_next", self_c.weight_error_next, true);
  de("mu_root", self_c.mu_root, true);
  de("evaluations", self_c.evaluations, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_13> = synth_constructor_9;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_13> = c2py::replace_tags(
   R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
status : {par_0}, default={}

filling : {par_1}, default=0

filling_next : {par_2}, default=0

weight_error : {par_3}, default=0

weight_error_next : {par_4}, default=0

mu_root : {par_5}, default={}

evaluations : {par_6}, default=0

)DOC",
   "par",
   {c2py::python_typename<std::string>(), c2py::python_typename<double>(), c2py::python_typename<double>(), c2py::python_typename<double>(),
    c2py::python_typename<double>(), c2py::python_typename<std::optional<double>>(), c2py::python_typename<long>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_13>[] = {
   {"__write_hdf5__", c2py::tpxx_write_h5<_c2py_cls_13>, METH_VARARGS, "  "},
   {"__getstate__", c2py::getstate_h5<_c2py_cls_13>, METH_NOARGS, ""},
   {"__setstate__", c2py::setstate_h5<_c2py_cls_13>, METH_O, ""},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_34 = R"DOC("fixed", "within_tolerance", "bracketed" or "limited")DOC";
constexpr auto _c2py_doc_member_35 = R"DOC(The filling at μ and at μ_next)DOC";
constexpr auto _c2py_doc_member_36 = R"DOC(The filling at μ and at μ_next)DOC";
constexpr auto _c2py_doc_member_37 = R"DOC(The largest |∫ A_ii dω − 1| at μ and at μ_next)DOC";
constexpr auto _c2py_doc_member_38 = R"DOC(The largest |∫ A_ii dω − 1| at μ and at μ_next)DOC";
constexpr auto _c2py_doc_member_39 = R"DOC(The root of n(μ) = goal, when bracketed)DOC";
constexpr auto _c2py_doc_member_40 = R"DOC(The number of lattice evaluations)DOC";
static PyObject *prop_get_dict_9(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_13> *)self)->_c);
  c2py::pydict dic;
  dic["status"]            = self_c.status;
  dic["filling"]           = self_c.filling;
  dic["filling_next"]      = self_c.filling_next;
  dic["weight_error"]      = self_c.weight_error;
  dic["weight_error_next"] = self_c.weight_error_next;
  dic["mu_root"]           = self_c.mu_root;
  dic["evaluations"]       = self_c.evaluations;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_13>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_13::status, _c2py_cls_13>("status", _c2py_doc_member_34),
   c2py::getsetdef_from_member<&_c2py_cls_13::filling, _c2py_cls_13>("filling", _c2py_doc_member_35),
   c2py::getsetdef_from_member<&_c2py_cls_13::filling_next, _c2py_cls_13>("filling_next", _c2py_doc_member_36),
   c2py::getsetdef_from_member<&_c2py_cls_13::weight_error, _c2py_cls_13>("weight_error", _c2py_doc_member_37),
   c2py::getsetdef_from_member<&_c2py_cls_13::weight_error_next, _c2py_cls_13>("weight_error_next", _c2py_doc_member_38),
   c2py::getsetdef_from_member<&_c2py_cls_13::mu_root, _c2py_cls_13>("mu_root", _c2py_doc_member_39),
   c2py::getsetdef_from_member<&_c2py_cls_13::evaluations, _c2py_cls_13>("evaluations", _c2py_doc_member_40),
   {"__dict__", (getter)prop_get_dict_9, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_13> =
   R"DOC(What the μ control did)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_13>;
// --------- class _c2py_cls_14 -----------
using _c2py_cls_14                                            = nrgljubljana_interface::dmft::mu_result_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_14>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_14> = "nrgljubljana_interface.dmft._dmft.MuResultT";

static int synth_constructor_10(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::mu_result_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_14> *)self)->_c = new _c2py_cls_14{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::mu_result_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_14> *)self)->_c);
  de("mu_next", self_c.mu_next, true);
  de("filling", self_c.filling, true);
  de("report", self_c.report, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_14> = synth_constructor_10;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_14> = c2py::replace_tags(
   R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
mu_next : {par_0}, default=0

filling : {par_1}, default=0

report : {par_2}, default={}

)DOC",
   "par", {c2py::python_typename<double>(), c2py::python_typename<double>(), c2py::python_typename<nrgljubljana_interface::dmft::mu_report_t>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_14>[] = {

   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_41 = R"DOC()DOC";
constexpr auto _c2py_doc_member_42 = R"DOC(The filling at μ)DOC";
constexpr auto _c2py_doc_member_43 = R"DOC()DOC";
static PyObject *prop_get_dict_10(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_14> *)self)->_c);
  c2py::pydict dic;
  dic["mu_next"] = self_c.mu_next;
  dic["filling"] = self_c.filling;
  dic["report"]  = self_c.report;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_14>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_14::mu_next, _c2py_cls_14>("mu_next", _c2py_doc_member_41),
   c2py::getsetdef_from_member<&_c2py_cls_14::filling, _c2py_cls_14>("filling", _c2py_doc_member_42),
   c2py::getsetdef_from_member<&_c2py_cls_14::report, _c2py_cls_14>("report", _c2py_doc_member_43),
   {"__dict__", (getter)prop_get_dict_10, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_14> =
   R"DOC(The result of adjusting μ)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_14>;
// --------- class _c2py_cls_15 -----------
using _c2py_cls_15                                            = nrgljubljana_interface::dmft::mu_control_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_15>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_15> = "nrgljubljana_interface.dmft._dmft.MuControlT";
static const auto _c2py_init_4 = c2py::dispatcher_c_kw_t{c2py::c_constructor<_c2py_cls_15, const nrgljubljana_interface::dmft::mu_control_params_t &>(
   "params"_a = nrgljubljana_interface::dmft::mu_control_params_t{})};
template <> constexpr initproc c2py::tp_init<_c2py_cls_15>    = c2py::pyfkw_constructor<_c2py_init_4>;
template <> const std::string c2py::tp_ctor_doc<_c2py_cls_15> = _c2py_init_4.doc(R"DOC()DOC");
// adjust
static auto const _c2py_fun_24 = c2py::dispatcher_f_kw_t{
   c2py::cmethod(c2py::cast<const nrgljubljana_interface::dmft::mu_control_t &, const nrgljubljana_interface::dmft::mu_input_t &,
                            const nrgljubljana_interface::dmft::lattice_variant_t &>(&nrgljubljana_interface::dmft::adjust),
                 "self", "input", "lattice")};

// params
static auto const _c2py_fun_25 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_15 const &self) -> decltype(auto) { return self.params(); }, "self")};

static const auto _c2py_doc_24 = _c2py_fun_24.doc(R"DOC(
μ for the next iteration, with Σ held fixed (see MuControlT)
)DOC");
static const auto _c2py_doc_25 = _c2py_fun_25.doc(R"DOC()DOC");

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_15>[] = {
   {"adjust", (PyCFunction)c2py::pyfkw<_c2py_fun_24>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_24.c_str()},
   {"params", (PyCFunction)c2py::pyfkw<_c2py_fun_25>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_25.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

template <>
const std::string c2py::tp_doc<_c2py_cls_15> =
   R"DOC(The chemical potential: fixed without a goal, otherwise adjusted towards the goal with Σ held fixed.

With a goal: the filling at μ; within root_tol of the goal μ is kept. Otherwise the filling at
μ ± maxdx / under in the direction of the error, which must not move the filling the wrong way. If the two
bracket the goal, Brent's method (GSL) to mu_tol, the residual checked against root_tol, and
μ_next = μ + clamp(under (μ* − μ), ±maxdx); otherwise μ_next = μ ± maxdx. At μ_next the filling
error must not have grown by more than root_tol. Evaluations are cached and limited to max_evaluations.

With a goal every evaluation checks that each ∫ A_ii dω is 1 within weight_tol and that the filling lies between
0 and the number of states; without one the weight error is only reported. Errors: std::invalid_argument for
invalid parameters or input, std::runtime_error for a failed search.)DOC"
   + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_15>;
// --------- class _c2py_cls_16 -----------
using _c2py_cls_16                                            = nrgljubljana_interface::dmft::convergence_params_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_16>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_16> = "nrgljubljana_interface.dmft._dmft.ConvergenceParamsT";

static int synth_constructor_11(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(
       PyExc_RuntimeError,
       ("Error in constructing nrgljubljana_interface::dmft::convergence_params_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_16> *)self)->_c = new _c2py_cls_16{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::convergence_params_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_16> *)self)->_c);
  de("eps", self_c.eps, true);
  de("window", self_c.window, true);
  de("min_iter", self_c.min_iter, true);
  de("max_iter", self_c.max_iter, true);
  de("filling_tol", self_c.filling_tol, true);
  de("mu_tol", self_c.mu_tol, true);
  de("eps_imp_loc", self_c.eps_imp_loc, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_16> = synth_constructor_11;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_16> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
eps : {par_0}, default=1e-6

window : {par_1}, default=3

min_iter : {par_2}, default=1

max_iter : {par_3}, default=50

filling_tol : {par_4}, default=1e-6

mu_tol : {par_5}, default=1e-7

eps_imp_loc : {par_6}, default={}

)DOC",
                      "par",
                      {c2py::python_typename<double>(), c2py::python_typename<long>(), c2py::python_typename<long>(), c2py::python_typename<long>(),
                       c2py::python_typename<double>(), c2py::python_typename<double>(), c2py::python_typename<std::optional<double>>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_16>[] = {
   {"__write_hdf5__", c2py::tpxx_write_h5<_c2py_cls_16>, METH_VARARGS, "  "},
   {"__getstate__", c2py::getstate_h5<_c2py_cls_16>, METH_NOARGS, ""},
   {"__setstate__", c2py::setstate_h5<_c2py_cls_16>, METH_O, ""},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_44 = R"DOC(Every distance of G_loc between consecutive iterations in the window must be below eps)DOC";
constexpr auto _c2py_doc_member_45 = R"DOC(The number of consecutive distances in the window)DOC";
constexpr auto _c2py_doc_member_46 = R"DOC(The first and the last iteration at which the loop may converge; at max_iter it stops)DOC";
constexpr auto _c2py_doc_member_47 = R"DOC(The first and the last iteration at which the loop may converge; at max_iter it stops)DOC";
constexpr auto _c2py_doc_member_48 = R"DOC(With a goal: the allowed |n(μ) − goal| and |μ_next − μ|)DOC";
constexpr auto _c2py_doc_member_49 = R"DOC(With a goal: the allowed |n(μ) − goal| and |μ_next − μ|)DOC";
constexpr auto _c2py_doc_member_50 = R"DOC(When set, the distance between G_imp and G_loc must also be at most this)DOC";
static PyObject *prop_get_dict_11(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_16> *)self)->_c);
  c2py::pydict dic;
  dic["eps"]         = self_c.eps;
  dic["window"]      = self_c.window;
  dic["min_iter"]    = self_c.min_iter;
  dic["max_iter"]    = self_c.max_iter;
  dic["filling_tol"] = self_c.filling_tol;
  dic["mu_tol"]      = self_c.mu_tol;
  dic["eps_imp_loc"] = self_c.eps_imp_loc;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_16>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_16::eps, _c2py_cls_16>("eps", _c2py_doc_member_44),
   c2py::getsetdef_from_member<&_c2py_cls_16::window, _c2py_cls_16>("window", _c2py_doc_member_45),
   c2py::getsetdef_from_member<&_c2py_cls_16::min_iter, _c2py_cls_16>("min_iter", _c2py_doc_member_46),
   c2py::getsetdef_from_member<&_c2py_cls_16::max_iter, _c2py_cls_16>("max_iter", _c2py_doc_member_47),
   c2py::getsetdef_from_member<&_c2py_cls_16::filling_tol, _c2py_cls_16>("filling_tol", _c2py_doc_member_48),
   c2py::getsetdef_from_member<&_c2py_cls_16::mu_tol, _c2py_cls_16>("mu_tol", _c2py_doc_member_49),
   c2py::getsetdef_from_member<&_c2py_cls_16::eps_imp_loc, _c2py_cls_16>("eps_imp_loc", _c2py_doc_member_50),
   {"__dict__", (getter)prop_get_dict_11, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_16> =
   R"DOC(Parameters of convergence_t)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_16>;
// --------- class _c2py_cls_17 -----------
using _c2py_cls_17                                            = nrgljubljana_interface::dmft::convergence_input_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_17>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_17> = "nrgljubljana_interface.dmft._dmft.ConvergenceInputT";

static int synth_constructor_12(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(
       PyExc_RuntimeError,
       ("Error in constructing nrgljubljana_interface::dmft::convergence_input_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_17> *)self)->_c = new _c2py_cls_17{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::convergence_input_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_17> *)self)->_c);
  de("iteration", self_c.iteration, true);
  de("G_loc", self_c.G_loc, false);
  de("G_loc_previous", self_c.G_loc_previous, true);
  de("G_imp", self_c.G_imp, false);
  de("filling", self_c.filling, true);
  de("goal", self_c.goal, true);
  de("mu_step", self_c.mu_step, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_17> = synth_constructor_12;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_17> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
G_loc : {par_0}

G_imp : {par_1}

iteration : {par_2}, default=0

G_loc_previous : {par_3}, default={}

filling : {par_4}, default=0

goal : {par_5}, default={}

mu_step : {par_6}, default=0

)DOC",
                      "par",
                      {c2py::python_typename<nrgljubljana_interface::g_w_t>(), c2py::python_typename<nrgljubljana_interface::g_w_t>(),
                       c2py::python_typename<long>(), c2py::python_typename<std::optional<nrgljubljana_interface::g_w_t>>(),
                       c2py::python_typename<double>(), c2py::python_typename<std::optional<double>>(), c2py::python_typename<double>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_17>[] = {

   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_51 = R"DOC(The iteration, from 1; consecutive)DOC";
constexpr auto _c2py_doc_member_52 = R"DOC(G_loc at μ, and that of the previous iteration (unset in the first))DOC";
constexpr auto _c2py_doc_member_53 = R"DOC()DOC";
constexpr auto _c2py_doc_member_54 = R"DOC(The impurity Green's function)DOC";
constexpr auto _c2py_doc_member_55 = R"DOC(The filling at μ; the goal when μ is controlled; μ_next − μ)DOC";
constexpr auto _c2py_doc_member_56 = R"DOC()DOC";
constexpr auto _c2py_doc_member_57 = R"DOC()DOC";
static PyObject *prop_get_dict_12(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_17> *)self)->_c);
  c2py::pydict dic;
  dic["iteration"]      = self_c.iteration;
  dic["G_loc"]          = self_c.G_loc;
  dic["G_loc_previous"] = self_c.G_loc_previous;
  dic["G_imp"]          = self_c.G_imp;
  dic["filling"]        = self_c.filling;
  dic["goal"]           = self_c.goal;
  dic["mu_step"]        = self_c.mu_step;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_17>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_17::iteration, _c2py_cls_17>("iteration", _c2py_doc_member_51),
   c2py::getsetdef_from_member<&_c2py_cls_17::G_loc, _c2py_cls_17>("G_loc", _c2py_doc_member_52),
   c2py::getsetdef_from_member<&_c2py_cls_17::G_loc_previous, _c2py_cls_17>("G_loc_previous", _c2py_doc_member_53),
   c2py::getsetdef_from_member<&_c2py_cls_17::G_imp, _c2py_cls_17>("G_imp", _c2py_doc_member_54),
   c2py::getsetdef_from_member<&_c2py_cls_17::filling, _c2py_cls_17>("filling", _c2py_doc_member_55),
   c2py::getsetdef_from_member<&_c2py_cls_17::goal, _c2py_cls_17>("goal", _c2py_doc_member_56),
   c2py::getsetdef_from_member<&_c2py_cls_17::mu_step, _c2py_cls_17>("mu_step", _c2py_doc_member_57),
   {"__dict__", (getter)prop_get_dict_12, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_17> =
   R"DOC(What the convergence test is given in one iteration)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_17>;
// --------- class _c2py_cls_18 -----------
using _c2py_cls_18                                            = nrgljubljana_interface::dmft::convergence_result_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_18>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_18> = "nrgljubljana_interface.dmft._dmft.ConvergenceResultT";

static int synth_constructor_13(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(
       PyExc_RuntimeError,
       ("Error in constructing nrgljubljana_interface::dmft::convergence_result_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_18> *)self)->_c = new _c2py_cls_18{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::convergence_result_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_18> *)self)->_c);
  de("decision", self_c.decision, true);
  de("distance", self_c.distance, true);
  de("window_max", self_c.window_max, true);
  de("window_samples", self_c.window_samples, true);
  de("distance_imp_loc", self_c.distance_imp_loc, true);
  de("window_ok", self_c.window_ok, true);
  de("iterations_ok", self_c.iterations_ok, true);
  de("filling_ok", self_c.filling_ok, true);
  de("mu_ok", self_c.mu_ok, true);
  de("imp_loc_ok", self_c.imp_loc_ok, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_18> = synth_constructor_13;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_18> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
decision : {par_0}, default={}

distance : {par_1}, default={}

window_max : {par_2}, default=0

window_samples : {par_3}, default=0

distance_imp_loc : {par_4}, default=0

window_ok : {par_5}, default=false

iterations_ok : {par_6}, default=false

filling_ok : {par_7}, default=false

mu_ok : {par_8}, default=false

imp_loc_ok : {par_9}, default=false

)DOC",
                      "par",
                      {c2py::python_typename<std::string>(), c2py::python_typename<std::optional<double>>(), c2py::python_typename<double>(),
                       c2py::python_typename<long>(), c2py::python_typename<double>(), c2py::python_typename<bool>(), c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(), c2py::python_typename<bool>(), c2py::python_typename<bool>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_18>[] = {
   {"__write_hdf5__", c2py::tpxx_write_h5<_c2py_cls_18>, METH_VARARGS, "  "},
   {"__getstate__", c2py::getstate_h5<_c2py_cls_18>, METH_NOARGS, ""},
   {"__setstate__", c2py::setstate_h5<_c2py_cls_18>, METH_O, ""},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_58 = R"DOC("continue", "converged" or "stop")DOC";
constexpr auto _c2py_doc_member_59 = R"DOC(The largest spectral distance of G_loc to the previous iteration (unset in the first))DOC";
constexpr auto _c2py_doc_member_60 = R"DOC(The largest distance in the window, and the number of distances in it)DOC";
constexpr auto _c2py_doc_member_61 = R"DOC()DOC";
constexpr auto _c2py_doc_member_62 = R"DOC(The largest spectral distance between G_imp and G_loc)DOC";
constexpr auto _c2py_doc_member_63 = R"DOC(Which criteria hold (an inactive criterion holds))DOC";
constexpr auto _c2py_doc_member_64 = R"DOC(Which criteria hold (an inactive criterion holds))DOC";
constexpr auto _c2py_doc_member_65 = R"DOC(Which criteria hold (an inactive criterion holds))DOC";
constexpr auto _c2py_doc_member_66 = R"DOC(Which criteria hold (an inactive criterion holds))DOC";
constexpr auto _c2py_doc_member_67 = R"DOC(Which criteria hold (an inactive criterion holds))DOC";
static PyObject *prop_get_dict_13(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_18> *)self)->_c);
  c2py::pydict dic;
  dic["decision"]         = self_c.decision;
  dic["distance"]         = self_c.distance;
  dic["window_max"]       = self_c.window_max;
  dic["window_samples"]   = self_c.window_samples;
  dic["distance_imp_loc"] = self_c.distance_imp_loc;
  dic["window_ok"]        = self_c.window_ok;
  dic["iterations_ok"]    = self_c.iterations_ok;
  dic["filling_ok"]       = self_c.filling_ok;
  dic["mu_ok"]            = self_c.mu_ok;
  dic["imp_loc_ok"]       = self_c.imp_loc_ok;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_18>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_18::decision, _c2py_cls_18>("decision", _c2py_doc_member_58),
   c2py::getsetdef_from_member<&_c2py_cls_18::distance, _c2py_cls_18>("distance", _c2py_doc_member_59),
   c2py::getsetdef_from_member<&_c2py_cls_18::window_max, _c2py_cls_18>("window_max", _c2py_doc_member_60),
   c2py::getsetdef_from_member<&_c2py_cls_18::window_samples, _c2py_cls_18>("window_samples", _c2py_doc_member_61),
   c2py::getsetdef_from_member<&_c2py_cls_18::distance_imp_loc, _c2py_cls_18>("distance_imp_loc", _c2py_doc_member_62),
   c2py::getsetdef_from_member<&_c2py_cls_18::window_ok, _c2py_cls_18>("window_ok", _c2py_doc_member_63),
   c2py::getsetdef_from_member<&_c2py_cls_18::iterations_ok, _c2py_cls_18>("iterations_ok", _c2py_doc_member_64),
   c2py::getsetdef_from_member<&_c2py_cls_18::filling_ok, _c2py_cls_18>("filling_ok", _c2py_doc_member_65),
   c2py::getsetdef_from_member<&_c2py_cls_18::mu_ok, _c2py_cls_18>("mu_ok", _c2py_doc_member_66),
   c2py::getsetdef_from_member<&_c2py_cls_18::imp_loc_ok, _c2py_cls_18>("imp_loc_ok", _c2py_doc_member_67),
   {"__dict__", (getter)prop_get_dict_13, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_18> =
   R"DOC(The decision and the measures behind it)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_18>;
// --------- class _c2py_cls_19 -----------
using _c2py_cls_19                                            = nrgljubljana_interface::dmft::convergence_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_19>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_19> = "nrgljubljana_interface.dmft._dmft.ConvergenceT";
static const auto _c2py_init_5 =
   c2py::dispatcher_c_kw_t{c2py::c_constructor<_c2py_cls_19, const nrgljubljana_interface::dmft::convergence_params_t &>(
      "params"_a = nrgljubljana_interface::dmft::convergence_params_t{})};
template <> constexpr initproc c2py::tp_init<_c2py_cls_19>    = c2py::pyfkw_constructor<_c2py_init_5>;
template <> const std::string c2py::tp_ctor_doc<_c2py_cls_19> = _c2py_init_5.doc(R"DOC()DOC");
// decide
static auto const _c2py_fun_26 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_19 &self, const nrgljubljana_interface::dmft::convergence_input_t &input) -> decltype(auto) { return self.decide(input); }, "self",
   "input")};

// history
static auto const _c2py_fun_27 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_19 const &self) -> decltype(auto) { return self.history(); }, "self")};

// params
static auto const _c2py_fun_28 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_19 const &self) -> decltype(auto) { return self.params(); }, "self")};

// reset
static auto const _c2py_fun_29 = c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_19 &self) -> decltype(auto) { return self.reset(); }, "self")};

static const auto _c2py_doc_26 = _c2py_fun_26.doc(R"DOC()DOC");
static const auto _c2py_doc_27 = _c2py_fun_27.doc(R"DOC(
(iteration, distance) of every iteration after the first
)DOC");
static const auto _c2py_doc_28 = _c2py_fun_28.doc(R"DOC()DOC");
static const auto _c2py_doc_29 = _c2py_fun_29.doc(R"DOC(
Forget the history, for a new run
)DOC");

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_19>[] = {
   {"decide", (PyCFunction)c2py::pyfkw<_c2py_fun_26>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_26.c_str()},
   {"history", (PyCFunction)c2py::pyfkw<_c2py_fun_27>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_27.c_str()},
   {"params", (PyCFunction)c2py::pyfkw<_c2py_fun_28>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_28.c_str()},
   {"reset", (PyCFunction)c2py::pyfkw<_c2py_fun_29>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_29.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

template <>
const std::string c2py::tp_doc<_c2py_cls_19> =
   R"DOC(The convergence test. converged: a full window of distances all strictly below eps, min_iter ≤ n ≤ max_iter,
with a goal |n − goal| ≤ filling_tol and |μ_next − μ| ≤ mu_tol, and with eps_imp_loc the G_imp–G_loc
distance at most that; stop: not converged at max_iter; continue otherwise. The state is the history of distances.

Throws std::invalid_argument for invalid parameters, a non-consecutive iteration, a missing previous G_loc after
the first iteration, and Green's functions that differ in structure or mesh.)DOC"
   + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_19>;
// --------- class _c2py_cls_20 -----------
using _c2py_cls_20                                            = nrgljubljana_interface::dmft::loop_params_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_20>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_20> = "nrgljubljana_interface.dmft._dmft.LoopParamsT";

static int synth_constructor_14(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::loop_params_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_20> *)self)->_c = new _c2py_cls_20{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::loop_params_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_20> *)self)->_c);
  de("solve_params", self_c.solve_params, true);
  de("options", self_c.options, true);
  de("Gamma_floor", self_c.Gamma_floor, true);
  de("sigma_floor", self_c.sigma_floor, true);
  de("run_dir", self_c.run_dir, true);
  de("keep_history", self_c.keep_history, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_20> = synth_constructor_14;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_20> = c2py::replace_tags(
   R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
solve_params : {par_0}, default={}

options : {par_1}, default={}

Gamma_floor : {par_2}, default=1e-6

sigma_floor : {par_3}, default=1e-12

run_dir : {par_4}, default={}

keep_history : {par_5}, default=true

)DOC",
   "par",
   {c2py::python_typename<nrgljubljana_interface::solve_params_t>(), c2py::python_typename<nrgljubljana_interface::options_t>(),
    c2py::python_typename<double>(), c2py::python_typename<double>(), c2py::python_typename<std::string>(), c2py::python_typename<bool>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_20>[] = {
   {"__write_hdf5__", c2py::tpxx_write_h5<_c2py_cls_20>, METH_VARARGS, "  "},
   {"__getstate__", c2py::getstate_h5<_c2py_cls_20>, METH_NOARGS, ""},
   {"__setstate__", c2py::setstate_h5<_c2py_cls_20>, METH_O, ""},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_68 = R"DOC(Solve parameters (T, Λ, Nz, the model parameters except the impurity levels))DOC";
constexpr auto _c2py_doc_member_69 = R"DOC(What each solve computes; must give G, F_l, F_r and I broadened, for the full set of z)DOC";
constexpr auto _c2py_doc_member_70 = R"DOC(Floors of the spectral parts of Δ and Σ)DOC";
constexpr auto _c2py_doc_member_71 = R"DOC(Floors of the spectral parts of Δ and Σ)DOC";
constexpr auto _c2py_doc_member_72 = R"DOC(Directory of the files of the run (rank 0): checkpoint.h5, history.h5, stats.dat; empty: no files)DOC";
constexpr auto _c2py_doc_member_73 = R"DOC(Append each iteration to history.h5)DOC";
static PyObject *prop_get_dict_14(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_20> *)self)->_c);
  c2py::pydict dic;
  dic["solve_params"] = self_c.solve_params;
  dic["options"]      = self_c.options;
  dic["Gamma_floor"]  = self_c.Gamma_floor;
  dic["sigma_floor"]  = self_c.sigma_floor;
  dic["run_dir"]      = self_c.run_dir;
  dic["keep_history"] = self_c.keep_history;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_20>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_20::solve_params, _c2py_cls_20>("solve_params", _c2py_doc_member_68),
   c2py::getsetdef_from_member<&_c2py_cls_20::options, _c2py_cls_20>("options", _c2py_doc_member_69),
   c2py::getsetdef_from_member<&_c2py_cls_20::Gamma_floor, _c2py_cls_20>("Gamma_floor", _c2py_doc_member_70),
   c2py::getsetdef_from_member<&_c2py_cls_20::sigma_floor, _c2py_cls_20>("sigma_floor", _c2py_doc_member_71),
   c2py::getsetdef_from_member<&_c2py_cls_20::run_dir, _c2py_cls_20>("run_dir", _c2py_doc_member_72),
   c2py::getsetdef_from_member<&_c2py_cls_20::keep_history, _c2py_cls_20>("keep_history", _c2py_doc_member_73),
   {"__dict__", (getter)prop_get_dict_14, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_20> =
   R"DOC(Parameters of the loop)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_20>;
// --------- class _c2py_cls_21 -----------
using _c2py_cls_21                                            = nrgljubljana_interface::dmft::step_result_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_21>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_21> = "nrgljubljana_interface.dmft._dmft.StepResultT";

static int synth_constructor_15(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::step_result_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_21> *)self)->_c = new _c2py_cls_21{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::dmft::step_result_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_21> *)self)->_c);
  de("iteration", self_c.iteration, true);
  de("mu", self_c.mu, true);
  de("mu_next", self_c.mu_next, true);
  de("filling", self_c.filling, true);
  de("sigma_projection", self_c.sigma_projection, true);
  de("delta_projection", self_c.delta_projection, true);
  de("mu_report", self_c.mu_report, true);
  de("mix_report", self_c.mix_report, true);
  de("convergence", self_c.convergence, true);
  de("expv", self_c.expv, true);
  de("kept_delta", self_c.kept_delta, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_21> = synth_constructor_15;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_21> = c2py::replace_tags(
   R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
iteration : {par_0}, default=0

mu : {par_1}, default=0

mu_next : {par_2}, default=0

filling : {par_3}, default=0

sigma_projection : {par_4}, default={}

delta_projection : {par_5}, default={}

mu_report : {par_6}, default={}

mix_report : {par_7}, default={}

convergence : {par_8}, default={}

expv : {par_9}, default={}

kept_delta : {par_10}, default=false

)DOC",
   "par",
   {c2py::python_typename<long>(), c2py::python_typename<double>(), c2py::python_typename<double>(), c2py::python_typename<double>(),
    c2py::python_typename<nrgljubljana_interface::dmft::causality_report_t>(),
    c2py::python_typename<nrgljubljana_interface::dmft::causality_report_t>(), c2py::python_typename<nrgljubljana_interface::dmft::mu_report_t>(),
    c2py::python_typename<nrgljubljana_interface::dmft::mix_report_t>(), c2py::python_typename<nrgljubljana_interface::dmft::convergence_result_t>(),
    c2py::python_typename<std::map<std::string, double>>(), c2py::python_typename<bool>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_21>[] = {
   {"__write_hdf5__", c2py::tpxx_write_h5<_c2py_cls_21>, METH_VARARGS, "  "},
   {"__getstate__", c2py::getstate_h5<_c2py_cls_21>, METH_NOARGS, ""},
   {"__setstate__", c2py::setstate_h5<_c2py_cls_21>, METH_O, ""},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_74 = R"DOC()DOC";
constexpr auto _c2py_doc_member_75 = R"DOC(μ of the solve, μ for the next iteration, and the filling at μ)DOC";
constexpr auto _c2py_doc_member_76 = R"DOC(μ of the solve, μ for the next iteration, and the filling at μ)DOC";
constexpr auto _c2py_doc_member_77 = R"DOC(μ of the solve, μ for the next iteration, and the filling at μ)DOC";
constexpr auto _c2py_doc_member_78 = R"DOC(The causal projections of Σ and of the next Δ)DOC";
constexpr auto _c2py_doc_member_79 = R"DOC()DOC";
constexpr auto _c2py_doc_member_80 = R"DOC()DOC";
constexpr auto _c2py_doc_member_81 = R"DOC()DOC";
constexpr auto _c2py_doc_member_82 = R"DOC()DOC";
constexpr auto _c2py_doc_member_83 = R"DOC(Expectation values of the solve)DOC";
constexpr auto _c2py_doc_member_84 = R"DOC(The next Δ could not be prepared after a terminal decision; Δ is kept)DOC";
static PyObject *prop_get_dict_15(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_21> *)self)->_c);
  c2py::pydict dic;
  dic["iteration"]        = self_c.iteration;
  dic["mu"]               = self_c.mu;
  dic["mu_next"]          = self_c.mu_next;
  dic["filling"]          = self_c.filling;
  dic["sigma_projection"] = self_c.sigma_projection;
  dic["delta_projection"] = self_c.delta_projection;
  dic["mu_report"]        = self_c.mu_report;
  dic["mix_report"]       = self_c.mix_report;
  dic["convergence"]      = self_c.convergence;
  dic["expv"]             = self_c.expv;
  dic["kept_delta"]       = self_c.kept_delta;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_21>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_21::iteration, _c2py_cls_21>("iteration", _c2py_doc_member_74),
   c2py::getsetdef_from_member<&_c2py_cls_21::mu, _c2py_cls_21>("mu", _c2py_doc_member_75),
   c2py::getsetdef_from_member<&_c2py_cls_21::mu_next, _c2py_cls_21>("mu_next", _c2py_doc_member_76),
   c2py::getsetdef_from_member<&_c2py_cls_21::filling, _c2py_cls_21>("filling", _c2py_doc_member_77),
   c2py::getsetdef_from_member<&_c2py_cls_21::sigma_projection, _c2py_cls_21>("sigma_projection", _c2py_doc_member_78),
   c2py::getsetdef_from_member<&_c2py_cls_21::delta_projection, _c2py_cls_21>("delta_projection", _c2py_doc_member_79),
   c2py::getsetdef_from_member<&_c2py_cls_21::mu_report, _c2py_cls_21>("mu_report", _c2py_doc_member_80),
   c2py::getsetdef_from_member<&_c2py_cls_21::mix_report, _c2py_cls_21>("mix_report", _c2py_doc_member_81),
   c2py::getsetdef_from_member<&_c2py_cls_21::convergence, _c2py_cls_21>("convergence", _c2py_doc_member_82),
   c2py::getsetdef_from_member<&_c2py_cls_21::expv, _c2py_cls_21>("expv", _c2py_doc_member_83),
   c2py::getsetdef_from_member<&_c2py_cls_21::kept_delta, _c2py_cls_21>("kept_delta", _c2py_doc_member_84),
   {"__dict__", (getter)prop_get_dict_15, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_21> =
   R"DOC(One iteration of the loop)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_21>;
// --------- class _c2py_cls_22 -----------
using _c2py_cls_22                                            = nrgljubljana_interface::dmft::py_dmft_loop;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_22>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_22> = "nrgljubljana_interface.dmft._dmft.PyDmftLoop";
static const auto _c2py_init_6                                = c2py::dispatcher_c_kw_t{
   c2py::c_constructor<_c2py_cls_22, nrgljubljana_interface::solver_core &, const nrgljubljana_interface::dmft::lattice_variant_t &,
                                                      const nrgljubljana_interface::dmft::mixer_variant_t &, nrgljubljana_interface::dmft::convergence_t &,
                                                      nrgljubljana_interface::dmft::loop_params_t, const std::optional<nrgljubljana_interface::dmft::mu_control_t> &>(
      "solver", "lattice", "mixer", "convergence", "params", "mu_control"_a = std::optional<nrgljubljana_interface::dmft::mu_control_t>{})};
template <> constexpr initproc c2py::tp_init<_c2py_cls_22>    = c2py::pyfkw_constructor<_c2py_init_6>;
template <> const std::string c2py::tp_ctor_doc<_c2py_cls_22> = _c2py_init_6.doc(R"DOC()DOC");
// Delta
static auto const _c2py_fun_30 = c2py::dispatcher_f_kw_t{
   c2py::cmethod(c2py::cast<const nrgljubljana_interface::dmft::py_dmft_loop &>(&nrgljubljana_interface::dmft::Delta), "self")};

// G_loc
static auto const _c2py_fun_31 = c2py::dispatcher_f_kw_t{
   c2py::cmethod(c2py::cast<const nrgljubljana_interface::dmft::py_dmft_loop &>(&nrgljubljana_interface::dmft::G_loc), "self")};

// Sigma
static auto const _c2py_fun_32 = c2py::dispatcher_f_kw_t{
   c2py::cmethod(c2py::cast<const nrgljubljana_interface::dmft::py_dmft_loop &>(&nrgljubljana_interface::dmft::Sigma), "self")};

// initialize
static auto const _c2py_fun_33 = c2py::dispatcher_f_kw_t{
   c2py::cmethod(
      [](_c2py_cls_22 &self, const nrgljubljana_interface::g_w_t &Sigma0, double mu0) -> decltype(auto) { return self.initialize(Sigma0, mu0); },
      "self", "Sigma0", "mu0"),
   c2py::cmethod([](_c2py_cls_22 &self, std::complex<double> Sigma0, double mu0) -> decltype(auto) { return self.initialize(Sigma0, mu0); }, "self",
                 "Sigma0", "mu0")};

// initialize_from
static auto const _c2py_fun_34 = c2py::dispatcher_f_kw_t{
   c2py::cmethod([](_c2py_cls_22 &self, const std::string &run_dir) -> decltype(auto) { return self.initialize_from(run_dir); }, "self", "run_dir")};

// iteration
static auto const _c2py_fun_35 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_22 const &self) -> decltype(auto) { return self.iteration(); }, "self")};

// last_step
static auto const _c2py_fun_36 = c2py::dispatcher_f_kw_t{
   c2py::cmethod(c2py::cast<const nrgljubljana_interface::dmft::py_dmft_loop &>(&nrgljubljana_interface::dmft::last_step), "self")};

// mu
static auto const _c2py_fun_37 = c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_22 const &self) -> decltype(auto) { return self.mu(); }, "self")};

// onsite_energy
static auto const _c2py_fun_38 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_22 const &self) -> decltype(auto) { return self.onsite_energy(); }, "self")};

// params
static auto const _c2py_fun_39 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_22 const &self) -> decltype(auto) { return self.params(); }, "self")};

// resume
static auto const _c2py_fun_40 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_22 &self, const std::string &run_dir, bool allow_changes) -> decltype(auto) { return self.resume(run_dir, allow_changes); }, "self",
   "run_dir", "allow_changes"_a = false)};

// run
static auto const _c2py_fun_41 = c2py::dispatcher_f_kw_t{
   c2py::cmethod([](_c2py_cls_22 &self, long max_steps) -> decltype(auto) { return self.run(max_steps); }, "self", "max_steps"_a = 0)};

// set_callback
static auto const _c2py_fun_42 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   c2py::cast<nrgljubljana_interface::dmft::py_dmft_loop &, std::function<c2py::pyref(const nrgljubljana_interface::dmft::step_result_t &)>>(
      &nrgljubljana_interface::dmft::set_callback),
   "self", "callback")};

// step
static auto const _c2py_fun_43 = c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_22 &self) -> decltype(auto) { return self.step(); }, "self")};

static const auto _c2py_doc_30 = _c2py_fun_30.doc(R"DOC(
Δ of the next step; throws before initialize
)DOC");
static const auto _c2py_doc_31 = _c2py_fun_31.doc(R"DOC(
G_loc of the last step; None before the first
)DOC");
static const auto _c2py_doc_32 = _c2py_fun_32.doc(R"DOC(
Σ of the last step; None before the first
)DOC");
static const auto _c2py_doc_33 = _c2py_fun_33.doc(R"DOC(
[1] Δ from the lattice at Σ0 and μ0, through hybridization_from_Gamma; Σ0 in the structure and mesh of G. The loop
controls are reset, and a history.h5 in the run directory is removed.

------

[2] As initialize, with Σ0 times the unit matrix in every block and at every frequency

------
)DOC");
static const auto _c2py_doc_34 = _c2py_fun_34.doc(R"DOC(
Start a new run from the Δ_mix and μ_next of the checkpoint in run_dir: iteration 0, the loop controls reset
)DOC");
static const auto _c2py_doc_35 = _c2py_fun_35.doc(R"DOC()DOC");
static const auto _c2py_doc_36 = _c2py_fun_36.doc(R"DOC(
The result of the last step; None before the first
)DOC");
static const auto _c2py_doc_37 = _c2py_fun_37.doc(R"DOC()DOC");
static const auto _c2py_doc_38 = _c2py_fun_38.doc(R"DOC()DOC");
static const auto _c2py_doc_39 = _c2py_fun_39.doc(R"DOC()DOC");
static const auto _c2py_doc_40 = _c2py_fun_40.doc(R"DOC(
Continue the run whose checkpoint.h5 is in run_dir: the next step is the iteration after it, with its Δ_mix and
μ_next, and the states of the mixer and the convergence test restored. The stored parameters must equal the
current ones (solve parameters, options, floors, μ-control and convergence parameters, the descriptions of the
mixer and the lattice, the manifest, the mesh); otherwise std::invalid_argument names each difference, unless
allow_changes, which continues with the current parameters.
)DOC");
static const auto _c2py_doc_41 = _c2py_fun_41.doc(R"DOC(
Steps until the decision is converged or stop, or max_steps (0: no limit) steps were taken; the last decision
)DOC");
static const auto _c2py_doc_42 = _c2py_fun_42.doc(R"DOC(
Call the function with the StepResult of each step, after it is committed; what it returns is ignored
)DOC");
static const auto _c2py_doc_43 = _c2py_fun_43.doc(R"DOC(
One iteration; throws std::logic_error before initialize
)DOC");

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_22>[] = {
   {"Delta", (PyCFunction)c2py::pyfkw<_c2py_fun_30>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_30.c_str()},
   {"G_loc", (PyCFunction)c2py::pyfkw<_c2py_fun_31>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_31.c_str()},
   {"Sigma", (PyCFunction)c2py::pyfkw<_c2py_fun_32>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_32.c_str()},
   {"initialize", (PyCFunction)c2py::pyfkw<_c2py_fun_33>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_33.c_str()},
   {"initialize_from", (PyCFunction)c2py::pyfkw<_c2py_fun_34>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_34.c_str()},
   {"iteration", (PyCFunction)c2py::pyfkw<_c2py_fun_35>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_35.c_str()},
   {"last_step", (PyCFunction)c2py::pyfkw<_c2py_fun_36>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_36.c_str()},
   {"mu", (PyCFunction)c2py::pyfkw<_c2py_fun_37>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_37.c_str()},
   {"onsite_energy", (PyCFunction)c2py::pyfkw<_c2py_fun_38>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_38.c_str()},
   {"params", (PyCFunction)c2py::pyfkw<_c2py_fun_39>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_39.c_str()},
   {"resume", (PyCFunction)c2py::pyfkw<_c2py_fun_40>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_40.c_str()},
   {"run", (PyCFunction)c2py::pyfkw<_c2py_fun_41>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_41.c_str()},
   {"set_callback", (PyCFunction)c2py::pyfkw<_c2py_fun_42>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_42.c_str()},
   {"step", (PyCFunction)c2py::pyfkw<_c2py_fun_43>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_43.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

template <>
const std::string c2py::tp_doc<_c2py_cls_22> =
   R"DOC(The DMFT loop with the solver of the package; the loop holds copies of the lattice and the mixer)DOC" + std::string{"\n\n----------\n\n"}
   + c2py::tp_ctor_doc<_c2py_cls_22>;

// ==================== module functions ====================

// causal_projection
static auto const _c2py_fun_44 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](const nrgljubljana_interface::g_w_t &F, double floor) { return nrgljubljana_interface::dmft::causal_projection(F, floor); }, "F", "floor")};

// high_frequency_limit
static auto const _c2py_fun_45 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](const nrgljubljana_interface::g_w_t &F) { return nrgljubljana_interface::dmft::high_frequency_limit(F); }, "F")};

// hilbert_transform
static auto const _c2py_fun_46 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](const nrgljubljana_interface::dmft::dos_t &dos, std::complex<double> z,
                 long moment) { return nrgljubljana_interface::dmft::hilbert_transform(dos, z, moment); },
              "dos", "z", "moment"_a = 0),
   c2py::cfun(
      [](const nrgljubljana_interface::dmft::dos_t &dos,
         const nda::basic_array<std::complex<double>, 1, nda::C_layout, 'A', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>> &z,
         long moment) { return nrgljubljana_interface::dmft::hilbert_transform(dos, z, moment); },
      "dos", "z", "moment"_a = 0),
   c2py::cfun(
      [](const nrgljubljana_interface::dmft::dos_t &dos,
         const nda::basic_array<std::complex<double>, 2, nda::C_layout, 'M', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>> &z,
         long moment) { return nrgljubljana_interface::dmft::hilbert_transform(dos, z, moment); },
      "dos", "z", "moment"_a = 0)};

// hybridization_from_Gamma
static auto const _c2py_fun_47 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](const nrgljubljana_interface::g_w_t &Gamma, const nrgljubljana_interface::dmft::block_matrices_t &onsite_energy,
                 double floor) { return nrgljubljana_interface::dmft::hybridization_from_Gamma(Gamma, onsite_energy, floor); },
              "Gamma", "onsite_energy", "floor")};

// impurity_gf
static auto const _c2py_fun_48 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](const nrgljubljana_interface::g_w_t &Delta, const nrgljubljana_interface::g_w_t &Sigma,
      const nrgljubljana_interface::dmft::block_matrices_t &h_imp) { return nrgljubljana_interface::dmft::impurity_gf(Delta, Sigma, h_imp); },
   "Delta", "Sigma", "h_imp")};

// occupancy
static auto const _c2py_fun_49 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](const nrgljubljana_interface::g_w_t &G, double T, const std::map<std::string, long> &degeneracy, double epsabs,
                 double epsrel) { return nrgljubljana_interface::dmft::occupancy(G, T, degeneracy, epsabs, epsrel); },
              "G", "T", "degeneracy", "epsabs"_a = 1e-10, "epsrel"_a = 1e-9)};

// rebuild_hermitian_part
static auto const _c2py_fun_50 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](const nrgljubljana_interface::g_w_t &F, const nrgljubljana_interface::dmft::block_matrices_t &high_frequency_limit) {
     return nrgljubljana_interface::dmft::rebuild_hermitian_part(F, high_frequency_limit);
   },
   "F", "high_frequency_limit")};

// sigma_FG
static auto const _c2py_fun_51 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](const nrgljubljana_interface::g_w_t &G, const nrgljubljana_interface::g_w_t &F_l) { return nrgljubljana_interface::dmft::sigma_FG(G, F_l); },
   "G", "F_l")};

// sigma_IFG
static auto const _c2py_fun_52 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](const nrgljubljana_interface::g_w_t &G, const nrgljubljana_interface::g_w_t &F_l, const nrgljubljana_interface::g_w_t &F_r,
      const nrgljubljana_interface::g_w_t &I, const nrgljubljana_interface::dmft::block_matrices_t &sigma_hartree) {
     return nrgljubljana_interface::dmft::sigma_IFG(G, F_l, F_r, I, sigma_hartree);
   },
   "G", "F_l", "F_r", "I", "sigma_hartree")};

// spectral_distance
static auto const _c2py_fun_53 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](const nrgljubljana_interface::g_w_t &a,
                                         const nrgljubljana_interface::g_w_t &b) { return nrgljubljana_interface::dmft::spectral_distance(a, b); },
                                      "a", "b")};

// spectral_part
static auto const _c2py_fun_54 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](const nrgljubljana_interface::g_w_t &F) { return nrgljubljana_interface::dmft::spectral_part(F); }, "F")};

// spectral_weight
static auto const _c2py_fun_55 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](const nrgljubljana_interface::g_w_t &G) { return nrgljubljana_interface::dmft::spectral_weight(G); }, "G")};

static const auto _c2py_doc_44 = _c2py_fun_44.doc(R"DOC(
Raise the eigenvalues of the spectral part below floor to it, keeping the eigenvectors (S' = S + δS), and
correct the Hermitian part by the Kramers-Kronig transform of the change: F' = H + KK(−δS) − iS'.

For an analytic F the result is the analytic function with spectral part S' and the same F_∞, which is not
needed. An eigenvalue within 1e-14 times the largest |F_ij(ω)| of the floor counts as on it, so a second
projection changes nothing. A block in which nothing is raised is returned unchanged.

Throws std::invalid_argument for a floor that is negative or not finite, and for a non-finite F (naming the
block, the element and the frequency).
)DOC");
static const auto _c2py_doc_45 = _c2py_fun_45.doc(R"DOC(
Estimate of the high-frequency limit F_∞: H − KK(−S) averaged over the two ends of the mesh, symmetrized.

Exact up to the accuracy of KK when ρ vanishes at and beyond the ends of the mesh. Throws std::invalid_argument
for a non-finite F.
)DOC");
static const auto _c2py_doc_46 = _c2py_fun_46.doc(R"DOC(
[1] The Hilbert transform G_n(z) = ∫ εⁿ ρ(ε) / (z − ε) dε of the Steffen interpolant of a tabulated DOS (NRG
Ljubljana numerics, as tabulated_lattice_t), for z off the real axis; a real z gives the limit z + i0. The DOS
needs finite values, strictly increasing energies and ρ ≥ 0, but no normalization.

For an array of z the transform is taken at each element; for a square matrix z through its eigenvectors,
z = U diag(λ) U⁻¹ and G_n(z) = U diag(G_n(λ)) U⁻¹, so z must be diagonalizable with an eigenvector matrix of
condition number at most 1e8. Errors are std::invalid_argument.

------

[2] The Hilbert transform at each element of z

------

[3] The Hilbert transform of a square matrix z

------
)DOC");
static const auto _c2py_doc_47 = _c2py_fun_47.doc(R"DOC(
The hybridization from a spectral part Γ, e.g. after mixing: the eigenvalues of Γ below floor raised to it (Γ'),
then Δ = ⟨h⟩ + KK(−Γ') − iΓ' with the on-site energy ⟨h⟩ (one matrix per block, symmetrized). The same as
rebuild_hermitian_part(causal_projection(−iΓ, floor).first, ⟨h⟩), without the correction that the rebuild
discards.

Throws std::invalid_argument for a negative or non-finite floor, a non-finite or non-Hermitian Γ (to 1e-12 of its
largest element), and on-site energies of the wrong number or size.
)DOC");
static const auto _c2py_doc_48 = _c2py_fun_48.doc(R"DOC(
The impurity Green's function from the impurity Dyson equation, G = (G0⁻¹ − Σ)⁻¹ with G0⁻¹ = ω − h_imp − Δ,
per block and mesh point, with the constant impurity levels h_imp (one matrix per block).

Σ must have the structure and mesh of Δ. Throws std::invalid_argument for a mismatch, a non-finite input, a
wrong number or size of h_imp, and a singular ω − h_imp − Δ − Σ (naming the block and the frequency).
)DOC");
static const auto _c2py_doc_49 = _c2py_fun_49.doc(R"DOC(
The occupancy of G at temperature T (T = 0: the step at ω = 0). degeneracy gives the number of equivalent copies
each block stands for (the manifest's [structure] degeneracy), for every block of G and no other.

Throws std::invalid_argument for a non-finite G, a negative or non-finite T, and a degeneracy that does not cover
exactly the blocks of G with positive counts.
)DOC");
static const auto _c2py_doc_50 = _c2py_fun_50.doc(R"DOC(
The Hermitian part rebuilt from the spectral part and the high-frequency limit: F' = F_∞ + KK(−S) − iS.

For a function whose Hermitian part is missing or not trusted, e.g. Δ after mixing Γ. F_∞ has one matrix per
block of the size of the block, Hermitian to 1e-12 relative; it is symmetrized. Throws std::invalid_argument
otherwise, and for a non-finite F.
)DOC");
static const auto _c2py_doc_51 = _c2py_fun_51.doc(R"DOC(
The self-energy Σ = F_l G⁻¹, the estimator without I and F_r. Errors as for sigma_IFG.
)DOC");
static const auto _c2py_doc_52 = _c2py_fun_52.doc(R"DOC(
The self-energy from G, F_l, F_r and I: Σ = Σ_H + I − F_l G⁻¹ F_r, per block and mesh point, with the constant
Hartree term Σ_H (one matrix per block).

F_l, F_r and I must have the structure and mesh of G. Throws std::invalid_argument for a mismatch, a non-finite
input, a wrong number or size of Hartree matrices, and a singular G (naming the block and the frequency).
)DOC");
static const auto _c2py_doc_53 = _c2py_fun_53.doc(R"DOC(
The distance between the spectral functions of a and b, which must have the same structure and mesh
)DOC");
static const auto _c2py_doc_54 = _c2py_fun_54.doc(R"DOC(
The spectral part S = (F† − F)/(2i) of each block: Γ for Δ, −Im Σ for a 1x1 Σ, π A for G. Exactly Hermitian.
)DOC");
static const auto _c2py_doc_55 = _c2py_fun_55.doc(R"DOC(
∫ A_ij(ω) dω per block, Hermitian: 1 on the diagonal for a normalized G, up to the weight beyond the mesh
)DOC");
//--------------------- module function table  -----------------------------

static PyMethodDef module_methods[] = {
   {"causal_projection", (PyCFunction)c2py::pyfkw<_c2py_fun_44>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_44.c_str()},
   {"high_frequency_limit", (PyCFunction)c2py::pyfkw<_c2py_fun_45>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_45.c_str()},
   {"hilbert_transform", (PyCFunction)c2py::pyfkw<_c2py_fun_46>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_46.c_str()},
   {"hybridization_from_Gamma", (PyCFunction)c2py::pyfkw<_c2py_fun_47>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_47.c_str()},
   {"impurity_gf", (PyCFunction)c2py::pyfkw<_c2py_fun_48>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_48.c_str()},
   {"occupancy", (PyCFunction)c2py::pyfkw<_c2py_fun_49>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_49.c_str()},
   {"rebuild_hermitian_part", (PyCFunction)c2py::pyfkw<_c2py_fun_50>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_50.c_str()},
   {"sigma_FG", (PyCFunction)c2py::pyfkw<_c2py_fun_51>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_51.c_str()},
   {"sigma_IFG", (PyCFunction)c2py::pyfkw<_c2py_fun_52>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_52.c_str()},
   {"spectral_distance", (PyCFunction)c2py::pyfkw<_c2py_fun_53>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_53.c_str()},
   {"spectral_part", (PyCFunction)c2py::pyfkw<_c2py_fun_54>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_54.c_str()},
   {"spectral_weight", (PyCFunction)c2py::pyfkw<_c2py_fun_55>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_55.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

//--------------------- module struct & init error definition ------------

//// module doc directly in the code or "" if not present...
/// Or mandatory ?
static struct PyModuleDef module_def = {
   PyModuleDef_HEAD_INIT,
   "_dmft",                                                                                                                      /* name of module */
   R"RAWDOC(The DMFT layer of nrgljubljana_interface: causality, self-energy, measurements, lattices and the DMFT loop.)RAWDOC", /* module documentation, may be NULL */
   -1, /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
   module_methods,
   NULL,
   NULL,
   NULL,
   NULL};

//--------------------- module init function -----------------------------

extern "C" __attribute__((visibility("default"))) PyObject *PyInit__dmft() {

  if (not c2py::check_python_version("_dmft")) return NULL;

  // import numpy iff 'numpy/arrayobject.h' included
#ifdef Py_ARRAYOBJECT_H
  import_array();
#endif

  PyObject *m;

  if (PyType_Ready(&c2py::wrap_pytype<c2py::py_range>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_0>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_1>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_2>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_3>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_4>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_5>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_6>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_7>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_8>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_9>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_10>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_11>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_12>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_13>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_14>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_15>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_16>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_17>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_18>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_19>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_20>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_21>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_22>) < 0) return NULL;

  m = PyModule_Create(&module_def);
  if (m == NULL) return NULL;

  auto &conv_table = *c2py::conv_table_sptr.get();

  conv_table[std::type_index(typeid(c2py::py_range)).name()] = &c2py::wrap_pytype<c2py::py_range>;
#define _add_type(T, N) c2py::add_type_object_to_main<T>(N, m, conv_table)
  _add_type(_c2py_cls_0, "CausalityReportT");
  _add_type(_c2py_cls_1, "OccupancyT");
  _add_type(_c2py_cls_2, "DistanceT");
  _add_type(_c2py_cls_3, "BetheLatticeT");
  _add_type(_c2py_cls_4, "DosT");
  _add_type(_c2py_cls_5, "TabulatedLatticeT");
  _add_type(_c2py_cls_6, "CallableLatticeT");
  _add_type(_c2py_cls_7, "MixInputT");
  _add_type(_c2py_cls_8, "MixReportT");
  _add_type(_c2py_cls_9, "MixResultT");
  _add_type(_c2py_cls_10, "LinearMixer");
  _add_type(_c2py_cls_11, "MuControlParamsT");
  _add_type(_c2py_cls_12, "MuInputT");
  _add_type(_c2py_cls_13, "MuReportT");
  _add_type(_c2py_cls_14, "MuResultT");
  _add_type(_c2py_cls_15, "MuControlT");
  _add_type(_c2py_cls_16, "ConvergenceParamsT");
  _add_type(_c2py_cls_17, "ConvergenceInputT");
  _add_type(_c2py_cls_18, "ConvergenceResultT");
  _add_type(_c2py_cls_19, "ConvergenceT");
  _add_type(_c2py_cls_20, "LoopParamsT");
  _add_type(_c2py_cls_21, "StepResultT");
  _add_type(_c2py_cls_22, "PyDmftLoop");
#undef _add_type

  c2py::pyref module = c2py::pyref::module("h5.formats");
  if (not module) return nullptr;
  c2py::pyref register_class = module.attr("register_class");

  register_h5_type<_c2py_cls_0>(register_class);
  register_h5_type<_c2py_cls_8>(register_class);
  register_h5_type<_c2py_cls_11>(register_class);
  register_h5_type<_c2py_cls_13>(register_class);
  register_h5_type<_c2py_cls_16>(register_class);
  register_h5_type<_c2py_cls_18>(register_class);
  register_h5_type<_c2py_cls_20>(register_class);
  register_h5_type<_c2py_cls_21>(register_class);

  return m;
}
#endif
// CLAIR_WRAP_GEN
