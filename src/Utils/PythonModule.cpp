#include <CepGen/Generator.h>

#include <boost/python.hpp>

#include "CepGenEPA/MatrixElements.h"
#include "CepGenEPA/TwoPartonFlux.h"
#include "CepGenEPA/TwoPartonFluxFactory.h"

BOOST_PYTHON_FUNCTION_OVERLOADS(sqme_sm, sm_aaaa::sqme, 2, 3)
BOOST_PYTHON_FUNCTION_OVERLOADS(sqme_eft, eft_aaaa::sqme, 2, 5)

namespace cepgen::epa::python {
  namespace py = boost::python;

  template <typename T>
  py::list to_python_list(const std::vector<T>& vec) {
    py::list list;
    std::for_each(vec.begin(), vec.end(), [&list](const auto& t) { list.append(t); });
    return list;
  }
}  // namespace cepgen::epa::python

BOOST_PYTHON_MODULE(libCepGenEPA) {
  namespace py = boost::python;

  cepgen::initialise();

  py::def("sqme_sm", sm_aaaa::sqme, sqme_sm((py::arg("s"), py::arg("t"), py::arg("exclude_loops") = false)));
  py::def(
      "sqme_eft",
      eft_aaaa::sqme,
      sqme_eft((
          py::arg("s"), py::arg("t"), py::arg("exclude_loops") = false, py::arg("zeta1") = 0., py::arg("zeta2") = 0.)));

  struct TwoPartonFluxWrap : cepgen::epa::TwoPartonFlux, py::wrapper<cepgen::epa::TwoPartonFlux> {
    std::pair<cepgen::spdgid_t, cepgen::spdgid_t> partons() const override {
      if (const py::override ov = this->get_override("partons"); ov)
        return ov();
      return TwoPartonFlux::partons();
    }
    double flux(double wpp) const override {
      if (const py::override ov = this->get_override("flux"); ov)
        return ov(wpp);
      return TwoPartonFlux::flux(wpp);
    }
  };

  py::class_<TwoPartonFluxWrap, boost::noncopyable>(
      "_TwoPartonFlux", "A modelling for the two-parton flux", py::no_init)
      .def("__call__", &cepgen::epa::TwoPartonFlux::flux);

  py::class_<cepgen::TwoPartonFluxFactory, boost::noncopyable>(
      "TwoPartonFlux", "Two-parton flux modelling retrieval tool", py::no_init)
      .def("build",
           py::make_function(
               +[](const std::string& name) {
                 return cepgen::TwoPartonFluxFactory::get().build(cepgen::ParametersList{}.setName(name)).release();
               },
               py::return_value_policy<py::manage_new_object>()))
      .add_static_property(
          "modules",
          +[]() { return cepgen::epa::python::to_python_list(cepgen::TwoPartonFluxFactory::get().modules()); });
}
