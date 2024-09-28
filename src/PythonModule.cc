#include <boost/python.hpp>

#include "ggMatrixElements/matrix_elements.h"

BOOST_PYTHON_MODULE(ggMatrixElements) {
  namespace py = boost::python;
  py::def("sqme_sm", sm_aaaa::sqme);
  py::def("sqme_eft", eft_aaaa::sqme);
}
