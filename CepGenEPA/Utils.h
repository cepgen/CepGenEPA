#ifndef ggMatrixElements_Utils_h
#define ggMatrixElements_Utils_h

#include <complex>

namespace cepgen::epa::utils {
  std::complex<double> B(double z);
  std::complex<double> T(double z);
  std::complex<double> F(double q, double a);
  std::complex<double> I(double z, double w);
}  // namespace cepgen::epa::utils

#endif
