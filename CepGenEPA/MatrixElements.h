#ifndef ggMatrixElements_MatrixElements_h
#define ggMatrixElements_MatrixElements_h

#include <complex>

namespace sm_aaaa {
  std::complex<double> me_SM(std::complex<double> (*me)(double, double, int),
                             double s,
                             double t,
                             bool exclude_loops = false);
  double sqme(double s, double t, bool exclude_loops = false);
}  // namespace sm_aaaa

namespace eft_aaaa {
  double sqme(double s, double t, bool exclude_loops_SM = false, double zeta1 = 0., double zeta2 = 0.);
}

#endif
