#ifndef ggMatrixElements_matrix_elements_h
#define ggMatrixElements_matrix_elements_h

static constexpr double alpha_em = 1. / 137.036;  // EM coupling at zero momentum (on shell scheme)
static constexpr double mW = 80.385;              // W mass in GeV

namespace sm_aaaa {
  void me_SM(void (*me)(double, double, double *, double *, int),
             double s,
             double t,
             double *re,
             double *im,
             bool exclude_loops = false);
  double sqme(double s, double t, bool exclude_loops = false);
}  // namespace sm_aaaa

namespace eft_aaaa {
  double sqme(double s, double t, bool exclude_loops_SM = false, double zeta1 = 0., double zeta2 = 0.);
}

#endif
