#include <gsl/gsl_sf.h>

#include <cmath>

#include "CepGenEPA/Utils.h"

namespace cepgen::epa::utils {
  std::complex<double> B(double z) {
    if (z == 0)  // b is infinite
      return 0;
    if (z < 0) {  // b > 1
      const auto b = sqrt(1 - 1 / z);
      //     return 0.5 * b * log( (1+b) / (b-1) ) - 1;
      return b * atanh(1 / b) - 1;
    }
    if (z >= 1) {  // 0 <= b < 1
      const auto b = sqrt(1 - 1 / z);
      //      return 0.5 * b * log( (1+b) / (1-b) ) - 1;
      return {b * atanh(b) - 1, -M_PI / 2 * sqrt(1 - 1 / z)};
    }
    // b is imaginary
    const auto b = sqrt(1 / z - 1);
    return b * atan(1 / b) - 1;
  }

  std::complex<double> T(double z) {
    if (z == 0)  // b is infinite
      return 0;
    if (z < 0) {  // b>1 is real
      const auto b = sqrt(1 - 1 / z);
      const auto temp = 0.5 * log((1 + b) / (b - 1));
      return temp * temp;
    }
    if (z > 1) {  // 0<b<1 is real
      const auto b = sqrt(1 - 1 / z);
      const auto temp = 0.5 * log((1 + b) / (1 - b));
      return {temp * temp - 0.25 * M_PI * M_PI, -M_PI * acosh(sqrt(z))};
    }
    if (z == 1)  // b=0
      return -0.25 * M_PI * M_PI;
    // b is imaginary
    const auto b = sqrt(1 / z - 1);
    const auto temp = atan(1 / b);
    return -temp * temp;
  }

  std::complex<double> F(double q, double a) {  // some auxiliary function used in ReI(z,w)
                                                // (sum of the 4 dilogs)
    if (q > 0 && q < 1) {                       // b(q) is imaginary, arguments of dilogs are complex
      const auto b = std::sqrt(1. / q - 1.);

      auto r = (a + 1) / std::hypot(a, b);
      const auto theta = -std::atan(b / a);  // r e^itheta = (a+1)/(a+b)

      gsl_sf_result result_re, result_im;
      gsl_sf_complex_dilog_e(r, theta, &result_re, &result_im);

      auto result = -2 * result_re.val;  // -Re[ Li2( a+1/a+b) + Li2( a+1/a-b)]

      r = (a - 1) / sqrt(a * a + b * b);  // r e^itheta = (a-1)/(a+b)

      gsl_sf_complex_dilog_e(r, theta, &result_re, &result_im);

      result += 2 * result_re.val;  // Re[ Li2( a-1/a+b) + Li2( a-1/a-b)]
      return result;
    }
    // b(q) real and so are all the arguments of the dilogs
    const auto b = std::sqrt(1. - 1. / q);
    return gsl_sf_dilog((a - 1) / (a + b)) + gsl_sf_dilog((a - 1) / (a - b)) - gsl_sf_dilog((a + 1) / (a + b)) -
           gsl_sf_dilog((a + 1) / (a - b));
  }

  std::complex<double> I(double z, double w) {  // allowed regions z,w<=0 || z>=0,-z<=w<=0
    if (z == 0 || w == 0)
      return 0;

    double i_real, i_imag;

    // because of exact cancellations Taylor-expand
    // the function near the origin to maintain double precision
    if (const auto lim = 1.e-4; z < lim && z > -lim && w < lim && w > -lim)
      //      return 2./3.*z*z;
      i_real = 2. / 3. * z * w + 4. / 15. * z * w * (z + w) + 16. / 105. * z * w * (z * z + w * w + z * w / 2);
    else {
      const auto a_real = std::sqrt(1 - 1 / w - 1 / z);
      i_real = std::real((F(z, a_real) + F(w, a_real)) / (2 * a_real));
    }

    {  // imaginary part
      // allowed regions z,w<=0 || z>=0,-z<=w<=0
      const auto a_imag = std::sqrt(1 - 1 / z - 1 / w);
      if (z > 1 && w < 0) {
        const auto b_imag = std::sqrt(1 - 1 / z);
        i_imag = M_PI / (2 * a_imag) * std::log((a_imag - b_imag) / (a_imag + b_imag));
      } else if (w > 1 && z < 0) {
        const auto b_imag = std::sqrt(1 - 1 / w);
        i_imag = M_PI / (2 * a_imag) * std::log((a_imag - b_imag) / (a_imag + b_imag));
      } else
        i_imag = 0.;
    }
    return {i_real, i_imag};
  }
}  // namespace cepgen::epa::utils
