// Computes different helicity amplitudes as defined in
// Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787

#include <CepGen/Physics/PDG.h>

#include <cmath>

#include "CepGenEPA/HelicityAmplitudes.h"
#include "CepGenEPA/Utils.h"

using namespace std::complex_literals;

enum struct Region { no_limits, low, high, forward, backward };

Region limits(double sred, double tred, double ured) {
  static constexpr double s_low = 1.e1, s_high = 1.e9, t_low = 1.e-4, u_low = 1.e-3;
  if (sred <= 1.e-3)
    return Region::low;  // EFT limit
  if ((sred <= s_low && -tred < t_low * sred) || (sred > s_low && sred <= s_high && -tred < 1.e-3) ||
      (sred > s_high && -tred < 1.))
    return Region::forward;  // forward limit
  if ((sred <= s_low && -ured < t_low * sred) || (sred > s_low && -ured < u_low) || (sred > s_high && -ured < 1.))
    return Region::backward;  // backward limit
  if (sred > s_high)
    return Region::high;     // high energy limit
  return Region::no_limits;  // no limit

  // explanation:
  // for sred>shigh, optimal value to switch from HE limit to forward limit is |tred|=1
  // (at these value both limits are somewhat bad but quickly converge at either side)
  // for sred<shigh, only switch from exact result to forward limit at |t| < 0.001 for better accuracy
}

std::complex<double> Mxxxx_fermion(double x, double y) {
  // some auxilliary function used in Mpppp, Mpmpm, Mpmmp.
  std::complex<double> output{1., 0.};

  const double z = -x - y;
  double temp;

  temp = 2 * (y * y + z * z) / (x * x) - 2 / x;
  output += temp * (cepgen::epa::utils::T(y) + cepgen::epa::utils::T(z));

  temp = 1 / (2 * x * y) - 1 / y;
  output += temp * cepgen::epa::utils::I(x, y);

  temp = 1 / (2 * x * z) - 1 / z;
  output += temp * cepgen::epa::utils::I(x, z);

  temp = 4 / x + 1 / y + 1 / z + 1 / (2 * z * y) - 2 * (y * y + z * z) / (x * x);
  output += temp * cepgen::epa::utils::I(y, z);

  temp = 2 * (y - z) / x;
  output += temp * (cepgen::epa::utils::B(y) - cepgen::epa::utils::B(z));

  return output;
}

std::complex<double> Mpppp_fermion(double sred, double tred, int exclude_loops) {
  // M++++ from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787

  const double ured = -sred - tred;

  if (exclude_loops == 1 || exclude_loops == 3)
    return 0.;
  switch (limits(sred, tred, ured)) {
    case Region::low:  // EFT limit
      return -4. * (4. * (-1. / 36.) + 3. * (7. / 90.)) * sred * sred;
    case Region::forward:
    case Region::backward:  // Forward and backward limit
      return 1. / (2. * sred * sred) *
             (2. * sred * sred + (-2. * sred + 4. * sred * sred) * cepgen::epa::utils::B(sred) +
              (2. * sred - 8. * sred * sred) * cepgen::epa::utils::B(-sred) +
              (-1. + 2. * sred) * cepgen::epa::utils::T(sred) +
              (-1. - 2. * sred + 4. * sred * sred) * cepgen::epa::utils::T(-sred));
    case Region::high:  // high energy limit
      return 1. + (tred - ured) / sred * log(tred / ured) +
             (tred * tred + ured * ured) / (2. * sred * sred) * (pow(log(tred / ured), 2) + M_PI * M_PI);
    case Region::no_limits:
    default:
      return Mxxxx_fermion(sred, tred);
  }
}

std::complex<double> Mpmmp_fermion(double sred, double tred, int exclude_loops) {
  // M+--+ from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787
  const double ured = -sred - tred;

  if (exclude_loops == 1 || exclude_loops == 3)
    return 0.;
  switch (limits(sred, tred, ured)) {
    case Region::low:  // EFT limit
      return -4. * (4. * (-1. / 36.) + 3. * (7. / 90.)) * tred * tred;
    case Region::forward:  // Forward limit
      return 0.;
    case Region::backward:  // Backward limit
      return 1. / (2. * sred * sred) *
             (2. * sred * sred + (2. * sred + 4. * sred * sred) * cepgen::epa::utils::B(-sred) +
              (-2. * sred - 8. * sred * sred) * cepgen::epa::utils::B(sred) +
              (-1. - 2. * sred) * cepgen::epa::utils::T(-sred) +
              (-1. + 2. * sred + 4. * sred * sred) * cepgen::epa::utils::T(sred));
    case Region::high:  // high energy limit
      return {1. + (sred - ured) / tred * log(-sred / ured) +
                  (sred * sred + ured * ured) / (2. * tred * tred) * pow(log(-sred / ured), 2),
              -M_PI * ((sred - ured) / tred + (sred * sred + ured * ured) / (tred * tred) * log(-sred / ured))};
    case Region::no_limits:
    default:
      return Mxxxx_fermion(tred, sred);
  }
}

std::complex<double> Mpmpm_fermion(double sred, double tred, int exclude_loops) {
  // M+-+- from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787
  const double ured = -sred - tred;

  if (exclude_loops == 1 || exclude_loops == 3)
    return 0.;
  switch (limits(sred, tred, ured)) {
    case Region::low:  // EFT limit
      return -4. * (4. * (-1. / 36.) + 3. * (7. / 90.)) * ured * ured;
    case Region::forward:  // Forward limit
      return 1. / (2. * sred * sred) *
             (2. * sred * sred + (2. * sred + 4. * sred * sred) * cepgen::epa::utils::B(-sred) +
              (-2. * sred - 8. * sred * sred) * cepgen::epa::utils::B(sred) +
              (-1. - 2. * sred) * cepgen::epa::utils::T(-sred) +
              (-1. + 2. * sred + 4. * sred * sred) * cepgen::epa::utils::T(sred));
    case Region::backward:  // Backward limit
      return 0.;
    case Region::high:  // high energy limit
      return {1. + (tred - sred) / ured * log(-tred / sred) +
                  (sred * sred + tred * tred) / (2. * ured * ured) * pow(log(-tred / sred), 2),
              M_PI * ((tred - sred) / ured + (sred * sred + tred * tred) / (ured * ured) * log(-tred / sred))};
    case Region::no_limits:
    default:
      return Mxxxx_fermion(ured, tred);
  }
}

std::complex<double> Mpppm_fermion(double sred, double tred, int exclude_loops) {
  // M+--- from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787
  if (exclude_loops == 1 || exclude_loops == 3)
    return 0.;
  const double ured = -sred - tred;
  switch (limits(sred, tred, ured)) {
    case Region::low:  // EFT limit
    case Region::forward:
    case Region::backward:  // Forward and backward limit
      return 0.;
    case Region::high:  // high energy limit
      return -1.;
    case Region::no_limits:
    default: {
      std::complex<double> output{-1., 0.};

      double temp = -1 / sred - 1 / tred - 1 / ured;
      output += temp * (cepgen::epa::utils::T(sred) + cepgen::epa::utils::T(tred) + cepgen::epa::utils::T(ured));

      temp = 1 / ured + 1 / (2 * sred * tred);
      output += temp * cepgen::epa::utils::I(sred, tred);

      temp = 1 / tred + 1 / (2 * sred * ured);
      output += temp * cepgen::epa::utils::I(sred, ured);

      temp = 1 / sred + 1 / (2 * tred * ured);
      output += temp * cepgen::epa::utils::I(tred, ured);

      return output;
    }
  }
}

std::complex<double> Mppmm_fermion(double sred, double tred, int exclude_loops) {
  // M++-- from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787

  const double ured = -sred - tred;

  if (exclude_loops == 1 || exclude_loops == 3)
    return 0.;
  switch (limits(sred, tred, ured)) {
    case Region::low:  // EFT limit
      return -4. * (4. * (-1. / 36.) + (7. / 90.)) * (sred * sred + tred * tred + ured * ured);
    case Region::forward:
    case Region::backward:  // Forward and backward limit
      return 1. / (2. * sred * sred) *
             (-2. * sred * sred + (-2. * sred * cepgen::epa::utils::B(sred) + 2. * sred * cepgen::epa::utils::B(-sred) -
                                   cepgen::epa::utils::T(sred) - cepgen::epa::utils::T(-sred)));
    case Region::high:  // high energy limit
      return -1.;
    case Region::no_limits:
    default: {
      std::complex<double> output{-1., 0.};

      double temp = 1 / (2 * sred * tred);
      output += temp * cepgen::epa::utils::I(sred, tred);

      temp = 1 / (2 * sred * ured);
      output += temp * cepgen::epa::utils::I(sred, ured);

      temp = 1 / (2 * tred * ured);
      output += temp * cepgen::epa::utils::I(tred, ured);

      return output;
    }
  }
}

std::complex<double> Mxxxx_vector(double x, double y) {
  // some auxilliary function used in Mpppp, Mpmpm, Mpmmp.
  std::complex<double> output{-1.5, 0.};

  const double z = -x - y;

  double temp = -3 * (y - z) / x;
  output += temp * (cepgen::epa::utils::B(y) - cepgen::epa::utils::B(z));

  temp = -1 / x * (8 * x - 3 - 6 * y * z / x);
  output += temp * (cepgen::epa::utils::T(y) + cepgen::epa::utils::T(z));

  temp = 1 / x * (8 * x - 6 - 6 * y * z / x) - 4 * (x - 0.25) * (x - 0.75) / (y * z);
  output += temp * cepgen::epa::utils::I(y, z);

  temp = -4 * (x - 0.25) * (x - 0.75) / (x * y);
  output += temp * cepgen::epa::utils::I(x, y);

  temp = -4 * (x - 0.25) * (x - 0.75) / (x * z);
  output += temp * cepgen::epa::utils::I(x, z);

  return output;
}

std::complex<double> Mpppp_vector(double sred, double tred, int exclude_loops) {
  if (exclude_loops == 2 || exclude_loops == 3)
    return 0.;

  double ured = -sred - tred;
  switch (limits(sred, tred, ured)) {
    case Region::low:  // EFT limit
      return -4. * (4. * (-5. / 32.) + 3. * (27. / 40.)) * sred * sred;
    case Region::forward:
    case Region::backward:  // Forward and backward limit
      return -1.5 + 8. * (sred - 0.25) * (sred - 0.75) / sred * cepgen::epa::utils::B(sred) +
             (-8. * (sred - 0.25) * (sred - 0.75) / sred + 3.) * cepgen::epa::utils::B(-sred) +
             4. * (sred - 0.25) * (sred - 0.75) / (sred * sred) * cepgen::epa::utils::T(sred) +
             (4. * (sred - 0.25) * (sred - 0.75) / (sred * sred) - (8. * sred - 3.) / sred) *
                 cepgen::epa::utils::T(-sred);
    case Region::high:  // high energy limit
      return {-1. * (1.5 + 1.5 * (ured - tred) / sred * log(ured / tred) +
                     2. * (1. - 0.75 * tred * ured / (sred * sred)) * (pow(log(ured / tred), 2) + M_PI * M_PI) +
                     2. * sred * sred *
                         (log(4. * sred) * log(-4. * tred) / (sred * tred) +
                          log(4. * sred) * log(-4. * ured) / (sred * ured) +
                          log(-4. * ured) * log(-4. * tred) / (ured * tred))),
              (2. * M_PI * sred * sred * (log(-4. * ured) / (sred * ured) + log(-4. * tred) / (sred * tred)))};
    case Region::no_limits:
    default:
      return Mxxxx_vector(sred, tred);
  }
}

std::complex<double> Mpmmp_vector(double sred, double tred, int exclude_loops) {
  if (exclude_loops == 2 || exclude_loops == 3)
    return 0.;

  double ured = -sred - tred;
  switch (limits(sred, tred, ured)) {
    case Region::low:  // EFT limit
      return -4. * (4. * (-5. / 32.) + 3. * (27. / 40.)) * tred * tred;
    case Region::forward:  // Forward limit
      return 0.;
    case Region::backward:  // Backward limit
      return -1.5 - 8. * (-sred - 0.25) * (-sred - 0.75) / sred * cepgen::epa::utils::B(-sred) +
             (8. * (-sred - 0.25) * (-sred - 0.75) / sred + 3.) * cepgen::epa::utils::B(sred) +
             4. * (-sred - 0.25) * (-sred - 0.75) / (sred * sred) * cepgen::epa::utils::T(-sred) +
             (4. * (-sred - 0.25) * (-sred - 0.75) / (sred * sred) + (-8. * sred - 3.) / sred) *
                 cepgen::epa::utils::T(sred);
    case Region::high:  // high energy limit
      return {
          -(1.5 + 1.5 * (ured - sred) / tred * log(-ured / sred) +
            2. * (1. - 0.75 * sred * ured / (tred * tred)) * pow(log(-ured / sred), 2) +
            2. * tred * tred *
                (log(4. * sred) * log(-4. * tred) / (sred * tred) + log(4. * sred) * log(-4. * ured) / (sred * ured) +
                 log(-4. * ured) * log(-4. * tred) / (ured * tred))),
          -(1.5 * (sred - ured) / tred * (-M_PI) +
            2. * (1. - 0.75 * sred * ured / (tred * tred)) * M_PI * 2. * log(-ured / sred) +
            2. * (-M_PI) * tred * tred * (log(-4. * ured) / (ured * sred) + log(-4. * tred) / (tred * sred)))};
    case Region::no_limits:
    default:
      return Mxxxx_vector(tred, sred);
  }
}

std::complex<double> Mpmpm_vector(double sred, double tred, int exclude_loops) {
  const double ured = -tred - sred;

  if (exclude_loops == 2 || exclude_loops == 3)
    return 0.;
  switch (limits(sred, tred, ured)) {
    case Region::low:  // EFT limit
      return -4. * (4. * (-5. / 32.) + 3. * (27. / 40.)) * ured * ured;
    case Region::forward:  // Forward limit
      return -1.5 - 8. * (-sred - 0.25) * (-sred - 0.75) / sred * cepgen::epa::utils::B(-sred) +
             (8. * (-sred - 0.25) * (-sred - 0.75) / sred + 3.) * cepgen::epa::utils::B(sred) +
             4. * (-sred - 0.25) * (-sred - 0.75) / (sred * sred) * cepgen::epa::utils::T(-sred) +
             (4. * (-sred - 0.25) * (-sred - 0.75) / (sred * sred) + (-8. * sred - 3.) / sred) *
                 cepgen::epa::utils::T(sred);
    case Region::backward:  // Backward limit
      return 0.;
    case Region::high:  // high energy limit
      return {-(1.5 + 1.5 * (tred - sred) / ured * std::log(-tred / sred) +
                2. * (1. - 0.75 * sred * tred / (ured * ured)) * std::pow(log(-tred / sred), 2) +
                2. * ured * ured *
                    (std::log(4. * sred) * std::log(-4. * tred) / (sred * tred) +
                     std::log(4. * sred) * std::log(-4. * ured) / (sred * ured) +
                     std::log(-4. * ured) * std::log(-4. * tred) / (ured * tred))),
              -(1.5 * (sred - tred) / ured * (-M_PI) +
                2. * (1. - 0.75 * sred * tred / (ured * ured)) * M_PI * 2. * std::log(-tred / sred) +
                2. * (-M_PI) * ured * ured *
                    (std::log(-4. * ured) / (ured * sred) + std::log(-4. * tred) / (tred * sred)))};
    case Region::no_limits:
    default:
      return Mxxxx_vector(ured, tred);
  }
}

std::complex<double> Mpppm_vector(double sred, double tred, int exclude_loops) {
  //double ured=-tred-sred;

  if (exclude_loops == 2 || exclude_loops == 3)
    return 0.;

  /*if (sred < 0.001)  // EFT limit
    return 0.;
  if (sred < 10000. && sred > 0.001 && (-tred < 0.0001 * sred || -ured < 0.0001 * sred))  // Forward and backward limit
    return 0.;*/

  return -1.5 * Mpppm_fermion(sred, tred, exclude_loops);
}

std::complex<double> Mppmm_vector(double sred, double tred, int exclude_loops) {
  //double ured = -sred-tred;

  if (exclude_loops == 2 || exclude_loops == 3)
    return 0.;

  /*if (sred < 0.001)  // EFT limit
    return -4. * (4. * (-5. / 32.) + (27. / 40.)) * (sred * sred + tred * tred + ured * ured);
  if (sred < 1.e4 && sred > 1.e-3 && (-tred < 1.e-4 * sred || -ured < 1.e-4 * sred))  // Forward and backward limit
    return -1.5 / (2. * sred * sred) *
           (-2. * sred * sred + (-2. * sred * cepgen::epa::utils::B(sred) + 2. * sred * cepgen::epa::utils::B(-sred) -
                                 cepgen::epa::utils::T(sred) - cepgen::epa::utils::T(-sred)));*/

  return -1.5 * Mppmm_fermion(sred, tred, exclude_loops);
}

std::complex<double> Mpppp_eft(double zeta1, double zeta2, double s, double t) {
  return -0.25 * (4. * zeta1 + 3 * zeta2) * s * s;
}

std::complex<double> Mpmmp_eft(double zeta1, double zeta2, double s, double t) {
  return -0.25 * (4. * zeta1 + 3 * zeta2) * t * t;
}

std::complex<double> Mpmpm_eft(double zeta1, double zeta2, double s, double t) {
  const double u = -s - t;
  return -0.25 * (4. * zeta1 + 3 * zeta2) * u * u;
}

std::complex<double> Mpppm_eft(double /*zeta1*/, double /*zeta2*/, double /*s*/, double /*t*/) { return 0.; }

std::complex<double> Mppmm_eft(double zeta1, double zeta2, double s, double t) {
  const double u = -s - t;
  return -0.25 * (4. * zeta1 + zeta2) * (s * s + t * t + u * u);
}
