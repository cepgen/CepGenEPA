/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2026  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <CepGen/Core/Exception.h>
#include <CepGen/Integration/Integrator.h>
#include <CepGen/Modules/IntegratorFactory.h>
#include <CepGen/Physics/Constants.h>
#include <CepGen/Physics/PDG.h>

#include "CepGenEPA/HelicityAmplitudes.h"
#include "CepGenEPA/MatrixElements.h"
#include "CepGenEPA/TwoPartonProcess.h"
#include "CepGenEPA/TwoPartonProcessFactory.h"

using namespace cepgen;

namespace eft_aaaa {
  // Computes the  squared matrix element and the SM interference from free zeta_1, zeta_2
  double sqme(double s, double t, bool exclude_loops_SM, double zeta1, double zeta2) {
    //NOTE: zeta1/zeta2 expressed in GeV^-4
    if (s < 0 || t > 0 || t < -s)
      throw CG_FATAL("eft_aaaa:sqme") << "Invalid domain. Valid range is s>=0 and -s<=t<=0.";

    double value = 0;

    {  // Mpppp:
      // factor 8 is needed because of the conventions in Costantini, DeTollis, Pistoni
      const auto me_ex = 8. * Mpppp_eft(zeta1, zeta2, s, t);                     // exotic matrix element
      const auto me_sm = sm_aaaa::me_SM(Mpppp_fermion, s, t, exclude_loops_SM);  //  SM matrix element
      value += std::real(me_ex) * (std::real(me_ex) + 2 * std::real(me_sm)) +
               std::imag(me_ex) * (std::imag(me_ex) + 2 * std::imag(me_sm));
    }
    {  // Mppmm:
      const auto me_ex = 8. * Mppmm_eft(zeta1, zeta2, s, t);
      const auto me_sm = sm_aaaa::me_SM(Mppmm_fermion, s, t, exclude_loops_SM);
      value += std::real(me_ex) * (std::real(me_ex) + 2 * std::real(me_sm)) +
               std::imag(me_ex) * (std::imag(me_ex) + 2 * std::imag(me_sm));
    }
    {  // Mpmmp:
      const auto me_ex = 8. * Mpmmp_eft(zeta1, zeta2, s, t);
      const auto me_sm = sm_aaaa::me_SM(Mpmmp_fermion, s, t, exclude_loops_SM);
      value += std::real(me_ex) * (std::real(me_ex) + 2 * std::real(me_sm)) +
               std::imag(me_ex) * (std::imag(me_ex) + 2 * std::imag(me_sm));
    }
    {  // Mpmpm:
      const auto me_ex = 8. * Mpmpm_eft(zeta1, zeta2, s, t);
      const auto me_sm = sm_aaaa::me_SM(Mpmpm_fermion, s, t, exclude_loops_SM);
      value += std::real(me_ex) * (std::real(me_ex) + 2 * std::real(me_sm)) +
               std::imag(me_ex) * (std::imag(me_ex) + 2 * std::imag(me_sm));
    }
    {  // Mpppm
      const auto me_ex = 8. * Mpppm_eft(zeta1, zeta2, s, t);
      const auto me_sm = sm_aaaa::me_SM(Mpppm_fermion, s, t, exclude_loops_SM);
      value += std::real(me_ex) * (std::real(me_ex) + 2 * std::real(me_sm)) +
               std::imag(me_ex) * (std::imag(me_ex) + 2 * std::imag(me_sm));
    }
    return 0.5 * value;
  }
}  //namespace eft_aaaa

class GammaGammaToGammaGammaEFT : public epa::TwoPartonProcess {
public:
  explicit GammaGammaToGammaGammaEFT(const ParametersList& params)
      : epa::TwoPartonProcess(params),
        integrator_(IntegratorFactory::get().build(steer<ParametersList>("integrator"))),
        exclude_loops_(steer<bool>("excludeLoops")),
        zeta1_(steer<double>("zeta1")),
        zeta2_(steer<double>("zeta2")) {}

  static ParametersDescription description() {
    auto desc = epa::TwoPartonProcess::description();
    desc.setDescription("Two-photon production of photon pair (EFT)");
    desc.add("integrator", IntegratorFactory::get().describeParameters("gsl"));
    desc.add("excludeLoops", false);
    desc.add("zeta1", 1.e-12);
    desc.add("zeta2", 1.e-12);
    return desc;
  }

  std::string processDescription() const override { return "$\\gamma\\gamma\\rightarrow\\gamma\\gamma$ (EFT)"; }
  double matrixElement(double w) const override {
    const auto s = w * w;
    return prefactor_ *
           integrator_->integrate(
               [this, &s](double t) { return eft_aaaa::sqme(s, t, exclude_loops_, zeta1_, zeta2_) / s / s; },
               Limits{-s, 0.});
  }

private:
  static constexpr double prefactor_ = constants::GEVM2_TO_PB / 16. * M_1_PI;
  const std::unique_ptr<Integrator> integrator_;
  const bool exclude_loops_;
  const double zeta1_;
  const double zeta2_;
};
REGISTER_TWOPARTON_PROCESS("gammagammatogammagamma:eft", GammaGammaToGammaGammaEFT);
