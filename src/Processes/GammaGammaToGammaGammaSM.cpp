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

namespace sm_aaaa {
  double prefac_W = -1.;
  const std::array<double, 9> SM_weight = {1, 1, 1, 16. / 27., 16. / 27., 16. / 27., 1. / 27., 1. / 27., 1. / 27.};
  std::array<double, 9> SM_masses;
  bool initialised = false;

  std::complex<double> me_SM(std::complex<double> (*me)(double, double, int), double s, double t, bool exclude_loops) {
    if (!initialised) {
      prefac_W = 0.25 / std::pow(PDG::get().mass(23 /*W*/), 2);
      SM_masses = {PDG::get().mass(11),
                   PDG::get().mass(13),
                   PDG::get().mass(15),
                   PDG::get().mass(2),
                   PDG::get().mass(4),
                   PDG::get().mass(6),
                   PDG::get().mass(1),
                   PDG::get().mass(3),
                   PDG::get().mass(5)};
      initialised = true;
    }
    // This routine computes the complex SM amplitude
    // The first argument can be any of the helicity amplitudes Mpppp,Mppmm,Mpmpm,Mpmmp,Mpppm
    std::complex<double> output;

    // SM fermion content: (e,mu,tau,u,c,t,d,s,b)
    // SM_weight equals (number of colors) * (el. charge)^4
    // SM masses in GeV
    for (size_t i = 0; i < SM_masses.size(); i++) {
      const auto prefac_f = 1. / (4 * SM_masses.at(i) * SM_masses.at(i));
      const auto d = me(s * prefac_f, t * prefac_f, exclude_loops);
      output += SM_weight.at(i) * d;
    }

    // Add also the W contribution
    if (me == Mpppp_fermion)
      output += Mpppp_vector(s * prefac_W, t * prefac_W, exclude_loops);
    else if (me == Mppmm_fermion)
      output += Mppmm_vector(s * prefac_W, t * prefac_W, exclude_loops);
    else if (me == Mpmpm_fermion)
      output += Mpmpm_vector(s * prefac_W, t * prefac_W, exclude_loops);
    else if (me == Mpmmp_fermion)
      output += Mpmmp_vector(s * prefac_W, t * prefac_W, exclude_loops);
    else if (me == Mpppm_fermion)
      output += Mpppm_vector(s * prefac_W, t * prefac_W, exclude_loops);

    // the factor of 8 is needed because of the conventions in
    // Costantini, DeTollis, Pistoni
    output *= 8 * constants::ALPHA_EM * constants::ALPHA_EM;

    return output;
  }

  // compute the SM squared matrix element, including leptons, quarks and the W boson
  double sqme(double s, double t, bool exclude_loops) {
    if (s < 0 || t > 0 || t < -s)
      throw CG_FATAL("sm_aaaa:sqme") << "Invalid domain. Valid range is s>=0 and -s<=t<=0.";

    return 0.5 * (4. * std::norm(me_SM(Mpppm_fermion, s, t, exclude_loops)) +
                  std::norm(me_SM(Mppmm_fermion, s, t, exclude_loops)) +
                  std::norm(me_SM(Mpppp_fermion, s, t, exclude_loops)) +
                  std::norm(me_SM(Mpmmp_fermion, s, t, exclude_loops)) +
                  std::norm(me_SM(Mpmpm_fermion, s, t, exclude_loops)));
  }

}  //namespace sm_aaaa

class GammaGammaToGammaGammaSM : public epa::TwoPartonProcess {
public:
  explicit GammaGammaToGammaGammaSM(const ParametersList& params)
      : epa::TwoPartonProcess(params),
        integrator_(IntegratorFactory::get().build(steer<ParametersList>("integrator"))),
        exclude_loops_(steer<bool>("excludeLoops")) {}

  static ParametersDescription description() {
    auto desc = epa::TwoPartonProcess::description();
    desc.setDescription("Two-photon production of photon pair (SM)");
    desc.add("integrator", IntegratorFactory::get().describeParameters("gsl"));
    desc.add("excludeLoops", false);
    return desc;
  }

  std::string processDescription() const override { return "$\\gamma\\gamma\\rightarrow\\gamma\\gamma$ (SM)"; }
  double matrixElement(double w) const override {
    const auto s = w * w;
    return prefactor_ *
           integrator_->integrate([this, &s](double t) { return sm_aaaa::sqme(s, t, exclude_loops_) / s / s; },
                                  Limits{-s, 0.});
  }

private:
  static constexpr double prefactor_ = constants::GEVM2_TO_PB / 16. * M_1_PI;
  const std::unique_ptr<Integrator> integrator_;
  const bool exclude_loops_;
};
REGISTER_TWOPARTON_PROCESS("gammagammatogammagamma:sm", GammaGammaToGammaGammaSM);
