/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Hamzeh Khanpour
 *                2024-2025  Laurent Forthomme
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
#include <CepGen/Integration/Integrand.h>
#include <CepGen/Integration/Integrator.h>
#include <CepGen/Modules/FormFactorsFactory.h>
#include <CepGen/Modules/IntegratorFactory.h>
#include <CepGen/Modules/PartonFluxFactory.h>
#include <CepGen/PartonFluxes/CollinearFlux.h>
#include <CepGen/Physics/PDG.h>

#include <cmath>

#include "CepGenEPA/TwoPartonFlux.h"
#include "CepGenEPA/TwoPartonFluxFactory.h"

using namespace cepgen;

class LeptonProtonTwoPhotonFlux final : public epa::TwoPartonFlux {
public:
  explicit LeptonProtonTwoPhotonFlux(const ParametersList& params)
      : epa::TwoPartonFlux(params),
        lepton_(steer<ParametersList>("lepton")),
        proton_(steer<ParametersList>("proton")),
        total_flux_integrator_(IntegratorFactory::get().build(steer<ParametersList>("fastIntegrator"))),
        lepton_flux_integrator_(IntegratorFactory::get().build(steer<ParametersList>("fastIntegrator"))),
        proton_flux_integrator_(IntegratorFactory::get().build(steer<ParametersList>("fastIntegrator"))) {}

  static ParametersDescription description() {
    auto desc = epa::TwoPartonFlux::description();

    // lepton beam properties
    auto lepton_desc = BeamProperties::description();
    auto flux_e = CollinearFluxFactory::get().describeParameters("EPAFlux");
    flux_e.add("formFactors", FormFactorsFactory::get().describeParameters("PointLikeFermion"));
    lepton_desc.add("flux", flux_e);
    lepton_desc.add("energy", 50.);
    lepton_desc.add("q2Range", Limits{0., 1.e5});
    desc.add("lepton", lepton_desc);

    // proton beam properties
    auto proton_desc = BeamProperties::description();
    auto flux_p_elastic = CollinearFluxFactory::get().describeParameters("EPAFlux");
    flux_p_elastic.add("formFactors", FormFactorsFactory::get().describeParameters("StandardDipole"));
    proton_desc.add("flux", flux_p_elastic);
    proton_desc.add("energy", 7000.);
    proton_desc.add("q2Range", Limits{0., 10.});
    desc.add("proton", proton_desc);

    desc.add("fastIntegrator", IntegratorFactory::get().describeParameters("root"));
    desc.add("preciseIntegrator", IntegratorFactory::get().describeParameters("Vegas"));
    return desc;
  }

  double flux(const std::vector<double>& arguments) const override {
    const auto& wgg = arguments.at(0);
    const auto flux_wgg = total_flux_integrator_->integrate(
        [this, &wgg](double x1) -> double {
          if (const auto q2min_e = q2min(x1, lepton_.flux->mass2()); q2min_e < lepton_.q2range.max()) {
            return lepton_flux_integrator_->integrate(
                [this, &wgg, &x1](double q2_1) {
                  const auto eb1sq = lepton_.energy * lepton_.energy;
                  // the tricky part is that x2 is constrained by s, and x1!
                  const auto x2 =
                      0.5 * (wgg * wgg + q2_1) /
                      (x1 * lepton_.energy * proton_.energy +
                       proton_.energy * std::sqrt(x1 * x1 * eb1sq + q2_1) * (1. - q2_1 / (2. * eb1sq * (1. - x1))));
                  if (const auto q2min_p = q2min(x2, proton_.flux->mass2()); q2min_p < proton_.q2range.max())
                    return lepton_.flux->fluxQ2(x1, q2_1) / q2_1 / x1 *
                           proton_flux_integrator_->integrate(
                               [this, &x2](double q2_2) { return proton_.flux->fluxQ2(x2, q2_2) / q2_2 / x2; },
                               proton_.q2range);
                  return 0.;
                },
                lepton_.q2range);
          }
          return 0.;
        },
        cepgen::Limits{0., 1.});
    CG_DEBUG("LeptonProtonTwoPhotonFlux:flux") << "Flux at w_gg=" << wgg << " GeV: " << flux_wgg << ".";
    return flux_wgg;
  }

  inline bool fragmenting() const override { return false; }
  inline spdgid_t partonPdgId() const override { return PDG::photon; }
  inline double mass2() const override { return 0.; }

private:
  static const double q2min(double x, double mass2) { return mass2 * x * x / (1. - x); }

  struct BeamProperties : SteeredObject<BeamProperties> {
    explicit BeamProperties(const ParametersList& params)
        : SteeredObject(params),
          flux(CollinearFluxFactory::get().build(steer<ParametersList>("flux"))),
          energy(steer<double>("energy")),
          q2range(steer<Limits>("q2Range")) {}

    static ParametersDescription description() {
      auto desc = ParametersDescription();
      desc.add("flux", ParametersDescription()).setDescription("parton-from-beam flux modelling");
      desc.add("energy", 0.).setDescription("beam particle energy, in GeV");
      desc.add("q2Range", Limits{0., 1.e5}).setDescription("parton virtuality range, in GeV^2");
      return desc;
    }

    const std::unique_ptr<CollinearFlux> flux;
    const double energy;
    const Limits q2range;
  };

  const BeamProperties lepton_;
  const BeamProperties proton_;
  const std::unique_ptr<Integrator> total_flux_integrator_;
  const std::unique_ptr<Integrator> lepton_flux_integrator_;
  const std::unique_ptr<Integrator> proton_flux_integrator_;
};
REGISTER_TWOPARTON_FLUX("gmgm:lp", LeptonProtonTwoPhotonFlux);
