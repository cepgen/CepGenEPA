/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2025  Laurent Forthomme
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

#include <CepGen/Modules/CouplingFactory.h>
#include <CepGen/Physics/Constants.h>
#include <CepGen/Physics/Coupling.h>
#include <CepGen/Physics/PDG.h>

#include "CepGenEPA/TwoPartonProcess.h"
#include "CepGenEPA/TwoPartonProcessFactory.h"

using namespace cepgen;

class GammaGammaToFF : public epa::TwoPartonProcess {
public:
  explicit GammaGammaToFF(const ParametersList& params)
      : epa::TwoPartonProcess(params),
        alpha_em_(AlphaEMFactory::get().build(steer<ParametersList>("alphaEM"))),
        fermion_properties_(steer<ParticleProperties>("fermion")),
        min_w2_(std::pow(2. * fermion_properties_.mass, 2)),
        prefactor_(4. * M_PI * constants::GEVM2_TO_PB * std::pow(fermion_properties_.integerCharge() / 3., 4) *
                   fermion_properties_.colours) {}

  static ParametersDescription description() {
    auto desc = epa::TwoPartonProcess::description();
    desc.setDescription("Two-photon production of fermion pair");
    desc.add("fermion", 13);
    desc.add("alphaEM", AlphaEMFactory::get().describeParameters("fixed"));
    return desc;
  }

  std::string processDescription() const override { return "$\\gamma\\gamma\\to f\\bar{f}$"; }
  double matrixElement(double w) const override {
    const auto w2 = w * w;
    if (w2 < min_w2_)
      return 0.;
    const auto beta2 = 1. - min_w2_ / w2;
    if (beta2 < 0.)
      return 0.;
    const auto alpha_em = alpha_em_->operator()(w);
    const auto beta = std::sqrt(beta2);
    return prefactor_ * alpha_em * alpha_em / w2 * beta * (3. - beta2 * beta2) / (2 * beta) *
               std::log((1. + beta) / (1. - beta)) -
           2 + beta2;
  }

private:
  const std::unique_ptr<Coupling> alpha_em_;
  const ParticleProperties fermion_properties_;
  const double min_w2_;
  const double prefactor_;
};
REGISTER_TWOPARTON_PROCESS("gammagammatoff", GammaGammaToFF);
