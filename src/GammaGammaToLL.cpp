/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#include <CepGen/Physics/Constants.h>
#include <CepGen/Physics/PDG.h>

#include "CepGenEPA/TwoPartonProcess.h"
#include "CepGenEPA/TwoPartonProcessFactory.h"

using namespace cepgen;

class GammaGammaToLL : public epa::TwoPartonProcess {
public:
  explicit GammaGammaToLL(const ParametersList& params)
      : epa::TwoPartonProcess(params), ml_(steer<ParticleProperties>("lepton").mass), ml2_(ml_ * ml_) {}

  static ParametersDescription description() {
    auto desc = epa::TwoPartonProcess::description();
    desc.setDescription("Two-photon production of lepton pair");
    desc.add("lepton", 13);
    return desc;
  }

  std::string processDescription() const override { return "$\\gamma\\gamma\\to l^{+}l^{-}$"; }
  double matrixElement(double w) const override {
    if (w < 2. * ml_)
      return 0.;
    const auto beta2 = 1. - 4. * ml2_ / w / w;
    if (beta2 < 0.)
      return 0.;
    const auto beta = std::sqrt(beta2);
    printf("%g->%g\n",
           w,
           prefactor_ / w / w * beta * (3. - beta2 * beta2) / (2 * beta) * std::log((1. + beta) / (1. - beta)) - 2 +
               beta2);
    return prefactor_ / w / w * beta * (3. - beta2 * beta2) / (2 * beta) * std::log((1. + beta) / (1. - beta)) - 2 +
           beta2;
  }

private:
  static constexpr double prefactor_ = 4. * M_PI * constants::GEVM2_TO_PB * constants::ALPHA_EM * constants::ALPHA_EM;
  const double ml_, ml2_;
};
REGISTER_TWOPARTON_PROCESS("gammagammatoll", GammaGammaToLL);
