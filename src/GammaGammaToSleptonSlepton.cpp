/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Hamzeh Khanpour
 *                2025  Laurent Forthomme
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

class GammaGammaToSleptonSlepton : public epa::TwoPartonProcess {
public:
  explicit GammaGammaToSleptonSlepton(const ParametersList& params)
      : epa::TwoPartonProcess(params), msl_(steer<double>("msl")) {}

  static ParametersDescription description() {
    auto desc = epa::TwoPartonProcess::description();
    desc.setDescription("Two-photon production of W boson pair");
    desc.add("msl", 100.).setDescription("slepton mass");
    return desc;
  }

  std::string processDescription() const override { return "$\\gamma\\gamma\\to\\tilde{l}^{+}\\tilde{l}^{-}$"; }
  double matrixElement(double wgg) const override {
    if (wgg <= msl_)
      return 0.;
    const auto alpha2 = 1. / 137. / 137.;
    const auto beta2 = 1. - 4. * msl_ * msl_ / wgg / wgg;
    if (beta2 < 0.)
      return 0.;
    const auto beta = std::sqrt(beta2);
    return 2. * constants::GEVM2_TO_PB * M_PI * alpha2 / wgg / wgg * beta *
           (2. - beta2 - (1. - beta2 * beta2) / (2. * beta) * std::log((1. + beta) / (1. - beta)));
  }

private:
  const double msl_;
};
REGISTER_TWOPARTON_PROCESS("gammagammatoslsl", GammaGammaToSleptonSlepton);
