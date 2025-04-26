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

class GammaGammaToZZ : public epa::TwoPartonProcess {
public:
  explicit GammaGammaToZZ(const ParametersList& params) : epa::TwoPartonProcess(params), mz_(PDG::get().mass(23)) {}

  static ParametersDescription description() {
    auto desc = epa::TwoPartonProcess::description();
    desc.setDescription("Two-photon production of Z boson pair");
    return desc;
  }

  std::string processDescription() const override { return "$\\gamma\\gamma\\to ZZ$"; }
  double matrixElement(double wgg) const override {
    if (wgg > 2. * mz_) {
      const auto inv_w2 = 1. / wgg / wgg;
      return 0.25786903395035327 / std::pow(1. + 5.749069613832837e11 * inv_w2 * inv_w2 * inv_w2 +
                                                6.914037195922673e7 * inv_w2 * inv_w2 + 23.264122861948383 * inv_w2,
                                            44.05927999125431);
    }
    return 0.;
  }

private:
  const double mz_;
};
REGISTER_TWOPARTON_PROCESS("gammagammatozz", GammaGammaToZZ);
