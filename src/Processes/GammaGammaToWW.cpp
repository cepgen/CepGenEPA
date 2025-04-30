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

class GammaGammaToWW : public epa::TwoPartonProcess {
public:
  explicit GammaGammaToWW(const ParametersList& params)
      : epa::TwoPartonProcess(params),
        me_(PDG::get().mass(PDG::electron)),
        mw_(PDG::get().mass(24)),
        inv_mw2_(1. / mw_ / mw_) {}

  static ParametersDescription description() {
    auto desc = epa::TwoPartonProcess::description();
    desc.setDescription("Two-photon production of W boson pair");
    return desc;
  }

  std::string processDescription() const override { return "$\\gamma\\gamma\\to W^{+}W^{-}$"; }
  double matrixElement(double wgg) const override {
    const auto alpha2 = 1. / 128. / 128.;
    if (wgg > 2. * mw_)
      return (19. / 2.) * M_PI * constants::GEVM2_TO_PB * alpha2 * inv_mw2_ * std::sqrt(wgg * wgg - 4. / inv_mw2_) /
             wgg;
    if (wgg > 300.)
      return 8. * M_PI * constants::GEVM2_TO_PB * alpha2 * inv_mw2_;
    return 0.;
  }

private:
  const double me_;
  const double mw_;
  const double inv_mw2_;
};
REGISTER_TWOPARTON_PROCESS("gammagammatoww", GammaGammaToWW);
