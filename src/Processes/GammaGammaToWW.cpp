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

#include <CepGen/Modules/CouplingFactory.h>
#include <CepGen/Physics/Constants.h>
#include <CepGen/Physics/Coupling.h>
#include <CepGen/Physics/PDG.h>

#include "CepGenEPA/TwoPartonProcess.h"
#include "CepGenEPA/TwoPartonProcessFactory.h"

using namespace cepgen;

class GammaGammaToWW : public epa::TwoPartonProcess {
public:
  explicit GammaGammaToWW(const ParametersList& params)
      : epa::TwoPartonProcess(params),
        alpha_em_(AlphaEMFactory::get().build(steer<ParametersList>("alphaEM"))),
        me_(PDG::get().mass(PDG::electron)),
        mw_(PDG::get().mass(24)),
        inv_mw2_(1. / mw_ / mw_) {}

  static ParametersDescription description() {
    auto desc = epa::TwoPartonProcess::description();
    desc.setDescription("Two-photon production of W boson pair");
    desc.add("alphaEM", AlphaEMFactory::get().describeParameters("burkhardt"));
    return desc;
  }

  std::string processDescription() const override { return "$\\gamma\\gamma\\to W^{+}W^{-}$"; }
  double matrixElement(double wgg) const override {
    const auto alpha_em = alpha_em_->operator()(wgg);
    if (wgg > 2. * mw_) {
      if (wgg > 300.)
        return 2. * prefactor_ * alpha_em * alpha_em * inv_mw2_;
      return (19. / 8.) * prefactor_ * alpha_em * alpha_em * inv_mw2_ * std::sqrt(wgg * wgg - 4. / inv_mw2_) / wgg;
    }
    return 0.;
  }

private:
  static constexpr double prefactor_ = 4. * M_PI * constants::GEVM2_TO_PB;
  const std::unique_ptr<Coupling> alpha_em_;
  const double me_;
  const double mw_;
  const double inv_mw2_;
};
REGISTER_TWOPARTON_PROCESS("gammagammatoww", GammaGammaToWW);
