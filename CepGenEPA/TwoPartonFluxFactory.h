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

#ifndef CepGenEPA_TwoPartonFluxFactory_h
#define CepGenEPA_TwoPartonFluxFactory_h

#include <CepGen/Modules/ModuleFactory.h>

/// Add a generic collinear, two-parton flux evaluator builder definition
#define REGISTER_TWOPARTON_FLUX(name, obj)                                           \
  namespace cepgen {                                                                 \
    struct BUILDER_NAME(obj) {                                                       \
      BUILDER_NAME(obj)() { TwoPartonFluxFactory::get().registerModule<obj>(name); } \
    };                                                                               \
    static const BUILDER_NAME(obj) gTwoPartonFlux##obj;                              \
  }                                                                                  \
  static_assert(true, "")

namespace cepgen::epa {
  class TwoPartonFlux;
}
namespace cepgen {
  /// A collinear, two-parton fluxes objects factory
  DEFINE_FACTORY(TwoPartonFluxFactory, epa::TwoPartonFlux, "Two-parton flux estimators factory");
}  // namespace cepgen

#endif
