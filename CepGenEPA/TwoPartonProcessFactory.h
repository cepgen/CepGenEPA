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

#ifndef CepGenEPA_TwoPartonProcessFactory_h
#define CepGenEPA_TwoPartonProcessFactory_h

#include <CepGen/Modules/ModuleFactory.h>

/// Add a generic collinear, two-parton matrix element builder implematation
#define REGISTER_TWOPARTON_PROCESS(name, obj)                                        \
  namespace cepgen {                                                                 \
    struct BUILDERNM(obj) {                                                          \
      BUILDERNM(obj)() { TwoPartonProcessFactory::get().registerModule<obj>(name); } \
    };                                                                               \
    static const BUILDERNM(obj) gTwoPartonProcess##obj;                              \
  }                                                                                  \
  static_assert(true, "")

namespace cepgen::epa {
  class TwoPartonProcess;
}
namespace cepgen {
  /// A collinear, two-parton fluxes objects factory
  DEFINE_FACTORY(TwoPartonProcessFactory, epa::TwoPartonProcess, "Two-parton-level process matrix elements factory");
}  // namespace cepgen

#endif
