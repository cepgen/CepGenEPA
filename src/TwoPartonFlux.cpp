/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Hamzeh Khanpour
 *                2024       Laurent Forthomme
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

#include "CepGenEPA/TwoPartonFlux.h"

namespace cepgen::epa {
  TwoPartonFlux::TwoPartonFlux(const ParametersList& params) : PartonFlux(params) {}

  ParametersDescription TwoPartonFlux::description() {
    auto desc = PartonFlux::description();
    desc.setDescription("Two-parton mass-dependent flux");
    return desc;
  }
}  // namespace cepgen::epa
