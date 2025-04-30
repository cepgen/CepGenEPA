/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2025  Laurent Forthomme
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

#ifndef CepGenEPA_BeamProperties_h
#define CepGenEPA_BeamProperties_h

#include <CepGen/Core/SteeredObject.h>

#include <memory>

namespace cepgen {
  class CollinearFlux;
}  // namespace cepgen

namespace cepgen::epa {
  struct BeamProperties : SteeredObject<BeamProperties> {
    explicit BeamProperties(const ParametersList&);

    static ParametersDescription description();

    std::unique_ptr<CollinearFlux> flux{};
    const double energy;
    const Limits q2range;
  };
}  // namespace cepgen::epa

#endif
