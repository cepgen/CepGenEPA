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

#include <CepGen/Core/Exception.h>
#include <CepGen/Modules/ModuleFactory.h>

#include "CepGenEPA/TwoPartonFlux.h"

namespace cepgen {
  template class ModuleFactory<epa::TwoPartonFlux>;

  template <>
  ModuleFactory<epa::TwoPartonFlux>::ModuleFactory(const std::string& desc) : description_(desc) {}

  template <>
  std::vector<std::string> ModuleFactory<epa::TwoPartonFlux>::modules() const {
    std::vector<std::string> out;
    std::transform(map_.begin(), map_.end(), std::back_inserter(out), [](const auto& val) { return val.first; });
    std::sort(out.begin(), out.end());
    return out;
  }

  template <>
  std::unique_ptr<epa::TwoPartonFlux> ModuleFactory<epa::TwoPartonFlux>::build(const ParametersList& params) const {
    CG_LOG << params;
    const auto mod_name = params.name();
    if (mod_name.empty())
      throw CG_FATAL("ModuleFactory")
          << "Failed to retrieve a flux name for the two-parton fluxes constructors lookup table.";
    return map_.at(mod_name)(params_map_.at(mod_name).validate(params));
  }
}  // namespace cepgen
