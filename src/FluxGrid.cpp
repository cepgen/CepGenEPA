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

#include <CepGen/Core/ParametersList.h>

#include <iostream>

#include "CepGenEPA/FluxGrid.h"

namespace cepgen::epa::grid {
  Header::Header(const ParametersList& params)
      : eb1(params.get<double>("eb1")),
        eb2(params.get<double>("eb2")),
        q2max1(params.get<double>("q2max1")),
        q2max2(params.get<double>("q2max2")),
        fragmenting(params.get<bool>("fragmenting")),
        parton_pdg_id(params.get<int>("partonPdgId")) {}

  bool Header::operator==(const Header& oth) const {
    // skip test of cepgen version
    return magic_number == oth.magic_number && eb1 == oth.eb1 && eb2 == oth.eb2 && q2max1 == oth.q2max1 &&
           q2max2 == oth.q2max2 && fragmenting == oth.fragmenting && parton_pdg_id == oth.parton_pdg_id;
  }

  bool Header::operator!=(const Header& oth) const { return !(*this == oth); }

  std::ostream& operator<<(std::ostream& os, const Header& header) {
    return os << "grid/Header{eb1:" << header.eb1 << ", eb2:" << header.eb2 << ", q2max1:" << header.q2max1
              << ", q2max2:" << header.q2max2 << ", fragmenting:" << std::boolalpha << header.fragmenting
              << ", parton PDGid:" << header.parton_pdg_id << ", CepGen version:" << header.cepgen_version << "}";
  }
}  // namespace cepgen::epa::grid
