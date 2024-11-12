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

#ifndef CepGenEPA_FluxGrid_h
#define CepGenEPA_FluxGrid_h

#include <iosfwd>
#include <memory>
#include <string>

namespace cepgen {
  class ParametersList;
}

namespace cepgen::epa::grid {
  struct Header {
    explicit Header(const ParametersList&);

    static int goodMagic() { return 0xdeadb33f; }
    bool operator==(const Header&) const;
    bool operator!=(const Header&) const;
    friend std::ostream& operator<<(std::ostream&, const Header&);

    int magic_number;
    char cepgen_version[10];
    double eb1, eb2;
    double q2max1, q2max2;
    bool fragmenting;
    int parton_pdg_id;
  };
  struct Value {
    double w, flux;
  };
}  // namespace cepgen::epa::grid

#endif
