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

#ifndef CepGenEPA_TwoPartonFlux_h
#define CepGenEPA_TwoPartonFlux_h

#include <CepGen/PartonFluxes/PartonFlux.h>

namespace cepgen::epa {
  /// Base object for a collinear parton flux parameterisation
  class TwoPartonFlux : public PartonFlux {
  public:
    explicit TwoPartonFlux(const ParametersList& params) : PartonFlux(params) {}

    static ParametersDescription description() {
      auto desc = PartonFlux::description();
      desc.setDescription("Two-parton mass-dependent flux");
      return desc;
    }

    virtual std::pair<spdgid_t, spdgid_t> partons() const = 0;  ///< List of partons emitted by the two-beam system
    virtual double flux(const std::vector<double>&) const = 0;  ///< Compute the collinear flux for this point

    // replace all PartonFlux pure virtual (and unused) attributes
    inline bool ktFactorised() const final { return false; }
    inline spdgid_t partonPdgId() const final { return 0; }
    inline double mass2() const final { return 0.; }
  };
}  // namespace cepgen::epa

#endif
