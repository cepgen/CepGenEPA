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

#ifndef CepGenEPA_TwoPartonProcess_h
#define CepGenEPA_TwoPartonProcess_h

#include <CepGen/Modules/NamedModule.h>

namespace cepgen::epa {
  /// Base object for a collinear two-parton-level process implementation
  class TwoPartonProcess : public NamedModule<TwoPartonProcess> {
  public:
    explicit TwoPartonProcess(const ParametersList& params)
        : NamedModule<TwoPartonProcess>(params), central_system_particles_(steer<std::vector<int> >("centralSystem")) {}

    static ParametersDescription description() {
      auto desc = ParametersDescription();
      desc.add("centralSystem", std::vector<int>{13, -13});
      return desc;
    }

    /// LaTeX-like description of the process
    virtual std::string processDescription() const = 0;
    /// Compute the collinear matrix element for this central mass w
    virtual double matrixElement(double w) const = 0;
    /// Retrieve the list of particles produced in the process
    virtual std::vector<int> centralParticles() const { return central_system_particles_; }

  protected:
    std::vector<int> central_system_particles_;
  };
}  // namespace cepgen::epa

#endif
