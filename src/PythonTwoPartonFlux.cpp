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

#include <CepGenPython/Functional.h>

#include "CepGenEPA/PythonUtils.h"
#include "CepGenEPA/TwoPartonFlux.h"
#include "CepGenEPA/TwoPartonFluxFactory.h"

using namespace cepgen;

class PythonTwoPartonFlux final : public epa::TwoPartonFlux {
public:
  explicit PythonTwoPartonFlux(const ParametersList& params)
      : epa::TwoPartonFlux(params),
        fragmenting_(steer<bool>("fragmenting")),
        parton_pdg_id_(steer<int>("partonPdgId")),
        functional_(python::functional(steer<std::string>("function"))) {}

  static ParametersDescription description() {
    auto desc = epa::TwoPartonFlux::description();
    desc.add<bool>("fragmenting", false);
    desc.add<int>("partonPdgId", 22);
    return desc;
  }

  double flux(double w) const override { return functional_->operator()(w); }

  inline bool fragmenting() const override { return fragmenting_; }

  inline pdgid_t partonPdgId() const override { return parton_pdg_id_; }

  inline double mass2() const override { return 0.; }

private:
  const bool fragmenting_;
  const int parton_pdg_id_;
  const std::unique_ptr<python::Functional> functional_;
};
REGISTER_TWOPARTON_FLUX("python", PythonTwoPartonFlux);
