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

#include <CepGen/Core/Exception.h>
#include <CepGenPython/Environment.h>
#include <CepGenPython/Functional.h>

#include "CepGenEPA/PythonUtils.h"
#include "CepGenEPA/TwoPartonFlux.h"
#include "CepGenEPA/TwoPartonFluxFactory.h"

using namespace cepgen;

class PythonTwoPartonFlux final : public epa::TwoPartonFlux {
public:
  explicit PythonTwoPartonFlux(const ParametersList& params)
      : epa::TwoPartonFlux(params),
        environment_(steer<ParametersList>("environment")),
        fragmenting_(steer<bool>("fragmenting")),
        parton_pdg_id_(steer<int>("partonPdgId")),
        functional_(python::functional(steer<std::string>("function"))),
        eb1_(steer<double>("eb1")),
        eb2_(steer<double>("eb2")),
        q2max1_(steer<double>("q2max1")),
        q2max2_(steer<double>("q2max2")) {
    if (!environment_.initialised())
      throw CG_ERROR("PythonTwoPartonFlux") << "Failed to initialise the Python environment.";
    if (!functional_)
      throw CG_ERROR("PythonTwoPartonFlux") << "Failed to retrieve the functional '" << steer<std::string>("function")
                                            << "' from the Python environment.";
  }

  static ParametersDescription description() {
    auto desc = epa::TwoPartonFlux::description();
    desc.add("environment", ParametersDescription()).setDescription("Python environment parameters");
    desc.add<bool>("fragmenting", false).setDescription("is the beam particle fragmenting after parton emission?");
    desc.add<int>("partonPdgId", 22).setDescription("PDG id of the emitted parton");
    desc.add<double>("eb1", 7000.).setDescription("positive-z beam particle energy, in GeV");
    desc.add<double>("q2max1", 1000.).setDescription("maximum positive-z parton virtuality, in GeV^2");
    desc.add<double>("eb2", 7000.).setDescription("negative-z beam particle energy, in GeV");
    desc.add<double>("q2max2", 1000.).setDescription("maximum negative-z parton virtuality, in GeV^2");
    return desc;
  }

  double flux(double w) const override {
    return functional_->operator()(std::vector<double>{w, eb1_, eb2_, q2max1_, q2max2_});
  }

  inline bool fragmenting() const override { return fragmenting_; }
  inline pdgid_t partonPdgId() const override { return parton_pdg_id_; }
  inline double mass2() const override { return 0.; }

private:
  const python::Environment environment_;
  const bool fragmenting_;
  const int parton_pdg_id_;
  const std::unique_ptr<python::Functional> functional_;
  const double eb1_, eb2_, q2max1_, q2max2_;
};
REGISTER_TWOPARTON_FLUX("python", PythonTwoPartonFlux);
