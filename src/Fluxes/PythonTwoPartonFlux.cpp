/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Hamzeh Khanpour
 *                2024-2025  Laurent Forthomme
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
#include <CepGen/PartonFluxes/CollinearFlux.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/Utils/String.h>
#include <CepGenPython/Environment.h>
#include <CepGenPython/Functional.h>

#include "CepGenEPA/BeamProperties.h"
#include "CepGenEPA/PythonUtils.h"
#include "CepGenEPA/TwoPartonFlux.h"
#include "CepGenEPA/TwoPartonFluxFactory.h"

using namespace cepgen;
using namespace std::string_literals;

class PythonTwoPartonFlux final : public epa::TwoPartonFlux {
public:
  explicit PythonTwoPartonFlux(const ParametersList& params)
      : epa::TwoPartonFlux(params),
        environment_(steer<ParametersList>("environment")),
        beam1_(steer<ParametersList>("beam1")),
        beam2_(steer<ParametersList>("beam2")),
        fragmenting_(steer<bool>("fragmenting")),
        functional_(python::make_functional(steer<std::string>("function"))) {
    if (!environment_.initialised())
      throw CG_ERROR("PythonTwoPartonFlux") << "Failed to initialise the Python environment.";
    if (!functional_)
      throw CG_ERROR("PythonTwoPartonFlux") << "Failed to retrieve the functional '" << steer<std::string>("function")
                                            << "' from the Python environment.";
    for (const auto& argument : functional_->arguments())
      arguments_names_.emplace_back(utils::toLower(argument));
  }

  static ParametersDescription description() {
    auto desc = epa::TwoPartonFlux::description();
    desc.add("environment", ParametersDescription()).setDescription("Python environment parameters");
    desc.add("beam1", epa::BeamProperties::description()).setDescription("positive-z beam properties");
    desc.add("beam2", epa::BeamProperties::description()).setDescription("negative-z beam properties");
    desc.add("fragmenting", false).setDescription("is the beam particle fragmenting after parton emission?");
    desc.add("function", ""s).setDescription("Python two-parton flux path (module.function)");
    return desc;
  }

  double flux(double w) const override {
    std::vector<double> arguments{w};
    if (arguments_names_.size() > 1) {
      for (size_t i = 1; i < arguments_names_.size(); ++i)
        if (const auto& argument = arguments_names_.at(i); argument == "eebeam"s)
          arguments.emplace_back(beam1_.energy);
        else if (argument == "pebeam"s)
          arguments.emplace_back(beam2_.energy);
    }
    const auto res = functional_->operator()(arguments);
    CG_DEBUG("PythonTwoPartonFlux:flux") << "Flux computed for arguments=" << arguments << ": " << res << ".";
    return res;
  }

  inline bool fragmenting() const override {
    if (beam1_.flux && beam2_.flux && (beam1_.flux->fragmenting() || beam2_.flux->fragmenting()))
      return true;
    return fragmenting_;
  }
  inline std::pair<spdgid_t, spdgid_t> partons() const override {
    return std::make_pair(PDG::photon, PDG::photon); /*FIXME*/
  }

private:
  const python::Environment environment_;
  const epa::BeamProperties beam1_;
  const epa::BeamProperties beam2_;
  const bool fragmenting_;
  const std::unique_ptr<python::Functional> functional_;
  std::vector<std::string> arguments_names_;
};
REGISTER_TWOPARTON_FLUX("python", PythonTwoPartonFlux);
