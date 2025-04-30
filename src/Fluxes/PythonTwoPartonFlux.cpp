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
        functional_(python::make_functional(steer<std::string>("function"))),
        eb1_(steer<double>("eb1")),
        eb2_(steer<double>("eb2")),
        q2range1_(steer<Limits>("q2Range1")),
        q2range2_(steer<Limits>("q2Range2")) {
    if (!environment_.initialised())
      throw CG_ERROR("PythonTwoPartonFlux") << "Failed to initialise the Python environment.";
    if (!functional_)
      throw CG_ERROR("PythonTwoPartonFlux") << "Failed to retrieve the functional '" << steer<std::string>("function")
                                            << "' from the Python environment.";
    arguments_names_ = functional_->arguments();
  }

  static ParametersDescription description() {
    auto desc = epa::TwoPartonFlux::description();
    desc.add("environment", ParametersDescription()).setDescription("Python environment parameters");
    desc.add<bool>("fragmenting", false).setDescription("is the beam particle fragmenting after parton emission?");
    desc.add<int>("partonPdgId", 22).setDescription("PDG id of the emitted parton");
    desc.add<double>("eb1", 7000.).setDescription("positive-z beam particle energy, in GeV");
    desc.add<double>("eb2", 7000.).setDescription("negative-z beam particle energy, in GeV");
    desc.add("q2Range1", Limits{0., 1000.}).setDescription("positive-z parton virtuality range, in GeV^2");
    desc.add("q2Range2", Limits{0., 1000.}).setDescription("negative-z parton virtuality range, in GeV^2");
    return desc;
  }

  double flux(const std::vector<double>& arguments) const override {
    return functional_->operator()(arguments);
    //return functional_->operator()(std::vector<double>{w, eb1_, eb2_, q2range1_.max(), q2range2_.max()});
  }

  inline bool fragmenting() const override { return fragmenting_; }
  inline spdgid_t partonPdgId() const override { return parton_pdg_id_; }
  inline double mass2() const override { return 0.; }

private:
  const python::Environment environment_;
  const bool fragmenting_;
  const spdgid_t parton_pdg_id_;
  const std::unique_ptr<python::Functional> functional_;
  const double eb1_, eb2_;
  const Limits q2range1_, q2range2_;
  std::vector<std::string> arguments_names_;
  mutable std::vector<double> arguments_;
};
REGISTER_TWOPARTON_FLUX("python", PythonTwoPartonFlux);
