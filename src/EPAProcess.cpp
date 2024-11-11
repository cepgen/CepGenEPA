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
#include <CepGen/Event/Event.h>
#include <CepGen/Modules/ProcessFactory.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/Process/Process.h>
#include <CepGen/Utils/Functional.h>
#include <CepGen/Utils/Math.h>
#include <CepGen/Utils/String.h>
#include <CepGenPython/Environment.h>
#include <CepGenPython/Functional.h>

#include "CepGenEPA/PythonUtils.h"
#include "CepGenEPA/TwoPartonFlux.h"
#include "CepGenEPA/TwoPartonFluxFactory.h"

using namespace std::string_literals;

namespace cepgen {
  class EPAProcess : public proc::Process {
  public:
    explicit EPAProcess(const ParametersList& params)
        : proc::Process(params), environment_(steer<ParametersList>("environment")) {}

    proc::ProcessPtr clone() const { return std::make_unique<EPAProcess>(parameters()); }

    static ParametersDescription description() {
      auto desc = proc::Process::description();
      desc.setDescription("Generic EPA process");
      auto me_description = ParametersDescription();
      me_description.add("function", ""s).setDescription("process internal name (or python functional)");
      me_description.add("centralSystem", std::vector<int>{13, -13});
      desc.add("matrixElement", me_description);
      return desc;
    }

    void addEventContent() override {
      proc::Process::setEventContent({{Particle::Role::IncomingBeam1, {PDG::electron}},
                                      {Particle::Role::IncomingBeam2, {PDG::proton}},
                                      {Particle::Role::Parton1, {PDG::photon}},
                                      {Particle::Role::Parton2, {PDG::photon}},
                                      {Particle::Role::OutgoingBeam1, {PDG::electron}},
                                      {Particle::Role::OutgoingBeam2, {PDG::proton}},
                                      {Particle::Role::CentralSystem, central_system_}});
    }
    void prepareKinematics() override {
      partons_flux_ = TwoPartonFluxFactory::get().build(ParametersList(steer<ParametersList>("partonsFlux"))
                                                            .set("eb1", pA().energy())
                                                            .set("eb2", pB().energy())
                                                            .set("q2max1", kinematics().cuts().initial.q2.at(0).max())
                                                            .set("q2max2", kinematics().cuts().initial.q2.at(1).max()));
      const auto matrix_element_definition = steer<ParametersList>("matrixElement");
      if (const auto name = matrix_element_definition.name(); name == "python") {
        central_function_ = python::functional(matrix_element_definition.get<std::string>("function"));
        const auto cs_particles = matrix_element_definition.get<std::vector<int> >("centralSystem");
        central_system_ = spdgids_t(cs_particles.begin(), cs_particles.end());
      } else
        throw CG_FATAL("EPAProcess") << "Invalid matrix element type requested: '" << name << "'.";
      defineVariable(
          m_w_central_, Mapping::linear, kinematics().cuts().central.mass_sum.truncate(Limits{0., 250.}), "w_central");
    }
    double computeWeight() override {
      const auto central_weight = central_function_->operator()({m_w_central_});
      if (!utils::positive(central_weight))
        return 0.;
      const auto fluxes_weight = partons_flux_->flux({m_w_central_});
      return central_weight * fluxes_weight;
    }
    void fillKinematics() override {
      /*pA() = Momentum(sp4vec_.vec[1].data());
      pB() = Momentum(sp4vec_.vec[0].data());
      pX() = Momentum(sp4vec_.vec[3].data());
      pY() = Momentum(sp4vec_.vec[2].data());
      pc(0) = Momentum(sp4vec_.vec[4].data());
      pc(1) = Momentum(sp4vec_.vec[5].data());*/
    }

  private:
    python::Environment environment_;
    std::unique_ptr<epa::TwoPartonFlux> partons_flux_;
    std::unique_ptr<utils::Functional> central_function_;

    spdgids_t central_system_;
    // mapped variables
    double m_w_central_{0.};
  };
}  // namespace cepgen
REGISTER_PROCESS("epa", EPAProcess);
