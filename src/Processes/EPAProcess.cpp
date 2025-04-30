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

#include <CepGen/Event/Event.h>
#include <CepGen/Modules/ProcessFactory.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/Process/Process.h>
#include <CepGen/Utils/Math.h>

#include "CepGenEPA/TwoPartonFlux.h"
#include "CepGenEPA/TwoPartonFluxFactory.h"
#include "CepGenEPA/TwoPartonProcess.h"
#include "CepGenEPA/TwoPartonProcessFactory.h"

using namespace std::string_literals;

namespace cepgen {
  class EPAProcess : public proc::Process {
  public:
    explicit EPAProcess(const ParametersList& params) : proc::Process(params) {}

    proc::ProcessPtr clone() const { return std::make_unique<EPAProcess>(parameters()); }

    static ParametersDescription description() {
      auto desc = proc::Process::description();
      desc.setDescription("Generic EPA process");
      desc.add("logW", true).setDescription("Use a logarithmic mapping of the w distribution?");
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
      const auto w_range = kinematics().cuts().central.mass_sum.truncate(Limits{1.e-9, (pA() + pB()).mass()});
      if (steer<bool>("logW"))
        defineVariable(m_w_central_, Mapping::exponential, w_range.compute(std::log), "w_central");
      else
        defineVariable(m_w_central_, Mapping::linear, w_range, "w_central");
      partons_flux_ = TwoPartonFluxFactory::get().build(ParametersList(steer<ParametersList>("partonsFlux"))
                                                            .set("eb1", pA().energy())
                                                            .set("eb2", pB().energy())
                                                            .set("wRange", w_range)
                                                            .set("q2Range1", kinematics().cuts().initial.q2.at(0))
                                                            .set("q2Range2", kinematics().cuts().initial.q2.at(1)));
      central_process_ = TwoPartonProcessFactory::get().build(steer<ParametersList>("matrixElement"));
      {  // register the list of particles in the central system
        central_system_.clear();
        for (const auto& central_particle : central_process_->centralParticles())
          central_system_.emplace_back(central_particle);
      }
    }
    double computeWeight() override {
      const auto central_weight = central_process_->matrixElement(m_w_central_);
      if (!utils::positive(central_weight))
        return 0.;
      const auto fluxes_weight = partons_flux_->flux({m_w_central_});
      if (!utils::positive(fluxes_weight))
        return 0.;
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
    std::unique_ptr<epa::TwoPartonFlux> partons_flux_;
    std::unique_ptr<epa::TwoPartonProcess> central_process_;

    spdgids_t central_system_;
    double m_w_central_{0.};  ///< central, two-parton invariant mass
  };
}  // namespace cepgen
REGISTER_PROCESS("epa", EPAProcess);
