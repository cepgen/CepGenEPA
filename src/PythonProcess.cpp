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
#include <CepGen/Event/Event.h>
#include <CepGen/Modules/ProcessFactory.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/Process/Process.h>
#include <CepGenPython/Environment.h>
#include <CepGenPython/ObjectPtr.h>

using namespace cepgen;

class PythonProcess final : public cepgen::proc::Process {
public:
  explicit PythonProcess(const ParametersList& params)
      : proc::Process(params),
        env_(steer<ParametersList>("environment")),
        //python_function_{},
        pair_(steer<ParticleProperties>("LPAIR")) {}

  proc::ProcessPtr clone() const override { return proc::ProcessPtr(new PythonProcess(*this)); }

  void addEventContent() override {
    proc::Process::setEventContent({{Particle::Role::IncomingBeam1, {PDG::electron}},
                                    {Particle::Role::IncomingBeam2, {PDG::proton}},
                                    {Particle::Role::Parton1, {PDG::photon}},
                                    {Particle::Role::Parton2, {PDG::photon}},
                                    {Particle::Role::OutgoingBeam1, {PDG::electron}},
                                    {Particle::Role::OutgoingBeam2, {PDG::proton}},
                                    {Particle::Role::CentralSystem, {+(spdgid_t)pair_.pdgid, -(spdgid_t)pair_.pdgid}}});
  }
  void prepareKinematics() override {}
  double computeWeight() override {
    return 0.;
    //return python_function_(m_x_.data());
  }
  void fillKinematics() override {
    /*pA() = Momentum(sp4vec_.vec[1].data());
    pB() = Momentum(sp4vec_.vec[0].data());
    pX() = Momentum(sp4vec_.vec[3].data());
    pY() = Momentum(sp4vec_.vec[2].data());
    pc(0) = Momentum(sp4vec_.vec[4].data());
    pc(1) = Momentum(sp4vec_.vec[5].data());*/
  }

  static ParametersDescription description() {
    auto desc = proc::Process::description();
    desc.setDescription("Python EPA function");
    return desc;
  }

private:
  const python::Environment env_;
  //const utils::Function1D python_function_;
  const ParticleProperties pair_;
  std::array<double, 10> m_x_;  ///< mapped variables
};
// register process
REGISTER_PROCESS("pythonEPA", PythonProcess);
