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
#include <CepGen/Utils/Environment.h>
#include <CepGen/Utils/Math.h>
#include <CepGen/Utils/String.h>
#include <CepGenPython/Error.h>
#include <CepGenPython/Functional.h>
#include <CepGenPython/ObjectPtr.h>

using namespace cepgen;

class PythonProcess final : public cepgen::proc::Process {
public:
  explicit PythonProcess(const ParametersList& params)
      : proc::Process(params), pair_(steer<ParticleProperties>("LPAIR")) {}

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
  void prepareKinematics() override {
    environment_.reset(new python::Environment(steer<ParametersList>("environment")));
    const auto get_functional = [this](const std::string& python_name, std::shared_ptr<python::Functional>& function) {
      const auto module_path = python_name.substr(0, python_name.rfind('.')),
                 function_path = python_name.substr(python_name.rfind('.') + 1);
      if (auto mod = python::ObjectPtr::importModule(module_path); mod) {
        CG_INFO("PythonProcess") << "Module '" << module_path << "' properly initialised. Will retrieve function '"
                                 << function_path << "'.";
        if (auto func = mod.attribute(function_path); func) {
          CG_INFO("PythonProcess") << "Function '" << function_path
                                   << "' was properly initialised. Attributes: " << func << ".";
          function.reset(new python::Functional(func));
        } else
          throw PY_ERROR << "Failed to retrieve a function '" << function_path << "' from Python module '"
                         << module_path << "'.";
      } else
        throw PY_ERROR << "Failed to import Python module '" << module_path << "'.";
    };
    get_functional(steer<std::string>("process"), central_function_);
    get_functional(steer<std::string>("fluxes"), fluxes_function_);
    defineVariable(
        m_w_central_, Mapping::linear, kinematics().cuts().central.mass_sum.truncate(Limits{0., 250.}), "w_central");
  }
  double computeWeight() override {
    CG_ASSERT(central_function_);
    const auto central_weight = central_function_->operator()(m_w_central_);
    if (!utils::positive(central_weight))
      return 0.;
    const auto fluxes_weight = fluxes_function_->operator()({m_w_central_,
                                                             pA().energy(),
                                                             pB().energy(),
                                                             kinematics().cuts().initial.q2.at(0).max(),
                                                             kinematics().cuts().initial.q2.at(1).max()});
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

  static ParametersDescription description() {
    auto desc = proc::Process::description();
    desc.setDescription("Python EPA function");
    return desc;
  }

private:
  std::shared_ptr<python::Environment> environment_;
  std::shared_ptr<python::Functional> central_function_{nullptr}, fluxes_function_{nullptr};
  const ParticleProperties pair_;
  // mapped variables
  double m_w_central_{0.};
};
// register process
REGISTER_PROCESS("pythonEPA", PythonProcess);
