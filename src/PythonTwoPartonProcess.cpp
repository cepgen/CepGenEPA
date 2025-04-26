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

#include <CepGenPython/Environment.h>
#include <CepGenPython/Functional.h>

#include "CepGenEPA/PythonUtils.h"
#include "CepGenEPA/TwoPartonProcess.h"
#include "CepGenEPA/TwoPartonProcessFactory.h"

using namespace cepgen;
using namespace std::string_literals;

class PythonTwoPartonProcess : public epa::TwoPartonProcess {
public:
  explicit PythonTwoPartonProcess(const ParametersList& params)
      : epa::TwoPartonProcess(params),
        environment_(steer<ParametersList>("environment")),
        central_function_(python::functional(steer<std::string>("function"))) {}

  static ParametersDescription description() {
    auto desc = epa::TwoPartonProcess::description();
    desc.setDescription("Python two-parton process");
    desc.add("function", ""s).setDescription("Python functional used for matrix element computation");
    return desc;
  }

  std::string processDescription() const override { return "Python process"; }  //FIXME
  double matrixElement(double w) const override { return central_function_->operator()({w}); }

private:
  python::Environment environment_;
  std::unique_ptr<python::Functional> central_function_;
};
REGISTER_TWOPARTON_PROCESS("python", PythonTwoPartonProcess);
