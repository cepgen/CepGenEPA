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

#include <CepGenPython/Error.h>
#include <CepGenPython/Functional.h>
#include <CepGenPython/ObjectPtr.h>

namespace cepgen::python {
  std::unique_ptr<Functional> functional(const std::string& python_name) {
    const auto module_path = python_name.substr(0, python_name.rfind('.')),
               function_path = python_name.substr(python_name.rfind('.') + 1);
    if (auto mod = ObjectPtr::importModule(module_path); mod) {
      CG_DEBUG("python::functional") << "Module '" << module_path << "' properly initialised. "
                                     << "Will retrieve function '" << function_path << "'.";
      if (auto func = mod.attribute(function_path); func) {
        CG_DEBUG("python::functional") << "Function '" << function_path << "' was properly initialised. "
                                       << "Attributes: " << func << ".";
        return std::make_unique<Functional>(func);
      }
      throw PY_ERROR << "Failed to retrieve a function '" << function_path << "' from Python module '" << module_path
                     << "'.";
    }
    throw PY_ERROR << "Failed to import Python module '" << module_path << "'.";
  }
}  // namespace cepgen::python
