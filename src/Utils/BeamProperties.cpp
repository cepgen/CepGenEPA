/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2025  Laurent Forthomme
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

#include <CepGen/Modules/PartonFluxFactory.h>
#include <CepGen/PartonFluxes/CollinearFlux.h>

#include "CepGenEPA/BeamProperties.h"

using namespace cepgen::epa;

BeamProperties::BeamProperties(const ParametersList& params)
    : SteeredObject(params), energy(steer<double>("energy")), q2range(steer<Limits>("q2Range")) {
  if (const auto flux_parameters = steer<ParametersList>("flux"); !flux_parameters.empty())
    flux = CollinearFluxFactory::get().build(flux_parameters);
}

cepgen::ParametersDescription BeamProperties::description() {
  auto desc = ParametersDescription();
  desc.add("flux", ParametersDescription()).setDescription("parton-from-beam flux modelling");
  desc.add("energy", 0.).setDescription("beam particle energy, in GeV");
  desc.add("q2Range", Limits{0., 1.e5}).setDescription("parton virtuality range, in GeV^2");
  return desc;
}
