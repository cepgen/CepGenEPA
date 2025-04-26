/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2025  Laurent Forthomme
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

#include <CepGen/Generator.h>
#include <CepGen/Modules/DrawerFactory.h>
#include <CepGen/Utils/ArgumentsParser.h>
#include <CepGen/Utils/Drawer.h>
#include <CepGen/Utils/Graph.h>

#include "CepGenEPA/TwoPartonFlux.h"
#include "CepGenEPA/TwoPartonFluxFactory.h"

using namespace std;

int main(int argc, char* argv[]) {
  vector<string> modellings;
  string plotter;
  cepgen::Limits range;
  int num_points;
  bool logx, logy, draw_grid;
  cepgen::initialise();
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument(
          "modellings,m", "flux modellings", &modellings, cepgen::TwoPartonFluxFactory::get().modules())
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .addOptionalArgument("range,r", "x-axis range", &range, cepgen::Limits{1.e-6, 1000.})
      .addOptionalArgument("num-points,n", "number of points to plot", &num_points, 100)
      .addOptionalArgument("logx", "logarithmic x-scale", &logx, false)
      .addOptionalArgument("logy,l", "logarithmic y-scale", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .parse();

  vector<cepgen::utils::Graph1D> graphs;
  for (const auto& modelling : modellings) {
    const auto partons_flux = cepgen::TwoPartonFluxFactory::get().build(
        cepgen::ParametersList().setName(modelling).set("checkHeader", false));
    auto& graph = graphs.emplace_back();
    for (const auto& wgg : range.generate(num_points, logx))
      graph.addPoint(wgg, partons_flux->flux({wgg}));
    graph.setTitle(modelling);
    graph.xAxis().setLabel("$w_{\\gamma\\gamma}$ (GeV)");
    graph.yAxis().setLabel("$S_{\\gamma\\gamma}$ (GeV${}^{-1}$)");
  }

  if (!plotter.empty()) {
    auto plot = cepgen::DrawerFactory::get().build(plotter);
    cepgen::utils::Drawer::Mode dm;
    if (logx)
      dm |= cepgen::utils::Drawer::Mode::logx;
    if (logy)
      dm |= cepgen::utils::Drawer::Mode::logy;
    if (draw_grid)
      dm |= cepgen::utils::Drawer::Mode::grid;

    cepgen::utils::DrawableColl collection;
    for (const auto& graph : graphs)
      collection.emplace_back(&graph);
    plot->draw(collection, "comparison_fluxes", "", dm);
  }
  return 0;
}
