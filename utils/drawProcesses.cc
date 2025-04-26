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

#include "CepGenEPA/TwoPartonProcess.h"
#include "CepGenEPA/TwoPartonProcessFactory.h"

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
          "modellings,m", "process modellings", &modellings, cepgen::TwoPartonProcessFactory::get().modules())
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .addOptionalArgument("range,r", "x-axis range", &range, cepgen::Limits{1.e-6, 1000.})
      .addOptionalArgument("num-points,n", "number of points to plot", &num_points, 100)
      .addOptionalArgument("logx", "logarithmic x-scale", &logx, false)
      .addOptionalArgument("logy,l", "logarithmic y-scale", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .parse();

  vector<cepgen::utils::Graph1D> graphs;
  for (const auto& modelling : modellings) {
    const auto process = cepgen::TwoPartonProcessFactory::get().build(cepgen::ParametersList().setName(modelling));
    auto& graph = graphs.emplace_back();
    for (const auto& wgg : range.generate(num_points, logx))
      graph.addPoint(wgg, process->matrixElement(wgg));
    graph.setTitle(process->processDescription());
    graph.xAxis().setLabel("$w_{\\gamma\\gamma}$ (GeV)");
    graph.yAxis().setLabel("$\\sigma_{\\gamma\\gamma}$ (pb)");
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
    plot->draw(collection, "comparison_processes", "", dm);
  }
  return 0;
}
