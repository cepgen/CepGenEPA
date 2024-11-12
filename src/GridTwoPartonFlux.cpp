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
#include <CepGen/Utils/GridHandler.h>
#include <CepGen/Utils/Timer.h>
#include <CepGen/Version.h>

#include <fstream>

#include "CepGenEPA/FluxGrid.h"
#include "CepGenEPA/TwoPartonFlux.h"
#include "CepGenEPA/TwoPartonFluxFactory.h"

using namespace cepgen;
using namespace std::string_literals;

class GridTwoPartonFlux final : public epa::TwoPartonFlux, private GridHandler<1, 1> {
public:
  explicit GridTwoPartonFlux(const ParametersList& params)
      : epa::TwoPartonFlux(params),
        GridHandler<1, 1>(GridType::linear),
        grid_path_(steerPath("path").empty() ? "flux.grid" : steerPath("path")),
        header_(params_) {
    if (steer<bool>("generateGrid") || steerPath("path").empty())  // grid is not provided by the user ; build it
      buildGrid();
    loadGrid();
  }

  static ParametersDescription description() {
    auto desc = epa::TwoPartonFlux::description();
    desc.setDescription("Grid interpolator for two-parton flux");
    desc.add("modelling", ParametersDescription()).setDescription("type of flux to use to build the grid");
    desc.add("path", ""s).setDescription("path to the interpolation grid");
    desc.add("wRange", Limits{1.e-9, 1.e3});
    desc.add("logW", true);
    desc.add<bool>("generateGrid", false).setDescription("(re-)generate the grid prior to run?");
    desc.add<int>("numPoints", 100).setDescription("number of points to compute for the grid construction");
    return desc;
  }

  double flux(double w) const override { return GridHandler<1, 1>::eval({w}).at(0); }

  inline bool fragmenting() const override { return header_.fragmenting; }
  inline pdgid_t partonPdgId() const override { return header_.parton_pdg_id; }
  inline double mass2() const override { return 0.; }

private:
  inline void buildGrid() {
    const auto flux_algorithm = TwoPartonFluxFactory::get().build(parameters() + steer<ParametersList>("modelling"));
    std::ofstream output_file(grid_path_, std::ios::out | std::ios::binary);
    epa::grid::Header header(params_);
    header.magic_number = header.goodMagic();
    cepgen::version::tag.copy(header.cepgen_version, 10);
    output_file.write(reinterpret_cast<char*>(&header), sizeof(epa::grid::Header));
    epa::grid::Value value;
    for (const auto& w : steer<Limits>("wRange").generate(steer<int>("numPoints"), steer<bool>("logW"))) {
      value.w = w;
      value.flux = flux_algorithm->flux(w);
      CG_LOG << "Adding a flux value f(" << value.w << ") = " << value.flux << ".";
      output_file.write(reinterpret_cast<char*>(&value), sizeof(epa::grid::Value));
    }
  }
  inline void loadGrid() {
    epa::grid::Header expected_header(params_);
    expected_header.magic_number = epa::grid::Header::goodMagic();
    cepgen::utils::Timer tmr;
    {  // file readout part
      std::ifstream file(grid_path_, std::ios::in);
      if (!file.is_open())
        throw CG_FATAL("GridTwoPartonFlux:loadGrid") << "Failed to load grid file \"" << grid_path_ << "\"!";
      file.read(reinterpret_cast<char*>(&header_), sizeof(epa::grid::Header));
      if (header_.magic_number != epa::grid::Header::goodMagic() || header_ != expected_header)
        throw CG_FATAL("GridTwoPartonFlux:loadGrid")
            << "Invalid grid read from file.\n"
            << "          Expected header: " << expected_header << ".\n"
            << "         Retrieved header: " << header_ << ",\n"
            << "             Magic number: 0x" << std::hex << header_.magic_number << std::dec << ",\n"
            << "  CepGen version for file: " << header_.cepgen_version << ".";
      epa::grid::Value value;
      while (!file.eof()) {
        file.read(reinterpret_cast<char*>(&value), sizeof(epa::grid::Value));
        insert({value.w}, {value.flux});
      }
      file.close();
      initialise();  // initialise the grid after filling its nodes
    }
    CG_INFO("GridTwoPartonFlux:loadGrid") << "Two-parton flux grid evaluator built in " << tmr.elapsed() << " s.\n\t"
                                          << " w in range " << boundaries().at(0) << ".";
  }

  const std::string grid_path_;
  epa::grid::Header header_;
};
REGISTER_TWOPARTON_FLUX("grid", GridTwoPartonFlux);
