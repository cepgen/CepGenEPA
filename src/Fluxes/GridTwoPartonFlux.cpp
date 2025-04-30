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
#include <CepGen/Core/ParametersList.h>
#include <CepGen/Utils/GridHandler.h>
#include <CepGen/Utils/Timer.h>
#include <CepGen/Version.h>

#include <fstream>
#include <iosfwd>
#include <iostream>
#include <memory>
#include <string>

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
        check_header_(steer<bool>("checkHeader")),
        header_(params_) {
    if (steer<bool>("generateGrid") || steerPath("path").empty())  // grid is not provided by the user ; build it
      buildGrid();
    loadGrid();
  }

  static ParametersDescription description() {
    auto desc = epa::TwoPartonFlux::description();
    desc.setDescription("Grid interpolator for two-parton flux");
    desc.add("modelling", ParametersDescription()).setDescription("type of flux to use to build the grid");
    desc.add("path", "flux.grid"s).setDescription("path to the interpolation grid");
    desc.add("checkHeader", true).setDescription("check the grid file header before parsing it?");
    desc.add("logW", true);
    desc.add<bool>("generateGrid", false).setDescription("(re-)generate the grid prior to run?");
    desc.add<int>("numPoints", 500).setDescription("number of points to compute for the grid construction");
    return desc;
  }

  double flux(const std::vector<double>& arguments) const override { return GridHandler<1, 1>::eval(arguments).at(0); }

  inline bool fragmenting() const override { return header_.fragmenting; }
  inline spdgid_t partonPdgId() const override { return header_.parton_pdg_id; }
  inline double mass2() const override { return 0.; }

private:
  inline void buildGrid() {
    const auto modelling = steer<ParametersList>("modelling");
    if (modelling.empty())
      throw CG_FATAL("GridTwoPartonFlux:buildGrid") << "A parton flux modelling should be provided using the "
                                                       "'modelling' parameter of this grid interpolator modelling.";
    if (modelling.name() == "grid")
      throw CG_FATAL("GridTwoPartonFlux:buildGrid") << "Cannot build a grid from a grid interpolator.";
    const auto flux_algorithm = TwoPartonFluxFactory::get().build(parameters() + modelling);
    std::ofstream output_file(grid_path_, std::ios::out | std::ios::binary);
    GridHeader header(params_);
    header.magic_number = header.goodMagic();
    cepgen::version::tag.copy(header.cepgen_version, 10);
    output_file.write(reinterpret_cast<char*>(&header), sizeof(GridHeader));
    GridValue value;
    for (const auto& w : steer<Limits>("wRange").generate(steer<int>("numPoints"), steer<bool>("logW"))) {
      value.w = w;
      value.flux = flux_algorithm->flux({w});
      CG_DEBUG("GridTwoPartonFlux") << "Adding a flux value f(" << value.w << ") = " << value.flux << ".";
      output_file.write(reinterpret_cast<char*>(&value), sizeof(GridValue));
    }
  }
  inline void loadGrid() {
    GridHeader expected_header(params_);
    expected_header.magic_number = GridHeader::goodMagic();
    cepgen::utils::Timer tmr;
    {  // file readout part
      std::ifstream file(grid_path_, std::ios::in);
      if (!file.is_open())
        throw CG_FATAL("GridTwoPartonFlux:loadGrid") << "Failed to load grid file \"" << grid_path_ << "\"!";
      file.read(reinterpret_cast<char*>(&header_), sizeof(GridHeader));
      if (header_.magic_number != GridHeader::goodMagic() || (check_header_ && header_ != expected_header))
        throw CG_FATAL("GridTwoPartonFlux:loadGrid")
            << "Invalid grid read from file.\n"
            << "   Expected header: " << expected_header << ".\n"
            << "  Retrieved header: " << header_ << ",\n"
            << "      Magic number: 0x" << std::hex << header_.magic_number << std::dec << ".";
      GridValue value;
      while (!file.eof()) {
        file.read(reinterpret_cast<char*>(&value), sizeof(GridValue));
        insert({value.w}, {value.flux});
      }
      file.close();
      initialise();  // initialise the grid after filling its nodes
    }
    CG_INFO("GridTwoPartonFlux:loadGrid") << "Two-parton flux grid evaluator built in " << tmr.elapsed() << " s.\n\t"
                                          << " w in range " << boundaries().at(0) << ".";
  }
  struct GridHeader {
    explicit GridHeader(const ParametersList& params)
        : eb1(params.get<double>("eb1")),
          eb2(params.get<double>("eb2")),
          q2max1(params.get<double>("q2max1")),
          q2max2(params.get<double>("q2max2")),
          fragmenting(params.get<bool>("fragmenting")),
          parton_pdg_id(params.get<int>("partonPdgId")) {}

    static int goodMagic() { return 0xdeadb33f; }
    bool operator==(const GridHeader& oth) const {
      // skip test of cepgen version
      return magic_number == oth.magic_number && eb1 == oth.eb1 && eb2 == oth.eb2 && q2max1 == oth.q2max1 &&
             q2max2 == oth.q2max2 && fragmenting == oth.fragmenting && parton_pdg_id == oth.parton_pdg_id;
    }
    bool operator!=(const GridHeader& oth) const { return !(*this == oth); }
    friend std::ostream& operator<<(std::ostream& os, const GridHeader& header) {
      return os << "GridHeader{eb1:" << header.eb1 << ", eb2:" << header.eb2 << ", q2max1:" << header.q2max1
                << ", q2max2:" << header.q2max2 << ", fragmenting:" << std::boolalpha << header.fragmenting
                << ", parton PDGid:" << header.parton_pdg_id << ", CepGen version:'" << header.cepgen_version << "'}";
    }

    int magic_number;
    char cepgen_version[10];
    double eb1, eb2;
    double q2max1, q2max2;
    bool fragmenting;
    int parton_pdg_id;
  } header_;
  struct GridValue {
    double w, flux;
  };

  const std::string grid_path_;
  const bool check_header_;
};
REGISTER_TWOPARTON_FLUX("grid", GridTwoPartonFlux);
