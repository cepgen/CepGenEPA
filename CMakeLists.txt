cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(CepGenEPA)
set(CEPGEN_PATH ${CEPGEN_PATH} $ENV{CEPGEN_PATH})
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CEPGEN_PATH}/cmake ${CEPGEN_PATH}/share/cmake ${CEPGEN_PATH}/build/cmake)

find_package(GSL REQUIRED)
find_package(Boost COMPONENTS system python REQUIRED)
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(CepGen COMPONENTS python REQUIRED)

file(GLOB sources src/Fluxes/*.cpp src/Modules/*.cpp src/Processes/*.cpp src/Utils/*.cpp)
file(GLOB utils_sources utils/*.cc)

add_library(CepGenEPA SHARED ${sources})
target_link_libraries(CepGenEPA PUBLIC CepGen::CepGen CepGen::python GSL::gsl ${Boost_LIBRARIES})
target_include_directories(CepGenEPA PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${Python_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
#set_target_properties(ggMatrixElements PROPERTIES PREFIX "")
#install(TARGETS ggMatrixElements DESTINATION ${CMAKE_INSTALL_LIBDIR})

foreach(_utils_src ${utils_sources})
  get_filename_component(utils_bin ${_utils_src} NAME_WE)
  add_executable(${utils_bin} ${_utils_src})
  target_link_libraries(${utils_bin} PRIVATE CepGen::CepGen CepGenEPA)
  target_include_directories(${utils_bin} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()

# copy the input cards
file(GLOB_RECURSE input_cards RELATIVE ${PROJECT_SOURCE_DIR} cards/*)
foreach(_files ${input_cards})
  configure_file(${_files} ${_files} COPYONLY)
endforeach()
