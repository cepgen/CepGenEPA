cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(CepGenEPA)
set(CEPGEN_PATH ${CEPGEN_PATH} $ENV{CEPGEN_PATH})
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CEPGEN_PATH}/cmake)

find_package(GSL REQUIRED)
find_package(Boost COMPONENTS system python REQUIRED)
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(CepGen REQUIRED)

file(GLOB SOURCES src/*.cpp)

add_library(CepGenEPA SHARED ${SOURCES})
target_link_libraries(CepGenEPA PRIVATE CepGen::CepGen GSL::gsl ${Boost_LIBRARIES})
target_include_directories(CepGenEPA PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${Python_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
#set_target_properties(ggMatrixElements PROPERTIES PREFIX "")
#install(TARGETS ggMatrixElements DESTINATION ${CMAKE_INSTALL_LIBDIR})
