cmake_minimum_required (VERSION 2.6)
project (growth-simulator)


set(CMAKE_C_FLAGS " -W -std=c++11 -O3")
set(CMAKE_CXX_FLAGS " -W -std=c++11 -O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )
set(SOURCE_FILES src/Experiment.cpp src/ExperimentConstantGrowth.cpp src/ExperimentIteratedGrowth.cpp
        src/GrowthModel.cpp src/GrowthModelCSC.cpp src/GrowthModelSimple.cpp src/Solver.h
        src/SolverGillespieOTL.cpp src/version.cpp)
find_package(Boost COMPONENTS system filesystem iostreams random REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})
add_executable(simulator src/simulator.cpp src/tinyxml2.cpp src/SolverGillespieOTL.cpp
        src/GrowthModelCSC.cpp src/GrowthModelSimple.cpp src/Experiment.cpp src/ExperimentIteratedGrowth.cpp
        src/ExperimentConstantGrowth.cpp)
target_link_libraries( simulator ${Boost_LIBRARIES} )


find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)
