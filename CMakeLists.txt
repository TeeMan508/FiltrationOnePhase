cmake_minimum_required(VERSION 3.28)
project(FiltrationOnePhase)
set(CMAKE_CXX_STANDARD 17)

find_package(inmost REQUIRED)

message(${INMOST_LIBRARY_DIRS})
link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})

#add_executable(FiltrationOnePhase main.cpp
#        utils/params/params.h
#        utils/params/params.cpp
#        model/darcy_one_phase.h
#        model/darcy_one_phase.cpp
#        utils/calculations/calculations.cpp
#        utils/calculations/calculations.h
#        solver/solver.cpp
#        solver/solver.h
#)

add_executable(FiltrationOnePhase main.cpp
        model/one_phase_model.h
        model/one_phase_model.cpp)

add_subdirectory(utils/grid_gen)
target_link_libraries(FiltrationOnePhase ${INMOST_LIBRARIES})

