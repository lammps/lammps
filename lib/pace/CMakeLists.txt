cmake_minimum_required(VERSION 3.7)  # CMake version check
project(aceevaluator)
set(CMAKE_CXX_STANDARD 11) # Enable c++11 standard


set(PACE_EVALUATOR_PATH ${CMAKE_CURRENT_LIST_DIR}/src/USER-PACE)
# message("CMakeLists.txt DEBUG: PACE_EVALUATOR_PATH=${PACE_EVALUATOR_PATH}")
set(PACE_EVALUATOR_SRC_PATH ${PACE_EVALUATOR_PATH})

FILE(GLOB PACE_EVALUATOR_SOURCE_FILES ${PACE_EVALUATOR_SRC_PATH}/*.cpp)
list(FILTER PACE_EVALUATOR_SOURCE_FILES EXCLUDE REGEX ".*pair_pace.*")
set(PACE_EVALUATOR_INCLUDE_DIR ${PACE_EVALUATOR_SRC_PATH})


##### aceevaluator #####
add_library(aceevaluator ${PACE_EVALUATOR_SOURCE_FILES})
target_include_directories(aceevaluator PUBLIC ${PACE_EVALUATOR_INCLUDE_DIR})
target_compile_options(aceevaluator PRIVATE -O3)
set_target_properties(aceevaluator PROPERTIES OUTPUT_NAME pace${LAMMPS_MACHINE})
