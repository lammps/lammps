set(PACE_EVALUATOR_PATH ${LAMMPS_LIB_SOURCE_DIR}/pace)
message("CMakeLists.txt DEBUG: PACE_EVALUATOR_PATH=${PACE_EVALUATOR_PATH}")
set(PACE_EVALUATOR_SRC_PATH ${PACE_EVALUATOR_PATH})

FILE(GLOB PACE_EVALUATOR_SOURCE_FILES ${PACE_EVALUATOR_SRC_PATH}/*.cpp)
set(PACE_EVALUATOR_INCLUDE_DIR ${PACE_EVALUATOR_SRC_PATH})


##### aceevaluator #####
add_library(aceevaluator ${PACE_EVALUATOR_SOURCE_FILES})
target_include_directories(aceevaluator PUBLIC ${PACE_EVALUATOR_INCLUDE_DIR})
target_compile_options(aceevaluator PRIVATE -O2)
set_target_properties(aceevaluator PROPERTIES OUTPUT_NAME lammps_pace${LAMMPS_MACHINE})
target_link_libraries(lammps PRIVATE aceevaluator)