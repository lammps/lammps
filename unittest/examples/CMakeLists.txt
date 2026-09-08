# -*- CMake -*- file for tests running abbreviated versions of example inputs
#
# The folders under examples/ that demonstrate complex property calculation
# workflows with production-like settings are excluded from the regression
# tests (see EXCLUDED_FOLDERS in tools/regression-tests/run_tests.py).  The
# tests here run abbreviated versions of those workflows instead: the run
# lengths and system sizes are defined as index-style variables in the
# example inputs, so the tests can preset shortened values and then include
# the unmodified example input.
#
# Each test runs in its own working directory, since the example inputs
# write output files with fixed names into the current directory.

function(add_example_test NAME TESTNAME FOLDER)
  add_executable(test_examples_${NAME} test_examples_${NAME}.cpp)
  target_compile_definitions(test_examples_${NAME} PRIVATE -DTEST_EXAMPLES_FOLDER=${LAMMPS_DIR}/examples/${FOLDER})
  target_link_libraries(test_examples_${NAME} PRIVATE lammps GTest::GMock)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/workdir-${NAME})
  add_test(NAME ${TESTNAME}
           COMMAND test_examples_${NAME}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/workdir-${NAME})
endfunction()

add_example_test(rerun     RerunExamples     rerun)
add_example_test(diffuse   DiffuseExamples   DIFFUSE)
add_example_test(elastic   ElasticExamples   ELASTIC)
add_example_test(elastic_t ElasticTExamples  ELASTIC_T)
add_example_test(heat      HeatExamples      HEAT)
add_example_test(kappa     KappaExamples     KAPPA)
add_example_test(viscosity ViscosityExamples VISCOSITY)
add_example_test(mc_loop   MCLoopExamples    MC-LOOP)

set_tests_properties(ElasticTExamples HeatExamples KappaExamples ViscosityExamples
                     PROPERTIES LABELS "slow")
