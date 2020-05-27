###############################################################################
# Coverage
#
# Requires latest gcovr (for GCC 8.1 support):#
# pip install git+https://github.com/gcovr/gcovr.git
###############################################################################
if(ENABLE_COVERAGE)
    find_program(GCOVR_BINARY gcovr)
    find_package_handle_standard_args(GCOVR DEFAULT_MSG GCOVR_BINARY)

    if(GCOVR_FOUND)
        get_filename_component(ABSOLUTE_LAMMPS_SOURCE_DIR ${LAMMPS_SOURCE_DIR} ABSOLUTE)

        add_custom_target(
            gen_coverage_xml
            COMMAND ${GCOVR_BINARY} -s -x -r ${ABSOLUTE_LAMMPS_SOURCE_DIR} --object-directory=${CMAKE_BINARY_DIR} -o coverage.xml
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMENT "Generating XML Coverage Report..."
        )

        add_custom_target(
            gen_coverage_html
            COMMAND ${GCOVR_BINARY} -s  --html --html-details -r ${ABSOLUTE_LAMMPS_SOURCE_DIR} --object-directory=${CMAKE_BINARY_DIR} -o coverage.html
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMENT "Generating HTML Coverage Report..."
        )
    endif()
endif()
