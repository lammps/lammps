###############################################################################
# Coverage
#
# Requires latest gcovr (for GCC 8.1 support):#
# pip install git+https://github.com/gcovr/gcovr.git
#
# For Python coverage the coverage package needs to be installed
###############################################################################
if(ENABLE_COVERAGE)
    find_program(GCOVR_BINARY gcovr)
    find_package_handle_standard_args(GCOVR DEFAULT_MSG GCOVR_BINARY)

    find_program(COVERAGE_BINARY coverage)
    find_package_handle_standard_args(COVERAGE DEFAULT_MSG COVERAGE_BINARY)

    if(GCOVR_FOUND)
        get_filename_component(ABSOLUTE_LAMMPS_SOURCE_DIR ${LAMMPS_SOURCE_DIR} ABSOLUTE)

        add_custom_target(
            gen_coverage_xml
            COMMAND ${GCOVR_BINARY} -s -x -r ${ABSOLUTE_LAMMPS_SOURCE_DIR} --object-directory=${CMAKE_BINARY_DIR} -o coverage.xml
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMENT "Generating XML coverage report..."
        )

        set(COVERAGE_HTML_DIR ${CMAKE_BINARY_DIR}/coverage_html)

        add_custom_target(coverage_html_folder
            COMMAND ${CMAKE_COMMAND} -E make_directory ${COVERAGE_HTML_DIR})

        add_custom_target(
            gen_coverage_html
            COMMAND ${GCOVR_BINARY} -s  --html --html-details -r ${ABSOLUTE_LAMMPS_SOURCE_DIR} --object-directory=${CMAKE_BINARY_DIR} -o ${COVERAGE_HTML_DIR}/index.html
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMENT "Generating HTML coverage report..."
        )
        add_dependencies(gen_coverage_html coverage_html_folder)

        add_custom_target(clean_coverage_html
            ${CMAKE_COMMAND} -E remove_directory ${COVERAGE_HTML_DIR}
            COMMENT "Deleting HTML coverage report..."
        )

        add_custom_target(reset_coverage
            ${CMAKE_COMMAND} -E remove -f */*.gcda */*/*.gcda */*/*/*.gcda
                              */*/*/*/*.gcda */*/*/*/*/*.gcda */*/*/*/*/*/*.gcda
                              */*/*/*/*/*/*/*.gcda */*/*/*/*/*/*/*/*.gcda
                              */*/*/*/*/*/*/*/*/*.gcda */*/*/*/*/*/*/*/*/*/*.gcda
            WORKIND_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMENT "Deleting coverage data files..."
        )
        add_dependencies(reset_coverage clean_coverage_html)
    endif()

    if(COVERAGE_FOUND)
        set(PYTHON_COVERAGE_HTML_DIR ${CMAKE_BINARY_DIR}/python_coverage_html)
        configure_file(.coveragerc.in ${CMAKE_BINARY_DIR}/.coveragerc @ONLY)

        add_custom_command(
            OUTPUT ${CMAKE_BINARY_DIR}/unittest/python/.coverage
            COMMAND ${COVERAGE_BINARY} combine
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/unittest/python
            COMMENT "Combine Python coverage files..."
        )

        add_custom_target(
            gen_python_coverage_html
            COMMAND ${COVERAGE_BINARY} html --rcfile=${CMAKE_BINARY_DIR}/.coveragerc -d ${PYTHON_COVERAGE_HTML_DIR}
            DEPENDS ${CMAKE_BINARY_DIR}/unittest/python/.coverage ${CMAKE_BINARY_DIR}/.coveragerc
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/unittest/python
            COMMENT "Generating HTML Python coverage report..."
        )

        add_custom_target(
            gen_python_coverage_xml
            COMMAND ${COVERAGE_BINARY} xml --rcfile=${CMAKE_BINARY_DIR}/.coveragerc -o ${CMAKE_BINARY_DIR}/python_coverage.xml
            DEPENDS ${CMAKE_BINARY_DIR}/unittest/python/.coverage ${CMAKE_BINARY_DIR}/.coveragerc
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/unittest/python
            COMMENT "Generating XML Python coverage report..."
        )
    endif()
endif()
