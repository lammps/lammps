###############################################################################
# Build documentation
###############################################################################
option(BUILD_DOC "Build LAMMPS documentation" OFF)
if(BUILD_DOC)
  include(ProcessorCount)
  ProcessorCount(NPROCS)
  find_package(PythonInterp 3 REQUIRED)

  set(VIRTUALENV ${PYTHON_EXECUTABLE} -m virtualenv)

  file(GLOB DOC_SOURCES ${LAMMPS_DOC_DIR}/src/[^.]*.rst)

  add_custom_command(
    OUTPUT docenv
    COMMAND ${VIRTUALENV} docenv
  )

  set(DOCENV_BINARY_DIR ${CMAKE_BINARY_DIR}/docenv/bin)

  add_custom_command(
    OUTPUT requirements.txt
    DEPENDS docenv
    COMMAND ${CMAKE_COMMAND} -E copy ${LAMMPS_DOC_DIR}/utils/requirements.txt requirements.txt
    COMMAND ${DOCENV_BINARY_DIR}/pip install -r requirements.txt --upgrade
    COMMAND ${DOCENV_BINARY_DIR}/pip install --upgrade ${LAMMPS_DOC_DIR}/utils/converters
  )

  add_custom_command(
    OUTPUT html
    DEPENDS ${DOC_SOURCES} docenv requirements.txt
    COMMAND ${DOCENV_BINARY_DIR}/sphinx-build -j ${NPROCS} -b html -c ${LAMMPS_DOC_DIR}/utils/sphinx-config -d ${CMAKE_BINARY_DIR}/doctrees ${LAMMPS_DOC_DIR}/src html
  )

  add_custom_target(
    doc ALL
    DEPENDS html
    SOURCES ${LAMMPS_DOC_DIR}/utils/requirements.txt ${DOC_SOURCES}
  )

  install(DIRECTORY ${CMAKE_BINARY_DIR}/html DESTINATION ${CMAKE_INSTALL_DOCDIR})
endif()
