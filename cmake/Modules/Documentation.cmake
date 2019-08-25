###############################################################################
# Build documentation
###############################################################################
option(BUILD_DOC "Build LAMMPS documentation" OFF)
if(BUILD_DOC)
  include(ProcessorCount)
  ProcessorCount(NPROCS)
  find_package(PythonInterp 3 REQUIRED)

  set(VIRTUALENV ${PYTHON_EXECUTABLE} -m virtualenv)

  file(GLOB DOC_SOURCES ${LAMMPS_DOC_DIR}/src/[^.]*.txt)
  file(GLOB PDF_EXTRA_SOURCES ${LAMMPS_DOC_DIR}/src/lammps_commands*.txt ${LAMMPS_DOC_DIR}/src/lammps_support.txt ${LAMMPS_DOC_DIR}/src/lammps_tutorials.txt)
  list(REMOVE_ITEM DOC_SOURCES ${PDF_EXTRA_SOURCES})

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

  set(RST_FILES "")
  set(RST_DIR ${CMAKE_BINARY_DIR}/rst)
  file(MAKE_DIRECTORY ${RST_DIR})
  foreach(TXT_FILE ${DOC_SOURCES})
    get_filename_component(FILENAME ${TXT_FILE} NAME_WE)
    set(RST_FILE ${RST_DIR}/${FILENAME}.rst)
    list(APPEND RST_FILES ${RST_FILE})
    add_custom_command(
      OUTPUT ${RST_FILE}
      DEPENDS requirements.txt docenv ${TXT_FILE}
      COMMAND ${DOCENV_BINARY_DIR}/txt2rst -o ${RST_DIR} ${TXT_FILE}
    )
  endforeach()

  add_custom_command(
    OUTPUT html
    DEPENDS ${RST_FILES}
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${LAMMPS_DOC_DIR}/src ${RST_DIR}
    COMMAND ${DOCENV_BINARY_DIR}/sphinx-build -j ${NPROCS} -b html -c ${LAMMPS_DOC_DIR}/utils/sphinx-config -d ${CMAKE_BINARY_DIR}/doctrees ${RST_DIR} html
  )

  add_custom_target(
    doc ALL
    DEPENDS html
    SOURCES ${LAMMPS_DOC_DIR}/utils/requirements.txt ${DOC_SOURCES}
  )

  install(DIRECTORY ${CMAKE_BINARY_DIR}/html DESTINATION ${CMAKE_INSTALL_DOCDIR})
endif()
