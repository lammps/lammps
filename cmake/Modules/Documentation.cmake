###############################################################################
# Build documentation
###############################################################################
option(BUILD_DOC "Build LAMMPS HTML documentation" OFF)

if(BUILD_DOC)
  # Sphinx 3.x requires at least Python 3.5
  if(CMAKE_VERSION VERSION_LESS 3.12)
    find_package(PythonInterp 3.5 REQUIRED)
    set(VIRTUALENV ${PYTHON_EXECUTABLE} -m venv)
  else()
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    if(Python3_VERSION VERSION_LESS 3.5)
      message(FATAL_ERROR "Python 3.5 and up is required to build the HTML documentation")
    endif()
    set(VIRTUALENV ${Python3_EXECUTABLE} -m venv)
  endif()
  find_package(Doxygen 1.8.10 REQUIRED)

  file(GLOB DOC_SOURCES ${CONFIGURE_DEPENDS} ${LAMMPS_DOC_DIR}/src/[^.]*.rst)


  add_custom_command(
    OUTPUT docenv
    COMMAND ${VIRTUALENV} docenv
  )

  set(DOCENV_BINARY_DIR ${CMAKE_BINARY_DIR}/docenv/bin)
  set(DOCENV_REQUIREMENTS_FILE ${LAMMPS_DOC_DIR}/utils/requirements.txt)

  set(SPHINX_CONFIG_DIR ${LAMMPS_DOC_DIR}/utils/sphinx-config)
  set(SPHINX_CONFIG_FILE_TEMPLATE ${SPHINX_CONFIG_DIR}/conf.py.in)
  set(SPHINX_STATIC_DIR  ${SPHINX_CONFIG_DIR}/_static)

  # configuration and static files are copied to binary dir to avoid collisions with parallel builds
  set(DOC_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/doc)
  set(DOC_BUILD_CONFIG_FILE ${DOC_BUILD_DIR}/conf.py)
  set(DOC_BUILD_STATIC_DIR ${DOC_BUILD_DIR}/_static)
  set(DOXYGEN_BUILD_DIR ${DOC_BUILD_DIR}/doxygen)
  set(DOXYGEN_XML_DIR ${DOXYGEN_BUILD_DIR}/xml)

  # copy entire configuration folder to doc build directory
  # files in _static are automatically copied during sphinx-build, so no need to copy them individually
  file(COPY ${SPHINX_CONFIG_DIR}/ DESTINATION ${DOC_BUILD_DIR})

  # configure paths in conf.py, since relative paths change when file is copied
  configure_file(${SPHINX_CONFIG_FILE_TEMPLATE} ${DOC_BUILD_CONFIG_FILE})

  add_custom_command(
    OUTPUT ${DOC_BUILD_DIR}/requirements.txt
    DEPENDS docenv ${DOCENV_REQUIREMENTS_FILE}
    COMMAND ${CMAKE_COMMAND} -E copy ${DOCENV_REQUIREMENTS_FILE} ${DOC_BUILD_DIR}/requirements.txt
    COMMAND ${DOCENV_BINARY_DIR}/pip $ENV{PIP_OPTIONS} install --upgrade pip
    COMMAND ${DOCENV_BINARY_DIR}/pip $ENV{PIP_OPTIONS} install --upgrade ${LAMMPS_DOC_DIR}/utils/converters
    COMMAND ${DOCENV_BINARY_DIR}/pip $ENV{PIP_OPTIONS} install -r ${DOC_BUILD_DIR}/requirements.txt --upgrade
  )

  set(MATHJAX_URL "https://github.com/mathjax/MathJax/archive/3.1.3.tar.gz" CACHE STRING "URL for MathJax tarball")
  set(MATHJAX_MD5 "d1c98c746888bfd52ca8ebc10704f92f" CACHE STRING "MD5 checksum of MathJax tarball")
  mark_as_advanced(MATHJAX_URL)

  # download mathjax distribution and unpack to folder "mathjax"
  if(NOT EXISTS ${DOC_BUILD_STATIC_DIR}/mathjax/es5)
    file(DOWNLOAD ${MATHJAX_URL}
      "${CMAKE_CURRENT_BINARY_DIR}/mathjax.tar.gz"
      EXPECTED_MD5 ${MATHJAX_MD5})
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf mathjax.tar.gz WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    file(GLOB MATHJAX_VERSION_DIR ${CONFIGURE_DEPENDS} ${CMAKE_CURRENT_BINARY_DIR}/MathJax-*)
    execute_process(COMMAND ${CMAKE_COMMAND} -E rename ${MATHJAX_VERSION_DIR} ${DOC_BUILD_STATIC_DIR}/mathjax)
  endif()

  # set up doxygen and add targets to run it
  file(MAKE_DIRECTORY ${DOXYGEN_BUILD_DIR})
  file(COPY ${LAMMPS_DOC_DIR}/doxygen/lammps-logo.png DESTINATION ${DOXYGEN_BUILD_DIR}/lammps-logo.png)
  configure_file(${LAMMPS_DOC_DIR}/doxygen/Doxyfile.in ${DOXYGEN_BUILD_DIR}/Doxyfile)
  get_target_property(LAMMPS_SOURCES lammps SOURCES)
  add_custom_command(
    OUTPUT ${DOXYGEN_XML_DIR}/index.xml
    DEPENDS ${DOC_SOURCES} ${LAMMPS_SOURCES}
    COMMAND Doxygen::doxygen ${DOXYGEN_BUILD_DIR}/Doxyfile WORKING_DIRECTORY ${DOXYGEN_BUILD_DIR}
    COMMAND ${CMAKE_COMMAND} -E touch ${DOXYGEN_XML_DIR}/run.stamp
  )

  if(EXISTS ${DOXYGEN_XML_DIR}/run.stamp)
    set(SPHINX_EXTRA_OPTS "-E")
  else()
    set(SPHINX_EXTRA_OPTS "")
  endif()
  add_custom_command(
    OUTPUT html
    DEPENDS ${DOC_SOURCES} docenv ${DOC_BUILD_DIR}/requirements.txt ${DOXYGEN_XML_DIR}/index.xml ${BUILD_DOC_CONFIG_FILE}
    COMMAND ${DOCENV_BINARY_DIR}/sphinx-build ${SPHINX_EXTRA_OPTS} -b html -c ${DOC_BUILD_DIR} -d ${DOC_BUILD_DIR}/doctrees ${LAMMPS_DOC_DIR}/src ${DOC_BUILD_DIR}/html
    COMMAND ${CMAKE_COMMAND} -E create_symlink Manual.html ${DOC_BUILD_DIR}/html/index.html
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${LAMMPS_DOC_DIR}/src/PDF ${DOC_BUILD_DIR}/html/PDF
    COMMAND ${CMAKE_COMMAND} -E remove -f ${DOXYGEN_XML_DIR}/run.stamp
  )

  add_custom_target(
    doc ALL
    DEPENDS html ${DOC_BUILD_STATIC_DIR}/mathjax/es5
    SOURCES ${LAMMPS_DOC_DIR}/utils/requirements.txt ${DOC_SOURCES}
  )

  install(DIRECTORY ${DOC_BUILD_DIR}/html DESTINATION ${CMAKE_INSTALL_DOCDIR})
endif()
