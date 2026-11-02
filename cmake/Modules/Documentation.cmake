###############################################################################
# Build documentation
###############################################################################
option(BUILD_DOC "Build LAMMPS HTML documentation" OFF)

if(BUILD_DOC)
  option(BUILD_DOC_VENV "Build LAMMPS documentation virtual environment" ON)
  mark_as_advanced(BUILD_DOC_VENV)
  # Current Sphinx versions require at least Python 3.8
  # use default (or custom) Python executable, if version is sufficient
  if(Python_VERSION VERSION_GREATER_EQUAL 3.8)
    set(Python3_EXECUTABLE ${Python_EXECUTABLE})
  endif()
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  if(Python3_VERSION VERSION_LESS 3.8)
    message(FATAL_ERROR "Python 3.8 and up is required to build the LAMMPS HTML documentation")
  endif()
  set(VIRTUALENV ${Python3_EXECUTABLE} -m venv)

  find_package(Doxygen 1.8.10 REQUIRED)
  file(GLOB DOC_SOURCES CONFIGURE_DEPENDS ${LAMMPS_DOC_DIR}/src/[^.]*.rst)

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
  # skip a MathJax checkout created by "make html" in the doc folder. it may be a different version
  file(COPY ${SPHINX_CONFIG_DIR}/ DESTINATION ${DOC_BUILD_DIR} REGEX "/_static/mathjax$" EXCLUDE)

  # configure paths in conf.py, since relative paths change when file is copied
  configure_file(${SPHINX_CONFIG_FILE_TEMPLATE} ${DOC_BUILD_CONFIG_FILE})

  if(BUILD_DOC_VENV)
    add_custom_command(
      OUTPUT docenv
      COMMAND ${VIRTUALENV} docenv
    )

    set(DOCENV_BINARY_DIR ${CMAKE_BINARY_DIR}/docenv/bin)
    set(DOCENV_REQUIREMENTS_FILE ${LAMMPS_DOC_DIR}/utils/requirements.txt)

    add_custom_command(
      OUTPUT ${DOC_BUILD_DIR}/requirements.txt
      DEPENDS docenv ${DOCENV_REQUIREMENTS_FILE}
      COMMAND ${CMAKE_COMMAND} -E copy ${DOCENV_REQUIREMENTS_FILE} ${DOC_BUILD_DIR}/requirements.txt
      COMMAND ${DOCENV_BINARY_DIR}/pip $ENV{PIP_OPTIONS} install --upgrade pip
      COMMAND ${DOCENV_BINARY_DIR}/pip $ENV{PIP_OPTIONS} install -r ${DOC_BUILD_DIR}/requirements.txt --upgrade
    )

    set(DOCENV_DEPS docenv ${DOC_BUILD_DIR}/requirements.txt)
    if(NOT TARGET Sphinx::sphinx-build)
      add_executable(Sphinx::sphinx-build IMPORTED GLOBAL)
      set_target_properties(Sphinx::sphinx-build PROPERTIES IMPORTED_LOCATION "${DOCENV_BINARY_DIR}/sphinx-build")
    endif()
  else()
    find_package(Sphinx)
  endif()

  # the MathJax version and checksum must be kept in sync with the MATHJAXTAG and MATHJAXSUM settings
  # in doc/Makefile. the version must be compatible with mathjax_path in doc/utils/sphinx-config/conf.py.in
  SetDownloadSettings(MATHJAX "MathJax"
    "https://github.com/mathjax/MathJax/archive/4.1.3.tar.gz"
    "f487c39d2913f371eb42dab078559a902da69acca38a9ebec7640a6581535ba3")
  GetFallbackURL(MATHJAX_URL MATHJAX_FALLBACK)

  # the MathJax fonts are distributed separately and their version must match the MathJax version.
  # the checksum must be kept in sync with the MATHJAXFONTSUM setting in doc/Makefile.
  # the archive has a unique name, so the fallback uses the same file name
  SetDownloadSettings(MATHJAX_FONT "MathJax fonts"
    "https://registry.npmjs.org/@mathjax/mathjax-newcm-font/-/mathjax-newcm-font-4.1.3.tgz"
    "87d7b869c6a2a6169d9a53acc4eab6c846a9cbe11752738226461bb5070c8b88")
  cmake_path(GET MATHJAX_FONT_URL FILENAME MATHJAX_FONT_FILE)
  set(MATHJAX_FONT_FALLBACK "${LAMMPS_THIRDPARTY_URL}/${MATHJAX_FONT_FILE}")

  # download an archive to the build folder unless it is already there with a matching checksum.
  # the result variable tells whether there is a new archive, e.g. after a MathJax version change
  function(FetchDocArchive url fallback sha256 archive result)
    set(DL_SHA256 "")
    if(EXISTS ${archive})
      file(SHA256 ${archive} DL_SHA256)
    endif()
    if("${DL_SHA256}" STREQUAL "${sha256}")
      set(${result} FALSE PARENT_SCOPE)
    else()
      file(DOWNLOAD ${url} ${archive} STATUS DL_STATUS SHOW_PROGRESS)
      file(SHA256 ${archive} DL_SHA256)
      if((NOT DL_STATUS EQUAL 0) OR (NOT "${DL_SHA256}" STREQUAL "${sha256}"))
        message(WARNING "Download from primary URL ${url} failed\nTrying fallback URL ${fallback}")
        file(DOWNLOAD ${fallback} ${archive} EXPECTED_HASH SHA256=${sha256} SHOW_PROGRESS)
      endif()
      set(${result} TRUE PARENT_SCOPE)
    endif()
  endfunction()

  # download mathjax distribution and unpack to folder "mathjax"
  set(MATHJAX_ARCHIVE ${CMAKE_CURRENT_BINARY_DIR}/mathjax.tar.gz)
  set(MATHJAX_DIR ${DOC_BUILD_STATIC_DIR}/mathjax)
  FetchDocArchive(${MATHJAX_URL} "${MATHJAX_FALLBACK}" ${MATHJAX_SHA256} ${MATHJAX_ARCHIVE} MATHJAX_NEW)
  if(MATHJAX_NEW OR (NOT EXISTS ${MATHJAX_DIR}/tex-mml-chtml.js))
    # remove the previously unpacked version and leftovers from incomplete previous attempts
    file(GLOB MATHJAX_VERSION_DIR ${CMAKE_CURRENT_BINARY_DIR}/MathJax-*)
    file(REMOVE_RECURSE ${MATHJAX_DIR} ${MATHJAX_VERSION_DIR})
    file(ARCHIVE_EXTRACT INPUT ${MATHJAX_ARCHIVE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(GLOB MATHJAX_VERSION_DIR ${CMAKE_CURRENT_BINARY_DIR}/MathJax-*)
    file(RENAME ${MATHJAX_VERSION_DIR} ${MATHJAX_DIR})
  endif()

  # download the fonts and unpack the files for HTML output to a subfolder of the "mathjax" folder,
  # so that they are not loaded from the internet. the location is set with mathjax_font_config
  # in doc/utils/sphinx-config/conf.py.in
  set(MATHJAX_FONT_ARCHIVE ${CMAKE_CURRENT_BINARY_DIR}/mathjax-font.tar.gz)
  set(MATHJAX_FONT_DIR ${MATHJAX_DIR}/mathjax-newcm-font)
  FetchDocArchive(${MATHJAX_FONT_URL} "${MATHJAX_FONT_FALLBACK}" ${MATHJAX_FONT_SHA256} ${MATHJAX_FONT_ARCHIVE} MATHJAX_FONT_NEW)
  if(MATHJAX_FONT_NEW OR (NOT EXISTS ${MATHJAX_FONT_DIR}/chtml/woff2))
    set(MATHJAX_FONT_UNPACK_DIR ${CMAKE_CURRENT_BINARY_DIR}/mathjax-font)
    file(REMOVE_RECURSE ${MATHJAX_FONT_DIR} ${MATHJAX_FONT_UNPACK_DIR})
    file(ARCHIVE_EXTRACT INPUT ${MATHJAX_FONT_ARCHIVE} DESTINATION ${MATHJAX_FONT_UNPACK_DIR}
      PATTERNS "package/chtml" "package/package.json")
    file(RENAME ${MATHJAX_FONT_UNPACK_DIR}/package ${MATHJAX_FONT_DIR})
    file(REMOVE_RECURSE ${MATHJAX_FONT_UNPACK_DIR})
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
    DEPENDS ${DOC_SOURCES} ${DOCENV_DEPS} ${DOXYGEN_XML_DIR}/index.xml ${BUILD_DOC_CONFIG_FILE}
    COMMAND ${Python3_EXECUTABLE} ${LAMMPS_DOC_DIR}/utils/make-globbed-tocs.py -d ${LAMMPS_DOC_DIR}/src
    COMMAND Sphinx::sphinx-build ${SPHINX_EXTRA_OPTS} -b html -c ${DOC_BUILD_DIR} -d ${DOC_BUILD_DIR}/doctrees ${LAMMPS_DOC_DIR}/src ${DOC_BUILD_DIR}/html
    COMMAND ${CMAKE_COMMAND} -E create_symlink Manual.html ${DOC_BUILD_DIR}/html/index.html
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${LAMMPS_DOC_DIR}/src/PDF ${DOC_BUILD_DIR}/html/PDF
    COMMAND ${CMAKE_COMMAND} -E remove -f ${DOXYGEN_XML_DIR}/run.stamp
  )

  add_custom_target(
    doc ALL
    DEPENDS html ${MATHJAX_DIR}/tex-mml-chtml.js ${MATHJAX_FONT_DIR}/package.json
    SOURCES ${LAMMPS_DOC_DIR}/utils/requirements.txt ${DOC_SOURCES}
  )

  install(DIRECTORY ${DOC_BUILD_DIR}/html DESTINATION ${CMAKE_INSTALL_DOCDIR})
endif()
