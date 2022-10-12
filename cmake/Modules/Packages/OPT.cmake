  set(OPT_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/OPT)
  set(OPT_SOURCES)
  set_property(GLOBAL PROPERTY "OPT_SOURCES" "${OPT_SOURCES}")

  # detects styles which have OPT version
  RegisterStylesExt(${OPT_SOURCES_DIR} opt OPT_SOURCES)

  get_property(OPT_SOURCES GLOBAL PROPERTY OPT_SOURCES)

  target_sources(lammps PRIVATE ${OPT_SOURCES})
  target_include_directories(lammps PRIVATE ${OPT_SOURCES_DIR})
