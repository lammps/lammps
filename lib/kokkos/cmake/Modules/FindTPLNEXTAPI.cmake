# Nextapi shared object and header both reside under NEXT_HOME directory
# Library under ${NEXT_HOME}/lib
# Headers under ${NEXT_HOME}/include

# FIXME_NEXTSILICON: See NextSilicon ticket CS-682

if(DEFINED ENV{NEXT_HOME})
  set(NEXT_HOME_PATH $ENV{NEXT_HOME})
else()
  # Give a default path in case env is not define.
  set(NEXT_HOME_PATH /opt/nextsilicon)
endif()

kokkos_find_imported(
  NEXTAPI
  LIBRARY
  nextapi
  LIBRARY_PATHS
  ${NEXT_HOME_PATH}
  HEADER_PATHS
  ${NEXT_HOME_PATH}
)
