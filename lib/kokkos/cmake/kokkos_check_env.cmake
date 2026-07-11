set(CRAYPE_VERSION $ENV{CRAYPE_VERSION})
if(CRAYPE_VERSION)
  set(KOKKOS_IS_CRAYPE TRUE)
  set(CRAYPE_LINK_TYPE $ENV{CRAYPE_LINK_TYPE})
  if(CRAYPE_LINK_TYPE)
    if(NOT CRAYPE_LINK_TYPE STREQUAL "dynamic")
      message(
        WARNING
          "CRAYPE_LINK_TYPE is set to ${CRAYPE_LINK_TYPE}. Linking is likely to fail unless this is set to 'dynamic'"
      )
    endif()
  else()
    message(WARNING "CRAYPE_LINK_TYPE is not set. Linking is likely to fail unless this is set to 'dynamic'")
  endif()
endif()
