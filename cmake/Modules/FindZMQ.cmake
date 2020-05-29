find_path(ZMQ_INCLUDE_DIR zmq.h)
find_library(ZMQ_LIBRARY NAMES zmq)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ZMQ DEFAULT_MSG ZMQ_LIBRARY ZMQ_INCLUDE_DIR)

# Copy the results to the output variables and target.
if(ZMQ_FOUND)
  set(ZMQ_LIBRARIES ${ZMQ_LIBRARY})
  set(ZMQ_INCLUDE_DIRS ${ZMQ_INCLUDE_DIR})

  if(NOT TARGET ZMQ::ZMQ)
    add_library(ZMQ::ZMQ UNKNOWN IMPORTED)
    set_target_properties(ZMQ::ZMQ PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      IMPORTED_LOCATION "${ZMQ_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${ZMQ_INCLUDE_DIR}")
  endif()
endif()
