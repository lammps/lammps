# utility script to call GenerateBinaryHeader function
include(${SOURCE_DIR}/Modules/LAMMPSUtils.cmake)
GenerateBinaryHeader(${VARNAME} ${HEADER_FILE} ${SOURCE_FILES})
