# utility script to call WriteOpenCLHeader function
include(${SOURCE_DIR}/Modules/OpenCLUtils.cmake)
WriteOpenCLHeader(${VARNAME} ${HEADER_FILE} ${SOURCE_FILES})
