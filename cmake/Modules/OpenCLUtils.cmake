function(WriteOpenCLHeader varname outfile files)
    file(WRITE ${outfile} "const char * ${varname} = \n")
    separate_arguments(files)

    foreach(filename ${files})
        file(READ ${filename} content)
        string(REGEX REPLACE "\\s*//[^\n]*\n" "\n" content "${content}")
        string(REGEX REPLACE "\\\\" "\\\\\\\\" content "${content}")
        string(REGEX REPLACE "\"" "\\\\\"" content "${content}")
        string(REGEX REPLACE "([^\n]+)\n" "\"\\1\\\\n\"\n" content "${content}")
        string(REGEX REPLACE "\n+" "\n" content "${content}")
        file(APPEND ${outfile} "${content}")
    endforeach()

    file(APPEND ${outfile} ";\n")
endfunction(WriteOpenCLHeader)

function(GenerateOpenCLHeader varname outfile files)
    list(REMOVE_AT ARGV 0 1)
    add_custom_command(OUTPUT ${outfile}
      COMMAND ${CMAKE_COMMAND} -D SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
                               -D VARNAME=${varname}
                               -D HEADER_FILE=${outfile}
                               -D SOURCE_FILES="${ARGV}"
                               -P ${CMAKE_CURRENT_SOURCE_DIR}/Modules/GenerateOpenCLHeader.cmake
      DEPENDS ${ARGV}
      COMMENT "Generating ${outfile}...")
endfunction(GenerateOpenCLHeader)
