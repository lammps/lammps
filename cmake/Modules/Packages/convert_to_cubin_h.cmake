foreach(CU_OBJ ${OBJS})
    get_filename_component(CU_NAME ${CU_OBJ} NAME_WE)
    string(REGEX REPLACE "^.*lal_" "" CU_NAME "${CU_NAME}")
    message(STATUS "Generating ${CU_NAME}_cubin.h")
    execute_process(
        COMMAND ${BIN2C} -c -n ${CU_NAME} ${CU_OBJ}
	OUTPUT_FILE ${OUTPUT_DIR}/${CU_NAME}_cubin.h
    )
endforeach()
