execute_process(COMMAND cythonize mliap_model_python_couple.pyx WORKING_DIRECTORY ../src/MLIAP)

target_compile_definitions(lammps PRIVATE -DLMP_MLIAPPY)
