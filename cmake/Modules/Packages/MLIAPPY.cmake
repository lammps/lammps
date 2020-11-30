execute_process(COMMAND cythonize mliap_model_python_couple.pyx WORKING_DIRECTORY ../src/MLIAPPY)

target_compile_definitions(lammps PRIVATE -DLMP_MLIAPPY)
