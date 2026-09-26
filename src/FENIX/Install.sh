# enforce using portable C locale
LC_ALL=C
export LC_ALL

echo "The FENIX package does not support the legacy build system. Please build LAMMPS with CMake instead."
exit 1
