# preset that turns on a wide range of packages, some of which require
# external libraries. Compared to all_on.cmake some more unusual packages
# are removed. The resulting binary should be able to run most inputs.

set(ALL_PACKAGES ASPHERE CLASS2 COLLOID CORESHELL DIPOLE
        GRANULAR KSPACE MANYBODY MC MISC MOLECULE OPT PERI
        PYTHON QEQ REPLICA RIGID SHOCK SNAP SRD VORONOI
        USER-CGDNA USER-CGSDK USER-COLVARS USER-DIFFRACTION USER-DPD
        USER-DRUDE USER-FEP USER-MEAMC USER-MESO
        USER-MISC USER-MOFFF USER-OMP USER-PHONON USER-REAXC
        USER-SPH USER-SMD USER-UEF USER-YAFF)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()
