# preset that turns on a wide range of packages none of which require
# external libraries. Some more unusual packages are removed as well.
# The resulting binary should be able to run most inputs.

set(ALL_PACKAGES ASPHERE CLASS2 COLLOID CORESHELL DIPOLE
        GRANULAR KSPACE MANYBODY MC MISC MOLECULE OPT PERI
        QEQ REPLICA RIGID SHOCK SRD
        USER-CGDNA USER-CGSDK USER-COLVARS USER-DIFFRACTION USER-DPD
        USER-DRUDE USER-FEP USER-MEAMC USER-MESO
        USER-MISC USER-MOFFF USER-OMP USER-PHONON USER-REAXC
        USER-SPH USER-UEF USER-YAFF)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()
