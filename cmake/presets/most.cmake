# preset that turns on a wide range of packages, some of which require
# external libraries. Compared to all_on.cmake some more unusual packages
# are removed. The resulting binary should be able to run most inputs.

set(ALL_PACKAGES ASPHERE BODY CLASS2 COLLOID COMPRESS CORESHELL DIPOLE
        GRANULAR KSPACE MANYBODY MC MISC ML-IAP MOLECULE OPT PERI
        PLUGIN POEMS PYTHON QEQ REPLICA RIGID SHOCK ML-SNAP SPIN SRD VORONOI
        BROWNIAN BOCS CG-DNA CG-SDK COLVARS
        DIFFRACTION DPD-REACT DRUDE EFF FEP MEAM
        DPD-MESO USER-MISC MOFFF OPENMP PHONON REACTION
        REAXFF DPD-SMOOTH SPH MACHDYN UEF YAFF DIELECTRIC)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()

set(BUILD_TOOLS ON CACHE BOOL "" FORCE)
