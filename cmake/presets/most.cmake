# preset that turns on a wide range of packages, some of which require
# external libraries. Compared to all_on.cmake some more unusual packages
# are removed. The resulting binary should be able to run most inputs.

set(ALL_PACKAGES
  AMOEBA
  ASPHERE
  BOCS
  BODY
  BPM
  BROWNIAN
  CG-DNA
  CG-SPICA
  CLASS2
  COLLOID
  COLVARS
  COMPRESS
  CORESHELL
  DIELECTRIC
  DIFFRACTION
  DIPOLE
  DPD-BASIC
  DPD-MESO
  DPD-REACT
  DPD-SMOOTH
  DRUDE
  ELECTRODE
  EFF
  EXTRA-COMPUTE
  EXTRA-DUMP
  EXTRA-FIX
  EXTRA-MOLECULE
  EXTRA-PAIR
  FEP
  GRANULAR
  INTERLAYER
  KSPACE
  MACHDYN
  MANYBODY
  MC
  MEAM
  MISC
  ML-IAP
  ML-SNAP
  MOFFF
  MOLECULE
  OPENMP
  OPT
  ORIENT
  PERI
  PHONON
  PLUGIN
  POEMS
  QEQ
  REACTION
  REAXFF
  REPLICA
  RIGID
  SHOCK
  SPH
  SPIN
  SRD
  TALLY
  UEF
  VORONOI
  YAFF)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()

set(BUILD_TOOLS ON CACHE BOOL "" FORCE)
