# Install/unInstall package classes in LAMMPS

# fix_shear_history.h must always be in src

if ($1 == 1) then

  cp style_granular.h ..

  cp atom_granular.cpp ..
  cp fix_freeze.cpp ..
  cp fix_gran_diag.cpp ..
  cp fix_nve_gran.cpp ..
  cp fix_pour.cpp ..
  cp fix_shear_history.cpp ..
  cp fix_wall_gran.cpp ..
  cp pair_gran_hertzian.cpp ..
  cp pair_gran_history.cpp ..
  cp pair_gran_no_history.cpp ..

  cp atom_granular.h ..
  cp fix_freeze.h ..
  cp fix_gran_diag.h ..
  cp fix_nve_gran.h ..
  cp fix_pour.h ..
#  cp fix_shear_history.h ..
  cp fix_wall_gran.h ..
  cp pair_gran_hertzian.h ..
  cp pair_gran_history.h ..
  cp pair_gran_no_history.h ..

else if ($1 == 0) then

  rm ../style_granular.h
  touch ../style_granular.h

  rm ../atom_granular.cpp
  rm ../fix_freeze.cpp
  rm ../fix_gran_diag.cpp
  rm ../fix_nve_gran.cpp
  rm ../fix_pour.cpp
  rm ../fix_shear_history.cpp
  rm ../fix_wall_gran.cpp
  rm ../pair_gran_hertzian.cpp
  rm ../pair_gran_history.cpp
  rm ../pair_gran_no_history.cpp

  rm ../atom_granular.h
  rm ../fix_freeze.h
  rm ../fix_gran_diag.h
  rm ../fix_nve_gran.h
  rm ../fix_pour.h
#  rm ../fix_shear_history.h
  rm ../fix_wall_gran.h
  rm ../pair_gran_hertzian.h
  rm ../pair_gran_history.h
  rm ../pair_gran_no_history.h

endif
