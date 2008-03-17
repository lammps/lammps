# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp style_granular.h ..

  cp atom_vec_granular.cpp ..
  cp fix_freeze.cpp ..
  cp fix_pour.cpp ..
  cp fix_wall_gran.cpp ..
  cp pair_gran_hertzian.cpp ..
  cp pair_gran_history.cpp ..
  cp pair_gran_no_history.cpp ..

  cp atom_vec_granular.h ..
  cp fix_freeze.h ..
  cp fix_pour.h ..
  cp fix_wall_gran.h ..
  cp pair_gran_hertzian.h ..
  cp pair_gran_history.h ..
  cp pair_gran_no_history.h ..

else if ($1 == 0) then

  rm ../style_granular.h
  touch ../style_granular.h

  rm ../atom_vec_granular.cpp
  rm ../fix_freeze.cpp
  rm ../fix_pour.cpp
  rm ../fix_wall_gran.cpp
  rm ../pair_gran_hertzian.cpp
  rm ../pair_gran_history.cpp
  rm ../pair_gran_no_history.cpp

  rm ../atom_vec_granular.h
  rm ../fix_freeze.h
  rm ../fix_pour.h
  rm ../fix_wall_gran.h
  rm ../pair_gran_hertzian.h
  rm ../pair_gran_history.h
  rm ../pair_gran_no_history.h

endif
