# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp style_granular.h ..

  cp atom_vec_granular.cpp ..
  cp fix_freeze.cpp ..
  cp fix_pour.cpp ..
  cp fix_wall_gran.cpp ..
  cp pair_gran_hertz_history.cpp ..
  cp pair_gran_hooke.cpp ..
  cp pair_gran_hooke_history.cpp ..

  cp atom_vec_granular.h ..
  cp fix_freeze.h ..
  cp fix_pour.h ..
  cp fix_wall_gran.h ..
  cp pair_gran_hertz_history.h ..
  cp pair_gran_hooke.h ..
  cp pair_gran_hooke_history.h ..

else if ($1 == 0) then

  rm ../style_granular.h
  touch ../style_granular.h

  rm ../atom_vec_granular.cpp
  rm ../fix_freeze.cpp
  rm ../fix_pour.cpp
  rm ../fix_wall_gran.cpp
  rm ../pair_gran_hertz_history.cpp
  rm ../pair_gran_hooke.cpp
  rm ../pair_gran_hooke_history.cpp

  rm ../atom_vec_granular.h
  rm ../fix_freeze.h
  rm ../fix_pour.h
  rm ../fix_wall_gran.h
  rm ../pair_gran_hertz_history.h
  rm ../pair_gran_hooke.h
  rm ../pair_gran_hooke_history.h

endif
