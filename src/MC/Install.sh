# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp fix_bond_break.cpp ..
  cp fix_bond_create.cpp ..
  cp fix_bond_swap.cpp ..
  cp fix_gcmc.cpp ..
  cp pair_dsmc.cpp ..

  cp fix_bond_break.h ..
  cp fix_bond_create.h ..
  cp fix_bond_swap.h ..
  cp fix_gcmc.h ..
  cp pair_dsmc.h ..

elif (test $1 = 0) then

  rm -f ../fix_bond_break.cpp
  rm -f ../fix_bond_create.cpp
  rm -f ../fix_bond_swap.cpp
  rm -f ../fix_gcmc.cpp
  rm -f ../pair_dsmc.cpp

  rm -f ../fix_bond_break.h
  rm -f ../fix_bond_create.h
  rm -f ../fix_bond_swap.h
  rm -f ../fix_gcmc.h
  rm -f ../pair_dsmc.h

fi
