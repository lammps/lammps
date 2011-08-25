# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp fix_append_atoms.cpp ..
  cp fix_msst.cpp ..
  cp fix_nphug.cpp ..
  cp fix_wall_piston.cpp ..

  cp fix_append_atoms.h ..
  cp fix_msst.h ..
  cp fix_nphug.h ..
  cp fix_wall_piston.h ..

elif (test $1 = 0) then

  rm -f ../fix_append_atoms.cpp
  rm -f ../fix_msst.cpp
  rm -f ../fix_nphug.cpp
  rm -f ../fix_wall_piston.cpp

  rm -f ../fix_append_atoms.h
  rm -f ../fix_msst.h
  rm -f ../fix_nphug.h
  rm -f ../fix_wall_piston.h

fi
