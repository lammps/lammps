# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp style_poems.h ..

  cp fix_poems.cpp ..

  cp fix_poems.h ..

else if ($1 == 0) then

  rm ../style_poems.h
  touch ../style_poems.h

  rm ../fix_poems.cpp

  rm ../fix_poems.h

endif
