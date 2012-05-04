# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp molfile_interface.cpp ..
  cp dump_molfile.cpp ..

  cp molfile_interface.h ..
  cp dump_molfile.h ..

  cp molfile_plugin.h ..
  cp vmdplugin.h ..

elif (test $1 = 0) then

  rm -f ../molfile_interface.cpp
  rm -f ../dump_molfile.cpp

  rm -f ../molfile_interface.h
  rm -f ../dump_molfile.h

  rm -f ../molfile_plugin.h
  rm -f ../vmdplugin.h
fi
