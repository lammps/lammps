# Install/unInstall package files in LAMMPS
# edit Makefile.package to include/exclude ATC library

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*awpmd[^ \t]* //g' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-I../../lib/awpmd/ivutils/include -I../../lib/awpmd/systems/interact |' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/awpmd |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-lawpmd |' ../Makefile.package
    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(user-awpmd_SYSPATH) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(user-awpmd_SYSLIB) |' ../Makefile.package
  fi

  cp atom_vec_wavepacket.cpp ..
  cp fix_nve_awpmd.cpp ..
  cp pair_awpmd_cut.cpp ..

  cp atom_vec_wavepacket.h ..
  cp fix_nve_awpmd.h ..
  cp pair_awpmd_cut.h ..

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*awpmd[^ \t]* //g' ../Makefile.package
  fi

  rm -f ../atom_vec_wavepacket.cpp
  rm -f ../fix_nve_awpmd.cpp
  rm -f ../pair_awpmd_cut.cpp

  rm -f ../atom_vec_wavepacket.h
  rm -f ../fix_nve_awpmd.h
  rm -f ../pair_awpmd_cut.h

fi
