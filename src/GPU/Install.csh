# Install/unInstall package classes in LAMMPS
# edit Makefile.package to include/exclude GPU library

if ($1 == 1) then

#  sed -i 's/\S*gpu //' ../Makefile.package
#  sed -i 's|^PKGPATH =\s*|&-L../../lib/gpu |' ../Makefile.package
#  sed -i 's|^PKGLIB =\s*|&-lpair_gpu_lib |' ../Makefile.package

  cp -f pair_lj_cut_gpu.h ../
  cp -f pair_lj_cut_gpu.cpp ../
  
  if (-e ../pair_gayberne.cpp) then
    cp -f style_gpu.h ../
    cp -f pair_gayberne_gpu.h ../
    cp -f pair_gayberne_gpu.cpp ../
  else
    grep -v erne style_gpu.h > ../style_gpu.h
  endif 
  
else if ($1 == 0) then

#  sed -i 's/\S*gpu //' ../Makefile.package

  rm -f ../style_gpu.h
  touch ../style_gpu.h

  rm -f ../pair_lj_cut_gpu.h
  rm -f ../pair_lj_cut_gpu.cpp

  rm -f ../pair_gayberne_gpu.h
  rm -f ../pair_gayberne_gpu.cpp

endif
