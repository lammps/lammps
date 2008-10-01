# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp style_manybody.h ..

  cp pair_airebo.cpp ..
  cp pair_eam.cpp ..
  cp pair_eam_alloy.cpp ..
  cp pair_eam_fs.cpp ..
  cp pair_sw.cpp ..
  cp pair_tersoff.cpp ..
  cp pair_tersoff_zbl.cpp ..

  cp pair_airebo.h ..
  cp pair_eam.h ..
  cp pair_eam_alloy.h ..
  cp pair_eam_fs.h ..
  cp pair_sw.h ..
  cp pair_tersoff.h ..
  cp pair_tersoff_zbl.h ..

else if ($1 == 0) then

  rm ../style_manybody.h
  touch ../style_manybody.h

  rm ../pair_airebo.cpp
  rm ../pair_eam.cpp
  rm ../pair_eam_alloy.cpp
  rm ../pair_eam_fs.cpp
  rm ../pair_sw.cpp
  rm ../pair_tersoff.cpp
  rm ../pair_tersoff_zbl.cpp

  rm ../pair_airebo.h
  rm ../pair_eam.h
  rm ../pair_eam_alloy.h
  rm ../pair_eam_fs.h
  rm ../pair_sw.h
  rm ../pair_tersoff.h
  rm ../pair_tersoff_zbl.h

endif
