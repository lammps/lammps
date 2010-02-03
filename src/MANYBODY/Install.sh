# Install/unInstall package files in LAMMPS
# for unInstall, also unInstall/Install OPT package if installed
#   so it will remove OPT files that depend on MANYBODY files,
#   then replace others

if (test $1 = 1) then

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

elif (test $1 = 0) then

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

  if (test -e ../pair_eam_opt.h) then
    cd ../OPT; sh Install.sh 0; sh Install.sh 1
  fi

fi
