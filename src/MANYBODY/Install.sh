# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp fix_qeq_comb.cpp ..
  cp pair_adp.cpp ..
  cp pair_airebo.cpp ..
  cp pair_comb.cpp ..
  cp pair_eam.cpp ..
  cp pair_eam_alloy.cpp ..
  cp pair_eam_fs.cpp ..
  cp pair_eim.cpp ..
  cp pair_lcbop.cpp ..
  cp pair_rebo.cpp ..
  cp pair_sw.cpp ..
  cp pair_tersoff.cpp ..
  cp pair_tersoff_zbl.cpp ..

  cp fix_qeq_comb.h ..
  cp pair_adp.h ..
  cp pair_airebo.h ..
  cp pair_comb.h ..
  cp pair_eam.h ..
  cp pair_eam_alloy.h ..
  cp pair_eam_fs.h ..
  cp pair_eim.h ..
  cp pair_lcbop.h ..
  cp pair_rebo.h ..
  cp pair_sw.h ..
  cp pair_tersoff.h ..
  cp pair_tersoff_zbl.h ..

elif (test $1 = 0) then

  rm -f ../fix_qeq_comb.cpp
  rm -f ../pair_adp.cpp
  rm -f ../pair_airebo.cpp
  rm -f ../pair_comb.cpp
  rm -f ../pair_eam.cpp
  rm -f ../pair_eam_alloy.cpp
  rm -f ../pair_eam_fs.cpp
  rm -f ../pair_eim.cpp
  rm -f ../pair_lcbop.cpp
  rm -f ../pair_rebo.cpp
  rm -f ../pair_sw.cpp
  rm -f ../pair_tersoff.cpp
  rm -f ../pair_tersoff_zbl.cpp

  rm -f ../fix_qeq_comb.h
  rm -f ../pair_adp.h
  rm -f ../pair_airebo.h
  rm -f ../pair_comb.h
  rm -f ../pair_eam.h
  rm -f ../pair_eam_alloy.h
  rm -f ../pair_eam_fs.h
  rm -f ../pair_eim.h
  rm -f ../pair_lcbop.h
  rm -f ../pair_rebo.h
  rm -f ../pair_sw.h
  rm -f ../pair_tersoff.h
  rm -f ../pair_tersoff_zbl.h

fi
