# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp compute_erotate_asphere.cpp ..
  cp compute_temp_asphere.cpp ..
  cp fix_nh_asphere.cpp ..
  cp fix_nph_asphere.cpp ..
  cp fix_npt_asphere.cpp ..
  cp fix_nve_asphere.cpp ..
  cp fix_nve_asphere_noforce.cpp ..
  cp fix_nve_line.cpp ..
  cp fix_nve_tri.cpp ..
  cp fix_nvt_asphere.cpp ..
  cp pair_gayberne.cpp ..
  cp pair_line_lj.cpp ..
  cp pair_resquared.cpp ..
  cp pair_tri_lj.cpp ..

  cp compute_erotate_asphere.h ..
  cp compute_temp_asphere.h ..
  cp fix_nh_asphere.h ..
  cp fix_nph_asphere.h ..
  cp fix_npt_asphere.h ..
  cp fix_nve_asphere.h ..
  cp fix_nve_asphere_noforce.h ..
  cp fix_nve_line.h ..
  cp fix_nve_tri.h ..
  cp fix_nvt_asphere.h ..
  cp pair_gayberne.h ..
  cp pair_line_lj.h ..
  cp pair_resquared.h ..
  cp pair_tri_lj.h ..

elif (test $1 = 0) then

  rm -f ../compute_erotate_asphere.cpp
  rm -f ../compute_temp_asphere.cpp
  rm -f ../fix_nh_asphere.cpp
  rm -f ../fix_nph_asphere.cpp
  rm -f ../fix_npt_asphere.cpp
  rm -f ../fix_nve_asphere.cpp
  rm -f ../fix_nve_asphere_noforce.cpp
  rm -f ../fix_nve_line.cpp
  rm -f ../fix_nve_tri.cpp
  rm -f ../fix_nvt_asphere.cpp
  rm -f ../pair_gayberne.cpp
  rm -f ../pair_line_lj.cpp
  rm -f ../pair_resquared.cpp
  rm -f ../pair_tri_lj.cpp

  rm -f ../compute_erotate_asphere.h
  rm -f ../compute_temp_asphere.h
  rm -f ../fix_nh_asphere.h
  rm -f ../fix_nph_asphere.h
  rm -f ../fix_npt_asphere.h
  rm -f ../fix_nve_asphere.h
  rm -f ../fix_nve_asphere_noforce.h
  rm -f ../fix_nve_line.h
  rm -f ../fix_nve_tri.h
  rm -f ../fix_nvt_asphere.h
  rm -f ../pair_gayberne.h
  rm -f ../pair_line_lj.h
  rm -f ../pair_resquared.h
  rm -f ../pair_tri_lj.h

fi
