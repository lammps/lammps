# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp style_asphere.h ..

  cp atom_vec_ellipsoid.cpp ..
  cp compute_erotate_asphere.cpp ..
  cp compute_temp_asphere.cpp ..
  cp fix_npt_asphere.cpp ..
  cp fix_nve_asphere.cpp ..
  cp fix_nvt_asphere.cpp ..
  cp pair_gayberne.cpp ..
  cp pair_resquared.cpp ..

  cp atom_vec_ellipsoid.h ..
  cp compute_erotate_asphere.h ..
  cp compute_temp_asphere.h ..
  cp fix_npt_asphere.h ..
  cp fix_nve_asphere.h ..
  cp fix_nvt_asphere.h ..
  cp pair_gayberne.h ..
  cp pair_resquared.h ..

else if ($1 == 0) then

  rm ../style_asphere.h
  touch ../style_asphere.h

  rm ../atom_vec_ellipsoid.cpp
  rm ../compute_erotate_asphere.cpp
  rm ../compute_temp_asphere.cpp
  rm ../fix_npt_asphere.cpp
  rm ../fix_nve_asphere.cpp
  rm ../fix_nvt_asphere.cpp
  rm ../pair_gayberne.cpp
  rm ../pair_resquared.cpp

  rm ../atom_vec_ellipsoid.h
  rm ../compute_erotate_asphere.h
  rm ../compute_temp_asphere.h
  rm ../fix_npt_asphere.h
  rm ../fix_nve_asphere.h
  rm ../fix_nvt_asphere.h
  rm ../pair_gayberne.h
  rm ../pair_resquared.h

endif
