# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp style_reax.h ..

  cp pair_reax.cpp ..

  cp pair_reax.h ..
  cp reax_cbkabo.h ..
  cp reax_cbkbo.h ..
  cp reax_cbkc.h ..
  cp reax_cbkch.h ..
  cp reax_cbkd.h ..
  cp reax_cbkia.h ..
  cp reax_cbklonpar.h ..
  cp reax_cbknubon2.h ..
  cp reax_cbkpairs.h ..
  cp reax_cbkqa.h ..
  cp reax_energies.h ..
  cp reax_fortran.h ..
  cp reax_functions.h ..
  cp reax_params.h ..
  cp reax_small.h ..

else if ($1 == 0) then

  rm ../style_reax.h
  touch ../style_reax.h

  rm ../pair_reax.cpp

  rm ../pair_reax.h
  rm ../reax_cbkabo.h
  rm ../reax_cbkbo.h
  rm ../reax_cbkc.h
  rm ../reax_cbkch.h
  rm ../reax_cbkd.h
  rm ../reax_cbkia.h
  rm ../reax_cbklonpar.h
  rm ../reax_cbknubon2.h
  rm ../reax_cbkpairs.h
  rm ../reax_cbkqa.h
  rm ../reax_energies.h
  rm ../reax_fortran.h
  rm ../reax_functions.h
  rm ../reax_params.h
  rm ../reax_small.h

endif
