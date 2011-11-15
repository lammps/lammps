# Install/unInstall package files in LAMMPS

if (test $1 == 1) then

  cp pair_brownian.cpp ..
  cp pair_brownian_poly.cpp ..
  cp pair_lubricate.cpp ..
  cp pair_lubricate_poly.cpp ..
  cp pair_lubricateU.cpp ..
  cp pair_lubricateU_poly.cpp ..

  cp pair_brownian.h ..
  cp pair_brownian_poly.h ..
  cp pair_lubricate.h ..
  cp pair_lubricate_poly.h ..
  cp pair_lubricateU.h ..
  cp pair_lubricateU_poly.h ..

elif (test $1 == 0) then

  rm ../pair_brownian.cpp
  rm ../pair_brownian_poly.cpp
  rm ../pair_lubricate.cpp
  rm ../pair_lubricate_poly.cpp
  rm ../pair_lubricateU.cpp
  rm ../pair_lubricateU_poly.cpp

  rm ../pair_brownian.h
  rm ../pair_brownian_poly.h
  rm ../pair_lubricate.h
  rm ../pair_lubricate_poly.h
  rm ../pair_lubricateU.h
  rm ../pair_lubricateU_poly.h

fi
