#! /bin/bash

REL_TOL_NVE=1e-8
REL_TOL_NVT=5e-3
UNITS=lj

LMPDIR=/Users/xwb17127/Work/code/lammps
SRCDIR=$LMPDIR/src
BUILDDIR=$LMPDIR/build

DAY=$(grep -e 'LAMMPS_VERSION' $SRCDIR/version.h | awk '{print $3}' | sed 's/"//g') 
MONTH=$(grep -e 'LAMMPS_VERSION' $SRCDIR/version.h | awk '{print $4}' | sed 's/"//g') 
YEAR=$(grep -e 'LAMMPS_VERSION' $SRCDIR/version.h | awk '{print $5}' | sed 's/"//g') 
DATE=$DAY$MONTH$YEAR

if [ $# -eq 1 ] && [ $1 = run ]; then

  if [ $UNITS = lj ]; then 
    EXDIR=$LMPDIR/examples/PACKAGES/cgdna/examples/lj_units
    echo '# Running tests with lj units' | tee -a $EXDIR/test.log

  elif [ $UNITS = real ]; then 
    EXDIR=$LMPDIR/examples/PACKAGES/cgdna/examples/real_units
    echo '# Running tests with real units' | tee -a $EXDIR/test.log

  else
    echo '# Choose lj or real units'
    exit 1

  fi

  echo '# Compiling executable in' $BUILDDIR | tee -a $EXDIR/test.log

  cd $LMPDIR
  rm -rf $BUILDDIR
  mkdir $BUILDDIR
  cd $BUILDDIR
  cmake ../cmake/ -DBUILD_MPI=yes -DLAMMPS_MACHINE=mpi -DPKG_ASPHERE=on -DPKG_MOLECULE=on -DPKG_CG-DNA=on -DCMAKE_CXX_STANDARD=23  | tee -a $EXDIR/test.log
  cmake --build . --target clean | tee -a $EXDIR/test.log
  cmake --build . -j14 | tee -a $EXDIR/test.log

  ######################################################
  printf '\n# Running oxDNA duplex1 NVE test\n' | tee -a $EXDIR/test.log
  cd $EXDIR/oxDNA/duplex1
  mkdir test
  cd test
  cp $BUILDDIR/lmp_mpi .
  cp ../in.duplex1 .
  cp ../data.duplex1 .

  ### 1 MPI-task ###
  mpirun -np 1 ./lmp_mpi -in in.duplex1 > /dev/null
  mv log.lammps log.$DATE.duplex1.g++.1
  grep -e '[0-9]  ekin' log.$DATE.duplex1.g++.1 > e_test.1.dat
  grep -e '[0-9]  ekin' ../log*1 > e_ref.1.dat

  paste e_ref.1.dat e_test.1.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 1 MPI-task passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ### 4 MPI-tasks ###
  mpirun -np 4 ./lmp_mpi -in in.duplex1 > /dev/null
  mv log.lammps log.$DATE.duplex1.g++.4
  grep -e '[0-9]  ekin' log.$DATE.duplex1.g++.4 > e_test.4.dat
  grep -e '[0-9]  ekin' ../log*4 > e_ref.4.dat

  paste e_ref.4.dat e_test.4.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 4 MPI-tasks passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ######################################################
  printf '\n# Running oxDNA duplex2 NVE test\n' | tee -a $EXDIR/test.log
  cd $EXDIR/oxDNA/duplex2
  mkdir test
  cd test
  cp $BUILDDIR/lmp_mpi .
  cp ../in.duplex2 .
  cp ../data.duplex2 .

  ### 1 MPI-task ###
  mpirun -np 1 ./lmp_mpi -in in.duplex2 > /dev/null
  mv log.lammps log.$DATE.duplex2.g++.1
  grep -e '[0-9]  ekin' log.$DATE.duplex2.g++.1 > e_test.1.dat
  grep -e '[0-9]  ekin' ../log*1 > e_ref.1.dat

  paste e_ref.1.dat e_test.1.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 1 MPI-task passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ### 4 MPI-tasks ###
  mpirun -np 4 ./lmp_mpi -in in.duplex2 > /dev/null
  mv log.lammps log.$DATE.duplex2.g++.4
  grep -e '[0-9]  ekin' log.$DATE.duplex2.g++.4 > e_test.4.dat
  grep -e '[0-9]  ekin' ../log*4 > e_ref.4.dat

  paste e_ref.4.dat e_test.4.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 4 MPI-tasks passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ######################################################
  printf '\n# Running oxDNA potential file NVE test\n' | tee -a $EXDIR/test.log
  cd $EXDIR/oxDNA/potential_file
  mkdir test
  cd test
  cp $BUILDDIR/lmp_mpi .
  cp ../in.duplex1 .
  cp ../data.duplex1 .
  if [ $UNITS = lj ]; then
    cp ../oxdna_lj.cgdna .
  elif [ $UNITS = real ]; then
    cp ../oxdna_real.cgdna .
  fi

  ### 1 MPI-task ###
  mpirun -np 1 ./lmp_mpi -in in.duplex1 > /dev/null
  mv log.lammps log.$DATE.duplex1.g++.1
  grep -e '[0-9]  ekin' log.$DATE.duplex1.g++.1 > e_test.1.dat
  grep -e '[0-9]  ekin' ../log*1 > e_ref.1.dat

  paste e_ref.1.dat e_test.1.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 1 MPI-task passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ### 4 MPI-tasks ###
  mpirun -np 4 ./lmp_mpi -in in.duplex1 > /dev/null
  mv log.lammps log.$DATE.duplex1.g++.4
  grep -e '[0-9]  ekin' log.$DATE.duplex1.g++.4 > e_test.4.dat
  grep -e '[0-9]  ekin' ../log*4 > e_ref.4.dat

  paste e_ref.4.dat e_test.4.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 4 MPI-tasks passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ######################################################
  printf '\n# Running oxDNA2 duplex1 NVE test\n' | tee -a $EXDIR/test.log
  cd $EXDIR/oxDNA2/duplex1
  mkdir test
  cd test
  cp $BUILDDIR/lmp_mpi .
  cp ../in.duplex1 .
  cp ../data.duplex1 .

  ### 1 MPI-task ###
  mpirun -np 1 ./lmp_mpi -in in.duplex1 > /dev/null
  mv log.lammps log.$DATE.duplex1.g++.1
  grep -e '[0-9]  ekin' log.$DATE.duplex1.g++.1 > e_test.1.dat
  grep -e '[0-9]  ekin' ../log*1 > e_ref.1.dat

  paste e_ref.1.dat e_test.1.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 1 MPI-task passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ### 4 MPI-tasks ###
  mpirun -np 4 ./lmp_mpi -in in.duplex1 > /dev/null
  mv log.lammps log.$DATE.duplex1.g++.4
  grep -e '[0-9]  ekin' log.$DATE.duplex1.g++.4 > e_test.4.dat
  grep -e '[0-9]  ekin' ../log*4 > e_ref.4.dat

  paste e_ref.4.dat e_test.4.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 4 MPI-tasks passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ######################################################
  printf '\n# Running oxDNA2 duplex2 NVE test\n' | tee -a $EXDIR/test.log
  cd $EXDIR/oxDNA2/duplex2
  mkdir test
  cd test
  cp $BUILDDIR/lmp_mpi .
  cp ../in.duplex2 .
  cp ../data.duplex2 .

  ### 1 MPI-task ###
  mpirun -np 1 ./lmp_mpi -in in.duplex2 > /dev/null
  mv log.lammps log.$DATE.duplex2.g++.1
  grep -e '[0-9]  ekin' log.$DATE.duplex2.g++.1 > e_test.1.dat
  grep -e '[0-9]  ekin' ../log*1 > e_ref.1.dat

  paste e_ref.1.dat e_test.1.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 1 MPI-task passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ### 4 MPI-tasks ###
  mpirun -np 4 ./lmp_mpi -in in.duplex2 > /dev/null
  mv log.lammps log.$DATE.duplex2.g++.4
  grep -e '[0-9]  ekin' log.$DATE.duplex2.g++.4 > e_test.4.dat
  grep -e '[0-9]  ekin' ../log*4 > e_ref.4.dat

  paste e_ref.4.dat e_test.4.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 4 MPI-tasks passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ######################################################
  printf '\n# Running oxDNA2 duplex3 NVE test\n' | tee -a $EXDIR/test.log
  cd $EXDIR/oxDNA2/duplex3
  mkdir test
  cd test
  cp $BUILDDIR/lmp_mpi .
  cp ../in.duplex3 .
  cp ../data.duplex3 .

  ### 1 MPI-task ###
  mpirun -np 1 ./lmp_mpi -in in.duplex3 > /dev/null
  mv log.lammps log.$DATE.duplex3.g++.1
  grep -e '[0-9]  ekin' log.$DATE.duplex3.g++.1 > e_test.1.dat
  grep -e '[0-9]  ekin' ../log*1 > e_ref.1.dat

  paste e_ref.1.dat e_test.1.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 1 MPI-task passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ### 4 MPI-tasks ###
  mpirun -np 4 ./lmp_mpi -in in.duplex3 > /dev/null
  mv log.lammps log.$DATE.duplex3.g++.4
  grep -e '[0-9]  ekin' log.$DATE.duplex3.g++.4 > e_test.4.dat
  grep -e '[0-9]  ekin' ../log*4 > e_ref.4.dat

  paste e_ref.4.dat e_test.4.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 4 MPI-tasks passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ######################################################
  printf '\n# Running oxDNA2 dsring NVE test\n' | tee -a $EXDIR/test.log
  cd $EXDIR/oxDNA2/dsring
  mkdir test
  cd test
  cp $BUILDDIR/lmp_mpi .
  cp ../in.dsring .
  cp ../data.dsring .

  ### 1 MPI-task ###
  mpirun -np 1 ./lmp_mpi -in in.dsring > /dev/null
  mv log.lammps log.$DATE.dsring.g++.1
  grep -e '[0-9]  ekin' log.$DATE.dsring.g++.1 > e_test.1.dat
  grep -e '[0-9]  ekin' ../log*1 > e_ref.1.dat

  paste e_ref.1.dat e_test.1.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 1 MPI-task passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ### 4 MPI-tasks ###
  mpirun -np 4 ./lmp_mpi -in in.dsring > /dev/null
  mv log.lammps log.$DATE.dsring.g++.4
  grep -e '[0-9]  ekin' log.$DATE.dsring.g++.4 > e_test.4.dat
  grep -e '[0-9]  ekin' ../log*4 > e_ref.4.dat

  paste e_ref.4.dat e_test.4.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 4 MPI-tasks passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ######################################################
  printf '\n# Running oxDNA2 potential file NVE test\n' | tee -a $EXDIR/test.log
  cd $EXDIR/oxDNA2/potential_file
  mkdir test
  cd test
  cp $BUILDDIR/lmp_mpi .
  cp ../in.duplex1 .
  cp ../data.duplex1 .
  if [ $UNITS = lj ]; then
    cp ../oxdna2_lj.cgdna .
  elif [ $UNITS = real ]; then
    cp ../oxdna2_real.cgdna .
  fi

  ### 1 MPI-task ###
  mpirun -np 1 ./lmp_mpi -in in.duplex1 > /dev/null
  mv log.lammps log.$DATE.duplex1.g++.1
  grep -e '[0-9]  ekin' log.$DATE.duplex1.g++.1 > e_test.1.dat
  grep -e '[0-9]  ekin' ../log*1 > e_ref.1.dat

  paste e_ref.1.dat e_test.1.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 1 MPI-task passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ### 4 MPI-tasks ###
  mpirun -np 4 ./lmp_mpi -in in.duplex1 > /dev/null
  mv log.lammps log.$DATE.duplex1.g++.4
  grep -e '[0-9]  ekin' log.$DATE.duplex1.g++.4 > e_test.4.dat
  grep -e '[0-9]  ekin' ../log*4 > e_ref.4.dat

  paste e_ref.4.dat e_test.4.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 4 MPI-tasks passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ######################################################
  printf '\n# Running oxDNA3 duplex2 / potential file NVE test\n' | tee -a $EXDIR/test.log
  cd $EXDIR/oxDNA3/duplex2
  mkdir test
  cd test
  cp $BUILDDIR/lmp_mpi .
  cp ../in.duplex2 .
  cp ../data.duplex2 .
  if [ $UNITS = lj ]; then
    cp ../oxdna3_lj.cgdna .
  elif [ $UNITS = real ]; then
    cp ../oxdna3_real.cgdna .
  fi

  ### 1 MPI-task ###
  mpirun -np 1 ./lmp_mpi -in in.duplex2 > /dev/null
  mv log.lammps log.$DATE.duplex2.g++.1
  grep -e '[0-9]  ekin' log.$DATE.duplex2.g++.1 > e_test.1.dat
  grep -e '[0-9]  ekin' ../log*1 > e_ref.1.dat

  paste e_ref.1.dat e_test.1.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = 0
      if($4!=0) diff = ($4-$52)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $52, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($8!=0) diff = ($8-$56)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $56, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($12!=0) diff = ($12-$60)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $60, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($16!=0) diff = ($16-$64)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $64, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($20!=0) diff = ($20-$68)/$20
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $20, $68, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($24!=0) diff = ($24-$72)/$24
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $24, $72, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($28!=0) diff = ($28-$76)/$28
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $28, $76, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($32!=0) diff = ($32-$80)/$32
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $32, $80, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($36!=0) diff = ($36-$84)/$36
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $36, $84, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($40!=0) diff = ($40-$88)/$40
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $40, $88, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($44!=0) diff = ($44-$92)/$44
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $44, $92, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($48!=0) diff = ($48-$96)/$48
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $48, $96, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }

    }
    END {
      if (failed == 0) print "# 1 MPI-task passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ### 4 MPI-tasks ###
  mpirun -np 4 ./lmp_mpi -in in.duplex2 > /dev/null
  mv log.lammps log.$DATE.duplex2.g++.4
  grep -e '[0-9]  ekin' log.$DATE.duplex2.g++.4 > e_test.4.dat
  grep -e '[0-9]  ekin' ../log*4 > e_ref.4.dat

  paste e_ref.4.dat e_test.4.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = 0
      if($4!=0) diff = ($4-$52)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $52, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($8!=0) diff = ($8-$56)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $56, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($12!=0) diff = ($12-$60)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $60, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($16!=0) diff = ($16-$64)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $64, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($20!=0) diff = ($20-$68)/$20
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $20, $68, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($24!=0) diff = ($24-$72)/$24
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $24, $72, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($28!=0) diff = ($28-$76)/$28
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $28, $76, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($32!=0) diff = ($32-$80)/$32
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $32, $80, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($36!=0) diff = ($36-$84)/$36
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $36, $84, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($40!=0) diff = ($40-$88)/$40
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $40, $88, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($44!=0) diff = ($44-$92)/$44
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $44, $92, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = 0
      if($48!=0) diff = ($48-$96)/$48
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $48, $96, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }

    }
    END {
      if (failed == 0) print "# 4 MPI-tasks passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log
  ######################################################
  printf '\n# Running oxDNA3 NVT and unique base pairing test\n' | tee -a $EXDIR/test.log
  cd $EXDIR/oxDNA3/unique_bp
  mkdir test
  cd test
  cp $BUILDDIR/lmp_mpi .
  cp ../in.dsring2 .
  cp ../data.dsring2 .

  ### 8 MPI-tasks ###
  mpirun -np 8 ./lmp_mpi -in in.dsring2 > /dev/null
  mv log.lammps log.$DATE.dsring2.g++.8
  grep -e '[0-9]  ekin' log.$DATE.dsring2.g++.8 | awk '{print $1, $12}' > edyn_test.8.dat
  grep -e '[0-9]  ekin' log.$DATE.dsring2.g++.8 | awk '{print $1, $28}' > ehbond_test.8.dat
  grep -e '[0-9]  ekin' ../log*dsring2.g++.8 | awk '{print $1, $12}' > edyn_ref.8.dat
  grep -e '[0-9]  ekin' ../log*dsring2.g++.8 | awk '{print $1, $28}' > ehbond_ref.8.dat

  avg_edyn_test=$(awk '{sum += $2; n++} END {if (n > 0) print sum / n}' edyn_test.8.dat)
  avg_edyn_ref=$(awk '{sum += $2; n++} END {if (n > 0) print sum / n}' edyn_ref.8.dat)
  avg_ehbond_test=$(awk 'NR > 2000 {sum += $2; n++} END {if (n > 0) print sum / n}' ehbond_test.8.dat) 
  avg_ehbond_ref=$(awk 'NR > 2000 {sum += $2; n++} END {if (n > 0) print sum / n}' ehbond_ref.8.dat)

  tol=$REL_TOL_NVT

  if [ $UNITS = lj ]; then
    ekin=44.4
  fi  

  if [ $UNITS = real ]; then
    ekin=264.6956072509588
  fi

  diff=$(echo "($avg_edyn_test - $ekin)/$ekin" | bc -l)
  diff=$(echo "if ($diff < 0) -1 * $diff else $diff" | bc -l)

  if (( $(echo "$diff > $REL_TOL_NVT" | bc -l) )); then
    printf "# Relative difference of kinetic energy %g > %g\n" "$diff" "$tol" | tee -a $EXDIR/test.log 
    echo "# 8 MPI-tasks NVT FAILED" | tee -a $EXDIR/test.log
  else
    echo "# 8 MPI-tasks NVT passed" | tee -a $EXDIR/test.log
  fi

  diff=$(echo "($avg_ehbond_test - $avg_ehbond_ref)/$avg_ehbond_ref" | bc -l)
  diff=$(echo "if ($diff < 0) -1 * $diff else $diff" | bc -l)

  if (( $(echo "$diff > $REL_TOL_NVT" | bc -l) )); then
    printf "# Relative difference of hydrogen bonding energy %g > %g\n" "$diff" "$tol" | tee -a $EXDIR/test.log
    echo "# 8 MPI-tasks unique base pairing FAILED" | tee -a $EXDIR/test.log
  else
    echo "# 8 MPI-tasks unique base pairing passed" | tee -a $EXDIR/test.log
  fi

  ######################################################
  printf '\n# Running oxRNA2 duplex2 NVE test\n' | tee -a $EXDIR/test.log
  cd $EXDIR/oxRNA2/duplex2
  mkdir test
  cd test
  cp $BUILDDIR/lmp_mpi .
  cp ../in.duplex2 .
  cp ../data.duplex2 .

  ### 1 MPI-task ###
  mpirun -np 1 ./lmp_mpi -in in.duplex2 > /dev/null
  mv log.lammps log.$DATE.duplex2.g++.1
  grep -e '[0-9]  ekin' log.$DATE.duplex2.g++.1 > e_test.1.dat
  grep -e '[0-9]  ekin' ../log*1 > e_ref.1.dat

  paste e_ref.1.dat e_test.1.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 1 MPI-task passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ### 4 MPI-tasks ###
  mpirun -np 4 ./lmp_mpi -in in.duplex2 > /dev/null
  mv log.lammps log.$DATE.duplex2.g++.4
  grep -e '[0-9]  ekin' log.$DATE.duplex2.g++.4 > e_test.4.dat
  grep -e '[0-9]  ekin' ../log*4 > e_ref.4.dat

  paste e_ref.4.dat e_test.4.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 4 MPI-tasks passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ######################################################
  printf '\n# Running oxRNA2 potential file NVE test\n' | tee -a $EXDIR/test.log
  cd $EXDIR/oxRNA2/potential_file
  mkdir test
  cd test
  cp $BUILDDIR/lmp_mpi .
  cp ../in.duplex2 .
  cp ../data.duplex2 .
  if [ $UNITS = lj ]; then
    cp ../oxrna2_lj.cgdna .
  elif [ $UNITS = real ]; then
    cp ../oxrna2_real.cgdna .
  fi

  ### 1 MPI-task ###
  mpirun -np 1 ./lmp_mpi -in in.duplex2 > /dev/null
  mv log.lammps log.$DATE.duplex2.g++.1
  grep -e '[0-9]  ekin' log.$DATE.duplex2.g++.1 > e_test.1.dat
  grep -e '[0-9]  ekin' ../log*1 > e_ref.1.dat

  paste e_ref.1.dat e_test.1.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 1 MPI-task FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 1 MPI-task passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

  ### 4 MPI-tasks ###
  mpirun -np 4 ./lmp_mpi -in in.duplex2 > /dev/null
  mv log.lammps log.$DATE.duplex2.g++.4
  grep -e '[0-9]  ekin' log.$DATE.duplex2.g++.4 > e_test.4.dat
  grep -e '[0-9]  ekin' ../log*4 > e_ref.4.dat

  paste e_ref.4.dat e_test.4.dat |

  awk -v tol="$REL_TOL_NVE" '
    failed == 0 {
      diff = ($4-$20)/$4
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $4, $20, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($8-$24)/$8
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $8, $24, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($12-$28)/$12
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $12, $28, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
      diff = ($16-$32)/$16
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $16, $32, diff, tol
        printf "# 4 MPI-tasks FAILED\n"
        failed = 1
        exit 1
      }
    }
    END {
      if (failed == 0) print "# 4 MPI-tasks passed"
    }
  ' 2>&1 | tee -a $EXDIR/test.log

 ######################################################

  echo | tee -a $EXDIR/test.log
  echo '# Done' | tee -a $EXDIR/test.log

elif [ $# -eq 1 ] && [ $1 = clean ]; then

  if [ $UNITS = lj ]; then 
    EXDIR=$LMPDIR/examples/PACKAGES/cgdna/examples/lj_units

  elif [ $UNITS = real ]; then 
    EXDIR=$LMPDIR/examples/PACKAGES/cgdna/examples/real_units

  else
    echo '# Choose lj or real units'
    exit 1

  fi

  echo '# Deleting test directories'
  rm -rf $EXDIR/oxDNA/duplex1/test
  rm -rf $EXDIR/oxDNA/duplex2/test
  rm -rf $EXDIR/oxDNA/potential_file/test
  rm -rf $EXDIR/oxDNA2/duplex1/test
  rm -rf $EXDIR/oxDNA2/duplex2/test
  rm -rf $EXDIR/oxDNA2/duplex3/test
  rm -rf $EXDIR/oxDNA2/unique_bp/test
  rm -rf $EXDIR/oxDNA2/dsring/test
  rm -rf $EXDIR/oxDNA2/potential_file/test
  rm -rf $EXDIR/oxDNA3/duplex2/test
  rm -rf $EXDIR/oxDNA3/unique_bp/test
  rm -rf $EXDIR/oxRNA2/duplex2/test
  rm -rf $EXDIR/oxRNA2/potential_file/test
  rm -rf $EXDIR/test.log

  echo '# Deleting build directory'
  rm -rf $LMPDIR/build

  echo '# Done'
  
else 
  echo '# Usage:'
  echo '# ./test.sh run  ... to run test suite'
  echo '# ./test.sh clean ... to delete test directories'
  
fi
