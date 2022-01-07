#! /bin/bash

DATE='14Dec21'
TOL=1e-8

LMPDIR=/Users/ohenrich/Work/code/lammps
SRCDIR=$LMPDIR/src
EXDIR=$LMPDIR/examples/PACKAGES/cgdna/examples

if [ $# -eq 1 ] && [ $1 = run ]; then
  echo '# Compiling executable in' $SRCDIR | tee -a $EXDIR/test.log

  cd $SRCDIR
  make clean-all | tee -a $EXDIR/test.log
  make purge | tee -a $EXDIR/test.log
  make pu | tee -a $EXDIR/test.log
  make ps | tee -a $EXDIR/test.log
  make -j8 mpi | tee -a $EXDIR/test.log

  ######################################################
  echo '# Running oxDNA duplex1 test' | tee -a $EXDIR/test.log
  cd $EXDIR/oxDNA/duplex1
  mkdir test
  cd test
  cp $SRCDIR/lmp_mpi .
  cp ../in.duplex1 .
  cp ../data.duplex1 .

  mpirun -np 1 ./lmp_mpi < in.duplex1 > /dev/null
  mv log.lammps log.$DATE.duplex1.g++.1
  grep etot log.$DATE.duplex1.g++.1 > e_test.1.dat
  grep etot ../log*1 > e_old.1.dat
  ndiff -relerr $TOL e_test.1.dat e_old.1.dat
  if [ $? -eq 0 ];
  then 
      echo "# 1 MPI-task passed" | tee -a $EXDIR/test.log
  else 
      echo "# 1 MPI-task unsuccessful" | tee -a $EXDIR/test.log
  fi

  mpirun -np 4 ./lmp_mpi < in.duplex1 > /dev/null
  mv log.lammps log.$DATE.duplex1.g++.4
  grep etot log.$DATE.duplex1.g++.4 > e_test.4.dat
  grep etot ../log*4 > e_old.4.dat
  ndiff -relerr $TOL e_test.4.dat e_old.4.dat
  if [ $? -eq 0 ];
  then 
      echo "# 4 MPI-tasks passed" | tee -a $EXDIR/test.log
  else 
      echo "# 4 MPI-tasks unsuccessful" | tee -a $EXDIR/test.log
  fi

  ######################################################

  ######################################################
  echo '# Running oxDNA duplex2 test' | tee -a $EXDIR/test.log
  cd $EXDIR/oxDNA/duplex2
  mkdir test
  cd test
  cp $SRCDIR/lmp_mpi .
  cp ../in.duplex2 .
  cp ../data.duplex2 .

  mpirun -np 1 ./lmp_mpi < in.duplex2 > /dev/null
  mv log.lammps log.$DATE.duplex2.g++.1
  grep etot log.$DATE.duplex2.g++.1 > e_test.1.dat
  grep etot ../log*1 > e_old.1.dat
  ndiff -relerr $TOL e_test.1.dat e_old.1.dat
  if [ $? -eq 0 ];
  then 
      echo "# 1 MPI-task passed" | tee -a $EXDIR/test.log
  else 
      echo "# 1 MPI-task unsuccessful" | tee -a $EXDIR/test.log
  fi

  mpirun -np 4 ./lmp_mpi < in.duplex2 > /dev/null
  mv log.lammps log.$DATE.duplex2.g++.4
  grep etot log.$DATE.duplex2.g++.4 > e_test.4.dat
  grep etot ../log*4 > e_old.4.dat
  ndiff -relerr $TOL e_test.4.dat e_old.4.dat
  if [ $? -eq 0 ];
  then 
      echo "# 4 MPI-tasks passed" | tee -a $EXDIR/test.log
  else 
      echo "# 4 MPI-tasks unsuccessful" | tee -a $EXDIR/test.log
  fi

  ######################################################

  ######################################################
  echo '# Running oxDNA2 duplex1 test' | tee -a $EXDIR/test.log
  cd $EXDIR/oxDNA2/duplex1
  mkdir test
  cd test
  cp $SRCDIR/lmp_mpi .
  cp ../in.duplex1 .
  cp ../data.duplex1 .

  mpirun -np 1 ./lmp_mpi < in.duplex1 > /dev/null
  mv log.lammps log.$DATE.duplex1.g++.1
  grep etot log.$DATE.duplex1.g++.1 > e_test.1.dat
  grep etot ../log*1 > e_old.1.dat
  ndiff -relerr $TOL e_test.1.dat e_old.1.dat
  if [ $? -eq 0 ];
  then 
      echo "# 1 MPI-task passed" | tee -a $EXDIR/test.log
  else 
      echo "# 1 MPI-task unsuccessful" | tee -a $EXDIR/test.log
  fi

  mpirun -np 4 ./lmp_mpi < in.duplex1 > /dev/null
  mv log.lammps log.$DATE.duplex1.g++.4
  grep etot log.$DATE.duplex1.g++.4 > e_test.4.dat
  grep etot ../log*4 > e_old.4.dat
  ndiff -relerr $TOL e_test.4.dat e_old.4.dat
  if [ $? -eq 0 ];
  then 
      echo "# 4 MPI-tasks passed" | tee -a $EXDIR/test.log
  else 
      echo "# 4 MPI-tasks unsuccessful" | tee -a $EXDIR/test.log
  fi

  ######################################################

  ######################################################
  echo '# Running oxDNA2 duplex2 test' | tee -a $EXDIR/test.log
  cd $EXDIR/oxDNA2/duplex2
  mkdir test
  cd test
  cp $SRCDIR/lmp_mpi .
  cp ../in.duplex2 .
  cp ../data.duplex2 .

  mpirun -np 1 ./lmp_mpi < in.duplex2 > /dev/null
  mv log.lammps log.$DATE.duplex2.g++.1
  grep etot log.$DATE.duplex2.g++.1 > e_test.1.dat
  grep etot ../log*1 > e_old.1.dat
  ndiff -relerr $TOL e_test.1.dat e_old.1.dat
  if [ $? -eq 0 ];
  then 
      echo "# 1 MPI-task passed" | tee -a $EXDIR/test.log
  else 
      echo "# 1 MPI-task unsuccessful" | tee -a $EXDIR/test.log
  fi

  mpirun -np 4 ./lmp_mpi < in.duplex2 > /dev/null
  mv log.lammps log.$DATE.duplex2.g++.4
  grep etot log.$DATE.duplex2.g++.4 > e_test.4.dat
  grep etot ../log*4 > e_old.4.dat
  ndiff -relerr $TOL e_test.4.dat e_old.4.dat
  if [ $? -eq 0 ];
  then 
      echo "# 4 MPI-tasks passed" | tee -a $EXDIR/test.log
  else 
      echo "# 4 MPI-tasks unsuccessful" | tee -a $EXDIR/test.log
  fi

  ######################################################

  ######################################################
  echo '# Running oxDNA2 duplex3 test' | tee -a $EXDIR/test.log
  cd $EXDIR/oxDNA2/duplex3
  mkdir test
  cd test
  cp $SRCDIR/lmp_mpi .
  cp ../in.duplex3 .
  cp ../data.duplex3 .

  mpirun -np 1 ./lmp_mpi < in.duplex3 > /dev/null
  mv log.lammps log.$DATE.duplex3.g++.1
  grep etot log.$DATE.duplex3.g++.1 > e_test.1.dat
  grep etot ../log*1 > e_old.1.dat
  ndiff -relerr $TOL e_test.1.dat e_old.1.dat
  if [ $? -eq 0 ];
  then 
      echo "# 1 MPI-task passed" | tee -a $EXDIR/test.log
  else 
      echo "# 1 MPI-task unsuccessful" | tee -a $EXDIR/test.log
  fi

  mpirun -np 4 ./lmp_mpi < in.duplex3 > /dev/null
  mv log.lammps log.$DATE.duplex3.g++.4
  grep etot  log.$DATE.duplex3.g++.4 > e_test.4.dat
  grep etot ../log*4 > e_old.4.dat
  ndiff -relerr $TOL e_test.4.dat e_old.4.dat
  if [ $? -eq 0 ];
  then 
      echo "# 4 MPI-tasks passed" | tee -a $EXDIR/test.log
  else 
      echo "# 4 MPI-tasks unsuccessful" | tee -a $EXDIR/test.log
  fi

  ######################################################

  ######################################################
  echo '# Running oxDNA2 unique_bp test' | tee -a $EXDIR/test.log
  cd $EXDIR/oxDNA2/unique_bp
  mkdir test
  cd test
  cp $SRCDIR/lmp_mpi .
  cp ../in.duplex4.4type .
  cp ../in.duplex4.8type .
  cp ../data.duplex4.4type .
  cp ../data.duplex4.8type .

  mpirun -np 1 ./lmp_mpi < in.duplex4.4type > /dev/null
  mv log.lammps log.$DATE.duplex4.4type.g++.1
  grep etot log.$DATE.duplex4.4type.g++.1 > e_test.4type.1.dat
  grep etot ../log*4type*1 > e_old.4type.1.dat
  ndiff -relerr $TOL e_test.4type.1.dat e_old.4type.1.dat
  if [ $? -eq 0 ];
  then 
      echo "# 1 MPI-task 4 types passed" | tee -a $EXDIR/test.log
  else 
      echo "# 1 MPI-task 4 types unsuccessful" | tee -a $EXDIR/test.log
  fi

  mpirun -np 4 ./lmp_mpi < in.duplex4.4type > /dev/null
  mv log.lammps log.$DATE.duplex4.4type.g++.4
  grep etot log.$DATE.duplex4.4type.g++.4 > e_test.4type.4.dat
  grep etot ../log*4type*4 > e_old.4type.4.dat
  ndiff -relerr $TOL e_test.4type.4.dat e_old.4type.4.dat
  if [ $? -eq 0 ];
  then 
      echo "# 4 MPI-tasks 4 types passed" | tee -a $EXDIR/test.log
  else 
      echo "# 4 MPI-tasks 4 types unsuccessful" | tee -a $EXDIR/test.log
  fi

  mpirun -np 1 ./lmp_mpi < in.duplex4.8type > /dev/null
  mv log.lammps log.$DATE.duplex4.8type.g++.1
  grep etot log.$DATE.duplex4.8type.g++.1 > e_test.8type.1.dat
  grep etot ../log*8type*1 > e_old.8type.1.dat
  ndiff -relerr $TOL e_test.8type.1.dat e_old.8type.1.dat
  if [ $? -eq 0 ];
  then 
      echo "# 1 MPI-task 8 types passed" | tee -a $EXDIR/test.log
  else 
      echo "# 1 MPI-task 8 types unsuccessful" | tee -a $EXDIR/test.log
  fi

  mpirun -np 4 ./lmp_mpi < in.duplex4.8type > /dev/null
  mv log.lammps log.$DATE.duplex4.8type.g++.4
  grep etot log.$DATE.duplex4.8type.g++.4 > e_test.8type.4.dat
  grep etot ../log*8type*4 > e_old.8type.4.dat
  ndiff -relerr $TOL e_test.8type.4.dat e_old.8type.4.dat
  if [ $? -eq 0 ];
  then 
      echo "# 4 MPI-tasks 8 types passed" | tee -a $EXDIR/test.log
  else 
      echo "# 4 MPI-tasks 8 types unsuccessful" | tee -a $EXDIR/test.log
  fi

  ######################################################

  ######################################################
  echo '# Running oxDNA2 dsring test' | tee -a $EXDIR/test.log
  cd $EXDIR/oxDNA2/dsring
  mkdir test
  cd test
  cp $SRCDIR/lmp_mpi .
  cp ../in.dsring .
  cp ../data.dsring .

  mpirun -np 1 ./lmp_mpi < in.dsring > /dev/null
  mv log.lammps log.$DATE.dsring.g++.1
  grep etot log.$DATE.dsring.g++.1 > e_test.1.dat
  grep etot ../log*1 > e_old.1.dat
  ndiff -relerr $TOL e_test.1.dat e_old.1.dat
  if [ $? -eq 0 ];
  then 
      echo "# 1 MPI-task passed" | tee -a $EXDIR/test.log
  else 
      echo "# 1 MPI-task unsuccessful" | tee -a $EXDIR/test.log
  fi

  mpirun -np 4 ./lmp_mpi < in.dsring > /dev/null
  mv log.lammps log.$DATE.dsring.g++.4
  grep etot log.$DATE.dsring.g++.4 > e_test.4.dat
  grep etot ../log*4 > e_old.4.dat
  ndiff -relerr $TOL e_test.4.dat e_old.4.dat
  if [ $? -eq 0 ];
  then 
      echo "# 4 MPI-tasks passed" | tee -a $EXDIR/test.log
  else 
      echo "# 4 MPI-tasks unsuccessful" | tee -a $EXDIR/test.log
  fi

  ######################################################

  ######################################################
  echo '# Running oxRNA2 duplex2 test' | tee -a $EXDIR/test.log
  cd $EXDIR/oxRNA2/duplex2
  mkdir test
  cd test
  cp $SRCDIR/lmp_mpi .
  cp ../in.duplex2 .
  cp ../data.duplex2 .

  mpirun -np 1 ./lmp_mpi < in.duplex2 > /dev/null
  mv log.lammps log.$DATE.duplex2.g++.1
  grep etot log.$DATE.duplex2.g++.1 > e_test.1.dat
  grep etot ../log*1 > e_old.1.dat
  ndiff -relerr $TOL e_test.1.dat e_old.1.dat
  if [ $? -eq 0 ];
  then 
      echo "# 1 MPI-task passed" | tee -a $EXDIR/test.log
  else 
      echo "# 1 MPI-task unsuccessful" | tee -a $EXDIR/test.log
  fi

  mpirun -np 4 ./lmp_mpi < in.duplex2 > /dev/null
  mv log.lammps log.$DATE.duplex2.g++.4
  grep etot log.$DATE.duplex2.g++.4 > e_test.4.dat
  grep etot ../log*4 > e_old.4.dat
  ndiff -relerr $TOL e_test.4.dat e_old.4.dat
  if [ $? -eq 0 ];
  then 
      echo "# 4 MPI-tasks passed" | tee -a $EXDIR/test.log
  else 
      echo "# 4 MPI-tasks unsuccessful" | tee -a $EXDIR/test.log
  fi

  ######################################################
  echo '# Done' | tee -a $EXDIR/test.log

elif [ $# -eq 1 ] && [ $1 = clean ]; then
  echo '# Deleting test directories'
  rm -rf $EXDIR/oxDNA/duplex1/test
  rm -rf $EXDIR/oxDNA/duplex2/test
  rm -rf $EXDIR/oxDNA2/duplex1/test
  rm -rf $EXDIR/oxDNA2/duplex2/test
  rm -rf $EXDIR/oxDNA2/duplex3/test
  rm -rf $EXDIR/oxDNA2/unique_bp/test
  rm -rf $EXDIR/oxDNA2/dsring/test
  rm -rf $EXDIR/oxRNA2/duplex2/test
  rm -rf $EXDIR/test.log
  echo '# Done'

else
  echo '# Usage:'
  echo '# ./test.sh run  ... to run test suite'
  echo '# ./test.sh clean ... to delete test directories'

fi
