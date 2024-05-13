/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// example of how a calling program can invoke
//   multiple instances of LAMMPS, each on a subset of procs,
//   and have each of them perform a different run, and collect output info
// Syntax: mpirun -np P multiple N in.lammps T Tdelta
//         P = # of total procs to run the driver program on
//         N = # of instances of LAMMPS to create, P must be divisible by N
//         in.lammps = LAMMPS input script,
//                     which takes "t" as a temperature variable
//         T = baseline temperature
//         Tdelta = incremental temperature for each of N runs
// See README for compilation instructions

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define LAMMPS_LIB_MPI  // to make lammps_open() visible
#include "library.h"

int main(int narg, char **arg)
{
  // setup MPI and various communicators
  // driver runs on all procs in MPI_COMM_WORLD
  // comm_lammps only has 1st P procs (could be all or any subset)

  MPI_Init(&narg,&arg);

  if (narg != 5) {
    printf("Syntax: multiple N in.lammps T Tdelta\n");
    exit(1);
  }

  int me,nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  int ninstance = atoi(arg[1]);
  if (nprocs % ninstance) {
    if (me == 0)
      printf("ERROR: Total procs must be divisble by N\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  char *infile = arg[2];
  double temperature = atof(arg[3]);
  double tdelta = atof(arg[4]);

  // create one communicator per instancem each with P/N procs

  int instance = me*ninstance / nprocs;
  MPI_Comm comm_lammps;
  MPI_Comm_split(MPI_COMM_WORLD,instance,0,&comm_lammps);

  // each instance: unique screen file, log file, temperature

  char str1[32],str2[32],str3[32];

  char **lmparg = new char*[8];
  lmparg[0] = (char *) "LAMMPS";              // required placeholder for program name
  lmparg[1] = (char *) "-screen";
  sprintf(str1,"screen.%d",instance);
  lmparg[2] = str1;
  lmparg[3] = (char *) "-log";
  sprintf(str2,"log.lammps.%d",instance);
  lmparg[4] = str2;
  lmparg[5] = (char *) "-var";
  lmparg[6] = (char *) "t";
  sprintf(str3,"%g",temperature + instance*tdelta);
  lmparg[7] = str3;

  // create N instances of LAMMPS

  void *lmp = lammps_open(8,lmparg,comm_lammps,NULL);

  delete [] lmparg;

  // run input script thru all instances of LAMMPS

  lammps_file(lmp,infile);

  // query final temperature and print result for each instance

  double *ptr = (double *)
    lammps_extract_compute(lmp,"thermo_temp",LMP_STYLE_GLOBAL,LMP_TYPE_SCALAR);
  double finaltemp = *ptr;

  double *temps = new double[ninstance];
  for (int i = 0; i < ninstance; i++) temps[i] = 0.0;

  int me_lammps;
  MPI_Comm_rank(comm_lammps,&me_lammps);
  if (me_lammps == 0) temps[instance] = finaltemp;

  double *alltemps = new double[ninstance];
  MPI_Allreduce(temps,alltemps,ninstance,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  if (me == 0)
    for (int i = 0; i < ninstance; i++)
      printf("Instance %d, final temp = %g\n",i+1,alltemps[i]);

  delete [] temps;
  delete [] alltemps;

  // delete LAMMPS instances

  lammps_close(lmp);

  // close down MPI

  MPI_Comm_free(&comm_lammps);
  MPI_Finalize();
}
