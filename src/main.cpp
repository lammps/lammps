/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "lammps.h"
#include "input.h"
#include "stdio.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   main program to drive LAMMPS
------------------------------------------------------------------------- */

int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);

  //try {
    LAMMPS *lammps = new LAMMPS(argc,argv,MPI_COMM_WORLD);
    lammps->input->file();
    delete lammps;
    //}
    //catch(...) {
    //int me;
    //MPI_Comm_rank(MPI_COMM_WORLD,&me);
    //fprintf(stderr,"Unknown exception caught on rank %d\n",me);
    //}

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}
