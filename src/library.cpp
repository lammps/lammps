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

// C or Fortran style library interface to LAMMPS
// new LAMMPS-specific functions can be added

#include "mpi.h"
#include "library.h"
#include "lammps.h"
#include "input.h"
#include "atom.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   create an instance of LAMMPS and return pointer to it
   pass in command-line args and MPI communicator to run on
------------------------------------------------------------------------- */

void lammps_open(int argc, char **argv, MPI_Comm communicator, void **ptr)
{
  LAMMPS *lammps = new LAMMPS(argc,argv,communicator);
  *ptr = (void *) lammps;
}

/* ----------------------------------------------------------------------
   destruct an instance of LAMMPS
------------------------------------------------------------------------- */

void lammps_close(void *ptr)
{
  LAMMPS *lammps = (LAMMPS *) ptr;
  delete lammps;
}

/* ----------------------------------------------------------------------
   process an input script in filename str
------------------------------------------------------------------------- */

void lammps_file(void *ptr, char *str)
{
  LAMMPS *lammps = (LAMMPS *) ptr;
  lammps->input->file(str);
}

/* ----------------------------------------------------------------------
   process a single input command in str
------------------------------------------------------------------------- */

char *lammps_command(void *ptr, char *str)
{
  LAMMPS *lammps = (LAMMPS *) ptr;
  return lammps->input->one(str);
}

/* ----------------------------------------------------------------------
   add LAMMPS-specific library functions
   all must receive LAMMPS pointer as argument
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

int lammps_get_natoms(void *ptr)
{
  LAMMPS *lammps = (LAMMPS *) ptr;
  int natoms = static_cast<int> (lammps->atom->natoms);
  return natoms;
}

/* ---------------------------------------------------------------------- */

void lammps_get_coords(void *ptr, double *coords)
{
  LAMMPS *lammps = (LAMMPS *) ptr;
  int natoms = static_cast<int> (lammps->atom->natoms);
  double *copy = new double[3*natoms];
  for (int i = 0; i < 3*natoms; i++) copy[i] = 0.0;

  double **x = lammps->atom->x;
  int *tag = lammps->atom->tag;
  int nlocal = lammps->atom->nlocal;

  int id,offset;
  for (int i = 0; i < nlocal; i++) {
    id = tag[i];
    offset = 3*(id-1);
    copy[offset+0] = x[i][0];
    copy[offset+1] = x[i][1];
    copy[offset+2] = x[i][2];
  }

  MPI_Allreduce(copy,coords,3*natoms,MPI_DOUBLE,MPI_SUM,lammps->world);
  delete [] copy;
}

/* ---------------------------------------------------------------------- */

void lammps_put_coords(void *ptr, double *coords)
{
  LAMMPS *lammps = (LAMMPS *) ptr;
  int natoms = static_cast<int> (lammps->atom->natoms);

  double **x = lammps->atom->x;

  int m,offset;
  for (int i = 0; i < natoms; i++) {
    if ((m = lammps->atom->map(i+1)) >= 0) {
      offset = 3*i;
      x[m][0] = coords[offset+0];
      x[m][1] = coords[offset+1];
      x[m][2] = coords[offset+2];
    }
  }
}
