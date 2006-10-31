/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// LAMMPS as a library that can be called from another program
// C-style interface

#include "mpi.h"
#include "library.h"
#include "lammps.h"
#include "input.h"
#include "atom.h"

// variable visible to all library functions

LAMMPS *lammps;

/* ---------------------------------------------------------------------- */

void lammps_open(int argc, char **argv, MPI_Comm communicator)
{
  lammps = new LAMMPS();
  lammps->open(argc,argv,communicator);
}

/* ---------------------------------------------------------------------- */

void lammps_close()
{
  lammps->close();
  delete lammps;
}

/* ---------------------------------------------------------------------- */

void lammps_file(char *str)
{
  lammps->input->file(str);
}

/* ---------------------------------------------------------------------- */

char *lammps_command(char *str)
{
  return lammps->input->one(str);
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

int lammps_get_natoms()
{
  int natoms = static_cast<int> (lammps->atom->natoms);
  return natoms;
}

/* ---------------------------------------------------------------------- */

void lammps_get_coords(double *coords)
{
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

void lammps_put_coords(double *coords)
{
  int natoms = static_cast<int> (lammps->atom->natoms);

  double **x = lammps->atom->x;
  int nlocal = lammps->atom->nlocal;

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
