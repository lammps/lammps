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
#include "string.h"
#include "library.h"
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   create an instance of LAMMPS and return pointer to it
   pass in command-line args and MPI communicator to run on
------------------------------------------------------------------------- */

void lammps_open(int argc, char **argv, MPI_Comm communicator, void **ptr)
{
  LAMMPS *lmp = new LAMMPS(argc,argv,communicator);
  *ptr = (void *) lmp;
}

/* ----------------------------------------------------------------------
   destruct an instance of LAMMPS
------------------------------------------------------------------------- */

void lammps_close(void *ptr)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  delete lmp;
}

/* ----------------------------------------------------------------------
   process an input script in filename str
------------------------------------------------------------------------- */

void lammps_file(void *ptr, char *str)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  lmp->input->file(str);
}

/* ----------------------------------------------------------------------
   process a single input command in str
------------------------------------------------------------------------- */

char *lammps_command(void *ptr, char *str)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  return lmp->input->one(str);
}

/* ----------------------------------------------------------------------
   add LAMMPS-specific library functions
   all must receive LAMMPS pointer as argument
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   extract a pointer to an internal LAMMPS value or data structure
   category: 0 = general, 1 = atom, 2 = fix, 3 = compute
   id = ID of fix or compute (else not used, just pass NULL)
   name = desired quantity, e.g. x or dt
   returns a void pointer which the caller can cast to the desired data type
   returns a NULL if LAMMPS does not recognize the name
------------------------------------------------------------------------- */

void *lammps_extract(void *ptr, int category, char *id, char *name)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  if (category == 0) {
    if (strcmp(name,"dt") == 0) return (void *) &lmp->update->dt;
    if (strcmp(name,"boxxlo") == 0) return (void *) &lmp->domain->boxlo[0];
    if (strcmp(name,"boxxhi") == 0) return (void *) &lmp->domain->boxhi[0];
    if (strcmp(name,"boxylo") == 0) return (void *) &lmp->domain->boxlo[1];
    if (strcmp(name,"boxyhi") == 0) return (void *) &lmp->domain->boxhi[1];
    if (strcmp(name,"boxzlo") == 0) return (void *) &lmp->domain->boxlo[2];
    if (strcmp(name,"boxzhi") == 0) return (void *) &lmp->domain->boxhi[2];
    return NULL;
  }

  if (category == 1) return lmp->atom->extract(name);
  if (category == 2) return lmp->modify->extract_fix(id,name);
  if (category == 3) return lmp->modify->extract_compute(id,name);

  return NULL;
}

/* ---------------------------------------------------------------------- */

int lammps_get_natoms(void *ptr)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  int natoms = static_cast<int> (lmp->atom->natoms);
  return natoms;
}

/* ---------------------------------------------------------------------- */

void lammps_get_coords(void *ptr, double *coords)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  int natoms = static_cast<int> (lmp->atom->natoms);
  double *copy = new double[3*natoms];
  for (int i = 0; i < 3*natoms; i++) copy[i] = 0.0;

  double **x = lmp->atom->x;
  int *tag = lmp->atom->tag;
  int nlocal = lmp->atom->nlocal;

  int id,offset;
  for (int i = 0; i < nlocal; i++) {
    id = tag[i];
    offset = 3*(id-1);
    copy[offset+0] = x[i][0];
    copy[offset+1] = x[i][1];
    copy[offset+2] = x[i][2];
  }

  MPI_Allreduce(copy,coords,3*natoms,MPI_DOUBLE,MPI_SUM,lmp->world);
  delete [] copy;
}

/* ---------------------------------------------------------------------- */

void lammps_put_coords(void *ptr, double *coords)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  int natoms = static_cast<int> (lmp->atom->natoms);

  double **x = lmp->atom->x;

  int m,offset;
  for (int i = 0; i < natoms; i++) {
    if ((m = lmp->atom->map(i+1)) >= 0) {
      offset = 3*i;
      x[m][0] = coords[offset+0];
      x[m][1] = coords[offset+1];
      x[m][2] = coords[offset+2];
    }
  }
}
