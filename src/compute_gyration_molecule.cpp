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

#include "math.h"
#include "string.h"
#include "compute_gyration_molecule.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeGyrationMolecule::ComputeGyrationMolecule(LAMMPS *lmp,
                                                 int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute gyration/molecule command");

  if (atom->molecular == 0)
    error->all(FLERR,"Compute gyration/molecule requires molecular atom style");

  tensor = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"tensor") == 0) {
      tensor = 1;
      iarg++;
    } else error->all(FLERR,"Illegal compute gyration/molecule command");
  }

  // setup molecule-based data

  nmolecules = molecules_in_group(idlo,idhi);

  memory->create(massproc,nmolecules,"gyration/molecule:massproc");
  memory->create(masstotal,nmolecules,"gyration/molecule:masstotal");
  memory->create(com,nmolecules,3,"gyration/molecule:com");
  memory->create(comall,nmolecules,3,"gyration/molecule:comall");

  rg = vector = NULL;
  rgt = array = NULL;

  if (tensor) {
    memory->create(rgt,nmolecules,6,"gyration/molecule:rgt");
    memory->create(array,nmolecules,6,"gyration/molecule:array");
    array_flag = 1;
    size_array_rows = nmolecules;
    size_array_cols = 6;
    extarray = 0;
  } else {
    memory->create(rg,nmolecules,"gyration/molecule:rg");
    memory->create(vector,nmolecules,"gyration/molecule:vector");
    vector_flag = 1;
    size_vector = nmolecules;
    extvector = 0;
  }

  // compute masstotal for each molecule

  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  tagint imol;
  double massone;

  for (int i = 0; i < nmolecules; i++) massproc[i] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      massproc[imol] += massone;
    }

  MPI_Allreduce(massproc,masstotal,nmolecules,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

ComputeGyrationMolecule::~ComputeGyrationMolecule()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(rg);
  memory->destroy(rgt);
}

/* ---------------------------------------------------------------------- */

void ComputeGyrationMolecule::init()
{
  int ntmp = molecules_in_group(idlo,idhi);
  if (ntmp != nmolecules)
    error->all(FLERR,"Molecule count changed in compute gyration/molecule");
}

/* ---------------------------------------------------------------------- */

void ComputeGyrationMolecule::compute_vector()
{
  tagint imol;
  double dx,dy,dz,massone;
  double unwrap[3];

  invoked_array = update->ntimestep;

  molcom();

  for (int i = 0; i < nmolecules; i++) rg[i] = 0.0;

  double **x = atom->x;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - comall[imol][0];
      dy = unwrap[1] - comall[imol][1];
      dz = unwrap[2] - comall[imol][2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      rg[imol] += (dx*dx + dy*dy + dz*dz) * massone;
    }

  MPI_Allreduce(rg,vector,nmolecules,MPI_DOUBLE,MPI_SUM,world);

  for (int i = 0; i < nmolecules; i++)
    vector[i] = sqrt(vector[i]/masstotal[i]);
}

/* ---------------------------------------------------------------------- */

void ComputeGyrationMolecule::compute_array()
{
  int i,j;
  tagint imol;
  double dx,dy,dz,massone;
  double unwrap[3];

  invoked_array = update->ntimestep;

  molcom();

  for (i = 0; i < nmolecules; i++)
    for (j = 0; j < 6; j++)
      rgt[i][j] = 0.0;

  double **x = atom->x;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - comall[imol][0];
      dy = unwrap[1] - comall[imol][1];
      dz = unwrap[2] - comall[imol][2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      rgt[imol][0] += dx*dx * massone;
      rgt[imol][1] += dy*dy * massone;
      rgt[imol][2] += dz*dz * massone;
      rgt[imol][3] += dx*dy * massone;
      rgt[imol][4] += dx*dz * massone;
      rgt[imol][5] += dy*dz * massone;
    }

  if (nmolecules)
    MPI_Allreduce(&rgt[0][0],&array[0][0],nmolecules*6,
                  MPI_DOUBLE,MPI_SUM,world);

  for (i = 0; i < nmolecules; i++)
    for (j = 0; j < 6; j++)
      array[i][j] /= masstotal[i];
}


/* ----------------------------------------------------------------------
   calculate per-molecule COM
------------------------------------------------------------------------- */

void ComputeGyrationMolecule::molcom()
{
  tagint imol;
  double massone;
  double unwrap[3];

  for (int i = 0; i < nmolecules; i++)
    com[i][0] = com[i][1] = com[i][2] = 0.0;

  double **x = atom->x;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      domain->unmap(x[i],image[i],unwrap);
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      com[imol][0] += unwrap[0] * massone;
      com[imol][1] += unwrap[1] * massone;
      com[imol][2] += unwrap[2] * massone;
    }

  MPI_Allreduce(&com[0][0],&comall[0][0],3*nmolecules,
                MPI_DOUBLE,MPI_SUM,world);
  for (int i = 0; i < nmolecules; i++) {
    comall[i][0] /= masstotal[i];
    comall[i][1] /= masstotal[i];
    comall[i][2] /= masstotal[i];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeGyrationMolecule::memory_usage()
{
  double bytes = (bigint) nmolecules * 2 * sizeof(double);
  if (molmap) bytes += (idhi-idlo+1) * sizeof(int);
  bytes += (bigint) nmolecules * 2*3 * sizeof(double);
  if (tensor) bytes += (bigint) nmolecules * 2*6 * sizeof(double);
  else bytes += (bigint) nmolecules * 2 * sizeof(double);
  return bytes;
}
