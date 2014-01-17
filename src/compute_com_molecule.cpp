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

#include "compute_com_molecule.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCOMMolecule::ComputeCOMMolecule(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute com/molecule command");

  if (atom->molecular == 0)
    error->all(FLERR,"Compute com/molecule requires molecular atom style");

  array_flag = 1;
  size_array_cols = 3;
  extarray = 0;

  // setup molecule-based data

  nmolecules = molecules_in_group(idlo,idhi);
  size_array_rows = nmolecules;

  memory->create(massproc,nmolecules,"com/molecule:massproc");
  memory->create(masstotal,nmolecules,"com/molecule:masstotal");
  memory->create(com,nmolecules,3,"com/molecule:com");
  memory->create(comall,nmolecules,3,"com/molecule:comall");
  array = comall;

  // compute masstotal for each molecule

  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  int i,imol;
  double massone;

  for (i = 0; i < nmolecules; i++) massproc[i] = 0.0;

  for (i = 0; i < nlocal; i++)
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

ComputeCOMMolecule::~ComputeCOMMolecule()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
}

/* ---------------------------------------------------------------------- */

void ComputeCOMMolecule::init()
{
  int ntmp = molecules_in_group(idlo,idhi);
  if (ntmp != nmolecules)
    error->all(FLERR,"Molecule count changed in compute com/molecule");
}

/* ---------------------------------------------------------------------- */

void ComputeCOMMolecule::compute_array()
{
  int i,imol;
  double massone;
  double unwrap[3];

  invoked_array = update->ntimestep;

  for (i = 0; i < nmolecules; i++)
    com[i][0] = com[i][1] = com[i][2] = 0.0;

  double **x = atom->x;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      domain->unmap(x[i],image[i],unwrap);
      com[imol][0] += unwrap[0] * massone;
      com[imol][1] += unwrap[1] * massone;
      com[imol][2] += unwrap[2] * massone;
    }

  MPI_Allreduce(&com[0][0],&comall[0][0],3*nmolecules,
                MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmolecules; i++) {
    comall[i][0] /= masstotal[i];
    comall[i][1] /= masstotal[i];
    comall[i][2] /= masstotal[i];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeCOMMolecule::memory_usage()
{
  double bytes = 2*nmolecules * sizeof(double);
  if (molmap) bytes += (idhi-idlo+1) * sizeof(int);
  bytes += 2*nmolecules*3 * sizeof(double);
  return bytes;
}
