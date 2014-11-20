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

#include "compute_vcm_molecule.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeVCMMolecule::ComputeVCMMolecule(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute vcm/molecule command");

  if (atom->molecular == 0)
    error->all(FLERR,"Compute vcm/molecule requires molecular atom style");

  array_flag = 1;
  size_array_cols = 3;
  extarray = 0;

  // setup molecule-based data

  nmolecules = molecules_in_group(idlo,idhi);
  size_array_rows = nmolecules;

  memory->create(massproc,nmolecules,"vcm/molecule:massproc");
  memory->create(masstotal,nmolecules,"vcm/molecule:masstotal");
  memory->create(vcm,nmolecules,3,"vcm/molecule:vcm");
  memory->create(vcmall,nmolecules,3,"vcm/molecule:vcmall");
  array = vcmall;

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

ComputeVCMMolecule::~ComputeVCMMolecule()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(vcm);
  memory->destroy(vcmall);
}

/* ---------------------------------------------------------------------- */

void ComputeVCMMolecule::init()
{
  int ntmp = molecules_in_group(idlo,idhi);
  if (ntmp != nmolecules)
    error->all(FLERR,"Molecule count changed in compute vcm/molecule");
}

/* ---------------------------------------------------------------------- */

void ComputeVCMMolecule::compute_array()
{
  tagint imol;
  double massone;

  invoked_array = update->ntimestep;

  for (int i = 0; i < nmolecules; i++)
    vcm[i][0] = vcm[i][1] = vcm[i][2] = 0.0;

  double **v = atom->v;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int *type = atom->type;
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
      vcm[imol][0] += v[i][0] * massone;
      vcm[imol][1] += v[i][1] * massone;
      vcm[imol][2] += v[i][2] * massone;
    }

  MPI_Allreduce(&vcm[0][0],&vcmall[0][0],3*nmolecules,
                MPI_DOUBLE,MPI_SUM,world);
  for (int i = 0; i < nmolecules; i++) {
    vcmall[i][0] /= masstotal[i];
    vcmall[i][1] /= masstotal[i];
    vcmall[i][2] /= masstotal[i];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeVCMMolecule::memory_usage()
{
  double bytes = (bigint) nmolecules * 2 * sizeof(double);
  if (molmap) bytes += (idhi-idlo+1) * sizeof(int);
  bytes += (bigint) nmolecules * 2*3 * sizeof(double);
  return bytes;
}
