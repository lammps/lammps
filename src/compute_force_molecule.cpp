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

#include "compute_force_molecule.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeForceMolecule::ComputeForceMolecule(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute force/molecule command");

  if (atom->molecular == 0)
    error->all(FLERR,"Compute force/molecule requires molecular atom style");

  array_flag = 1;
  size_array_cols = 3;
  extarray = 0;

  // setup molecule-based data

  nmolecules = molecules_in_group(idlo,idhi);
  size_array_rows = nmolecules;

  memory->create(force,nmolecules,3,"force/molecule:force");
  memory->create(forceall,nmolecules,3,"force/molecule:forceall");
  array = forceall;
}

/* ---------------------------------------------------------------------- */

ComputeForceMolecule::~ComputeForceMolecule()
{
  memory->destroy(force);
  memory->destroy(forceall);
}

/* ---------------------------------------------------------------------- */

void ComputeForceMolecule::init()
{
  int ntmp = molecules_in_group(idlo,idhi);
  if (ntmp != nmolecules)
    error->all(FLERR,"Molecule count changed in compute force/molecule");
}

/* ---------------------------------------------------------------------- */

void ComputeForceMolecule::compute_array()
{
  tagint imol;
  double massone;
  double unwrap[3];

  invoked_array = update->ntimestep;

  for (int i = 0; i < nmolecules; i++)
    force[i][0] = force[i][1] = force[i][2] = 0.0;

  double **f = atom->f;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      force[imol][0] += f[i][0];
      force[imol][1] += f[i][1];
      force[imol][2] += f[i][2];
    }

  MPI_Allreduce(&force[0][0],&forceall[0][0],3*nmolecules,
                MPI_DOUBLE,MPI_SUM,world);
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeForceMolecule::memory_usage()
{
  double bytes = 0;
  if (molmap) bytes += (idhi-idlo+1) * sizeof(int);
  bytes += (bigint) nmolecules * 2*3 * sizeof(double);
  return bytes;
}
