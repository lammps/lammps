// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ----------------------------------------------------------------------- */


#include "compute_momentum.h"

#include "atom.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;

ComputeMomentum::ComputeMomentum(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute momentum command");

  vector_flag = 1;
  size_vector = 3;
  extvector = 1;
  vector = new double[size_vector];
}

ComputeMomentum::~ComputeMomentum() {
  delete[] vector;
}

void ComputeMomentum::init()
{
}

void ComputeMomentum::compute_vector()
{
  invoked_vector = update->ntimestep;

  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double mom[3] = {0.0, 0.0, 0.0};

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        for (int j = 0; j < 3; ++j)
          mom[j] += rmass[i] * v[i][j];
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        for (int j = 0; j < 3; ++j)
          mom[j] += mass[type[i]] * v[i][j];
  }

  MPI_Allreduce(&mom, vector, 3, MPI_DOUBLE, MPI_SUM, world);
}
