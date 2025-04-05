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

/* ----------------------------------------------------------------------
   Contributing author: Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#include "compute_quadrupole.h"

#include "atom.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeQuadrupole::ComputeQuadrupole(LAMMPS *lmp, int narg, char **arg)
  : Compute(lmp, narg, arg) {
  if (narg != 3) error->all(FLERR, "Illegal compute quadrupole command");

  vector_flag = 1;
  size_vector = 6;  // Q_xx, Q_yy, Q_zz, Q_xy, Q_xz, Q_yz
  extvector = 1;

  vector = new double[size_vector];

}

/* ---------------------------------------------------------------------- */

ComputeQuadrupole::~ComputeQuadrupole()
{
  delete[] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeQuadrupole::compute_vector() {
  invoked_vector = update->ntimestep;

  double **x = atom->x;
  double *q = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < 6; i++) vector[i] = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double xi = x[i][0];
      double yi = x[i][1];
      double zi = x[i][2];
      double qi = q[i];

      // Populate LAMMPS vector (only 6 independent symmetric components)
      vector[0] += qi * xi * xi;  // Q_xx
      vector[1] += qi * yi * yi;  // Q_yy
      vector[2] += qi * zi * zi;  // Q_zz
      vector[3] += qi * xi * yi;  // Q_xy
      vector[4] += qi * xi * zi;  // Q_xz
      vector[5] += qi * yi * zi;  // Q_yz
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, vector, 6, MPI_DOUBLE, MPI_SUM, world);

}

