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

/* ------------------------------------------------------------------------
   Contributing author: Zhengtao Huang (The University of Hong Kong)
                        hzt990224@gmail.com
------------------------------------------------------------------------- */

#include "compute_temp_tspin.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "tspin.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempTSpin::ComputeTempTSpin(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), index_vs(-1), index_sm(-1), tfactor(0.0)
{
  if (narg != 3) error->all(FLERR, "Incorrect number of arguments for compute temp/tspin");
  if (!atom->sp_flag) error->all(FLERR, "Compute temp/tspin requires atom style spin");

  scalar_flag = 1;
  extscalar = 0;

  // A spin-temperature diagnostic must not supply translational kinetic
  // energy or a kinetic pressure tensor to the ordinary temperature users.

  tempflag = 0;
}

/* ---------------------------------------------------------------------- */

void ComputeTempTSpin::init()
{
  tfactor = force->mvv2e / force->boltz;
  tspin_bind_state(atom, error, "Compute temp/tspin", index_vs, index_sm);
}

/* ---------------------------------------------------------------------- */

double ComputeTempTSpin::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **v_s = atom->darray[index_vs];
  double *s_mass = atom->dvector[index_sm];
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  double local[2] = {0.0, 0.0};
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit) || (s_mass[i] <= 0.0)) continue;
    local[0] += s_mass[i] *
        (v_s[i][0] * v_s[i][0] + v_s[i][1] * v_s[i][1] + v_s[i][2] * v_s[i][2]);
    local[1] += 3.0;
  }

  double total[2];
  MPI_Allreduce(local, total, 2, MPI_DOUBLE, MPI_SUM, world);
  dof = total[1];
  scalar = (dof > 0.0) ? tfactor * total[0] / dof : 0.0;
  return scalar;
}
