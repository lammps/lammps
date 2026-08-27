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

   Kinetic energy stored in the spin degrees of freedom of inertial spin
   dynamics, sum over the group of 1/2 s_mass |v_s|^2.  Because this energy
   is not part of the thermodynamic keyword ke, add it explicitly to the
   thermo output with a c_ID reference when running fix nve/tspin,
   fix nvt/tspin or fix langevin/tspin.
------------------------------------------------------------------------- */

#include "compute_ke_tspin.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "tspin.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeKETSpin::ComputeKETSpin(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), index_vs(-1), index_sm(-1), pfactor(0.0)
{
  if (narg != 3) error->all(FLERR, "Illegal compute ke/tspin command");

  if (!atom->sp_flag)
    error->all(FLERR, "Compute ke/tspin requires atom style spin");

  scalar_flag = 1;
  extscalar = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeKETSpin::init()
{
  pfactor = 0.5 * force->mvv2e;
  tspin_bind_state(atom, error, "Compute ke/tspin", index_vs, index_sm);
}

/* ---------------------------------------------------------------------- */

double ComputeKETSpin::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **v_s = atom->darray[index_vs];
  double *s_mass = atom->dvector[index_sm];
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  double ke = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      ke += s_mass[i] * (v_s[i][0] * v_s[i][0] + v_s[i][1] * v_s[i][1] + v_s[i][2] * v_s[i][2]);

  MPI_Allreduce(&ke, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
  scalar *= pfactor;
  return scalar;
}
