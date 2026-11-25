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

   Kinetic temperature and kinetic energy of the spin degrees of freedom of
   inertial spin dynamics.  The scalar is the temperature and the first vector
   element the energy, sum over the group of 1/2 s_mass |v_s|^2.  Because this
   energy is not part of the thermodynamic keyword ke, add it explicitly to
   the thermo output with a c_ID[1] reference when running fix nve/tspin,
   fix nvt/tspin or fix langevin/tspin.
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
    Compute(lmp, narg, arg), index_vs(-1), index_sm(-1)
{
  if (narg != 3) error->all(FLERR, "Incorrect number of arguments for compute temp/tspin");
  if (!atom->sp_flag) error->all(FLERR, "Compute temp/tspin requires atom style spin");

  scalar_flag = vector_flag = 1;
  size_vector = 1;
  extscalar = 0;
  extvector = 1;

  // A spin-temperature diagnostic must not supply translational kinetic
  // energy or a kinetic pressure tensor to the ordinary temperature users.

  tempflag = 0;

  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputeTempTSpin::~ComputeTempTSpin()
{
  delete[] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempTSpin::init()
{
  tspin_bind_state(atom, error, "Compute temp/tspin", index_vs, index_sm);
}

/* ---------------------------------------------------------------------- */

double ComputeTempTSpin::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  compute_vector();
  scalar = (dof > 0.0) ? 2.0 * vector[0] / (dof * force->boltz) : 0.0;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempTSpin::compute_vector()
{
  invoked_vector = update->ntimestep;

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
  vector[0] = 0.5 * force->mvv2e * total[0];
  dof = total[1];
}
