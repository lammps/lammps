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

#include "compute_pressure_alchemy.h"

#include "domain.h"
#include "error.h"
#include "fix.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePressureAlchemy::ComputePressureAlchemy(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR, "Illegal compute pressure/alchemy command");
  if (igroup) error->all(FLERR, "Compute pressure/alchemy must use group all");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 0;
  pressflag = 1;
  timeflag = 1;

  id_fix = arg[3];
  if (!modify->get_fix_by_id(id_fix))
    error->all(FLERR, "Could not find compute pressure/alchemy fix ID {} for fix alchemy", id_fix);

  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputePressureAlchemy::~ComputePressureAlchemy()
{
  delete[] vector;
}

/* ---------------------------------------------------------------------- */

void ComputePressureAlchemy::init()
{

  fix = modify->get_fix_by_id(id_fix);
  if (!fix)
    error->all(FLERR, "Could not find compute pressure/alchemy fix ID {} for fix alchemy", id_fix);

  int dim = 0;
  void *ptr = fix->extract("pressure", dim);
  if (!ptr || (dim != 1)) error->all(FLERR, "Could not extract pressure from fix alchemy");
}

/* ----------------------------------------------------------------------
   compute total pressure from tensor, averaged over Pxx, Pyy, Pzz
------------------------------------------------------------------------- */

double ComputePressureAlchemy::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if (update->vflag_global != invoked_scalar)
    error->all(FLERR, "Virial was not tallied on needed timestep");

  compute_vector();

  if (domain->dimension == 3) {
    scalar = (vector[0] + vector[1] + vector[2]) / 3.0;
  } else {
    scalar = (vector[0] + vector[1]) / 2.0;
  }
  return scalar;
}

/* ----------------------------------------------------------------------
   extract compute combined system pressure tensor from alchemy fix
------------------------------------------------------------------------- */

void ComputePressureAlchemy::compute_vector()
{
  invoked_vector = update->ntimestep;
  if (update->vflag_global != invoked_vector)
    error->all(FLERR, "Virial was not tallied on needed timestep");

  int dim = 0;
  double *pressure = (double *) fix->extract("pressure", dim);
  if (!pressure || (dim != 1)) error->all(FLERR, "Could not extract pressure from fix alchemy");

  for (int i = 0; i < 6; i++) vector[i] = pressure[i];
}
