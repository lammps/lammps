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

#include "compute_temp_uvt.h"

#include "error.h"
#include "fix_uvt.h"
#include "force.h"
#include "modify.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempUVT::ComputeTempUVT(LAMMPS *lmp, int narg, char **arg) :
    ComputeTemp(lmp, 3, arg), id_fix(nullptr), ne_dot(nullptr), ne_mass(nullptr)
{
  if (narg != 4) error->all(FLERR, "Illegal compute temp/uvt command");
  id_fix = utils::strdup(arg[3]);
}

/* ---------------------------------------------------------------------- */

ComputeTempUVT::~ComputeTempUVT()
{
  if (!copymode) delete[] id_fix;
}

/* ---------------------------------------------------------------------- */

void ComputeTempUVT::init()
{
  // Resolve the fix after construction, so FixUVT may create this compute itself.
  auto *fix = dynamic_cast<FixUVT *>(modify->get_fix_by_id(id_fix));
  if (!fix) error->all(FLERR, "Compute temp/uvt requires a fix uvt with ID {}", id_fix);
  if (fix->igroup != igroup)
    error->all(FLERR, "Compute temp/uvt and fix {} must use the same group", id_fix);

  int dim;
  ne_dot = static_cast<double *>(fix->extract("ne_dot", dim));
  ne_mass = static_cast<double *>(fix->extract("ne_mass", dim));
  if (!ne_dot || !ne_mass)
    error->all(FLERR, "Fix {} does not expose electronic state for compute temp/uvt", id_fix);
}

/* ---------------------------------------------------------------------- */

void ComputeTempUVT::dof_compute()
{
  ComputeTemp::dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR, "Temperature compute degrees of freedom < 0");

  // Normalize the nuclear sum by the combined DOF; compute_scalar adds Ne's energy.
  dof += 1.0;
  if (dof <= 0.0)
    error->all(FLERR, "Compute temp/uvt requires positive combined degrees of freedom");
  tfactor = force->mvv2e / (dof * force->boltz);
}

/* ---------------------------------------------------------------------- */

double ComputeTempUVT::compute_scalar()
{
  ComputeTemp::compute_scalar();
  // Ne is replicated on all ranks. Add its contribution once, after the nuclear sum.
  scalar += (*ne_mass) * (*ne_dot) * (*ne_dot) / (dof * force->boltz);
  return scalar;
}
