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
------------------------------------------------------------------------- */

#include "compute_erotate_rigid.h"

#include "error.h"
#include "fix.h"
#include "fix_rigid.h"
#include "fix_rigid_small.h"
#include "force.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeERotateRigid::ComputeERotateRigid(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), rfix(nullptr)
{
  if (narg != 4) error->all(FLERR,"Illegal compute erotate/rigid command");

  scalar_flag = 1;
  extscalar = 1;

  rfix = utils::strdup(arg[3]);
}

/* ---------------------------------------------------------------------- */

ComputeERotateRigid::~ComputeERotateRigid()
{
  delete[] rfix;
}

/* ---------------------------------------------------------------------- */

void ComputeERotateRigid::init()
{
  irfix = modify->get_fix_by_id(rfix);
  if (!irfix) error->all(FLERR,"Fix ID {} for compute erotate/rigid does not exist", rfix);

  if (!utils::strmatch(irfix->style,"^rigid"))
    error->all(FLERR,"Compute erotate/rigid with non-rigid fix {}", rfix);
}

/* ---------------------------------------------------------------------- */

double ComputeERotateRigid::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  if (utils::strmatch(irfix->style,"^rigid")) {
    if (auto *smallfix = dynamic_cast<FixRigidSmall *>(irfix)) {
      scalar = smallfix->extract_erotational();
    } else if (auto *plainfix =dynamic_cast<FixRigid *>(irfix)) {
      scalar = plainfix->extract_erotational();
    } else scalar = 0.0;
  }
  scalar *= force->mvv2e;
  return scalar;
}
