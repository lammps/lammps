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
  delete [] rfix;
}

/* ---------------------------------------------------------------------- */

void ComputeERotateRigid::init()
{
  irfix = modify->find_fix(rfix);
  if (irfix < 0)
    error->all(FLERR,"Fix ID for compute erotate/rigid does not exist");

  if (strncmp(modify->fix[irfix]->style,"rigid",5))
    error->all(FLERR,"Compute erotate/rigid with non-rigid fix-ID");
}

/* ---------------------------------------------------------------------- */

double ComputeERotateRigid::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  if (strncmp(modify->fix[irfix]->style,"rigid",5) == 0) {
    if (strstr(modify->fix[irfix]->style,"/small")) {
      scalar = (dynamic_cast<FixRigidSmall *>(modify->fix[irfix]))->extract_erotational();
    } else scalar = (dynamic_cast<FixRigid *>(modify->fix[irfix]))->extract_erotational();
  }
  scalar *= force->mvv2e;
  return scalar;
}
