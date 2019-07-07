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

#include "compute_ke_rigid.h"
#include <cstring>
#include "update.h"
#include "force.h"
#include "modify.h"
#include "fix.h"
#include "fix_rigid.h"
#include "fix_rigid_small.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeKERigid::ComputeKERigid(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), rfix(NULL)
{
  if (narg != 4) error->all(FLERR,"Illegal compute ke/rigid command");

  scalar_flag = 1;
  extscalar = 1;

  int n = strlen(arg[3]) + 1;
  rfix = new char[n];
  strcpy(rfix,arg[3]);
}

/* ---------------------------------------------------------------------- */

ComputeKERigid::~ComputeKERigid()
{
  delete [] rfix;
}

/* ---------------------------------------------------------------------- */

void ComputeKERigid::init()
{
  irfix = modify->find_fix(rfix);
  if (irfix < 0) error->all(FLERR,"Fix ID for compute ke/rigid does not exist");

  if (strncmp(modify->fix[irfix]->style,"rigid",5))
    error->all(FLERR,"Compute ke/rigid with non-rigid fix-ID");
}

/* ---------------------------------------------------------------------- */

double ComputeKERigid::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  if (strncmp(modify->fix[irfix]->style,"rigid",5) == 0) {
    if (strstr(modify->fix[irfix]->style,"/small")) {
      scalar = ((FixRigidSmall *) modify->fix[irfix])->extract_ke();
    } else scalar = ((FixRigid *) modify->fix[irfix])->extract_ke();
  }
  scalar *= force->mvv2e;
  return scalar;
}
