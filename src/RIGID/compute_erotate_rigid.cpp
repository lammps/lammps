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

#include "mpi.h"
#include "string.h"
#include "compute_erotate_rigid.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_rigid.h"
#include "fix_rigid_small.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeERotateRigid::ComputeERotateRigid(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal compute erotate/rigid command");

  scalar_flag = 1;
  extscalar = 1;

  int n = strlen(arg[3]) + 1;
  rfix = new char[n];
  strcpy(rfix,arg[3]);
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

  int flag = 1;
  if (strcmp(modify->fix[irfix]->style,"rigid/small") == 0) flag = 0;
  else if (strstr(modify->fix[irfix]->style,"rigid")) flag = 0;
  if (flag) error->all(FLERR,"Compute ke/rigid with non-rigid fix-ID");
}

/* ---------------------------------------------------------------------- */

double ComputeERotateRigid::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  if (strcmp(modify->fix[irfix]->style,"rigid/small") == 0)
    scalar = ((FixRigidSmall *) modify->fix[irfix])->extract_erotational();
  else if (strstr(modify->fix[irfix]->style,"rigid"))
    scalar = ((FixRigid *) modify->fix[irfix])->extract_erotational();

  return scalar;
}
