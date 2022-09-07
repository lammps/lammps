/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Akhlak Mahmood

   Contact:
     Department of Materials Science and Engineering,
     North Carolina State University,
     Raleigh, NC, USA

     amahmoo3@ncsu.edu; mahmoodakhlak@gmail.com
------------------------------------------------------------------------- */

#include "compute_mspin_distance.h"
#include "error.h"
#include "fix_mspin_nh.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include <cstdio>
#include <cstring>

using namespace LAMMPS_NS;
using namespace std;

ComputeMSDistance::ComputeMSDistance(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), rfix(NULL)
{
  if (narg != 6) error->all(FLERR, "Illegal compute mspin/distance command.");

  vector_flag = 1;
  size_vector = 1;
  extvector = 0;

  // get the fix id
  int n = strlen(arg[3]) + 1;
  rfix = new char[n];
  strcpy(rfix, arg[3]);

  ibody = utils::inumeric(FLERR, arg[4], false, lmp) - 1;
  jbody = utils::inumeric(FLERR, arg[5], false, lmp) - 1;

  memory->create(vector, size_vector, "compute/mspin:distance");
}

ComputeMSDistance::~ComputeMSDistance()
{
  delete[] rfix;
  memory->destroy(vector);
}

void ComputeMSDistance::init()
{
  irfix = modify->find_fix(rfix);
  if (irfix < 0) error->all(FLERR, "Fix ID for compute mspin/distance does not exist");

  if (!utils::strmatch(modify->fix[irfix]->style,"^rigid/n.t/mspin"))
    error->all(FLERR, "Compute mspin with non-mspin fix-ID");
}

void ComputeMSDistance::compute_vector()
{
  double dist;
  invoked_vector = update->ntimestep;

  if (utils::strmatch(modify->fix[irfix]->style,"^rigid/n.t/mspin")) {
    dist = ((FixMspinNH *) modify->fix[irfix])->extract_distance(ibody, jbody);
  }

  vector[0] = dist;
}
