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

/* ----------------------------------------------------------------------
   Contributing author: Richard Berger (Temple U)
------------------------------------------------------------------------- */

#include "fix_pour_opt.h"
#include "update.h"
#include "comm.h"
#include "domain.h"
#include "region.h"
#include "error.h"
#include "cell_list.h"

using namespace LAMMPS_NS;
using namespace FixConst;


/* ---------------------------------------------------------------------- */

FixPourOpt::FixPourOpt(LAMMPS *lmp, int narg, char **arg) :
    FixPour(lmp, narg, arg), clist(lmp), dist_clist(lmp)
{
}

/* ---------------------------------------------------------------------- */

FixPourOpt::~FixPourOpt()
{
}


/* ----------------------------------------------------------------------
   perform particle insertion
------------------------------------------------------------------------- */

void FixPourOpt::init_near_lists(int, INearList*&  local_nlist, IDistributedNearList*& dist_nlist)
{
  double x_lo[3], x_hi[3];
  domain->regions[iregion]->get_bounding_box(x_lo, x_hi);

  if (domain->dimension == 3) {
    double lo_current = zlo + (update->ntimestep - nfirst) * update->dt * rate;
    double hi_current = zhi + (update->ntimestep - nfirst) * update->dt * rate;
    x_lo[2] = lo_current;
    x_hi[2] = hi_current;
  } else {
    double lo_current = ylo + (update->ntimestep - nfirst) * update->dt * rate;
    double hi_current = yhi + (update->ntimestep - nfirst) * update->dt * rate;
    x_lo[1] = lo_current;
    x_hi[1] = hi_current;
  }

  clist.clear();
  clist.setup(x_lo, x_hi, 4*radius_max);
  local_nlist = &clist;

  if(comm->nprocs > 1) {
    dist_clist.clear();
    dist_clist.setup(x_lo, x_hi, 4*radius_max);
    dist_nlist = &dist_clist;
  }
}

void FixPourOpt::cleanup_near_lists(INearList*&, IDistributedNearList*&)
{
  // no cleanup, we'll reuse existing allocations for next insertion
}

