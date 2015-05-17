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

#include "string.h"
#include "compute_tally_stress.h"
#include "atom.h"
#include "pair.h"
#include "update.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTallyStress::ComputeTallyStress(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute tally/stress command");
  scalar_flag = 1;
  extscalar = 0;
  peflag = 1;                   // we need Pair::ev_tally() to be run
  timeflag = 1;
}

/* ---------------------------------------------------------------------- */

ComputeTallyStress::~ComputeTallyStress()
{
  if (force->pair) force->pair->del_tally_callback(this);
}

/* ---------------------------------------------------------------------- */

void ComputeTallyStress::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Trying to use compute tally/stress with no pair style");
  else
    force->pair->add_tally_callback(this);

  did_compute = 0;
}


/* ---------------------------------------------------------------------- */
void ComputeTallyStress::pair_tally_callback(int i, int j, int, int,
                                         double, double, double,
                                         double, double, double)
{
  did_compute = 1;
}
