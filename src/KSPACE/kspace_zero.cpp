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

#include "kspace_zero.h"
#include <mpi.h>
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

KSpaceZero::KSpaceZero(LAMMPS *lmp) : KSpace(lmp)
{
  ewaldflag = 1;
}

/* ---------------------------------------------------------------------- */

void KSpaceZero::init()
{
  if (comm->me == 0) {
    if (screen) fprintf(screen,"KSpaceZero initialization ...\n");
    if (logfile) fprintf(logfile,"KSpaceZero initialization ...\n");
  }

  if (force->pair == NULL)
    error->all(FLERR,"KSpace solver requires a pair style");

  if (!force->pair->ewaldflag && !force->pair->pppmflag && !force->pair->msmflag)
    error->all(FLERR,"KSpace style is incompatible with Pair style");

  int itmp;
  double *p_cutoff = (double *) force->pair->extract("cut_coul",itmp);
  if (p_cutoff == NULL)
    error->all(FLERR,"KSpace style is incompatible with Pair style");

  // compute qsum & qsqsum and warn if not charge-neutral

  scale = 1.0;
  qqrd2e = force->qqrd2e;
  qsum_qsq();

  // set accuracy (force units) from accuracy_relative or accuracy_absolute

  accuracy = 0.0;
  setup();
}

/* ----------------------------------------------------------------------
   compute the KSpaceZero long-range force, energy, virial
------------------------------------------------------------------------- */

void KSpaceZero::compute(int eflag, int vflag)
{
  // set energy/virial flags

  ev_init(eflag,vflag);
}

/* ----------------------------------------------------------------------
   memory usage of local arrays
------------------------------------------------------------------------- */

double KSpaceZero::memory_usage()
{
  return (double)sizeof(KSpaceZero);
}



