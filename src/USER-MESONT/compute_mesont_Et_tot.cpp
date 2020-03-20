/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: Maxim Shugaev (UVA), mvs9t@virginia.edu
------------------------------------------------------------------------- */

#include "compute_mesont_Et_tot.h"
#include <mpi.h>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "force.h"
#include "modify.h"
#include "domain.h"
#include "error.h"
#include "pair.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeMESONT_Et_tot::ComputeMESONT_Et_tot(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute mesont/Et_tot command");
  if (igroup) error->all(FLERR,"Compute mesont/Et_tot must use group all");
  timeflag = 1;
  extscalar = 1;
  scalar_flag = 1;
}

/* ---------------------------------------------------------------------- */

double ComputeMESONT_Et_tot::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if (update->eflag_global != invoked_scalar)
    error->all(FLERR,"Energy was not tallied on needed timestep");

  int i;
  double* ptr = static_cast<double*>(force->pair->extract("mesonttpm_Et_tot",i));
  if(!ptr) error->all(FLERR,
   "mesont/Et_tot is allowed only with mesont/tpm pair style");
  MPI_Allreduce(ptr,&scalar,1,MPI_DOUBLE,MPI_SUM,world);

  return scalar;
}
