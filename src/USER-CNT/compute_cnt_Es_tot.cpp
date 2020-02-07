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

#include "compute_cnt_Es_tot.h"
#include "pair_cnt_tpm.h"
#include <mpi.h>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "force.h"
#include "bond.h"
#include "modify.h"
#include "domain.h"
#include "error.h"
#include <typeinfo>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCNT_Es_tot::ComputeCNT_Es_tot(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute cnt/Es_tot command");
  if (igroup) error->all(FLERR,"Compute cnt/Es_tot must use group all");
  timeflag = 1;
  extscalar = 1;
  scalar_flag = 1;
}

/* ---------------------------------------------------------------------- */

double ComputeCNT_Es_tot::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if (update->eflag_global != invoked_scalar)
    error->all(FLERR,"Energy was not tallied on needed timestep");

  double one = 0.0;
  if (force->pair){
    try {
      PairCNTTPM* pair = dynamic_cast<PairCNTTPM*>(force->pair);
      one += pair->energy_s;
    }
    catch (std::bad_cast& bc){
      error->all(FLERR,"cnt/Es_tot is allowed only with cnt pair style");
    }
  }
  MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_SUM,world);

  return scalar;
}
