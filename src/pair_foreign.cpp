/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
   See the README file in the top-level LAMMPS directory.
   Pair zero is a dummy pair interaction useful for requiring a
   force cutoff distance in the absence of pair-interactions or
   with hybrid/overlay if a larger force cutoff distance is required.
   This can be used in conjunction with bond/create to create bonds
   that are longer than the cutoff of a given force field, or to
   calculate radial distribution functions for models without
   pair interactions.
------------------------------------------------------------------------- */

#include "pair_foreign.h"

#include "atom.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;

PairForeign::PairForeign(LAMMPS *lmp, void* ctx, ForeignCompute compute) : Pair(lmp), ctx(ctx), compute_fptr(compute){
}

/* ---------------------------------------------------------------------- */

PairForeign::~PairForeign()
{
   if (allocated) {
      memory->destroy(setflag);
      memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairForeign::compute(int eflag, int vflag)
{
   compute_fptr(ctx, eflag, vflag);
}

/* ---------------------------------------------------------------------- */

void PairForeign::compute_outer(int eflag, int vflag)
{
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairForeign::settings(int narg, char **arg)
{
   if (narg != 1)
      error->all(FLERR,"Illegal pair_style command");

   cut_global = utils::numeric(FLERR,arg[0],false,lmp);
}

void PairForeign::allocate() {
   allocated = 1;
   int n = atom->ntypes;

   memory->create(setflag,n+1,n+1,"pair:setflag");
   for (int i = 1; i <= n; i++)
      for (int j = i; j <= n; j++)
         setflag[i][j] = 0;

   memory->create(cutsq,n+1,n+1,"pair:cutsq");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairForeign::coeff(int narg, char **arg)
{  
   const int ntypes = atom->ntypes;

   if (!allocated) allocate();

   for (int i = 1; i <= ntypes ; i++) {
      for (int j = i; j <= ntypes ; j++) {
         setflag[i][j] = 1;
         cutsq[i][j] = cut_global*cut_global;
      }
   }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairForeign::init_one(int i, int j)
{
   return cut_global;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairForeign::write_restart(FILE *fp)
{
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairForeign::read_restart(FILE *fp)
{
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairForeign::write_restart_settings(FILE *fp)
{
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairForeign::read_restart_settings(FILE *fp)
{
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairForeign::write_data(FILE *fp)
{
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairForeign::write_data_all(FILE *fp)
{
}

/* ---------------------------------------------------------------------- */

double PairForeign::single(int /*i*/, int /*j*/, int /* itype */, int /* jtype */,
                        double /* rsq */, double /*factor_coul*/,
                        double /* factor_lj */, double &fforce)
{
}

