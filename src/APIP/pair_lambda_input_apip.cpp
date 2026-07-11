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
/* ----------------------------------------------------------------------
   Contributing author: David Immel (d.immel@fz-juelich.de, FZJ, Germany)
------------------------------------------------------------------------- */

#include "pair_lambda_input_apip.h"

#include "atom.h"
#include "error.h"
#include "fix_lambda_apip.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLambdaInputAPIP::PairLambdaInputAPIP(LAMMPS *lmp) : Pair(lmp), fix_lambda(nullptr), cut(nullptr)
{

  cut_global = -1;
  ignore_group_bit = 0;
  timer = 0;
  n_calculations = 0;
  time_per_atom = -1;
}

/* ---------------------------------------------------------------------- */

PairLambdaInputAPIP::~PairLambdaInputAPIP()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
  }
}

/* ---------------------------------------------------------------------- */
void PairLambdaInputAPIP::coeff(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      setflag[i][j] = 1;
      cut[i][j] = cut_global;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLambdaInputAPIP::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;

  memory->create(setflag, n, n, "pair:setflag");
  for (int i = 1; i < n; i++)
    for (int j = i; j < n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n, n, "pair:cutsq");
  memory->create(cut, n, n, "pair:cut");
}

/* ---------------------------------------------------------------------- */

void PairLambdaInputAPIP::compute(int eflag, int vflag)
{
  // basic stuff (see pair_zero)
  ev_init(eflag, vflag);
  if (vflag_fdotr) virial_fdotr_compute();

  double timer_start = platform::walltime();

  n_calculations += calculate_lambda_input();

  timer += platform::walltime() - timer_start;

  fix_lambda->update_lambda_input_history();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLambdaInputAPIP::settings(int narg, char **arg)
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "pair_style lambda_input", error);

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);

  // reset cutoffs that have been explicitly set
  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLambdaInputAPIP::init_style()
{
  if (!atom->apip_lambda_input_flag)
    error->all(FLERR, "pair_lambda input requires an atom style with lambda_input");

  // find fix lambda
  int count = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style, "lambda/apip") == 0) {
      fix_lambda = (FixLambdaAPIP *) modify->fix[i];
      count++;
    }
  }
  if (count != 1) error->all(FLERR, "Exact one fix lambda required");

  // get group whose input is ignored from fix lambda
  ignore_group_bit = fix_lambda->group_bit_ignore_lambda_input;

  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLambdaInputAPIP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) { cut[i][j] = mix_distance(cut[i][i], cut[j][j]); }
  return cut[i][j];
}

/**
  * Compute lambda_input and write it to atom->apip_lambda_input.
  * Count the number of computations and measure the compute time for
  * fix atom_weight/apip.
  */

int PairLambdaInputAPIP::calculate_lambda_input()
{
  int i, ii, inum;
  int *ilist;

  inum = list->inum;
  ilist = list->ilist;

  double *lambda_input = atom->apip_lambda_input;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    // "compute" and set lambda input
    lambda_input[i] = 0;
  }
  // return number of calculations
  return inum;
}

/* ----------------------------------------------------------------------
   set return values for timers and counted particles
------------------------------------------------------------------------- */

void PairLambdaInputAPIP::calculate_time_per_atom()
{
  if (n_calculations > 0)
    time_per_atom = timer / n_calculations;
  else
    time_per_atom = -1;

  // reset
  timer = 0;
  n_calculations = 0;
}

/* ---------------------------------------------------------------------- */

void *PairLambdaInputAPIP::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "lambda/input/apip:time_per_atom") == 0) {
    calculate_time_per_atom();
    return (void *) &time_per_atom;
  }
  return nullptr;
}
