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

#include "pair_lambda_zone_apip.h"

#include "atom.h"
#include "error.h"
#include "fix.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLambdaZoneAPIP::PairLambdaZoneAPIP(LAMMPS *lmp) : Pair(lmp), cut(nullptr), lambda_ta(nullptr)
{
  // set defaults

  cut_global = cut_lo = cut_hi = cut_hi_sq = cut_width = lambda_non_group = -1;
  groupbit = -1;
  nmax_ta = 0;

  timer = 0;
  n_calculations = 0;
  time_per_atom = -1;
}

/* ---------------------------------------------------------------------- */

PairLambdaZoneAPIP::~PairLambdaZoneAPIP()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
  }

  if (nmax_ta > 0) memory->destroy(lambda_ta);
}

/* ---------------------------------------------------------------------- */
void PairLambdaZoneAPIP::coeff(int narg, char **arg)
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

void PairLambdaZoneAPIP::allocate()
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

void PairLambdaZoneAPIP::compute(int eflag, int vflag)
{
  // basic stuff (see pair_zero)
  ev_init(eflag, vflag);
  if (vflag_fdotr) virial_fdotr_compute();

  calculate_lambda();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLambdaZoneAPIP::settings(int narg, char **arg)
{
  // parse arguments
  if (narg != 1) error->all(FLERR, "pair_lambda_zone: expected 1 instead of {} arguments", narg);

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);

  if (cut_global <= 0) error->all(FLERR, "pair_lambda_zone: cut_global = {} <= 0", cut_global);

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

void PairLambdaZoneAPIP::init_style()
{
  if (!atom->apip_lambda_input_ta_flag)
    error->all(FLERR, "pair_lambda_zone requires an atom style with lambda_input_ta");

  // find fix lambda/apip
  class Fix *fix_lambda = nullptr;
  int count = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style, "lambda/apip") == 0) {
      fix_lambda = modify->fix[i];
      count++;
    }
  }
  if (count != 1) error->all(FLERR, "Exact one fix lambda required");

  int dim = 0;
  cut_lo = *((double *) fix_lambda->extract("fix_lambda:cut_lo", dim));
  cut_hi = *((double *) fix_lambda->extract("fix_lambda:cut_hi", dim));
  cut_hi_sq = *((double *) fix_lambda->extract("fix_lambda:cut_hi_sq", dim));
  cut_width = *((double *) fix_lambda->extract("fix_lambda:cut_width", dim));
  lambda_non_group = *((double *) fix_lambda->extract("fix_lambda:lambda_non_group", dim));
  groupbit = fix_lambda->groupbit;

  if (cut_hi > cut_global) error->all(FLERR, "The r_lambda_hi > neighbour_list_cutoff.");

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_GHOST);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLambdaZoneAPIP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) { cut[i][j] = mix_distance(cut[i][i], cut[j][j]); }
  return cut[i][j];
}

/**
  * calculate new lambda with lambda_input_ta of own and ghost atoms.
  * Search local maximum of neighbouring atoms for each atom.
  * input: lambda_input_ta of own + ghost particles
  * output: lambda_input_ta of own particles
  */

void PairLambdaZoneAPIP::calculate_lambda()
{
  double timer_start = platform::walltime();

  double xtmp, ytmp, ztmp, delx, dely, delz, lambda_tmp, rsq, lambda_new;
  double **x, *lambda_input_ta;
  int allnum, ii, i, nlocal, jnum, j, jj, loop_iterations;
  int *ilist, *mask, *jlist, *numneigh, **firstneigh;

  mask = atom->mask;
  x = atom->x;
  lambda_input_ta = atom->apip_lambda_input_ta;

  nlocal = atom->nlocal;

  if (nlocal > nmax_ta) {
    memory->destroy(lambda_ta);
    nmax_ta = nlocal;
    memory->create(lambda_ta, nmax_ta, "pair/lambda:lambda_ta");
  }

  // 1 set lambda for own particles
  for (i = 0; i < nlocal; i++) {
    lambda_ta[i] = (mask[i] & groupbit) ? lambda_input_ta[i] : lambda_non_group;
  }

  // 2 loop over all atoms with non-simple lambda
  loop_iterations = 0;

  if (cut_hi_sq > 0) {
    ilist = list->ilist;
    allnum = list->inum + list->gnum;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    for (ii = 0; ii < allnum; ii++) {
      i = ilist[ii];

      // skip simple atoms and non-group atoms
      // which do not influence the lambda_ta of neighbouring atoms
      if (lambda_input_ta[i] == 1 || (!(mask[i] & groupbit))) { continue; }

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      lambda_tmp = 1 - lambda_input_ta[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      loop_iterations++;

      // 3 loop over neighbours to set their lambda_ta
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        // it is not required to set lambda_ta of ghosts or non-group atoms
        if (j >= nlocal || (!(mask[j] & groupbit))) { continue; }

        // the neighbour j is already complex -> skip
        if (lambda_ta[j] == 0) { continue; }

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;

        if (rsq >= cut_hi_sq) { continue; }

        lambda_new = 1 - lambda_tmp * switching_function_poly_distance(sqrt(rsq));

        // more complex lambda_ta found ? set lambda_ta for ngh
        if (lambda_new < lambda_ta[j]) lambda_ta[j] = lambda_new;
      }
    }
  }

  // copy calculated lambda max back to lambda_ta
  for (i = 0; i < nlocal; i++) lambda_input_ta[i] = lambda_ta[i];

  timer += platform::walltime() - timer_start;
  n_calculations += loop_iterations;
}

// helper function
// similar to cutoff_func_poly in ace_radial.cpp
// compare Phys Rev Mat 6, 013804 (2022) APPENDIX C: RADIAL AND CUTOFF FUNCTIONS 2. Cutoff function
// the first two derivatives of the switching function lambda vanishes at the boundaries of the switching region
double PairLambdaZoneAPIP::switching_function_poly_distance(double input)
{
  // calculate lambda
  if (input <= cut_lo) {
    return 1;
  } else if (input >= cut_hi) {
    return 0;
  } else {
    double deltatmp = 1 - 2 * (1 + (input - cut_hi) / (cut_width));
    return 0.5 + 7.5 / 2. * (deltatmp / 4. - pow(deltatmp, 3) / 6. + pow(deltatmp, 5) / 20.);
  }
}

/* ----------------------------------------------------------------------
   set return values for timers and counted particles
------------------------------------------------------------------------- */

void PairLambdaZoneAPIP::calculate_time_per_atom()
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

void *PairLambdaZoneAPIP::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "lambda/zone/apip:time_per_atom") == 0) {
    calculate_time_per_atom();
    return (void *) &time_per_atom;
  }
  return nullptr;
}
