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

#include "pair_lambda_input_csp_apip.h"

#include "atom.h"
#include "error.h"
#include "memory.h"
#include "neigh_list.h"

#include <cstring>
#include <utility>
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLambdaInputCSPAPIP::PairLambdaInputCSPAPIP(LAMMPS *lmp) :
    PairLambdaInputAPIP(lmp), distsq(nullptr), nearest(nullptr)
{
  // set defaults
  nnn = 0;
  nnn_buffer = 0;
  maxneigh = 0;
  cut_csp_sq = -1;
}

/* ---------------------------------------------------------------------- */

PairLambdaInputCSPAPIP::~PairLambdaInputCSPAPIP()
{
  if (copymode) return;

  memory->destroy(distsq);
  memory->destroy(nearest);
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLambdaInputCSPAPIP::settings(int narg, char **arg)
{
  double cut_csp = 5;

  // parse arguments
  if (narg < 1) error->all(FLERR, "pair lambda_input/csp: lattice requires one argument");

  if (strcmp(arg[0], "fcc") == 0)
    nnn = 12;
  else if (strcmp(arg[0], "bcc") == 0)
    nnn = 8;
  else
    nnn = utils::inumeric(FLERR, arg[0], false, lmp);

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "pair lambda_input/csp: threshold requires an argument");
      cut_csp = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "N_buffer") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "pair lambda_input/csp: N_buffer requires an argument");
      nnn_buffer = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "pair_lambda_input_csp: unknown argument {}", arg[iarg]);
  }

  if (nnn <= 1 || nnn % 2)
    error->all(FLERR,
               "pair_lambda_input_csp: even number of neighbours > 1 for csp calculation required");
  if (cut_csp <= 0) error->all(FLERR, "pair_lambda_input_csp: cut_csp <= 0");
  if (nnn_buffer < 0) error->all(FLERR, "pair_lambda_input_csp: N_buffer negative");

  cut_global = cut_csp;
  cut_csp_sq = cut_csp * cut_csp;

  // reset cutoffs that have been explicitly set
  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/**
  * Compute CSP and write it to atom->apip_lambda_input.
  * Count the number of computations and measure the compute time for
  * fix atom_weight/apip.
  */

int PairLambdaInputCSPAPIP::calculate_lambda_input()
{
  int i, j, k, ii, jj, kk, n, n_cutoff, inum, jnum;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq, value;
  int *ilist, *jlist, *numneigh, **firstneigh, *mask;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  mask = atom->mask;

  // npairs = number of unique pairs

  int nhalf = nnn / 2;
  int nnn_all = nnn + nnn_buffer;
  int npairs = nnn_all * (nnn_all - 1) / 2;
  auto *pairs = new double[npairs];

  double **x = atom->x;
  double *lambda_input = atom->apip_lambda_input;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    if (mask[i] & ignore_group_bit) {
      // do not calculate the input as it is not used later
      lambda_input[i] = 0;
      continue;
    }

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // ensure distsq and nearest arrays are long enough

    if (jnum > maxneigh) {
      memory->destroy(distsq);
      memory->destroy(nearest);
      maxneigh = jnum;
      memory->create(distsq, maxneigh, "pair lambda_input/csp:distsq");
      memory->create(nearest, maxneigh, "pair lambda_input/csp:nearest");
    }

    // loop over list of all neighbors within force cutoff
    // distsq[] = distance sq to each
    // nearest[] = atom indices of neighbors

    n_cutoff = 0;
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      // do not use cutsq since cutsq may be quite large due to the maximum search of lambda
      if (rsq < cut_csp_sq) {
        distsq[n_cutoff] = rsq;
        nearest[n_cutoff++] = j;
      }
    }

    if (n_cutoff >= nnn_all) {
      // calculate the values of the centro symmetry parameter for this atom

      // store nnn_all nearest neighs in 1st nnn_all locations of distsq and nearest

      select2(nnn_all, n_cutoff, distsq, nearest);

      // R = Ri + Rj for each of npairs i,j pairs among nnn_all neighbors
      // pairs = squared length of each R

      n = 0;
      for (j = 0; j < nnn_all; j++) {
        jj = nearest[j];
        for (k = j + 1; k < nnn_all; k++) {
          kk = nearest[k];
          delx = x[jj][0] + x[kk][0] - 2.0 * xtmp;
          dely = x[jj][1] + x[kk][1] - 2.0 * ytmp;
          delz = x[jj][2] + x[kk][2] - 2.0 * ztmp;
          pairs[n++] = delx * delx + dely * dely + delz * delz;
        }
      }

      // store nhalf smallest pair distances in 1st nhalf locations of pairs

      select(nhalf, npairs, pairs);

      // centrosymmetry = sum of nhalf smallest squared values

      // calculate centro symmetry parameter of this atom
      value = 0.0;
      for (j = 0; j < nhalf; j++) value += pairs[j];

    } else {
      // cannot calculate nnn/2 neighbour pairs
      // -> just set a high value
      value = 1000.0;
    }

    // store cs for this atom
    lambda_input[i] = value;
  }

  delete[] pairs;

  // return number of calculations
  return inum;
}

/* ----------------------------------------------------------------------
   2 select routines from Numerical Recipes (slightly modified)
   find k smallest values in array of length n
   2nd routine sorts auxiliary array at same time
------------------------------------------------------------------------- */

void PairLambdaInputCSPAPIP::select(int k, int n, double *arr)
{
  int i, ir, j, l, mid;
  double a;

  arr--;
  l = 1;
  ir = n;
  for (;;) {
    if (ir <= l + 1) {
      if (ir == l + 1 && arr[ir] < arr[l]) { std::swap(arr[l], arr[ir]); }
      return;
    } else {
      mid = (l + ir) >> 1;
      std::swap(arr[mid], arr[l + 1]);
      if (arr[l] > arr[ir]) { std::swap(arr[l], arr[ir]); }
      if (arr[l + 1] > arr[ir]) { std::swap(arr[l + 1], arr[ir]); }
      if (arr[l] > arr[l + 1]) { std::swap(arr[l], arr[l + 1]); }
      i = l + 1;
      j = ir;
      a = arr[l + 1];
      for (;;) {
        do i++;
        while (arr[i] < a);
        do j--;
        while (arr[j] > a);
        if (j < i) break;
        std::swap(arr[i], arr[j]);
      }
      arr[l + 1] = arr[j];
      arr[j] = a;
      if (j >= k) ir = j - 1;
      if (j <= k) l = i;
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLambdaInputCSPAPIP::select2(int k, int n, double *arr, int *iarr)
{
  int i, ir, j, l, mid, ia;
  double a;

  arr--;
  iarr--;
  l = 1;
  ir = n;
  for (;;) {
    if (ir <= l + 1) {
      if (ir == l + 1 && arr[ir] < arr[l]) {
        std::swap(arr[l], arr[ir]);
        std::swap(iarr[l], iarr[ir]);
      }
      return;
    } else {
      mid = (l + ir) >> 1;
      std::swap(arr[mid], arr[l + 1]);
      std::swap(iarr[mid], iarr[l + 1]);
      if (arr[l] > arr[ir]) {
        std::swap(arr[l], arr[ir]);
        std::swap(iarr[l], iarr[ir]);
      }
      if (arr[l + 1] > arr[ir]) {
        std::swap(arr[l + 1], arr[ir]);
        std::swap(iarr[l + 1], iarr[ir]);
      }
      if (arr[l] > arr[l + 1]) {
        std::swap(arr[l], arr[l + 1]);
        std::swap(iarr[l], iarr[l + 1]);
      }
      i = l + 1;
      j = ir;
      a = arr[l + 1];
      ia = iarr[l + 1];
      for (;;) {
        do i++;
        while (arr[i] < a);
        do j--;
        while (arr[j] > a);
        if (j < i) break;
        std::swap(arr[i], arr[j]);
        std::swap(iarr[i], iarr[j]);
      }
      arr[l + 1] = arr[j];
      arr[j] = a;
      iarr[l + 1] = iarr[j];
      iarr[j] = ia;
      if (j >= k) ir = j - 1;
      if (j <= k) l = i;
    }
  }
}
