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
   Contributing author: Luca Ferraro (CASPUR)
   email: luca.ferraro@caspur.it

   Environment Dependent Interatomic Potential
   References:
    1) J. F. Justo, M. Z. Bazant, E. Kaxiras, V. V. Bulatov, S. Yip
       Phys. Rev. B 58, 2539 (1998)
------------------------------------------------------------------------- */

#include "pair_edip.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>
#include <limits>

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define DELTA 4

#define GRIDDENSITY 8000
#define GRIDSTART 0.1

// max number of interaction per atom for f(Z) environment potential

static constexpr int leadDimInteractionList = 64;

/* ---------------------------------------------------------------------- */

PairEDIP::PairEDIP(LAMMPS *lmp) :
    Pair(lmp), preInvR_ij(nullptr), preExp3B_ij(nullptr), preExp3BDerived_ij(nullptr),
    preExp2B_ij(nullptr), preExp2BDerived_ij(nullptr), prePow2B_ij(nullptr), preForceCoord(nullptr),
    cutoffFunction(nullptr), cutoffFunctionDerived(nullptr), pow2B(nullptr), exp2B(nullptr),
    exp3B(nullptr), qFunctionGrid(nullptr), expMinusBetaZeta_iZeta_iGrid(nullptr),
    tauFunctionGrid(nullptr), tauFunctionDerivedGrid(nullptr)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;

  params = nullptr;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairEDIP::~PairEDIP()
{
  memory->destroy(params);
  memory->destroy(elem3param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    deallocateGrids();
    deallocatePreLoops();
  }
}

/* ---------------------------------------------------------------------- */

void PairEDIP::compute(int eflag, int vflag)
{
  int i, j, k, ii, inum, jnum;
  int itype, jtype, ktype, ijparam, ikparam;
  double xtmp, ytmp, ztmp, evdwl;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int preForceCoord_counter;

  double invR_ij;
  double invR_ik;
  double directorCos_ij_x;
  double directorCos_ij_y;
  double directorCos_ij_z;
  double directorCos_ik_x;
  double directorCos_ik_y;
  double directorCos_ik_z;
  double cosTeta;

  int interpolIDX;
  double interpolTMP;
  double interpolDeltaX;
  double interpolY1;
  double interpolY2;

  double invRMinusCutoffA;
  double sigmaInvRMinusCutoffA;
  double gammInvRMinusCutoffA;
  double cosTetaDiff;
  double cosTetaDiffCosTetaDiff;
  double cutoffFunction_ij;
  double exp2B_ij;
  double exp2BDerived_ij;
  double pow2B_ij;
  double pow2BDerived_ij;
  double exp3B_ij;
  double exp3BDerived_ij;
  double exp3B_ik;
  double exp3BDerived_ik;
  double qFunction;
  double tauFunction;
  double tauFunctionDerived;
  double expMinusBetaZeta_iZeta_i;
  double qFunctionCosTetaDiffCosTetaDiff;
  double expMinusQFunctionCosTetaDiffCosTetaDiff;
  double zeta_i;
  double zeta_iDerived;
  double zeta_iDerivedInvR_ij;

  double forceModCoord_factor;
  double forceModCoord;
  double forceModCoord_ij;
  double forceMod2B;
  double forceMod3B_factor1_ij;
  double forceMod3B_factor2_ij;
  double forceMod3B_factor2;
  double forceMod3B_factor1_ik;
  double forceMod3B_factor2_ik;
  double potentia3B_factor;
  double potential2B_factor;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    zeta_i = 0.0;
    int numForceCoordPairs = 0;

    i = ilist[ii];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    // pre-loop to compute environment coordination f(Z)

    for (int neighbor_j = 0; neighbor_j < jnum; neighbor_j++) {
      j = jlist[neighbor_j];
      j &= NEIGHMASK;

      double dr_ij[3], r_ij;

      dr_ij[0] = xtmp - x[j][0];
      dr_ij[1] = ytmp - x[j][1];
      dr_ij[2] = ztmp - x[j][2];
      r_ij = dr_ij[0] * dr_ij[0] + dr_ij[1] * dr_ij[1] + dr_ij[2] * dr_ij[2];

      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];
      if (r_ij > params[ijparam].cutsq) continue;

      r_ij = sqrt(r_ij);

      invR_ij = 1.0 / r_ij;
      preInvR_ij[neighbor_j] = invR_ij;

      invRMinusCutoffA = 1.0 / (r_ij - cutoffA);
      sigmaInvRMinusCutoffA = sigma * invRMinusCutoffA;
      gammInvRMinusCutoffA = gamm * invRMinusCutoffA;

      interpolDeltaX = r_ij - GRIDSTART;
      interpolTMP = (interpolDeltaX * GRIDDENSITY);
      interpolIDX = (int) interpolTMP;

      interpolY1 = exp3B[interpolIDX];
      interpolY2 = exp3B[interpolIDX + 1];
      exp3B_ij = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

      exp3BDerived_ij = -exp3B_ij * gammInvRMinusCutoffA * invRMinusCutoffA;

      preExp3B_ij[neighbor_j] = exp3B_ij;
      preExp3BDerived_ij[neighbor_j] = exp3BDerived_ij;

      interpolY1 = exp2B[interpolIDX];
      interpolY2 = exp2B[interpolIDX + 1];
      exp2B_ij = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

      exp2BDerived_ij = -exp2B_ij * sigmaInvRMinusCutoffA * invRMinusCutoffA;

      preExp2B_ij[neighbor_j] = exp2B_ij;
      preExp2BDerived_ij[neighbor_j] = exp2BDerived_ij;

      interpolY1 = pow2B[interpolIDX];
      interpolY2 = pow2B[interpolIDX + 1];
      pow2B_ij = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

      prePow2B_ij[neighbor_j] = pow2B_ij;

      // zeta and its derivative

      if (r_ij < cutoffC)
        zeta_i += 1.0;
      else {
        interpolY1 = cutoffFunction[interpolIDX];
        interpolY2 = cutoffFunction[interpolIDX + 1];
        cutoffFunction_ij = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

        zeta_i += cutoffFunction_ij;

        interpolY1 = cutoffFunctionDerived[interpolIDX];
        interpolY2 = cutoffFunctionDerived[interpolIDX + 1];
        zeta_iDerived = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

        zeta_iDerivedInvR_ij = zeta_iDerived * invR_ij;

        preForceCoord_counter = numForceCoordPairs * 5;
        preForceCoord[preForceCoord_counter + 0] = zeta_iDerivedInvR_ij;
        preForceCoord[preForceCoord_counter + 1] = dr_ij[0];
        preForceCoord[preForceCoord_counter + 2] = dr_ij[1];
        preForceCoord[preForceCoord_counter + 3] = dr_ij[2];
        preForceCoord[preForceCoord_counter + 4] = j;
        numForceCoordPairs++;
      }
    }

    // quantities depending on zeta_i

    interpolDeltaX = zeta_i;
    interpolTMP = (interpolDeltaX * GRIDDENSITY);
    interpolIDX = (int) interpolTMP;

    interpolY1 = expMinusBetaZeta_iZeta_iGrid[interpolIDX];
    interpolY2 = expMinusBetaZeta_iZeta_iGrid[interpolIDX + 1];
    expMinusBetaZeta_iZeta_i = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

    interpolY1 = qFunctionGrid[interpolIDX];
    interpolY2 = qFunctionGrid[interpolIDX + 1];
    qFunction = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

    interpolY1 = tauFunctionGrid[interpolIDX];
    interpolY2 = tauFunctionGrid[interpolIDX + 1];
    tauFunction = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

    interpolY1 = tauFunctionDerivedGrid[interpolIDX];
    interpolY2 = tauFunctionDerivedGrid[interpolIDX + 1];
    tauFunctionDerived = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

    forceModCoord_factor = 2.0 * beta * zeta_i * expMinusBetaZeta_iZeta_i;

    forceModCoord = 0.0;

    // two-body interactions, skip half of them

    for (int neighbor_j = 0; neighbor_j < jnum; neighbor_j++) {
      double dr_ij[3], r_ij, f_ij[3];

      j = jlist[neighbor_j];
      j &= NEIGHMASK;

      dr_ij[0] = x[j][0] - xtmp;
      dr_ij[1] = x[j][1] - ytmp;
      dr_ij[2] = x[j][2] - ztmp;
      r_ij = dr_ij[0] * dr_ij[0] + dr_ij[1] * dr_ij[1] + dr_ij[2] * dr_ij[2];

      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];
      if (r_ij > params[ijparam].cutsq) continue;

      r_ij = sqrt(r_ij);

      invR_ij = preInvR_ij[neighbor_j];
      pow2B_ij = prePow2B_ij[neighbor_j];

      potential2B_factor = pow2B_ij - expMinusBetaZeta_iZeta_i;

      exp2B_ij = preExp2B_ij[neighbor_j];

      pow2BDerived_ij = -rho * invR_ij * pow2B_ij;

      forceModCoord += (forceModCoord_factor * exp2B_ij);

      exp2BDerived_ij = preExp2BDerived_ij[neighbor_j];
      forceMod2B = exp2BDerived_ij * potential2B_factor + exp2B_ij * pow2BDerived_ij;

      directorCos_ij_x = invR_ij * dr_ij[0];
      directorCos_ij_y = invR_ij * dr_ij[1];
      directorCos_ij_z = invR_ij * dr_ij[2];

      exp3B_ij = preExp3B_ij[neighbor_j];
      exp3BDerived_ij = preExp3BDerived_ij[neighbor_j];

      f_ij[0] = forceMod2B * directorCos_ij_x;
      f_ij[1] = forceMod2B * directorCos_ij_y;
      f_ij[2] = forceMod2B * directorCos_ij_z;

      f[i][0] += f_ij[0];
      f[i][1] += f_ij[1];
      f[i][2] += f_ij[2];

      f[j][0] -= f_ij[0];
      f[j][1] -= f_ij[1];
      f[j][2] -= f_ij[2];

      // potential energy

      evdwl = (exp2B_ij * potential2B_factor);

      if (evflag)
        ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, -forceMod2B * invR_ij, dr_ij[0], dr_ij[1],
                 dr_ij[2]);

      // three-body Forces

      for (int neighbor_k = neighbor_j + 1; neighbor_k < jnum; neighbor_k++) {
        double dr_ik[3], r_ik, f_ik[3];

        k = jlist[neighbor_k];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        ikparam = elem3param[itype][ktype][ktype];

        dr_ik[0] = x[k][0] - xtmp;
        dr_ik[1] = x[k][1] - ytmp;
        dr_ik[2] = x[k][2] - ztmp;
        r_ik = dr_ik[0] * dr_ik[0] + dr_ik[1] * dr_ik[1] + dr_ik[2] * dr_ik[2];

        if (r_ik > params[ikparam].cutsq) continue;

        r_ik = sqrt(r_ik);

        invR_ik = preInvR_ij[neighbor_k];

        directorCos_ik_x = invR_ik * dr_ik[0];
        directorCos_ik_y = invR_ik * dr_ik[1];
        directorCos_ik_z = invR_ik * dr_ik[2];

        cosTeta = directorCos_ij_x * directorCos_ik_x + directorCos_ij_y * directorCos_ik_y +
            directorCos_ij_z * directorCos_ik_z;

        cosTetaDiff = cosTeta + tauFunction;
        cosTetaDiffCosTetaDiff = cosTetaDiff * cosTetaDiff;
        qFunctionCosTetaDiffCosTetaDiff = cosTetaDiffCosTetaDiff * qFunction;
        expMinusQFunctionCosTetaDiffCosTetaDiff = exp(-qFunctionCosTetaDiffCosTetaDiff);

        potentia3B_factor = lambda *
            ((1.0 - expMinusQFunctionCosTetaDiffCosTetaDiff) +
             eta * qFunctionCosTetaDiffCosTetaDiff);

        exp3B_ik = preExp3B_ij[neighbor_k];
        exp3BDerived_ik = preExp3BDerived_ij[neighbor_k];

        forceMod3B_factor1_ij = -exp3BDerived_ij * exp3B_ik * potentia3B_factor;
        forceMod3B_factor2 = 2.0 * lambda * exp3B_ij * exp3B_ik * qFunction * cosTetaDiff *
            (eta + expMinusQFunctionCosTetaDiffCosTetaDiff);
        forceMod3B_factor2_ij = forceMod3B_factor2 * invR_ij;

        f_ij[0] = forceMod3B_factor1_ij * directorCos_ij_x +
            forceMod3B_factor2_ij * (cosTeta * directorCos_ij_x - directorCos_ik_x);
        f_ij[1] = forceMod3B_factor1_ij * directorCos_ij_y +
            forceMod3B_factor2_ij * (cosTeta * directorCos_ij_y - directorCos_ik_y);
        f_ij[2] = forceMod3B_factor1_ij * directorCos_ij_z +
            forceMod3B_factor2_ij * (cosTeta * directorCos_ij_z - directorCos_ik_z);

        forceMod3B_factor1_ik = -exp3BDerived_ik * exp3B_ij * potentia3B_factor;
        forceMod3B_factor2_ik = forceMod3B_factor2 * invR_ik;

        f_ik[0] = forceMod3B_factor1_ik * directorCos_ik_x +
            forceMod3B_factor2_ik * (cosTeta * directorCos_ik_x - directorCos_ij_x);
        f_ik[1] = forceMod3B_factor1_ik * directorCos_ik_y +
            forceMod3B_factor2_ik * (cosTeta * directorCos_ik_y - directorCos_ij_y);
        f_ik[2] = forceMod3B_factor1_ik * directorCos_ik_z +
            forceMod3B_factor2_ik * (cosTeta * directorCos_ik_z - directorCos_ij_z);

        forceModCoord += (forceMod3B_factor2 * (tauFunctionDerived - 0.5 * mu * cosTetaDiff));

        f[j][0] += f_ij[0];
        f[j][1] += f_ij[1];
        f[j][2] += f_ij[2];

        f[k][0] += f_ik[0];
        f[k][1] += f_ik[1];
        f[k][2] += f_ik[2];

        f[i][0] -= f_ij[0] + f_ik[0];
        f[i][1] -= f_ij[1] + f_ik[1];
        f[i][2] -= f_ij[2] + f_ik[2];

        // potential energy

        evdwl = (exp3B_ij * exp3B_ik * potentia3B_factor);

        if (evflag) ev_tally3(i, j, k, evdwl, 0.0, f_ij, f_ik, dr_ij, dr_ik);
      }
    }

    // forces due to environment coordination f(Z)

    for (int idx = 0; idx < numForceCoordPairs; idx++) {
      double dr_ij[3], f_ij[3];

      preForceCoord_counter = idx * 5;
      zeta_iDerivedInvR_ij = preForceCoord[preForceCoord_counter + 0];
      dr_ij[0] = preForceCoord[preForceCoord_counter + 1];
      dr_ij[1] = preForceCoord[preForceCoord_counter + 2];
      dr_ij[2] = preForceCoord[preForceCoord_counter + 3];
      j = static_cast<int>(preForceCoord[preForceCoord_counter + 4]);

      forceModCoord_ij = forceModCoord * zeta_iDerivedInvR_ij;

      f_ij[0] = forceModCoord_ij * dr_ij[0];
      f_ij[1] = forceModCoord_ij * dr_ij[1];
      f_ij[2] = forceModCoord_ij * dr_ij[2];

      f[i][0] -= f_ij[0];
      f[i][1] -= f_ij[1];
      f[i][2] -= f_ij[2];

      f[j][0] += f_ij[0];
      f[j][1] += f_ij[1];
      f[j][2] += f_ij[2];

      if (evflag)
        ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, -forceModCoord_ij, dr_ij[0], dr_ij[1],
                 dr_ij[2]);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairEDIP::allocateGrids(void)
{
  int numGridPointsOneCutoffFunction;
  int numGridPointsNotOneCutoffFunction;
  int numGridPointsCutoffFunction;
  int numGridPointsR;
  int numGridPointsRTotal;
  int numGridPointsQFunctionGrid;
  int numGridPointsExpMinusBetaZeta_iZeta_i;
  int numGridPointsTauFunctionGrid;
  double maxArgumentTauFunctionGrid;
  double maxArgumentQFunctionGrid;
  double maxArgumentExpMinusBetaZeta_iZeta_i;
  double const leftLimitToZero = -std::numeric_limits<double>::min() * 1000.0;

  deallocateGrids();

  // tauFunctionGrid

  maxArgumentTauFunctionGrid = leadDimInteractionList;
  numGridPointsTauFunctionGrid = (int) ((maxArgumentTauFunctionGrid) *GRIDDENSITY) + 2;

  memory->create(tauFunctionGrid, numGridPointsTauFunctionGrid, "edip:tauFunctionGrid");
  memory->create(tauFunctionDerivedGrid, numGridPointsTauFunctionGrid,
                 "edip:tauFunctionDerivedGrid");

  // expMinusBetaZeta_iZeta_iGrid

  maxArgumentExpMinusBetaZeta_iZeta_i = leadDimInteractionList;
  numGridPointsExpMinusBetaZeta_iZeta_i =
      (int) ((maxArgumentExpMinusBetaZeta_iZeta_i) *GRIDDENSITY) + 2;
  memory->create(expMinusBetaZeta_iZeta_iGrid, numGridPointsExpMinusBetaZeta_iZeta_i,
                 "edip:expMinusBetaZeta_iZeta_iGrid");

  // qFunctionGrid

  maxArgumentQFunctionGrid = leadDimInteractionList;
  numGridPointsQFunctionGrid = (int) ((maxArgumentQFunctionGrid) *GRIDDENSITY) + 2;
  memory->create(qFunctionGrid, numGridPointsQFunctionGrid, "edip:qFunctionGrid");

  // cutoffFunction

  numGridPointsOneCutoffFunction = (int) ((cutoffC - GRIDSTART) * GRIDDENSITY);
  numGridPointsNotOneCutoffFunction = (int) ((cutoffA - cutoffC) * GRIDDENSITY);
  numGridPointsCutoffFunction =
      numGridPointsOneCutoffFunction + numGridPointsNotOneCutoffFunction + 2;

  memory->create(cutoffFunction, numGridPointsCutoffFunction, "edip:cutoffFunction");
  memory->create(cutoffFunctionDerived, numGridPointsCutoffFunction, "edip:cutoffFunctionDerived");

  // pow2B

  numGridPointsR = (int) ((cutoffA + leftLimitToZero - GRIDSTART) * GRIDDENSITY);
  numGridPointsRTotal = numGridPointsR + 2;

  memory->create(pow2B, numGridPointsRTotal, "edip:pow2B");
  memory->create(exp2B, numGridPointsRTotal, "edip:exp2B");
  memory->create(exp3B, numGridPointsRTotal, "edip:exp3B");
}

/* ----------------------------------------------------------------------
   pre-calculated structures
------------------------------------------------------------------------- */

void PairEDIP::allocatePreLoops(void)
{
  int nthreads = comm->nthreads;

  deallocatePreLoops();
  memory->create(preInvR_ij, nthreads * leadDimInteractionList, "edip:preInvR_ij");
  memory->create(preExp3B_ij, nthreads * leadDimInteractionList, "edip:preExp3B_ij");
  memory->create(preExp3BDerived_ij, nthreads * leadDimInteractionList, "edip:preExp3BDerived_ij");
  memory->create(preExp2B_ij, nthreads * leadDimInteractionList, "edip:preExp2B_ij");
  memory->create(preExp2BDerived_ij, nthreads * leadDimInteractionList, "edip:preExp2BDerived_ij");
  memory->create(prePow2B_ij, nthreads * leadDimInteractionList, "edip:prePow2B_ij");
  memory->create(preForceCoord, 5 * nthreads * leadDimInteractionList, "edip:preForceCoord");
}

/* ----------------------------------------------------------------------
   deallocate grids
------------------------------------------------------------------------- */

void PairEDIP::deallocateGrids(void)
{
  memory->destroy(cutoffFunction);
  memory->destroy(cutoffFunctionDerived);
  memory->destroy(pow2B);
  memory->destroy(exp2B);
  memory->destroy(exp3B);
  memory->destroy(qFunctionGrid);
  memory->destroy(expMinusBetaZeta_iZeta_iGrid);
  memory->destroy(tauFunctionGrid);
  memory->destroy(tauFunctionDerivedGrid);
}

/* ----------------------------------------------------------------------
   deallocate preLoops
------------------------------------------------------------------------- */

void PairEDIP::deallocatePreLoops(void)
{
  memory->destroy(preInvR_ij);
  memory->destroy(preExp3B_ij);
  memory->destroy(preExp3BDerived_ij);
  memory->destroy(preExp2B_ij);
  memory->destroy(preExp2BDerived_ij);
  memory->destroy(prePow2B_ij);
  memory->destroy(preForceCoord);
}

/* ---------------------------------------------------------------------- */

void PairEDIP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  map = new int[n + 1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairEDIP::settings(int narg, char ** /*arg*/)
{
  if (narg != 0) error->all(FLERR, "Illegal pair_style command");
}

/* ---------------------------------------------------------------------- */

void PairEDIP::initGrids(void)
{
  int l;
  int numGridPointsOneCutoffFunction;
  int numGridPointsNotOneCutoffFunction;
  int numGridPointsCutoffFunction;
  int numGridPointsR;
  int numGridPointsQFunctionGrid;
  int numGridPointsExpMinusBetaZeta_iZeta_i;
  int numGridPointsTauFunctionGrid;
  double maxArgumentTauFunctionGrid;
  double maxArgumentQFunctionGrid;
  double maxArgumentExpMinusBetaZeta_iZeta_i;
  double r;
  double temp;
  double temp3;
  double temp4;
  double deltaArgumentR;
  double deltaArgumentCutoffFunction;
  double deltaArgumentQFunctionGrid;
  double deltaArgumentTauFunctionGrid;
  double deltaArgumentExpMinusBetaZeta_iZeta_i;
  double const leftLimitToZero = -std::numeric_limits<double>::min() * 1000.0;

  // tauFunctionGrid

  maxArgumentTauFunctionGrid = leadDimInteractionList;

  numGridPointsTauFunctionGrid = (int) ((maxArgumentTauFunctionGrid) *GRIDDENSITY) + 2;

  r = 0.0;
  deltaArgumentTauFunctionGrid = 1.0 / GRIDDENSITY;

  for (l = 0; l < numGridPointsTauFunctionGrid; l++) {
    tauFunctionGrid[l] = u1 + u2 * u3 * exp(-u4 * r) - u2 * exp(-2.0 * u4 * r);
    tauFunctionDerivedGrid[l] = -u2 * u3 * u4 * exp(-u4 * r) + 2.0 * u2 * u4 * exp(-2.0 * u4 * r);
    r += deltaArgumentTauFunctionGrid;
  }

  // expMinusBetaZeta_iZeta_iGrid

  maxArgumentExpMinusBetaZeta_iZeta_i = leadDimInteractionList;

  numGridPointsExpMinusBetaZeta_iZeta_i =
      (int) ((maxArgumentExpMinusBetaZeta_iZeta_i) *GRIDDENSITY) + 2;

  r = 0.0;
  deltaArgumentExpMinusBetaZeta_iZeta_i = 1.0 / GRIDDENSITY;

  for (l = 0; l < numGridPointsExpMinusBetaZeta_iZeta_i; l++) {
    expMinusBetaZeta_iZeta_iGrid[l] = exp(-beta * r * r);
    r += deltaArgumentExpMinusBetaZeta_iZeta_i;
  }

  // qFunctionGrid

  maxArgumentQFunctionGrid = leadDimInteractionList;
  numGridPointsQFunctionGrid = (int) ((maxArgumentQFunctionGrid) *GRIDDENSITY) + 2;

  r = 0.0;
  deltaArgumentQFunctionGrid = 1.0 / GRIDDENSITY;

  for (l = 0; l < numGridPointsQFunctionGrid; l++) {
    qFunctionGrid[l] = Q0 * exp(-mu * r);
    r += deltaArgumentQFunctionGrid;
  }

  // cutoffFunction

  numGridPointsOneCutoffFunction = (int) ((cutoffC - GRIDSTART) * GRIDDENSITY);
  numGridPointsNotOneCutoffFunction = (int) ((cutoffA - cutoffC) * GRIDDENSITY);
  numGridPointsCutoffFunction =
      numGridPointsOneCutoffFunction + numGridPointsNotOneCutoffFunction + 2;

  r = GRIDSTART;
  deltaArgumentCutoffFunction = 1.0 / GRIDDENSITY;

  for (l = 0; l < numGridPointsOneCutoffFunction; l++) {
    cutoffFunction[l] = 1.0;
    cutoffFunctionDerived[l] = 0.0;
    r += deltaArgumentCutoffFunction;
  }

  for (l = numGridPointsOneCutoffFunction; l < numGridPointsCutoffFunction; l++) {
    temp = (cutoffA - cutoffC) / (r - cutoffC);
    temp3 = temp * temp * temp;
    temp4 = temp3 * temp;
    cutoffFunction[l] = exp(alpha / (1.0 - temp3));
    cutoffFunctionDerived[l] = (-3 * alpha / (cutoffA - cutoffC)) *
        (temp4 / ((1 - temp3) * (1 - temp3))) * exp(alpha / (1.0 - temp3));
    r += deltaArgumentCutoffFunction;
  }

  // pow2B

  numGridPointsR = (int) ((cutoffA + leftLimitToZero - GRIDSTART) * GRIDDENSITY);

  r = GRIDSTART;
  deltaArgumentR = 1.0 / GRIDDENSITY;
  for (l = 0; l < numGridPointsR; l++) {
    pow2B[l] = pow((B / r), rho);
    exp2B[l] = A * exp(sigma / (r - cutoffA));
    exp3B[l] = exp(gamm / (r - cutoffA));
    r += deltaArgumentR;
  }

  pow2B[numGridPointsR] = pow((B / r), rho);
  exp2B[numGridPointsR] = 0;
  exp3B[numGridPointsR] = 0;
  r += deltaArgumentR;
  pow2B[numGridPointsR + 1] = pow((B / r), rho);
  exp2B[numGridPointsR + 1] = 0;
  exp3B[numGridPointsR + 1] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairEDIP::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  map_element2type(narg - 3, arg + 3);
  if (nelements != 1) error->all(FLERR, "Pair style edip only supports single element potentials");

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup_params();

  // allocate tables and internal structures

  allocatePreLoops();
  allocateGrids();
  initGrids();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEDIP::init_style()
{
  if (force->newton_pair == 0) error->all(FLERR, "Pair style edip requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairEDIP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairEDIP::read_file(char *file)
{
  int params_per_line = 20;
  char **words = new char *[params_per_line + 1];

  memory->sfree(params);
  params = nullptr;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = utils::open_potential(file, lmp, nullptr);
    if (fp == nullptr)
      error->one(FLERR,"Cannot open EDIP potential file {}: {}", file,utils::getsyserror());
  }

  // read each set of params from potential file
  // one set of params can span multiple lines
  // store params if all 3 element tags are in element list

  int n, nwords, ielement, jelement, kelement;
  char line[MAXLINE], *ptr;
  int eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line, MAXLINE, fp);
      if (ptr == nullptr) {
        eof = 1;
        fclose(fp);
      } else
        n = strlen(line) + 1;
    }
    MPI_Bcast(&eof, 1, MPI_INT, 0, world);
    if (eof) break;
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    MPI_Bcast(line, n, MPI_CHAR, 0, world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line, '#'))) *ptr = '\0';
    nwords = utils::count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n], MAXLINE - n, fp);
        if (ptr == nullptr) {
          eof = 1;
          fclose(fp);
        } else
          n = strlen(line) + 1;
      }
      MPI_Bcast(&eof, 1, MPI_INT, 0, world);
      if (eof) break;
      MPI_Bcast(&n, 1, MPI_INT, 0, world);
      MPI_Bcast(line, n, MPI_CHAR, 0, world);
      if ((ptr = strchr(line, '#'))) *ptr = '\0';
      nwords = utils::count_words(line);
    }

    if (nwords != params_per_line) error->all(FLERR, "Incorrect format in EDIP potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line, " \t\n\r\f");
    while ((words[nwords++] = strtok(nullptr, " \t\n\r\f"))) continue;

    // ielement,jelement,kelement = 1st args
    // if all 3 args are in element list, then parse this line
    // else skip to next entry in file

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0], elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1], elements[jelement]) == 0) break;
    if (jelement == nelements) continue;
    for (kelement = 0; kelement < nelements; kelement++)
      if (strcmp(words[2], elements[kelement]) == 0) break;
    if (kelement == nelements) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      maxparam += DELTA;
      params = (Param *) memory->srealloc(params, maxparam * sizeof(Param), "pair:params");

      // make certain all addional allocated storage is initialized
      // to avoid false positives when checking with valgrind

      memset(params + nparams, 0, DELTA * sizeof(Param));
    }

    params[nparams].ielement = ielement;
    params[nparams].jelement = jelement;
    params[nparams].kelement = kelement;
    params[nparams].A = atof(words[3]);
    params[nparams].B = atof(words[4]);
    params[nparams].cutoffA = atof(words[5]);
    params[nparams].cutoffC = atof(words[6]);
    params[nparams].alpha = atof(words[7]);
    params[nparams].beta = atof(words[8]);
    params[nparams].eta = atof(words[9]);
    params[nparams].gamm = atof(words[10]);
    params[nparams].lambda = atof(words[11]);
    params[nparams].mu = atof(words[12]);
    params[nparams].rho = atof(words[13]);
    params[nparams].sigma = atof(words[14]);
    params[nparams].Q0 = atof(words[15]);
    params[nparams].u1 = atof(words[16]);
    params[nparams].u2 = atof(words[17]);
    params[nparams].u3 = atof(words[18]);
    params[nparams].u4 = atof(words[19]);

    if (params[nparams].A < 0.0 || params[nparams].B < 0.0 || params[nparams].cutoffA < 0.0 ||
        params[nparams].cutoffC < 0.0 || params[nparams].alpha < 0.0 ||
        params[nparams].beta < 0.0 || params[nparams].eta < 0.0 || params[nparams].gamm < 0.0 ||
        params[nparams].lambda < 0.0 || params[nparams].mu < 0.0 || params[nparams].rho < 0.0 ||
        params[nparams].sigma < 0.0)
      error->all(FLERR, "Illegal EDIP parameter");

    nparams++;
  }

  delete[] words;
}

/* ---------------------------------------------------------------------- */

void PairEDIP::setup_params()
{
  int i, j, k, m, n;
  double rtmp;

  // set elem3param for all triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem3param);
  memory->create(elem3param, nelements, nelements, nelements, "pair:elem3param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement && k == params[m].kelement) {
            if (n >= 0) error->all(FLERR, "Potential file has duplicate entry");
            n = m;
          }
        }
        if (n < 0) error->all(FLERR, "Potential file is missing an entry");
        elem3param[i][j][k] = n;
      }

  // set cutoff square

  for (m = 0; m < nparams; m++) { params[m].cutsq = params[m].cutoffA * params[m].cutoffA; }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++) {
    rtmp = sqrt(params[m].cutsq);
    if (rtmp > cutmax) cutmax = rtmp;
  }

  // this should be removed for multi species parameterization

  A = params[0].A;
  B = params[0].B;
  rho = params[0].rho;
  cutoffA = params[0].cutoffA;
  cutoffC = params[0].cutoffC;
  sigma = params[0].sigma;
  lambda = params[0].lambda;
  gamm = params[0].gamm;
  eta = params[0].eta;
  Q0 = params[0].Q0;
  mu = params[0].mu;
  beta = params[0].beta;
  alpha = params[0].alpha;
  u1 = params[0].u1;
  u2 = params[0].u2;
  u3 = params[0].u3;
  u4 = params[0].u4;
}
