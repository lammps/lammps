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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meißner (TUHH)
------------------------------------------------------------------------- */

#include "electrode_matrix.h"

#include "atom.h"
#include "comm.h"
#include "electrode_kspace.h"
#include "electrode_math.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "math_const.h"
#include "neigh_list.h"
#include "pair.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ElectrodeMatrix::ElectrodeMatrix(LAMMPS *lmp, int electrode_group, double eta) : Pointers(lmp)
{
  igroup = electrode_group;    // group of all electrode atoms
  groupbit = group->bitmask[igroup];
  ngroup = group->count(igroup);
  this->eta = eta;
  tfflag = false;
}

/* ---------------------------------------------------------------------- */

void ElectrodeMatrix::setup(const std::unordered_map<tagint, int> &tag_ids, class Pair *fix_pair,
                            class NeighList *fix_neighlist)
{
  pair = fix_pair;
  cutsq = pair->cutsq;
  list = fix_neighlist;

  electrode_kspace = dynamic_cast<ElectrodeKSpace *>(force->kspace);
  if (electrode_kspace == nullptr) error->all(FLERR, "KSpace does not implement ElectrodeKSpace");
  g_ewald = force->kspace->g_ewald;

  tag_to_iele = tag_ids;
}

/* ---------------------------------------------------------------------- */

void ElectrodeMatrix::setup_tf(const std::map<int, double> &tf_types)
{
  tfflag = true;
  this->tf_types = tf_types;
}

/* ---------------------------------------------------------------------- */

void ElectrodeMatrix::compute_array(double **array, bool timer_flag)
{
  // setting all entries of coulomb matrix to zero
  size_t nbytes = sizeof(double) * ngroup * ngroup;
  if (nbytes) memset(&array[0][0], 0, nbytes);

  MPI_Barrier(world);
  double kspace_time = MPI_Wtime();
  update_mpos();
  electrode_kspace->compute_matrix(&mpos[0], array, timer_flag);
  MPI_Barrier(world);
  if (timer_flag && (comm->me == 0))
    utils::logmesg(lmp, fmt::format("KSpace time: {:.4g} s\n", MPI_Wtime() - kspace_time));
  pair_contribution(array);
  self_contribution(array);
  electrode_kspace->compute_matrix_corr(&mpos[0], array);
  if (tfflag) tf_contribution(array);

  // reduce coulomb matrix with contributions from all procs
  // all procs need to know full matrix for matrix inversion
  for (int i = 0; i < ngroup; i++) {
    MPI_Allreduce(MPI_IN_PLACE, &array[i][0], ngroup, MPI_DOUBLE, MPI_SUM, world);
  }
}

/* ---------------------------------------------------------------------- */

void ElectrodeMatrix::pair_contribution(double **array)
{
  int inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double r, rinv, rsq, aij;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double **x = atom->x;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  double etaij = eta * eta / sqrt(2.0 * eta * eta);    // see mw ewald theory eq. (29)-(30)

  // neighbor list will be ready because called from post_neighbor
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // skip if I,J are not in 2 groups

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    // skip if atom I is not in either group
    if (!(mask[i] & groupbit)) continue;

    bigint const ipos = mpos[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // real-space part of matrix is symmetric
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      if (!(mask[j] & groupbit)) continue;

      delx = xtmp - x[j][0];    // neighlists take care of pbc
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        rinv = 1.0 / r;
        aij = rinv;
        aij *= ElectrodeMath::safe_erfc(g_ewald * r);
        aij -= ElectrodeMath::safe_erfc(etaij * r) * rinv;
        // newton on or off?
        if (!(newton_pair || j < nlocal)) aij *= 0.5;
        bigint jpos = tag_to_iele[tag[j]];
        array[ipos][jpos] += aij;
        array[jpos][ipos] += aij;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void ElectrodeMatrix::self_contribution(double **array)
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  const double selfint = 2.0 / MY_PIS * g_ewald;
  const double preta = MY_SQRT2 / MY_PIS;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) { array[mpos[i]][mpos[i]] += preta * eta - selfint; }
}

/* ---------------------------------------------------------------------- */

void ElectrodeMatrix::tf_contribution(double **array)
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  int *mask = atom->mask;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) array[mpos[i]][mpos[i]] += tf_types[type[i]];
}

/* ---------------------------------------------------------------------- */

void ElectrodeMatrix::update_mpos()
{
  int const nall = atom->nlocal + atom->nghost;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  mpos = std::vector<bigint>(nall, -1);

  for (int i = 0; i < nall; i++) {
    if (mask[i] & groupbit)
      mpos[i] = tag_to_iele[tag[i]];
    else
      mpos[i] = -1;
  }
}
