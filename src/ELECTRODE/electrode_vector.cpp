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

#include "electrode_vector.h"

#include "atom.h"
#include "comm.h"
#include "electrode_kspace.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "neigh_list.h"
#include "pair.h"

using namespace LAMMPS_NS;

#define EWALD_P 0.3275911
#define A1 0.254829592
#define A2 -0.284496736
#define A3 1.421413741
#define A4 -1.453152027
#define A5 1.061405429

ElectrodeVector::ElectrodeVector(LAMMPS *lmp, int electrode_group, double eta) : Pointers(lmp)
{
  igroup = electrode_group;    // group of all electrode atoms
  groupbit = group->bitmask[igroup];
  ngroup = group->count(igroup);
  vector = new double[ngroup]();    // init to zero
  this->eta = eta;

  setup_time_total = 0;
  reduce_time_total = 0;
  kspace_time_total = 0;
  pair_time_total = 0;
  boundary_time_total = 0;
  b_time_total = 0;
  alloc_time_total = 0;
  mpos_time_total = 0;
}

/* ---------------------------------------------------------------------- */

ElectrodeVector::~ElectrodeVector()
{
  if (comm->me == 0) {
    utils::logmesg(lmp, fmt::format("B time: {}\n", b_time_total));
    utils::logmesg(lmp, fmt::format("B kspace time: {}\n", kspace_time_total));
    utils::logmesg(lmp, fmt::format("B pair time: {}\n", pair_time_total));
    utils::logmesg(lmp, fmt::format("B boundary time: {}\n", boundary_time_total));
    utils::logmesg(lmp, fmt::format("B setup time: {}\n", setup_time_total));
    utils::logmesg(lmp, fmt::format("B reduce time: {}\n", reduce_time_total));
    utils::logmesg(lmp, fmt::format("B alloc time: {}\n", alloc_time_total));
    utils::logmesg(lmp, fmt::format("B mpos time: {}\n", mpos_time_total));
  }
  delete[] vector;
}

/* ---------------------------------------------------------------------- */

void ElectrodeVector::setup(const std::map<tagint, int> &tag_ids, class Pair *fix_pair,
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

void ElectrodeVector::compute_vector()
{
  MPI_Barrier(world);
  double start_time = MPI_Wtime();
  // setup
  double setup_start_time = MPI_Wtime();
  update_mpos();
  for (int i = 0; i < ngroup; i++) vector[i] = 0.;
  MPI_Barrier(world);
  setup_time_total += MPI_Wtime() - setup_start_time;
  // pair
  double pair_start_time = MPI_Wtime();
  pair_contribution();
  MPI_Barrier(world);
  pair_time_total += MPI_Wtime() - pair_start_time;
  // kspace
  double kspace_start_time = MPI_Wtime();
  electrode_kspace->compute_vector(&mpos[0], vector);
  MPI_Barrier(world);
  kspace_time_total += MPI_Wtime() - kspace_start_time;
  // boundary
  double boundary_start_time = MPI_Wtime();
  electrode_kspace->compute_vector_corr(&mpos[0], vector);
  MPI_Barrier(world);
  boundary_time_total += MPI_Wtime() - boundary_start_time;
  // reduce
  double reduce_start_time = MPI_Wtime();
  MPI_Allreduce(MPI_IN_PLACE, vector, ngroup, MPI_DOUBLE, MPI_SUM, world);
  MPI_Barrier(world);
  reduce_time_total += MPI_Wtime() - reduce_start_time;
  b_time_total += MPI_Wtime() - start_time;
}

/* ---------------------------------------------------------------------- */

void ElectrodeVector::pair_contribution()
{
  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  int *mask = atom->mask;
  // neighbor list will be ready because called from post_neighbor
  int const nlocal = atom->nlocal;
  int const inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int newton_pair = force->newton_pair;

  for (int ii = 0; ii < inum; ii++) {
    int const i = ilist[ii];
    bool const i_in_electrode = (mask[i] & groupbit);
    double const xtmp = x[i][0];
    double const ytmp = x[i][1];
    double const ztmp = x[i][2];
    int itype = type[i];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    for (int jj = 0; jj < jnum; jj++) {
      int const j = jlist[jj] & NEIGHMASK;
      bool const j_in_electrode = (mask[j] & groupbit);
      if (i_in_electrode == j_in_electrode) continue;
      double const delx = xtmp - x[j][0];    // neighlists take care of pbc
      double const dely = ytmp - x[j][1];
      double const delz = ztmp - x[j][2];
      double const rsq = delx * delx + dely * dely + delz * delz;
      int jtype = type[j];
      if (rsq >= cutsq[itype][jtype]) continue;
      double const r = sqrt(rsq);
      double const rinv = 1.0 / r;
      double aij = rinv;
      aij *= calc_erfc(g_ewald * r);
      aij -= calc_erfc(eta * r) * rinv;
      if (!(newton_pair || j < nlocal)) aij *= 0.5;
      if (i_in_electrode) {
        vector[mpos[i]] += aij * q[j];
      } else if (j_in_electrode) {
        vector[mpos[j]] += aij * q[i];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void ElectrodeVector::update_mpos()
{
  MPI_Barrier(world);
  double alloc_start = MPI_Wtime();
  int const nall = atom->nlocal + atom->nghost;
  int *tag = atom->tag;
  int *mask = atom->mask;
  mpos = std::vector<bigint>(nall, -1);

  MPI_Barrier(world);
  alloc_time_total += MPI_Wtime() - alloc_start;
  double mpos_start = MPI_Wtime();
  for (int i = 0; i < nall; i++) {
    if (mask[i] & groupbit)
      mpos[i] = tag_to_iele[tag[i]];
    else
      mpos[i] = -1;
  }
  MPI_Barrier(world);
  mpos_time_total += MPI_Wtime() - mpos_start;
}

/* ---------------------------------------------------------------------- */

double ElectrodeVector::calc_erfc(double x)
{
  double expm2 = exp(-x * x);
  double t = 1.0 / (1.0 + EWALD_P * x);
  return t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
}

