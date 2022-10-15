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
#include "electrode_math.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "neigh_list.h"
#include "pair.h"

using namespace LAMMPS_NS;

ElectrodeVector::ElectrodeVector(LAMMPS *lmp, int sensor_group, int source_group, double eta,
                                 bool invert_source) :
    Pointers(lmp)
{
  igroup = sensor_group;                // group of all atoms at which we calculate potential
  this->source_group = source_group;    // group of all atoms influencing potential
  this->invert_source = invert_source;
  groupbit = group->bitmask[igroup];
  ngroup = group->count(igroup);
  source_grpbit = group->bitmask[source_group];
  this->eta = eta;

  kspace_time_total = 0;
  pair_time_total = 0;
  boundary_time_total = 0;
  b_time_total = 0;
}

/* ---------------------------------------------------------------------- */

ElectrodeVector::~ElectrodeVector()
{
  if (timer_flag && (comm->me == 0)) {
    try {
      utils::logmesg(lmp, fmt::format("B time: {:.4g} s\n", b_time_total));
      utils::logmesg(lmp, fmt::format("B kspace time: {:.4g} s\n", kspace_time_total));
      utils::logmesg(lmp, fmt::format("B pair time: {:.4g} s\n", pair_time_total));
      utils::logmesg(lmp, fmt::format("B boundary time: {:.4g} s\n", boundary_time_total));
    } catch (std::exception &) {
    }
  }
}

/* ---------------------------------------------------------------------- */

void ElectrodeVector::setup(class Pair *fix_pair, class NeighList *fix_neighlist, bool timer_flag)
{
  pair = fix_pair;
  cutsq = pair->cutsq;
  list = fix_neighlist;
  this->timer_flag = timer_flag;

  electrode_kspace = dynamic_cast<ElectrodeKSpace *>(force->kspace);
  if (electrode_kspace == nullptr) error->all(FLERR, "KSpace does not implement ElectrodeKSpace");
  g_ewald = force->kspace->g_ewald;
}

/* ---------------------------------------------------------------------- */

void ElectrodeVector::compute_vector(double *vector)
{
  MPI_Barrier(world);
  double start_time = MPI_Wtime();
  // pair
  double pair_start_time = MPI_Wtime();
  pair_contribution(vector);
  MPI_Barrier(world);
  pair_time_total += MPI_Wtime() - pair_start_time;
  // kspace
  double kspace_start_time = MPI_Wtime();
  electrode_kspace->compute_vector(vector, groupbit, source_grpbit, invert_source);
  MPI_Barrier(world);
  kspace_time_total += MPI_Wtime() - kspace_start_time;
  // boundary
  double boundary_start_time = MPI_Wtime();
  electrode_kspace->compute_vector_corr(vector, groupbit, source_grpbit, invert_source);
  MPI_Barrier(world);
  boundary_time_total += MPI_Wtime() - boundary_start_time;
  b_time_total += MPI_Wtime() - start_time;
}

/* ---------------------------------------------------------------------- */

void ElectrodeVector::pair_contribution(double *vector)
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
    bool const i_in_sensor = (mask[i] & groupbit);
    bool const i_in_source = !!(mask[i] & source_grpbit) != invert_source;
    double const xtmp = x[i][0];
    double const ytmp = x[i][1];
    double const ztmp = x[i][2];
    int itype = type[i];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    for (int jj = 0; jj < jnum; jj++) {
      int const j = jlist[jj] & NEIGHMASK;
      bool const j_in_sensor = (mask[j] & groupbit);
      bool const j_in_source = !!(mask[j] & source_grpbit) != invert_source;
      bool const compute_ij = i_in_sensor && j_in_source;
      bool const compute_ji = (newton_pair || j < nlocal) && (j_in_sensor && i_in_source);
      if (!(compute_ij || compute_ji)) continue;
      double const delx = xtmp - x[j][0];    // neighlists take care of pbc
      double const dely = ytmp - x[j][1];
      double const delz = ztmp - x[j][2];
      double const rsq = delx * delx + dely * dely + delz * delz;
      int jtype = type[j];
      if (rsq >= cutsq[itype][jtype]) continue;
      double const r = sqrt(rsq);
      double const rinv = 1.0 / r;
      double aij = rinv;
      aij *= ElectrodeMath::safe_erfc(g_ewald * r);
      aij -= ElectrodeMath::safe_erfc(eta * r) * rinv;
      if (i_in_sensor) {
        vector[i] += aij * q[j];
      } else if (j_in_sensor) {
        vector[j] += aij * q[i];
      }
    }
  }
}
