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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meissner (TUHH)
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
#include "math_const.h"
#include "neigh_list.h"
#include "pair.h"

#include <cassert>
#include <cmath>
#include <exception>

using namespace LAMMPS_NS;
using namespace MathConst;

ElectrodeVector::ElectrodeVector(LAMMPS *lmp, int sensor_group, int source_group, double eta,
                                 bool invert_source) :
    Pointers(lmp), cutsq(nullptr), pair(nullptr), list(nullptr), electrode_kspace(nullptr)
{
  igroup = sensor_group;                // group of all atoms at which we calculate potential
  this->source_group = source_group;    // group of all atoms influencing potential
  this->invert_source = invert_source;
  groupbit = group->bitmask[igroup];
  ngroup = group->count(igroup);
  source_grpbit = group->bitmask[source_group];
  this->eta = eta;
  tfflag = false;
  etaflag = false;

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
      utils::logmesg(lmp, "B time: {:.4g} s\n", b_time_total);
      utils::logmesg(lmp, "B kspace time: {:.4g} s\n", kspace_time_total);
      utils::logmesg(lmp, "B pair time: {:.4g} s\n", pair_time_total);
      utils::logmesg(lmp, "B boundary time: {:.4g} s\n", boundary_time_total);
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

void ElectrodeVector::setup_tf(const std::map<int, double> &tf_types)
{
  tfflag = true;
  this->tf_types = tf_types;
}

/* ---------------------------------------------------------------------- */

void ElectrodeVector::setup_eta(int index)
{
  etaflag = true;
  eta_index = index;
}

/* ---------------------------------------------------------------------- */

void ElectrodeVector::compute_vector(double *vector)
{
  MPI_Barrier(world);
  double start_time = MPI_Wtime();
  // pair
  double pair_start_time = MPI_Wtime();
  pair_contribution(vector);
  self_contribution(vector);
  if (tfflag) tf_contribution(vector);
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
  const int nlocal = atom->nlocal;
  const int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int newton_pair = force->newton_pair;

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    bool const i_in_sensor = (mask[i] & groupbit);
    bool const i_in_source = !!(mask[i] & source_grpbit) != invert_source;
    if (!(i_in_sensor || i_in_source)) continue;
    double xi[3];
    charge_position(i, xi);
    double const eta_i = etaflag ? atom->dvector[eta_index][i] : eta;
    int itype = type[i];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj] & NEIGHMASK;
      bool const j_in_sensor = (mask[j] & groupbit);
      bool const j_in_source = !!(mask[j] & source_grpbit) != invert_source;
      bool const compute_ij = i_in_sensor && j_in_source;
      bool const compute_ji = (newton_pair || j < nlocal) && (j_in_sensor && i_in_source);
      if (!(compute_ij || compute_ji)) continue;
      double xj[3];
      charge_position(j, xj);
      double const delx = xi[0] - xj[0];
      double const dely = xi[1] - xj[1];
      double const delz = xi[2] - xj[2];
      double const rsq = delx * delx + dely * dely + delz * delz;
      int jtype = type[j];
      if (rsq >= pair_cutsq(itype, jtype)) continue;
      double const eta_j = etaflag ? atom->dvector[eta_index][j] : eta;
      double etaij;
      if (i_in_sensor && j_in_sensor) {
        etaij = eta_i * eta_j / sqrt(eta_i * eta_i + eta_j * eta_j);
      } else if (i_in_sensor) {
        etaij = eta_i;
      } else {
        assert(j_in_sensor);
        etaij = eta_j;
      }
      double const r = sqrt(rsq);
      double const rinv = 1.0 / r;
      double aij = rinv;
      aij *= ElectrodeMath::safe_erfc(g_ewald * r);
      aij -= ElectrodeMath::safe_erfc(etaij * r) * rinv;
      if (i_in_sensor) { vector[i] += aij * q[j]; }
      if (j_in_sensor && (!invert_source || !i_in_sensor)) { vector[j] += aij * q[i]; }
    }
  }
}

/* ---------------------------------------------------------------------- */

void ElectrodeVector::self_contribution(double *vector)
{
  const int inum = list->inum;
  int *mask = atom->mask;
  int *ilist = list->ilist;
  double *q = atom->q;

  const double selfint = 2.0 / MY_PIS * g_ewald;
  const double preta = MY_SQRT2 / MY_PIS;

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    double const eta_i = etaflag ? atom->dvector[eta_index][i] : eta;
    bool const i_in_sensor = (mask[i] & groupbit);
    bool const i_in_source = !!(mask[i] & source_grpbit) != invert_source;
    if (i_in_sensor && i_in_source) vector[i] += (preta * eta_i - selfint) * q[i];
  }
}

/* ---------------------------------------------------------------------- */

void ElectrodeVector::tf_contribution(double *vector)
{
  const int inum = list->inum;
  int *mask = atom->mask;
  int *type = atom->type;
  int *ilist = list->ilist;
  double *q = atom->q;

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    bool const i_in_sensor = (mask[i] & groupbit);
    bool const i_in_source = !!(mask[i] & source_grpbit) != invert_source;
    if (i_in_sensor && i_in_source) vector[i] += tf_types[type[i]] * q[i];
  }
}

/* ---------------------------------------------------------------------- */

void ElectrodeVector::get_charge_position(int i, double *xsite)
{
  charge_position(i, xsite);
}

/* ---------------------------------------------------------------------- */

void ElectrodeVector::add_charge_force(int i, const double *fcharge)
{
  charge_force(i, fcharge);
}

/* ---------------------------------------------------------------------- */

double ElectrodeVector::get_pair_cutsq(int itype, int jtype) const
{
  return pair_cutsq(itype, jtype);
}

/* ---------------------------------------------------------------------- */

double ElectrodeVector::get_charge_force_alpha() const
{
  return charge_force_alpha();
}

/* ---------------------------------------------------------------------- */

int ElectrodeVector::get_charge_force_virial(int i, const double *fcharge,
                                             double *v, int *vlist)
{
  return charge_force_virial(i, fcharge, v, vlist);
}

/* ---------------------------------------------------------------------- */
void ElectrodeVector::charge_position(int i, double *xsite)
{
  xsite[0] = atom->x[i][0];
  xsite[1] = atom->x[i][1];
  xsite[2] = atom->x[i][2];
}

/* ---------------------------------------------------------------------- */

void ElectrodeVector::charge_force(int i, const double *fcharge)
{
  atom->f[i][0] += fcharge[0];
  atom->f[i][1] += fcharge[1];
  atom->f[i][2] += fcharge[2];
}

/* ---------------------------------------------------------------------- */

int ElectrodeVector::charge_force_virial(int i, const double *fcharge,
                                         double *v, int *vlist)
{
  vlist[0] = i;
  const double *xi = atom->x[i];

  v[0] = xi[0] * fcharge[0];
  v[1] = xi[1] * fcharge[1];
  v[2] = xi[2] * fcharge[2];
  v[3] = xi[0] * fcharge[1];
  v[4] = xi[0] * fcharge[2];
  v[5] = xi[1] * fcharge[2];

  return 1;
}

/* ---------------------------------------------------------------------- */
double ElectrodeVector::pair_cutsq(int itype, int jtype) const
{
  return cutsq[itype][jtype];
}
