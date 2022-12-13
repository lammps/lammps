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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meißner (TUHH)
------------------------------------------------------------------------- */

#include "slab_2d.h"

#include "atom.h"
#include "domain.h"
#include "force.h"
#include "kspace.h"
#include "math_const.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ----------------------------------------------------------------------
   Slab-geometry correction term (k=0) of EW2D. See Hu, JCTC 10:12 (2014)
   pp. 5254-5264 or metalwalls ewald and parallelization documentation.
------------------------------------------------------------------------- */
Slab2d::Slab2d(LAMMPS *lmp) : BoundaryCorrection(lmp){};

void Slab2d::compute_corr(double /*qsum*/, int eflag_atom, int eflag_global, double &energy,
                          double *eatom)
{
  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  double const g_ewald = force->kspace->g_ewald;
  bigint natoms = atom->natoms;

  std::vector<double> z = std::vector<double>(nlocal);
  for (int i = 0; i < nlocal; i++) z[i] = x[i][2];
  std::vector<double> z_all = std::vector<double>(natoms);
  std::vector<double> q_all = std::vector<double>(natoms);
  std::vector<int> recvcounts = gather_recvcounts(nlocal);
  std::vector<int> displs = gather_displs(recvcounts);
  MPI_Allgatherv(q, nlocal, MPI_DOUBLE, &q_all.front(), &recvcounts.front(), &displs.front(),
                 MPI_DOUBLE, world);
  MPI_Allgatherv(&z.front(), nlocal, MPI_DOUBLE, &z_all.front(), &recvcounts.front(),
                 &displs.front(), MPI_DOUBLE, world);

  const double g_ewald_inv = 1.0 / g_ewald;
  double const scale = 1.0;
  const double qscale = force->qqrd2e * scale;
  double const area = domain->xprd * domain->yprd;
  const double ffact = qscale * MY_2PI / area;
  const double efact = qscale * MY_PIS / area;
  double e_keq0 = 0;
  for (int i = 0; i < nlocal; i++) {
    double pot_ij = 0.0;

    for (bigint j = 0; j < natoms; j++) {
      double const zij = z_all[j] - x[i][2];
      double const g_zij = g_ewald * zij;

      // coulomb potential; see eq. (4) in metalwalls parallelization doc
      pot_ij += q_all[j] * (exp(-g_zij * g_zij) * g_ewald_inv + MY_PIS * zij * erf(g_zij));
      f[i][2] -= ffact * q[i] * q_all[j] * erf(g_zij);
    }

    // per-atom energy; see eq. (20) in metalwalls ewald doc
    if (eflag_atom) eatom[i] -= efact * q[i] * pot_ij;
    if (eflag_global) e_keq0 -= q[i] * pot_ij;
  }
  if (eflag_global) {
    MPI_Allreduce(MPI_IN_PLACE, &e_keq0, 1, MPI_DOUBLE, MPI_SUM, world);
    energy += efact * e_keq0;
  }
}

void Slab2d::vector_corr(double *vec, int sensor_grpbit, int source_grpbit, bool invert_source)
{
  int const nlocal = atom->nlocal;
  double **x = atom->x;
  double *q = atom->q;
  int *mask = atom->mask;
  std::vector<double> z_local;    // z coordinates of electrolyte atoms
  std::vector<double> q_local;    // charges of electrolyte atoms
  for (int i = 0; i < nlocal; i++) {
    if (!!(mask[i] & source_grpbit) != invert_source) {
      z_local.push_back(x[i][2]);
      q_local.push_back(q[i]);
    }
  }

  int n_electrolyte_local = z_local.size();
  int n_electrolyte;
  MPI_Allreduce(&n_electrolyte_local, &n_electrolyte, 1, MPI_INT, MPI_SUM, world);
  std::vector<double> z_all = std::vector<double>(n_electrolyte);
  std::vector<double> q_all = std::vector<double>(n_electrolyte);
  std::vector<int> recvcounts = gather_recvcounts(n_electrolyte_local);
  std::vector<int> displs = gather_displs(recvcounts);
  MPI_Allgatherv(&z_local.front(), n_electrolyte_local, MPI_DOUBLE, &z_all.front(),
                 &recvcounts.front(), &displs.front(), MPI_DOUBLE, world);
  MPI_Allgatherv(&q_local.front(), n_electrolyte_local, MPI_DOUBLE, &q_all.front(),
                 &recvcounts.front(), &displs.front(), MPI_DOUBLE, world);
  double const g_ewald = force->kspace->g_ewald;
  double const area = domain->xprd * domain->yprd;
  double const prefac = 2 * MY_PIS / area;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & sensor_grpbit)) continue;
    double b = 0;
    double zi = x[i][2];
    for (size_t j = 0; j < z_all.size(); j++) {
      double zij = z_all[j] - zi;
      double gzij = g_ewald * zij;
      double zfac = (gzij > 7)
          ? MY_PIS * zij
          : exp(-(gzij * gzij)) / g_ewald + MY_PIS * zij * erf(gzij);    // avoid exp underflow
      b += q_all[j] * zfac;
    }
    vec[i] -= prefac * b;
  }
}

void Slab2d::matrix_corr(bigint *imat, double **matrix)
{
  int nlocal = atom->nlocal;
  double **x = atom->x;

  // how many local and total group atoms?
  int ngrouplocal = 0;
  for (int i = 0; i < nlocal; i++)
    if (imat[i] > -1) ngrouplocal++;
  bigint ngroup = 0;
  MPI_Allreduce(&ngrouplocal, &ngroup, 1, MPI_INT, MPI_SUM, world);

  // gather non-periodic positions of groups
  std::vector<double> nprd_local = std::vector<double>(ngrouplocal);
  for (int i = 0, n = 0; i < nlocal; i++) {
    if (imat[i] < 0) continue;
    nprd_local[n++] = x[i][2];
  }

  // gather subsets nprd positions
  std::vector<int> recvcounts = gather_recvcounts(ngrouplocal);
  std::vector<int> displs = gather_displs(recvcounts);
  std::vector<double> nprd_all = std::vector<double>(ngroup);
  MPI_Allgatherv(&nprd_local.front(), ngrouplocal, MPI_DOUBLE, &nprd_all.front(),
                 &recvcounts.front(), &displs.front(), MPI_DOUBLE, world);

  double const g_ewald = force->kspace->g_ewald;
  const double g_ewald_inv = 1.0 / g_ewald;
  const double g_ewald_sq = g_ewald * g_ewald;
  double const area = domain->xprd * domain->yprd;
  const double prefac = 2.0 * MY_PIS / area;
  std::vector<bigint> jmat = gather_jmat(imat);
  for (int i = 0; i < nlocal; i++) {
    if (imat[i] < 0) continue;
    for (bigint j = 0; j < ngroup; j++) {
      // matrix is symmetric
      if (jmat[j] > imat[i]) continue;
      double dij = nprd_all[j] - x[i][2];
      double aij =
          prefac * (exp(-dij * dij * g_ewald_sq) * g_ewald_inv + MY_PIS * dij * erf(dij * g_ewald));
      matrix[imat[i]][jmat[j]] -= aij;
      if (imat[i] != jmat[j]) matrix[jmat[j]][imat[i]] -= aij;
    }
  }
}
