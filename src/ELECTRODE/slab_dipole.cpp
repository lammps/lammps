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

#include "slab_dipole.h"

#include "atom.h"
#include "domain.h"
#include "force.h"
#include "kspace.h"
#include "math_const.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr double SMALL = 0.00001;

/* ----------------------------------------------------------------------
   Slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2D Ewald if
   adequate empty space is left between repeating slabs (J. Chem. Phys.
   111, 3155).  Slabs defined here to be parallel to the xy plane. Also
   extended to non-neutral systems (J. Chem. Phys. 131, 094107).
-------------------------------------------------------------------------
*/
SlabDipole::SlabDipole(LAMMPS *lmp) : BoundaryCorrection(lmp){};

void SlabDipole::compute_corr(double qsum, int eflag_atom, int eflag_global, double &energy,
                              double *eatom)
{
  // compute local contribution to global dipole moment
  double const volume = get_volume();
  double *q = atom->q;
  double **x = atom->x;
  double zprd_slab = domain->zprd * force->kspace->slab_volfactor;
  int nlocal = atom->nlocal;
  double dipole = 0.0;
  for (int i = 0; i < nlocal; i++) dipole += q[i] * x[i][2];

  // sum local contributions to get global dipole moment
  double dipole_all;
  MPI_Allreduce(&dipole, &dipole_all, 1, MPI_DOUBLE, MPI_SUM, world);

  // need to make non-neutral systems and/or
  //  per-atom energy translationally invariant
  double dipole_r2 = 0.0;
  if (eflag_atom || fabs(qsum) > SMALL) {
    for (int i = 0; i < nlocal; i++) dipole_r2 += q[i] * x[i][2] * x[i][2];
    // sum local contributions
    double tmp;
    MPI_Allreduce(&dipole_r2, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
    dipole_r2 = tmp;
  }

  // compute corrections
  double const e_slabcorr = MY_2PI *
      (dipole_all * dipole_all - qsum * dipole_r2 - qsum * qsum * zprd_slab * zprd_slab / 12.0) /
      volume;
  double const scale = 1.0;
  double const qscale = force->qqrd2e * scale;
  if (eflag_global) energy += qscale * e_slabcorr;

  // per-atom energy
  if (eflag_atom) {
    double efact = qscale * MY_2PI / volume;
    for (int i = 0; i < nlocal; i++)
      eatom[i] += efact * q[i] *
          (x[i][2] * dipole_all - 0.5 * (dipole_r2 + qsum * x[i][2] * x[i][2]) -
           qsum * zprd_slab * zprd_slab / 12.0);
  }

  // add on force corrections
  double const ffact = qscale * (-4.0 * MY_PI / volume);
  double **f = atom->f;
  for (int i = 0; i < nlocal; i++) f[i][2] += ffact * q[i] * (dipole_all - qsum * x[i][2]);
}

void SlabDipole::vector_corr(double *vec, int sensor_grpbit, int source_grpbit, bool invert_source)
{
  double const volume = get_volume();
  int const nlocal = atom->nlocal;
  double **x = atom->x;
  double *q = atom->q;
  int *mask = atom->mask;
  double dipole = 0.;
  for (int i = 0; i < nlocal; i++) {
    if (!!(mask[i] & source_grpbit) != invert_source) dipole += q[i] * x[i][2];
  }
  MPI_Allreduce(MPI_IN_PLACE, &dipole, 1, MPI_DOUBLE, MPI_SUM, world);
  dipole *= 4.0 * MY_PI / volume;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & sensor_grpbit) vec[i] += x[i][2] * dipole;
  }
}

void SlabDipole::matrix_corr(bigint *imat, double **matrix)
{
  double const volume = get_volume();
  int nlocal = atom->nlocal;
  double **x = atom->x;

  // how many local and total group atoms?
  int ngrouplocal = 0;
  for (int i = 0; i < nlocal; i++)
    if (imat[i] > -1) ngrouplocal++;
  bigint ngroup = 0;
  MPI_Allreduce(&ngrouplocal, &ngroup, 1, MPI_INT, MPI_SUM, world);

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

  std::vector<bigint> jmat = gather_jmat(imat);
  const double prefac = MY_4PI / volume;
  for (int i = 0; i < nlocal; i++) {
    if (imat[i] < 0) continue;
    for (bigint j = 0; j < ngroup; j++) {
      if (jmat[j] > imat[i]) continue;    // matrix is symmetric
      double aij = prefac * x[i][2] * nprd_all[j];
      matrix[imat[i]][jmat[j]] += aij;
      if (imat[i] != jmat[j]) matrix[jmat[j]][imat[i]] += aij;
    }
  }
}
