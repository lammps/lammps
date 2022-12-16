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

#include "wire_dipole.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ----------------------------------------------------------------------
   Wire-geometry correction term to dampen inter-wire interactions between
   periodically repeating wires.  Yields good approximation to 1D Ewald if
   adequate empty space is left between repeating wires (J. Mol. Struct.
   704, 101). x and y are non-periodic.
-------------------------------------------------------------------------
*/
WireDipole::WireDipole(LAMMPS *lmp) : BoundaryCorrection(lmp){};

void WireDipole::compute_corr(double /*qsum*/, int eflag_atom, int eflag_global, double &energy,
                              double *eatom)
{
  double const volume = get_volume();
  double *q = atom->q;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  double xdipole = 0.0;
  double ydipole = 0.0;
  for (int i = 0; i < nlocal; i++) {
    xdipole += q[i] * x[i][0];
    ydipole += q[i] * x[i][1];
  }

  // sum local contributions to get global dipole moment
  double xdipole_all, ydipole_all;
  MPI_Allreduce(&xdipole, &xdipole_all, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&ydipole, &ydipole_all, 1, MPI_DOUBLE, MPI_SUM, world);

  // need to make per-atom energy translationally invariant
  double xdipole_r2 = 0.0;
  double ydipole_r2 = 0.0;
  if (eflag_atom) {
    for (int i = 0; i < nlocal; i++) {
      xdipole_r2 += q[i] * x[i][0] * x[i][0];
      ydipole_r2 += q[i] * x[i][1] * x[i][1];
    }

    // sum local contributions
    double tmp;
    MPI_Allreduce(&xdipole_r2, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
    xdipole_r2 = tmp;
    MPI_Allreduce(&ydipole_r2, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
    ydipole_r2 = tmp;
  }

  // compute corrections
  const double e_wirecorr =
      MY_PI * (xdipole_all * xdipole_all + ydipole_all * ydipole_all) / volume;
  double const scale = 1.0;
  const double qscale = force->qqrd2e * scale;
  if (eflag_global) energy += qscale * e_wirecorr;

  // per-atom energy
  if (eflag_atom) {
    double efact = qscale * MY_PI / volume;
    for (int i = 0; i < nlocal; i++)
      eatom[i] += efact * q[i] *
          (x[i][0] * xdipole_all + x[i][1] * ydipole_all - 0.5 * (xdipole_r2 + ydipole_r2));
  }

  // add on force corrections
  double ffact = qscale * (-MY_2PI / volume);
  double **f = atom->f;
  for (int i = 0; i < nlocal; i++) {
    f[i][0] += ffact * q[i] * xdipole_all;
    f[i][1] += ffact * q[i] * ydipole_all;
  }
}

void WireDipole::vector_corr(double *vec, int sensor_grpbit, int source_grpbit, bool invert_source)
{
  double const volume = get_volume();
  int const nlocal = atom->nlocal;
  double **x = atom->x;
  double *q = atom->q;
  int *mask = atom->mask;
  double dipole[2] = {0., 0.};    // dipole in x and y direction
  for (int i = 0; i < nlocal; i++) {
    if (!!(mask[i] & source_grpbit) != invert_source) {
      for (int dim : {0, 1}) dipole[dim] += q[i] * x[i][dim];
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &dipole, 2, MPI_DOUBLE, MPI_SUM, world);
  for (double &d : dipole) d *= MY_2PI / volume;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & sensor_grpbit) vec[i] += x[i][0] * dipole[0] + x[i][1] * dipole[1];
  }
}

void WireDipole::matrix_corr(bigint *imat, double **matrix)
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

  // gather non-periodic positions of groups
  std::vector<double> xprd_local = std::vector<double>(ngrouplocal);
  std::vector<double> yprd_local = std::vector<double>(ngrouplocal);
  for (int i = 0, n = 0; i < nlocal; i++) {
    if (imat[i] < 0) continue;
    xprd_local[n] = x[i][0];
    yprd_local[n] = x[i][1];
    n++;
  }

  // gather subsets nprd positions
  std::vector<int> recvcounts = gather_recvcounts(ngrouplocal);
  std::vector<int> displs = gather_displs(recvcounts);
  std::vector<double> xprd_all = std::vector<double>(ngroup);
  std::vector<double> yprd_all = std::vector<double>(ngroup);
  MPI_Allgatherv(&xprd_local.front(), ngrouplocal, MPI_DOUBLE, &xprd_all.front(),
                 &recvcounts.front(), &displs.front(), MPI_DOUBLE, world);
  MPI_Allgatherv(&yprd_local.front(), ngrouplocal, MPI_DOUBLE, &yprd_all.front(),
                 &recvcounts.front(), &displs.front(), MPI_DOUBLE, world);

  std::vector<bigint> jmat = gather_jmat(imat);
  const double prefac = MY_2PI / volume;
  for (int i = 0; i < nlocal; i++) {
    if (imat[i] < 0) continue;
    for (bigint j = 0; j < ngroup; j++) {
      // matrix is symmetric
      if (jmat[j] > imat[i]) continue;
      double aij = prefac * (x[i][0] * xprd_all[j] + x[i][1] * yprd_all[j]);
      matrix[imat[i]][jmat[j]] += aij;
      if (imat[i] != jmat[j]) matrix[jmat[j]][imat[i]] += aij;
    }
  }
}
