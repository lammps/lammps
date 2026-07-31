/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This file is distributed under the GNU General Public License.

   Direct Ewald response for periodic QM images.  The formulas follow the
   half-reciprocal-space convention used by Amber sander's qm_ewald module.
------------------------------------------------------------------------- */

#include "qmmm_xtb_ewald.h"

#include "math_const.h"

#include <cmath>
#include <stdexcept>

using namespace LAMMPS_NS;
using namespace MathConst;

void QMMMXTBEwald::setup(const CellMatrix &cell, double alpha, const std::array<int, 3> &kmax,
                         int ksqmax)
{
  alpha_ = alpha;
  kterms_.clear();

  // H stores the direct lattice vectors as columns.  Its inverse transpose
  // maps an integer reciprocal vector n to H^-T n; multiplying by 2*pi then
  // gives the Cartesian k vector used in the Ewald phase and force.
  const double determinant = cell[0] * (cell[4] * cell[8] - cell[5] * cell[7]) -
      cell[1] * (cell[3] * cell[8] - cell[5] * cell[6]) +
      cell[2] * (cell[3] * cell[7] - cell[4] * cell[6]);
  volume_ = determinant;

  if (alpha <= 0.0 || volume_ <= 0.0 || ksqmax <= 0)
    throw std::invalid_argument("Invalid direct-Ewald parameters");

  const double inverse[9] = {(cell[4] * cell[8] - cell[5] * cell[7]) / determinant,
                             (cell[2] * cell[7] - cell[1] * cell[8]) / determinant,
                             (cell[1] * cell[5] - cell[2] * cell[4]) / determinant,
                             (cell[5] * cell[6] - cell[3] * cell[8]) / determinant,
                             (cell[0] * cell[8] - cell[2] * cell[6]) / determinant,
                             (cell[2] * cell[3] - cell[0] * cell[5]) / determinant,
                             (cell[3] * cell[7] - cell[4] * cell[6]) / determinant,
                             (cell[1] * cell[6] - cell[0] * cell[7]) / determinant,
                             (cell[0] * cell[4] - cell[1] * cell[3]) / determinant};

  for (int ix = 0; ix <= kmax[0]; ++ix) {
    const int iylo = (ix == 0) ? 0 : -kmax[1];
    for (int iy = iylo; iy <= kmax[1]; ++iy) {
      const int izlo = (ix == 0 && iy == 0) ? 1 : -kmax[2];
      for (int iz = izlo; iz <= kmax[2]; ++iz) {
        const int isq = ix * ix + iy * iy + iz * iz;
        if (isq == 0 || isq > ksqmax) continue;

        KTerm term;
        const double integer_k[3] = {static_cast<double>(ix), static_cast<double>(iy),
                                     static_cast<double>(iz)};
        for (int dim = 0; dim < 3; ++dim)
          term.k[dim] = MY_2PI *
              (inverse[dim] * integer_k[0] + inverse[3 + dim] * integer_k[1] +
               inverse[6 + dim] * integer_k[2]);
        const double ksq = term.k[0] * term.k[0] + term.k[1] * term.k[1] + term.k[2] * term.k[2];
        term.coefficient = 4.0 * MY_PI / volume_ * std::exp(-ksq / (4.0 * alpha * alpha)) / ksq;
        kterms_.push_back(term);
      }
    }
  }

  if (kterms_.empty()) throw std::invalid_argument("Direct-Ewald K-vector set is empty");
}

void QMMMXTBEwald::response(const std::vector<double> &x, std::vector<double> &matrix) const
{
  const int n = x.size() / 3;
  matrix.assign(static_cast<std::size_t>(n) * n, -MY_PI / (alpha_ * alpha_ * volume_));

  for (const auto &term : kterms_) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        const double phase = term.k[0] * (x[3 * i] - x[3 * j]) +
            term.k[1] * (x[3 * i + 1] - x[3 * j + 1]) + term.k[2] * (x[3 * i + 2] - x[3 * j + 2]);
        matrix[static_cast<std::size_t>(i) * n + j] += 2.0 * term.coefficient * std::cos(phase);
      }
    }
  }

  for (int i = 0; i < n; ++i) {
    matrix[static_cast<std::size_t>(i) * n + i] -= 2.0 * alpha_ / MY_PIS;
    for (int j = 0; j < i; ++j) {
      const double dx = x[3 * i] - x[3 * j];
      const double dy = x[3 * i + 1] - x[3 * j + 1];
      const double dz = x[3 * i + 2] - x[3 * j + 2];
      const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (r == 0.0) throw std::invalid_argument("Coincident QM atoms in direct Ewald");
      const double correction = std::erf(alpha_ * r) / r;
      matrix[static_cast<std::size_t>(i) * n + j] -= correction;
      matrix[static_cast<std::size_t>(j) * n + i] -= correction;
    }
  }
}

double QMMMXTBEwald::energy(const std::vector<double> &x, const std::vector<double> &q) const
{
  std::vector<double> matrix;
  response(x, matrix);
  const int n = q.size();
  double value = 0.0;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      value += 0.5 * q[i] * matrix[static_cast<std::size_t>(i) * n + j] * q[j];
  return value;
}

void QMMMXTBEwald::add_force(const std::vector<double> &x, const std::vector<double> &q,
                             std::vector<double> &force) const
{
  const int n = q.size();
  if (force.size() != x.size()) force.assign(x.size(), 0.0);

  // Derivative of the subtracted primary-cell erf(alpha*r)/r term.
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      add_erf_pair_force(&x[3 * i], &x[3 * j], q[i], q[j], alpha_, &force[3 * i], &force[3 * j]);
    }
  }

  // Reciprocal force from U_k = coefficient * |sum_i q_i exp(i k.r_i)|^2
  // for one representative of every +/-k pair.
  for (const auto &term : kterms_) {
    double structure_cos = 0.0;
    double structure_sin = 0.0;
    for (int i = 0; i < n; ++i) {
      const double phase =
          term.k[0] * x[3 * i] + term.k[1] * x[3 * i + 1] + term.k[2] * x[3 * i + 2];
      structure_cos += q[i] * std::cos(phase);
      structure_sin += q[i] * std::sin(phase);
    }
    for (int i = 0; i < n; ++i) {
      const double phase =
          term.k[0] * x[3 * i] + term.k[1] * x[3 * i + 1] + term.k[2] * x[3 * i + 2];
      const double scalar = 2.0 * term.coefficient * q[i] *
          (std::sin(phase) * structure_cos - std::cos(phase) * structure_sin);
      force[3 * i] += scalar * term.k[0];
      force[3 * i + 1] += scalar * term.k[1];
      force[3 * i + 2] += scalar * term.k[2];
    }
  }
}

void QMMMXTBEwald::add_erf_pair_force(const double *xi, const double *xj, double qi, double qj,
                                      double alpha, double *fi, double *fj)
{
  const double delta[3] = {xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]};
  const double rsq = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];
  const double r = std::sqrt(rsq);
  if (r == 0.0) throw std::invalid_argument("Coincident charges in erf correction");

  const double ar = alpha * r;
  const double prefactor =
      qi * qj * (2.0 * alpha / MY_PIS * std::exp(-ar * ar) / rsq - std::erf(ar) / (rsq * r));
  for (int dim = 0; dim < 3; ++dim) {
    const double fij = prefactor * delta[dim];
    fi[dim] += fij;
    fj[dim] -= fij;
  }
}
