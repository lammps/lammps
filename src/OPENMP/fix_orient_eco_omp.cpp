// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_orient_eco_omp.h"

#include "atom.h"
#include "comm.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"

#include <cmath>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "omp_compat.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixOrientECOOMP::FixOrientECOOMP(LAMMPS *lmp, int narg, char **arg) :
  FixOrientECO(lmp, narg, arg)
{
}

/* ----------------------------------------------------------------------
   threaded version of FixOrientECO::post_force().

   Phase 1 (order parameter + wave-function terms) writes only the owning
   atom's nbr[i]/order[i] slots and reduces added_energy across threads.
   Phase 2 (force) writes ONLY f[i] -- it reads neighbor data nbr[j]
   (acquired for ghosts via the forward_comm) but never writes f[j].  All
   per-atom scratch is declared inside the loop bodies so it is thread
   private.  Bit-for-bit identical to the serial result per atom.
------------------------------------------------------------------------- */

void FixOrientECOOMP::post_force(int /* vflag */)
{
  const double omega_pre = MY_PI2 * inv_eta;
  const double duchi_pre = half_u * MY_PI * inv_eta * inv_norm_fac;

  added_energy = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  const int nall = atom->nlocal + atom->nghost;

  const int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // ensure nbr and order data structures are adequate size

  if (nall > nmax) {
    nmax = nall;
    memory->destroy(nbr);
    memory->destroy(order);
    nbr = (Nbr *) memory->smalloc(nmax * sizeof(Nbr), "orient/eco:nbr");
    memory->create(order, nmax, 2, "orient/eco:order");
    array_atom = order;
  }

  // loop over owned atoms and compute order parameter

  double e_acc = 0.0;

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic) default(shared) reduction(+:e_acc)
#endif
  for (int ii = 0; ii < inum; ++ii) {
    const int i = ilist[ii];
    int *jlist = firstneigh[i];
    const int jnum = numneigh[i];

    // initializations
    double chi = 0.0;
    for (int k = 0; k < 3; ++k) {
      nbr[i].real_phi[0][k] = nbr[i].real_phi[1][k] = 0.0;
      nbr[i].imag_phi[0][k] = nbr[i].imag_phi[1][k] = 0.0;
    }

    // loop over all neighbors of atom i
    for (int jj = 0; jj < jnum; ++jj) {
      int j = jlist[jj];
      j &= NEIGHMASK;

      const double dx = x[i][0] - x[j][0];
      const double dy = x[i][1] - x[j][1];
      const double dz = x[i][2] - x[j][2];
      double squared_distance = dx * dx + dy * dy + dz * dz;

      if (squared_distance < squared_cutoff) {
        squared_distance *= inv_squared_cutoff;
        const double weight = squared_distance * (squared_distance - 2.0) + 1.0;

        for (int lambda = 0; lambda < 2; ++lambda) {
          for (int k = 0; k < 3; ++k) {
            const double scalar_product = reciprocal_vectors[lambda][k][0] * dx +
              reciprocal_vectors[lambda][k][1] * dy + reciprocal_vectors[lambda][k][2] * dz;
            nbr[i].real_phi[lambda][k] += weight * cos(scalar_product);
            nbr[i].imag_phi[lambda][k] += weight * sin(scalar_product);
          }
        }
      }
    }

    // collect contributions
    for (int k = 0; k < 3; ++k) {
      chi += (nbr[i].real_phi[0][k] * nbr[i].real_phi[0][k] + nbr[i].imag_phi[0][k] * nbr[i].imag_phi[0][k] -
                nbr[i].real_phi[1][k] * nbr[i].real_phi[1][k] - nbr[i].imag_phi[1][k] * nbr[i].imag_phi[1][k]);
    }
    chi *= inv_norm_fac;
    order[i][0] = chi;

    // compute normalized order parameter
    // and potential energy
    if (chi > eta) {
      e_acc += half_u;
      nbr[i].duchi = 0.0;
      order[i][1] = sign;
    } else if (chi < -eta) {
      e_acc -= half_u;
      nbr[i].duchi = 0.0;
      order[i][1] = -sign;
    } else {
      const double omega = omega_pre * chi;
      const double sin_om = sin(omega);

      e_acc += half_u * sin_om;
      nbr[i].duchi = duchi_pre * cos(omega);
      order[i][1] = sign * sin_om;
    }

    // compute product with potential derivative
    for (int k = 0; k < 3; ++k) {
      for (int lambda = 0; lambda < 2; ++lambda) {
        nbr[i].real_phi[lambda][k] *= nbr[i].duchi;
        nbr[i].imag_phi[lambda][k] *= nbr[i].duchi;
      }
    }
  }

  added_energy = e_acc;

  // compute force only if synthetic potential is not zero

  if (u_0 != 0.0) {
    // communicate to acquire nbr data for ghost atoms
    comm->forward_comm(this);

    // loop over all atoms; only f[i] is written -> race-free partition

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic) default(shared)
#endif
    for (int ii = 0; ii < inum; ++ii) {
      const int i = ilist[ii];

      // skip atoms not in group
      if (!(mask[i] & groupbit)) continue;

      int *jlist = firstneigh[i];
      const int jnum = numneigh[i];
      const bool no_boundary_atom = (nbr[i].duchi == 0.0);

      // per-atom scratch (thread private)
      double gradient_ii_cos[2][3][3];
      double gradient_ii_sin[2][3][3];
      double gradient_ij_vec[2][3][3];
      double gradient_ij_sca[2][3];

      for (int k = 0; k < 3; ++k) {
        for (int lambda = 0; lambda < 2; ++lambda) {
          for (int dim = 0; dim < 3; ++dim) {
            gradient_ii_cos[lambda][k][dim] = 0.0;
            gradient_ii_sin[lambda][k][dim] = 0.0;
            gradient_ij_vec[lambda][k][dim] = 0.0;
          }
          gradient_ij_sca[lambda][k] = 0.0;
        }
      }

      // loop over all neighbors of atom i
      // for those within squared_cutoff, compute force
      for (int jj = 0; jj < jnum; ++jj) {
        int j = jlist[jj];
        j &= NEIGHMASK;

        // do not compute force on atom i if it is far from boundary
        if ((nbr[j].duchi == 0.0) && no_boundary_atom) continue;

        const double dx = x[i][0] - x[j][0];
        const double dy = x[i][1] - x[j][1];
        const double dz = x[i][2] - x[j][2];
        double squared_distance = dx * dx + dy * dy + dz * dz;

        if (squared_distance < squared_cutoff) {
          squared_distance *= inv_squared_cutoff;
          const double weight = squared_distance * (squared_distance - 2.0) + 1.0;
          const double weight_gradient_prefactor = 4.0 * (squared_distance - 1.0) * inv_squared_cutoff;
          double weight_gradient[3];
          weight_gradient[0] = weight_gradient_prefactor * dx;
          weight_gradient[1] = weight_gradient_prefactor * dy;
          weight_gradient[2] = weight_gradient_prefactor * dz;

          for (int lambda = 0; lambda < 2; ++lambda) {
            for (int k = 0; k < 3; ++k) {
              const double scalar_product = reciprocal_vectors[lambda][k][0] * dx +
                reciprocal_vectors[lambda][k][1] * dy + reciprocal_vectors[lambda][k][2] * dz;
              const double cos_scalar_product = cos(scalar_product);
              const double sin_scalar_product = sin(scalar_product);
              for (int dim = 0; dim < 3; ++dim) {
                const double gcos_scalar_product = weight_gradient[dim] * cos_scalar_product;
                const double gsin_scalar_product = weight_gradient[dim] * sin_scalar_product;
                gradient_ii_cos[lambda][k][dim] += gcos_scalar_product;
                gradient_ii_sin[lambda][k][dim] += gsin_scalar_product;
                gradient_ij_vec[lambda][k][dim] += (nbr[j].real_phi[lambda][k] * gcos_scalar_product -
                                                    nbr[j].imag_phi[lambda][k] * gsin_scalar_product);
              }
              gradient_ij_sca[lambda][k] += weight * (nbr[j].real_phi[lambda][k] * sin_scalar_product +
                                                      nbr[j].imag_phi[lambda][k] * cos_scalar_product);
            }
          }
        }
      }

      // sum contributions
      for (int k = 0; k < 3; ++k) {
        for (int dim = 0; dim < 3; ++dim) {
          f[i][dim] -= (nbr[i].real_phi[0][k] * gradient_ii_cos[0][k][dim] + nbr[i].imag_phi[0][k] * gradient_ii_sin[0][k][dim] + gradient_ij_vec[0][k][dim] + reciprocal_vectors[1][k][dim] * gradient_ij_sca[1][k]);
          f[i][dim] += (nbr[i].real_phi[1][k] * gradient_ii_cos[1][k][dim] + nbr[i].imag_phi[1][k] * gradient_ii_sin[1][k][dim] + gradient_ij_vec[1][k][dim] + reciprocal_vectors[0][k][dim] * gradient_ij_sca[0][k]);
        }
      }
    }
  }
}
