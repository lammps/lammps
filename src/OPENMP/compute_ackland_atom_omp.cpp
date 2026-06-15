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

#include "compute_ackland_atom_omp.h"

#include "atom.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

// this mirrors the file-static enum in compute_ackland_atom.cpp

enum { UNKNOWN, BCC, FCC, HCP, ICO };

/* ---------------------------------------------------------------------- */

ComputeAcklandAtomOMP::ComputeAcklandAtomOMP(LAMMPS *lmp, int narg, char **arg) :
    ComputeAcklandAtom(lmp, narg, arg)
{
}

/* ----------------------------------------------------------------------
   threaded variant of ComputeAcklandAtom::compute_peratom().  The per-atom
   neighbor scratch (distsq/nearest/nearest_n0/nearest_n1) is allocated
   thread-locally inside the parallel region, sized to the largest neighbor
   count over the looped atoms; the chi[] histogram is per-iteration local.
   select2() is reentrant given thread-private scratch, and each atom writes
   only its own output structure[i], so the result is bit-identical to the
   serial compute.
------------------------------------------------------------------------- */

void ComputeAcklandAtomOMP::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow structure array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(structure);
    nmax = atom->nmax;
    memory->create(structure, nmax, "compute/ackland/atom:ackland");
    vector_atom = structure;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  // compute structure parameter for each atom in group
  // use full neighbor list

  const double *const *const x = atom->x;
  const int *const mask = atom->mask;
  const double cutsq = force->pair->cutforce * force->pair->cutforce;

  // largest neighbor count over the looped atoms, for sizing thread scratch

  int maxn = 1;
  for (int ii = 0; ii < inum; ii++) {
    const int jnum = numneigh[ilist[ii]];
    if (jnum > maxn) maxn = jnum;
  }

#if defined(_OPENMP)
#pragma omp parallel
#endif
  {
    auto *const t_distsq = new double[maxn];
    auto *const t_nearest = new int[maxn];
    auto *const t_nearest_n0 = new int[maxn];
    auto *const t_nearest_n1 = new int[maxn];

#if defined(_OPENMP)
#pragma omp for schedule(dynamic)
#endif
    for (int ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      if (mask[i] & groupbit) {
        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];
        const int *const jlist = firstneigh[i];
        const int jnum = numneigh[i];

        // loop over list of all neighbors within force cutoff
        // t_distsq[] = distance sq to each
        // t_nearest[] = atom indices of neighbors

        int n = 0;
        for (int jj = 0; jj < jnum; jj++) {
          const int j = jlist[jj] & NEIGHMASK;
          const double delx = xtmp - x[j][0];
          const double dely = ytmp - x[j][1];
          const double delz = ztmp - x[j][2];
          const double rsq = delx * delx + dely * dely + delz * delz;
          if (rsq < cutsq) {
            t_distsq[n] = rsq;
            t_nearest[n++] = j;
          }
        }

        // Select 6 nearest neighbors

        select2(6, n, t_distsq, t_nearest);

        // Mean squared separation

        double r0_sq = 0.;
        for (int j = 0; j < 6; j++) r0_sq += t_distsq[j];
        r0_sq /= 6.;

        // n0 near neighbors with: distsq<1.45*r0_sq
        // n1 near neighbors with: distsq<1.55*r0_sq

        double n0_dist_sq = 1.45 * r0_sq, n1_dist_sq = 1.55 * r0_sq;
        int n0 = 0, n1 = 0;
        for (int j = 0; j < n; j++) {
          if (t_distsq[j] < n1_dist_sq) {
            t_nearest_n1[n1++] = t_nearest[j];
            if (t_distsq[j] < n0_dist_sq) { t_nearest_n0[n0++] = t_nearest[j]; }
          }
        }

        // Evaluate all angles <(r_ij,rik) forall n0 particles with:
        // distsq < 1.45*r0_sq

        int chi[8];
        chi[0] = chi[1] = chi[2] = chi[3] = chi[4] = chi[5] = chi[6] = chi[7] = 0;
        for (int j = 0; j < n0; j++) {
          const double x_ij = x[i][0] - x[t_nearest_n0[j]][0];
          const double y_ij = x[i][1] - x[t_nearest_n0[j]][1];
          const double z_ij = x[i][2] - x[t_nearest_n0[j]][2];
          const double norm_j = sqrt(x_ij * x_ij + y_ij * y_ij + z_ij * z_ij);
          if (norm_j <= 0.) continue;
          for (int k = j + 1; k < n0; k++) {
            const double x_ik = x[i][0] - x[t_nearest_n0[k]][0];
            const double y_ik = x[i][1] - x[t_nearest_n0[k]][1];
            const double z_ik = x[i][2] - x[t_nearest_n0[k]][2];
            const double norm_k = sqrt(x_ik * x_ik + y_ik * y_ik + z_ik * z_ik);
            if (norm_k <= 0.) continue;

            const double bond_angle =
                (x_ij * x_ik + y_ij * y_ik + z_ij * z_ik) / (norm_j * norm_k);

            // Histogram for identifying the relevant peaks

            if (bond_angle < -0.945)
              chi[0]++;
            else if (bond_angle < -0.915)
              chi[1]++;
            else if (bond_angle < -0.755)
              chi[2]++;
            else if (bond_angle < -0.195)
              chi[3]++;
            else if (bond_angle < 0.195)
              chi[4]++;
            else if (bond_angle < 0.245)
              chi[5]++;
            else if (bond_angle < 0.795)
              chi[6]++;
            else
              chi[7]++;
          }
        }

        if (legacy) {

          // This is the original implementation by Gerolf Ziegenhain
          // Deviations from the different lattice structures

          double delta_bcc = 0.35 * chi[4] / (double) (chi[5] + chi[6] - chi[4]);
          double delta_cp = fabs(1. - (double) chi[6] / 24.);
          double delta_fcc = 0.61 * (fabs((double) (chi[0] + chi[1] - 6.)) + (double) chi[2]) / 6.0;
          double delta_hcp = (fabs((double) chi[0] - 3.) +
                              fabs((double) chi[0] + (double) chi[1] + (double) chi[2] +
                                   (double) chi[3] - 9.0)) /
              12.0;

          // Identification of the local structure according to the reference

          if (chi[0] == 7) {
            delta_bcc = 0.;
          } else if (chi[0] == 6) {
            delta_fcc = 0.;
          } else if (chi[0] <= 3) {
            delta_hcp = 0.;
          }

          if (chi[7] > 0.)
            structure[i] = UNKNOWN;
          else if (chi[4] < 3.) {
            if (n1 > 13 || n1 < 11)
              structure[i] = UNKNOWN;
            else
              structure[i] = ICO;
          } else if (delta_bcc <= delta_cp) {
            if (n1 < 11)
              structure[i] = UNKNOWN;
            else
              structure[i] = BCC;
          } else if (n1 > 12 || n1 < 11)
            structure[i] = UNKNOWN;
          else if (delta_fcc < delta_hcp)
            structure[i] = FCC;
          else
            structure[i] = HCP;

        } else {

          // This is the updated implementation by Brian Barnes

          if (chi[7] > 0 || n0 < 11)
            structure[i] = UNKNOWN;
          else if (chi[0] == 7)
            structure[i] = BCC;
          else if (chi[0] == 6)
            structure[i] = FCC;
          else if (chi[0] == 3)
            structure[i] = HCP;
          else {
            // Deviations from the different lattice structures

            double delta_cp = fabs(1. - (double) chi[6] / 24.);

            // ensure we do not get divide by zero
            // and if we will, make delta_bcc irrelevant
            double delta_bcc = delta_cp + 1.0;
            int chi56m4 = chi[5] + chi[6] - chi[4];

            // note that chi[7] presumed zero
            if (chi56m4 != 0) delta_bcc = 0.35 * chi[4] / (double) chi56m4;

            double delta_fcc = 0.61 * (fabs((double) (chi[0] + chi[1] - 6)) + (double) chi[2]) / 6.0;

            double delta_hcp = (fabs((double) chi[0] - 3.) +
                                fabs((double) chi[0] + (double) chi[1] + (double) chi[2] +
                                     (double) chi[3] - 9.0)) /
                12.0;

            // Identification of the local structure according to the reference

            if (delta_bcc >= 0.1 && delta_cp >= 0.1 && delta_fcc >= 0.1 && delta_hcp >= 0.1)
              structure[i] = UNKNOWN;

            // not part of Ackland-Jones 2006; included for backward compatibility
            if (chi[4] < 3. && n1 == 12)
              structure[i] = ICO;

            else {
              if (delta_bcc <= delta_cp && n1 > 10 && n1 < 13)
                structure[i] = BCC;
              else {
                if (n0 > 12)
                  structure[i] = UNKNOWN;
                else {
                  if (delta_fcc < delta_hcp)
                    structure[i] = FCC;
                  else
                    structure[i] = HCP;
                }
              }
            }
          }
        }
      } else
        structure[i] = 0.0;
    }

    delete[] t_distsq;
    delete[] t_nearest;
    delete[] t_nearest_n0;
    delete[] t_nearest_n1;
  }
}
