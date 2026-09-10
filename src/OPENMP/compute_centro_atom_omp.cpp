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

#include "compute_centro_atom_omp.h"

#include "atom.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCentroAtomOMP::ComputeCentroAtomOMP(LAMMPS *lmp, int narg, char **arg) :
    ComputeCentroAtom(lmp, narg, arg)
{
}

/* ----------------------------------------------------------------------
   threaded variant of ComputeCentroAtom::compute_peratom().  The per-atom
   neighbor scratch (distsq/nearest) and the pair-distance scratch (pairs)
   are allocated thread-locally inside the parallel region, sized to the
   largest neighbor count over the looped atoms.  The base select()/select2()
   helpers are reentrant given thread-private scratch, and each atom writes
   only its own output (centro[i] / array_atom[i]), so the result is
   bit-identical to the serial compute.
------------------------------------------------------------------------- */

void ComputeCentroAtomOMP::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow centro array if necessary
  // grow array_atom array if axes_flag set

  if (atom->nmax > nmax) {
    if (!axes_flag) {
      memory->destroy(centro);
      nmax = atom->nmax;
      memory->create(centro, nmax, "centro/atom:centro");
      vector_atom = centro;
    } else {
      memory->destroy(centro);
      memory->destroy(array_atom);
      nmax = atom->nmax;
      memory->create(centro, nmax, "centro/atom:centro");
      memory->create(array_atom, nmax, size_peratom_cols, "centro/atom:array_atom");
    }
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  // nhalf = number of pairs to sum, npairs = number of unique pairs

  const int nhalf = nnn / 2;
  const int npairs = nnn * (nnn - 1) / 2;

  // compute centro-symmetry parameter for each atom in group
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
    auto *const t_pairs = new double[npairs > 0 ? npairs : 1];

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

        // check whether to include local crystal symmetry axes

        if (!axes_flag) {

          // if not nnn neighbors, centro = 0.0

          if (n < nnn) {
            centro[i] = 0.0;
            continue;
          }

          // store nnn nearest neighs in 1st nnn locations of distsq and nearest

          select2(nnn, n, t_distsq, t_nearest);

          // R = Ri + Rj for each of npairs i,j pairs among nnn neighbors
          // pairs = squared length of each R

          int m = 0;
          for (int j = 0; j < nnn; j++) {
            const int jj = t_nearest[j];
            for (int k = j + 1; k < nnn; k++) {
              const int kk = t_nearest[k];
              const double delx = x[jj][0] + x[kk][0] - 2.0 * xtmp;
              const double dely = x[jj][1] + x[kk][1] - 2.0 * ytmp;
              const double delz = x[jj][2] + x[kk][2] - 2.0 * ztmp;
              t_pairs[m++] = delx * delx + dely * dely + delz * delz;
            }
          }

        } else {

          // calculate local crystal symmetry axes
          // rsq1, rsq2 are two smallest values of R^2
          // R1, R2 are corresponding vectors Ri - Rj
          // R3 is normal to R1, R2

          double *r1 = &array_atom[i][1];
          double *r2 = &array_atom[i][4];
          double *r3 = &array_atom[i][7];

          if (n < nnn) {
            centro[i] = 0.0;
            MathExtra::zero3(r1);
            MathExtra::zero3(r2);
            MathExtra::zero3(r3);
            continue;
          }

          // store nnn nearest neighs in 1st nnn locations of distsq and nearest

          select2(nnn, n, t_distsq, t_nearest);

          int m = 0;
          double rsq1, rsq2;
          rsq1 = rsq2 = cutsq;
          for (int j = 0; j < nnn; j++) {
            const int jj = t_nearest[j];
            for (int k = j + 1; k < nnn; k++) {
              const int kk = t_nearest[k];
              const double delx = x[jj][0] + x[kk][0] - 2.0 * xtmp;
              const double dely = x[jj][1] + x[kk][1] - 2.0 * ytmp;
              const double delz = x[jj][2] + x[kk][2] - 2.0 * ztmp;
              const double rsq = delx * delx + dely * dely + delz * delz;
              t_pairs[m++] = rsq;

              if (rsq < rsq2) {
                if (rsq < rsq1) {
                  rsq2 = rsq1;
                  MathExtra::copy3(r1, r2);
                  rsq1 = rsq;
                  MathExtra::sub3(x[jj], x[kk], r1);
                } else {
                  rsq2 = rsq;
                  MathExtra::sub3(x[jj], x[kk], r2);
                }
              }
            }
          }

          MathExtra::cross3(r1, r2, r3);
          MathExtra::norm3(r1);
          MathExtra::norm3(r2);
          MathExtra::norm3(r3);
        }

        // store nhalf smallest pair distances in 1st nhalf locations of pairs

        select(nhalf, npairs, t_pairs);

        // centrosymmetry = sum of nhalf smallest squared values

        double value = 0.0;
        for (int j = 0; j < nhalf; j++) value += t_pairs[j];
        centro[i] = value;

      } else {
        centro[i] = 0.0;
        if (axes_flag) {
          MathExtra::zero3(&array_atom[i][1]);
          MathExtra::zero3(&array_atom[i][4]);
          MathExtra::zero3(&array_atom[i][7]);
        }
      }
    }

    delete[] t_distsq;
    delete[] t_nearest;
    delete[] t_pairs;
  }

  if (axes_flag)
    for (int ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      if (mask[i] & groupbit) array_atom[i][0] = centro[i];
    }
}
