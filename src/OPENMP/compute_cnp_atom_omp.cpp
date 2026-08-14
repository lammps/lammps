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

#include "compute_cnp_atom_omp.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;

// these mirror the file-static constants in compute_cnp_atom.cpp

static constexpr int MAXNEAR = 24;
static constexpr int MAXCOMMON = 12;

/* ---------------------------------------------------------------------- */

ComputeCNPAtomOMP::ComputeCNPAtomOMP(LAMMPS *lmp, int narg, char **arg) :
    ComputeCNPAtom(lmp, narg, arg)
{
}

/* ----------------------------------------------------------------------
   threaded variant of ComputeCNPAtom::compute_peratom().  Two threaded
   phases with an implicit barrier between them: phase 1 fills the per-atom
   nearest[]/nnearest[] lists, phase 2 computes cnpv[i] reading neighbors'
   rows with thread-private common[]/onenearest[] scratch.  Each atom writes
   only its own cnpv[i], so the result is bit-identical to the serial
   compute; the "too many neighbors" counters use omp reductions and the
   MPI_Allreduce calls stay collective outside the parallel regions.
------------------------------------------------------------------------- */

void ComputeCNPAtomOMP::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow arrays if necessary

  if (atom->nmax > nmax) {
    memory->destroy(nearest);
    memory->destroy(nnearest);
    memory->destroy(cnpv);
    nmax = atom->nmax;
    memory->create(nearest, nmax, MAXNEAR, "cnp:nearest");
    memory->create(nnearest, nmax, "cnp:nnearest");
    memory->create(cnpv, nmax, "cnp:cnp_cnpv");
    vector_atom = cnpv;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  // find the neighbors of each atom within cutoff using full neighbor list
  // nearest[] = atom indices of nearest neighbors, up to MAXNEAR
  // do this for all atoms, not just compute group
  // since CNP calculation requires neighbors of neighbors

  const double *const *const x = atom->x;
  const int *const mask = atom->mask;
  const int nlocal = atom->nlocal;

  int nerror = 0;

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic) reduction(+ : nerror)
#endif
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int *const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    int n = 0;
    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj] & NEIGHMASK;
      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq = delx * delx + dely * dely + delz * delz;
      if (rsq < cutsq) {
        if (n < MAXNEAR)
          nearest[i][n++] = j;
        else {
          nerror++;
          break;
        }
      }
    }
    nnearest[i] = n;
  }

  // warning message

  int nerrorall;
  MPI_Allreduce(&nerror, &nerrorall, 1, MPI_INT, MPI_SUM, world);
  if (nerrorall && comm->me == 0)
    error->warning(FLERR, "Too many neighbors in CNP for {} atoms", nerrorall);

  // compute CNP value for each atom in group

  nerror = 0;

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic) reduction(+ : nerror)
#endif
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    // reset cnpv
    cnpv[i] = 0.0;

    // skip computation of cnpv for atoms outside the compute group

    if (!(mask[i] & groupbit)) continue;

    int onenearest[MAXNEAR];
    int common[MAXCOMMON];

    // loop over nearest neighbors of I to build cnp data structure
    //  cnp[k][NCOMMON] = # of common neighbors of I with each of its neighbors

    for (int m = 0; m < nnearest[i]; m++) {
      const int jatom = nearest[i][m];
      const double xjtmp = x[jatom][0];
      const double yjtmp = x[jatom][1];
      const double zjtmp = x[jatom][2];

      // common = list of neighbors common to atom I and atom J
      // if J is an owned atom, use its near neighbor list to find them
      // if J is a ghost atom, use full neighbor list of I to find them
      // in latter case, must exclude J from I's neighbor list

      int firstflag, ncommon;

      // find common neighbors of i and j using near neighbor list
      if (jatom < nlocal) {
        firstflag = 1;
        ncommon = 0;
        for (int inear = 0; inear < nnearest[i]; inear++)
          for (int jnear = 0; jnear < nnearest[jatom]; jnear++)
            if (nearest[i][inear] == nearest[jatom][jnear]) {
              if (ncommon < MAXCOMMON)
                common[ncommon++] = nearest[i][inear];
              else if (firstflag) {
                nerror++;
                firstflag = 0;
              }
            }

        // find common neighbors of i and j using full neighbor list
      } else {
        const int *const jlist = firstneigh[i];
        const int jnum = numneigh[i];

        int n = 0;
        for (int kk = 0; kk < jnum; kk++) {
          const int k = jlist[kk] & NEIGHMASK;
          if (k == jatom) continue;

          const double delx = xjtmp - x[k][0];
          const double dely = yjtmp - x[k][1];
          const double delz = zjtmp - x[k][2];
          const double rsq = delx * delx + dely * dely + delz * delz;
          if (rsq < cutsq) {
            if (n < MAXNEAR)
              onenearest[n++] = k;
            else
              break;
          }
        }

        firstflag = 1;
        ncommon = 0;
        for (int inear = 0; inear < nnearest[i]; inear++)
          for (int jnear = 0; (jnear < n) && (n < MAXNEAR); jnear++)
            if (nearest[i][inear] == onenearest[jnear]) {
              if (ncommon < MAXCOMMON)
                common[ncommon++] = nearest[i][inear];
              else if (firstflag) {
                nerror++;
                firstflag = 0;
              }
            }
      }

      // Calculate and update sum |Rik+Rjk|^2
      double rjkx = 0.0;
      double rjky = 0.0;
      double rjkz = 0.0;
      for (int kk = 0; kk < ncommon; kk++) {
        const int k = common[kk];
        rjkx += 2.0 * x[k][0] - xjtmp - xtmp;
        rjky += 2.0 * x[k][1] - yjtmp - ytmp;
        rjkz += 2.0 * x[k][2] - zjtmp - ztmp;
      }
      // update cnpv with summed (valuejk)
      cnpv[i] += rjkx * rjkx + rjky * rjky + rjkz * rjkz;

      // end of loop over j atoms
    }

    // normalize cnp by the number of nearest neighbors
    cnpv[i] = cnpv[i] / nnearest[i];

    // end of loop over i atoms
  }

  // warning message
  MPI_Allreduce(&nerror, &nerrorall, 1, MPI_INT, MPI_SUM, world);
  if (nerrorall && comm->me == 0)
    error->warning(FLERR, "Too many common neighbors in CNP {} times", nerrorall);
}
