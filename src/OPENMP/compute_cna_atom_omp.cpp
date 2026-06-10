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

#include "compute_cna_atom_omp.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;

// these mirror the file-static constants in compute_cna_atom.cpp

static constexpr int MAXNEAR = 16;
static constexpr int MAXCOMMON = 8;

enum { UNKNOWN, FCC, HCP, BCC, ICOS, OTHER };
enum { NCOMMON, NBOND, MAXBOND, MINBOND };

/* ---------------------------------------------------------------------- */

ComputeCNAAtomOMP::ComputeCNAAtomOMP(LAMMPS *lmp, int narg, char **arg) :
    ComputeCNAAtom(lmp, narg, arg)
{
}

/* ----------------------------------------------------------------------
   threaded variant of ComputeCNAAtom::compute_peratom().  Two threaded
   phases: phase 1 fills the per-atom nearest[]/nnearest[] lists (each atom
   writes only its own row), phase 2 classifies each atom reading its
   neighbors' rows.  The implicit barrier between the two omp-for regions
   guarantees all rows are built before classification.  The classification
   scratch (cna/onenearest/common/bonds) is thread-private.  Each atom writes
   only its own pattern[i], so the result is bit-identical to the serial
   compute.  The two "too many neighbors" counters are summed across threads
   and then MPI-reduced (collective calls stay outside the parallel region).
------------------------------------------------------------------------- */

void ComputeCNAAtomOMP::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow arrays if necessary

  if (atom->nmax > nmax) {
    memory->destroy(nearest);
    memory->destroy(nnearest);
    memory->destroy(pattern);
    nmax = atom->nmax;

    memory->create(nearest, nmax, MAXNEAR, "cna:nearest");
    memory->create(nnearest, nmax, "cna:nnearest");
    memory->create(pattern, nmax, "cna:cna_pattern");
    vector_atom = pattern;
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
  // since CNA calculation requires neighbors of neighbors

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
    error->warning(FLERR, "Too many neighbors in CNA for {} atoms", nerrorall);

  // compute CNA for each atom in group
  // only performed if # of nearest neighbors = 12 or 14 (fcc,hcp)

  nerror = 0;

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic) reduction(+ : nerror)
#endif
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];

    if (!(mask[i] & groupbit)) {
      pattern[i] = UNKNOWN;
      continue;
    }

    if (nnearest[i] != 12 && nnearest[i] != 14) {
      pattern[i] = OTHER;
      continue;
    }

    // per-thread scratch

    int cna[MAXNEAR][4], onenearest[MAXNEAR];
    int common[MAXCOMMON], bonds[MAXCOMMON];

    // loop over near neighbors of I to build cna data structure
    // cna[k][NCOMMON] = # of common neighbors of I with each of its neighs
    // cna[k][NBONDS] = # of bonds between those common neighbors
    // cna[k][MAXBOND] = max # of bonds of any common neighbor
    // cna[k][MINBOND] = min # of bonds of any common neighbor

    for (int m = 0; m < nnearest[i]; m++) {
      const int jatom = nearest[i][m];

      // common = list of neighbors common to atom I and atom J
      // if J is an owned atom, use its near neighbor list to find them
      // if J is a ghost atom, use full neighbor list of I to find them
      // in latter case, must exclude J from I's neighbor list

      int firstflag, ncommon;
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

      } else {
        const double xtmp = x[jatom][0];
        const double ytmp = x[jatom][1];
        const double ztmp = x[jatom][2];
        const int *const jlist = firstneigh[i];
        const int jnum = numneigh[i];

        int n = 0;
        for (int kk = 0; kk < jnum; kk++) {
          const int k = jlist[kk] & NEIGHMASK;
          if (k == jatom) continue;

          const double delx = xtmp - x[k][0];
          const double dely = ytmp - x[k][1];
          const double delz = ztmp - x[k][2];
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

      cna[m][NCOMMON] = ncommon;

      // calculate total # of bonds between common neighbor atoms
      // also max and min # of common atoms any common atom is bonded to
      // bond = pair of atoms within cutoff

      for (int n = 0; n < ncommon; n++) bonds[n] = 0;

      int nbonds = 0;
      for (int jj = 0; jj < ncommon - 1; jj++) {
        const int j = common[jj];
        const double xtmp = x[j][0];
        const double ytmp = x[j][1];
        const double ztmp = x[j][2];
        for (int kk = jj + 1; kk < ncommon; kk++) {
          const int k = common[kk];
          const double delx = xtmp - x[k][0];
          const double dely = ytmp - x[k][1];
          const double delz = ztmp - x[k][2];
          const double rsq = delx * delx + dely * dely + delz * delz;
          if (rsq < cutsq) {
            nbonds++;
            bonds[jj]++;
            bonds[kk]++;
          }
        }
      }

      cna[m][NBOND] = nbonds;

      int maxbonds = 0;
      int minbonds = MAXCOMMON;
      for (int n = 0; n < ncommon; n++) {
        maxbonds = MAX(bonds[n], maxbonds);
        minbonds = MIN(bonds[n], minbonds);
      }
      cna[m][MAXBOND] = maxbonds;
      cna[m][MINBOND] = minbonds;
    }

    // detect CNA pattern of the atom

    int nfcc, nhcp, nbcc4, nbcc6, nico;
    nfcc = nhcp = nbcc4 = nbcc6 = nico = 0;
    pattern[i] = OTHER;

    if (nnearest[i] == 12) {
      for (int inear = 0; inear < 12; inear++) {
        const int cj = cna[inear][NCOMMON];
        const int ck = cna[inear][NBOND];
        const int cl = cna[inear][MAXBOND];
        const int cm = cna[inear][MINBOND];
        if (cj == 4 && ck == 2 && cl == 1 && cm == 1)
          nfcc++;
        else if (cj == 4 && ck == 2 && cl == 2 && cm == 0)
          nhcp++;
        else if (cj == 5 && ck == 5 && cl == 2 && cm == 2)
          nico++;
      }
      if (nfcc == 12)
        pattern[i] = FCC;
      else if (nfcc == 6 && nhcp == 6)
        pattern[i] = HCP;
      else if (nico == 12)
        pattern[i] = ICOS;

    } else if (nnearest[i] == 14) {
      for (int inear = 0; inear < 14; inear++) {
        const int cj = cna[inear][NCOMMON];
        const int ck = cna[inear][NBOND];
        const int cl = cna[inear][MAXBOND];
        const int cm = cna[inear][MINBOND];
        if (cj == 4 && ck == 4 && cl == 2 && cm == 2)
          nbcc4++;
        else if (cj == 6 && ck == 6 && cl == 2 && cm == 2)
          nbcc6++;
      }
      if (nbcc4 == 6 && nbcc6 == 8) pattern[i] = BCC;
    }
  }

  // warning message

  MPI_Allreduce(&nerror, &nerrorall, 1, MPI_INT, MPI_SUM, world);
  if (nerrorall && comm->me == 0)
    error->warning(FLERR, "Too many common neighbors in CNA: {}x", nerrorall);
}
