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
   Contributing author: Wan Liang (Chinese Academy of Sciences)
------------------------------------------------------------------------- */

#include "compute_cna_atom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

static constexpr int MAXNEAR = 16;
static constexpr int MAXCOMMON = 8;

enum { UNKNOWN, FCC, HCP, BCC, ICOS, OTHER };
enum { NCOMMON, NBOND, MAXBOND, MINBOND };

/* ---------------------------------------------------------------------- */

ComputeCNAAtom::ComputeCNAAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), list(nullptr), nearest(nullptr), nnearest(nullptr), pattern(nullptr)
{
  if (narg != 4) error->all(FLERR, "Illegal compute cna/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  double cutoff = utils::numeric(FLERR, arg[3], false, lmp);
  if (cutoff < 0.0) error->all(FLERR, "Illegal compute cna/atom command");
  cutsq = cutoff * cutoff;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeCNAAtom::~ComputeCNAAtom()
{
  memory->destroy(nearest);
  memory->destroy(nnearest);
  memory->destroy(pattern);
}

/* ---------------------------------------------------------------------- */

void ComputeCNAAtom::init()
{
  if (force->pair == nullptr)
    error->all(FLERR, "Compute cna/atom requires a pair style be defined");
  if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR, "Compute cna/atom cutoff is longer than pairwise cutoff");

  // cannot use neighbor->cutneighmax b/c neighbor has not yet been init

  if ((2.0 * sqrt(cutsq)) > (force->pair->cutforce + neighbor->skin) && (comm->me == 0))
    error->warning(FLERR, "Compute cna/atom cutoff may be too large to find ghost atom neighbors");

  // need an occasional full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);

  if (modify->get_compute_by_style(style).size() > 1)
    if (comm->me == 0) error->warning(FLERR, "More than one compute {}", style);
}

/* ---------------------------------------------------------------------- */

void ComputeCNAAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeCNAAtom::compute_peratom()
{
  int i, j, k, ii, jj, kk, m, n, inum, jnum, inear, jnear;
  int firstflag, ncommon, nbonds, maxbonds, minbonds;
  int nfcc, nhcp, nbcc4, nbcc6, nico, cj, ck, cl, cm;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int cna[MAXNEAR][4], onenearest[MAXNEAR];
  int common[MAXCOMMON], bonds[MAXCOMMON];
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;

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

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // find the neighbors of each atom within cutoff using full neighbor list
  // nearest[] = atom indices of nearest neighbors, up to MAXNEAR
  // do this for all atoms, not just compute group
  // since CNA calculation requires neighbors of neighbors

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int nerror = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    n = 0;
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
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
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    if (!(mask[i] & groupbit)) {
      pattern[i] = UNKNOWN;
      continue;
    }

    if (nnearest[i] != 12 && nnearest[i] != 14) {
      pattern[i] = OTHER;
      continue;
    }

    // loop over near neighbors of I to build cna data structure
    // cna[k][NCOMMON] = # of common neighbors of I with each of its neighs
    // cna[k][NBONDS] = # of bonds between those common neighbors
    // cna[k][MAXBOND] = max # of bonds of any common neighbor
    // cna[k][MINBOND] = min # of bonds of any common neighbor

    for (m = 0; m < nnearest[i]; m++) {
      j = nearest[i][m];

      // common = list of neighbors common to atom I and atom J
      // if J is an owned atom, use its near neighbor list to find them
      // if J is a ghost atom, use full neighbor list of I to find them
      // in latter case, must exclude J from I's neighbor list

      if (j < nlocal) {
        firstflag = 1;
        ncommon = 0;
        for (inear = 0; inear < nnearest[i]; inear++)
          for (jnear = 0; jnear < nnearest[j]; jnear++)
            if (nearest[i][inear] == nearest[j][jnear]) {
              if (ncommon < MAXCOMMON)
                common[ncommon++] = nearest[i][inear];
              else if (firstflag) {
                nerror++;
                firstflag = 0;
              }
            }

      } else {
        xtmp = x[j][0];
        ytmp = x[j][1];
        ztmp = x[j][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        n = 0;
        for (kk = 0; kk < jnum; kk++) {
          k = jlist[kk];
          k &= NEIGHMASK;
          if (k == j) continue;

          delx = xtmp - x[k][0];
          dely = ytmp - x[k][1];
          delz = ztmp - x[k][2];
          rsq = delx * delx + dely * dely + delz * delz;
          if (rsq < cutsq) {
            if (n < MAXNEAR)
              onenearest[n++] = k;
            else
              break;
          }
        }

        firstflag = 1;
        ncommon = 0;
        for (inear = 0; inear < nnearest[i]; inear++)
          for (jnear = 0; (jnear < n) && (n < MAXNEAR); jnear++)
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

      for (n = 0; n < ncommon; n++) bonds[n] = 0;

      nbonds = 0;
      for (jj = 0; jj < ncommon - 1; jj++) {
        j = common[jj];
        xtmp = x[j][0];
        ytmp = x[j][1];
        ztmp = x[j][2];
        for (kk = jj + 1; kk < ncommon; kk++) {
          k = common[kk];
          delx = xtmp - x[k][0];
          dely = ytmp - x[k][1];
          delz = ztmp - x[k][2];
          rsq = delx * delx + dely * dely + delz * delz;
          if (rsq < cutsq) {
            nbonds++;
            bonds[jj]++;
            bonds[kk]++;
          }
        }
      }

      cna[m][NBOND] = nbonds;

      maxbonds = 0;
      minbonds = MAXCOMMON;
      for (n = 0; n < ncommon; n++) {
        maxbonds = MAX(bonds[n], maxbonds);
        minbonds = MIN(bonds[n], minbonds);
      }
      cna[m][MAXBOND] = maxbonds;
      cna[m][MINBOND] = minbonds;
    }

    // detect CNA pattern of the atom

    nfcc = nhcp = nbcc4 = nbcc6 = nico = 0;
    pattern[i] = OTHER;

    if (nnearest[i] == 12) {
      for (inear = 0; inear < 12; inear++) {
        cj = cna[inear][NCOMMON];
        ck = cna[inear][NBOND];
        cl = cna[inear][MAXBOND];
        cm = cna[inear][MINBOND];
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
      for (inear = 0; inear < 14; inear++) {
        cj = cna[inear][NCOMMON];
        ck = cna[inear][NBOND];
        cl = cna[inear][MAXBOND];
        cm = cna[inear][MINBOND];
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

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeCNAAtom::memory_usage()
{
  double bytes = (double) nmax * sizeof(int);
  bytes += (double) nmax * MAXNEAR * sizeof(int);
  bytes += (double) nmax * sizeof(double);
  return bytes;
}
