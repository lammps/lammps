// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Common Neighbor Parameter as proposed in:
   Tsuzuki, Branicio, Rino, Comput Phys Comm, 177, 518 (2007)
   Cite: https://doi.org/10.1063/1.2197987

------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Paulo Branicio (University of Southern California)
   branicio@usc.edu
------------------------------------------------------------------------- */

#include "compute_cnp_atom.h"

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
#include <cstring>

using namespace LAMMPS_NS;

//define maximum values
static constexpr int MAXNEAR = 24;
static constexpr int MAXCOMMON = 12;

enum{NCOMMON};

/* ---------------------------------------------------------------------- */

ComputeCNPAtom::ComputeCNPAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  list(nullptr), nearest(nullptr), nnearest(nullptr), cnpv(nullptr)
{
  if (narg != 4) error->all(FLERR,"Illegal compute cnp/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  double cutoff = utils::numeric(FLERR,arg[3],false,lmp);
  if (cutoff < 0.0) error->all(FLERR,"Illegal compute cnp/atom command");
  cutsq = cutoff*cutoff;

  // apply check for single type atoms in compute group
  int lasttype = -1;
  int n = -1;
  for (int i=0; i < atom->nlocal; ++i) {
    if (atom->mask[i] & groupbit) {
      if (lasttype != atom->type[i]) {
        lasttype = atom->type[i];
        ++n;
      }
    }
  }
  int all_n = 0;
  MPI_Allreduce(&n,&all_n,1,MPI_INT,MPI_MAX,world);
  if (all_n > 0)
    error->warning(FLERR,"Compute cnp/atom requested on multi-type system");

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeCNPAtom::~ComputeCNPAtom()
{
  memory->destroy(nearest);
  memory->destroy(nnearest);
  memory->destroy(cnpv);
}

/* ---------------------------------------------------------------------- */

void ComputeCNPAtom::init()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Compute cnp/atom requires a pair style be defined");

  if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR,"Compute cnp/atom cutoff is longer than pairwise cutoff");

  if (2.0*sqrt(cutsq) > force->pair->cutforce + neighbor->skin &&
      comm->me == 0)
    error->warning(FLERR,"Compute cnp/atom cutoff may be too large to find "
                   "ghost atom neighbors");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"cnp/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute cnp/atom defined");

  // need an occasional full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);
}

/* ---------------------------------------------------------------------- */

void ComputeCNPAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeCNPAtom::compute_peratom()
{
  int i,j,k,ii,jj,kk,m,n,inum,jnum,inear,jnear;
  int firstflag,ncommon;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int onenearest[MAXNEAR];
  int common[MAXCOMMON];
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double xjtmp,yjtmp,zjtmp,rjkx,rjky,rjkz;

  invoked_peratom = update->ntimestep;

  // grow arrays if necessary

  if (atom->nmax > nmax) {
    memory->destroy(nearest);
    memory->destroy(nnearest);
    memory->destroy(cnpv);
    nmax = atom->nmax;
    memory->create(nearest,nmax,MAXNEAR,"cnp:nearest");
    memory->create(nnearest,nmax,"cnp:nnearest");
    memory->create(cnpv,nmax,"cnp:cnp_cnpv");
    vector_atom = cnpv;
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
  // since CNP calculation requires neighbors of neighbors

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
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < cutsq) {
        if (n < MAXNEAR) nearest[i][n++] = j;
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
  MPI_Allreduce(&nerror,&nerrorall,1,MPI_INT,MPI_SUM,world);
  if (nerrorall && comm->me == 0)
    error->warning(FLERR,"Too many neighbors in CNP for {} atoms",nerrorall);

  // compute CNP value for each atom in group
  // only performed if # of nearest neighbors = 12 or 14 (fcc,hcp)

  nerror = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    // reset cnpv
    cnpv[i] = 0.0;

    // skip computation of cnpv for atoms outside the compute group

    if (!(mask[i] & groupbit)) continue;

    // loop over nearest neighbors of I to build cnp data structure
    //  cnp[k][NCOMMON] = # of common neighbors of I with each of its neighbors
    for (m = 0; m < nnearest[i]; m++) {
      j = nearest[i][m];
      xjtmp = x[j][0];
      yjtmp = x[j][1];
      zjtmp = x[j][2];

      // common = list of neighbors common to atom I and atom J
      // if J is an owned atom, use its near neighbor list to find them
      // if J is a ghost atom, use full neighbor list of I to find them
      // in latter case, must exclude J from I's neighbor list

      // find common neighbors of i and j using near neighbor list
      if (j < nlocal) {
        firstflag = 1;
        ncommon = 0;
        for (inear = 0; inear < nnearest[i]; inear++)
          for (jnear = 0; jnear < nnearest[j]; jnear++)
            if (nearest[i][inear] == nearest[j][jnear]) {
              if (ncommon < MAXCOMMON) common[ncommon++] = nearest[i][inear];
              else if (firstflag) {
                nerror++;
                firstflag = 0;
              }
            }

      // find common neighbors of i and j using full neighbor list
      } else {
        jlist = firstneigh[i];
        jnum = numneigh[i];

        n = 0;
        for (kk = 0; kk < jnum; kk++) {
          k = jlist[kk];
          k &= NEIGHMASK;
          if (k == j) continue;

          delx = xjtmp - x[k][0];
          dely = yjtmp - x[k][1];
          delz = zjtmp - x[k][2];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutsq) {
            if (n < MAXNEAR) onenearest[n++] = k;
            else break;
          }
        }

        firstflag = 1;
        ncommon = 0;
        for (inear = 0; inear < nnearest[i]; inear++)
          for (jnear = 0; (jnear < n) && (n < MAXNEAR); jnear++)
            if (nearest[i][inear] == onenearest[jnear]) {
              if (ncommon < MAXCOMMON) common[ncommon++] = nearest[i][inear];
              else if (firstflag) {
                nerror++;
                firstflag = 0;
              }
            }
      }

      // Calculate and update sum |Rik+Rjk|ˆ2
      rjkx = 0.0;
      rjky = 0.0;
      rjkz = 0.0;
      for (kk = 0; kk < ncommon; kk++) {
        k = common[kk];
        rjkx += 2.0*x[k][0] - xjtmp - xtmp;
        rjky += 2.0*x[k][1] - yjtmp - ytmp;
        rjkz += 2.0*x[k][2] - zjtmp - ztmp;
      }
      // update cnpv with summed (valuejk)
      cnpv[i] += rjkx*rjkx + rjky*rjky + rjkz*rjkz;

    // end of loop over j atoms
    }

    // normalize cnp by the number of nearest neighbors
    cnpv[i] = cnpv[i] / nnearest[i];

  // end of loop over i atoms
  }

  // warning message
  MPI_Allreduce(&nerror,&nerrorall,1,MPI_INT,MPI_SUM,world);
  if (nerrorall && comm->me == 0)
    error->warning(FLERR,"Too many common neighbors in CNP {} times",nerrorall);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeCNPAtom::memory_usage()
{
  double bytes = (double)nmax * sizeof(int);
  bytes += (double)nmax * MAXNEAR * sizeof(int);
  bytes += (double)nmax * sizeof(double);
  return bytes;
}
