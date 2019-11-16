/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_snap.h"
#include <cstring>
#include <cstdlib>
#include "sna.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

ComputeSnap::ComputeSnap(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), cutsq(NULL), list(NULL), snap(NULL),
  radelem(NULL), wjelem(NULL)
{
  double rfac0, rmin0;
  int twojmax, switchflag, bzeroflag;
  radelem = NULL;
  wjelem = NULL;

  int ntypes = atom->ntypes;
  int nargmin = 6+2*ntypes;

  if (narg < nargmin) error->all(FLERR,"Illegal compute snap command");

  // default values

  rmin0 = 0.0;
  switchflag = 1;
  bzeroflag = 1;
  quadraticflag = 0;

  // process required arguments

  memory->create(radelem,ntypes+1,"snap:radelem"); // offset by 1 to match up with types
  memory->create(wjelem,ntypes+1,"snap:wjelem");
  rcutfac = atof(arg[3]);
  rfac0 = atof(arg[4]);
  twojmax = atoi(arg[5]);
  for(int i = 0; i < ntypes; i++)
    radelem[i+1] = atof(arg[6+i]);
  for(int i = 0; i < ntypes; i++)
    wjelem[i+1] = atof(arg[6+ntypes+i]);

  // construct cutsq

  double cut;
  cutmax = 0.0;
  memory->create(cutsq,ntypes+1,ntypes+1,"snap:cutsq");
  for(int i = 1; i <= ntypes; i++) {
    cut = 2.0*radelem[i]*rcutfac;
    if (cut > cutmax) cutmax = cut;
    cutsq[i][i] = cut*cut;
    for(int j = i+1; j <= ntypes; j++) {
      cut = (radelem[i]+radelem[j])*rcutfac;
      cutsq[i][j] = cutsq[j][i] = cut*cut;
    }
  }

  // process optional args

  int iarg = nargmin;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"rmin0") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute snap command");
      rmin0 = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"bzeroflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute snap command");
      bzeroflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"switchflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute snap command");
      switchflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"quadraticflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute snap command");
      quadraticflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal compute snap command");
  }

  snaptr = new SNA(lmp,rfac0,twojmax,
                   rmin0,switchflag,bzeroflag);

  ncoeff = snaptr->ncoeff;
  nperdim = ncoeff;
  if (quadraticflag) nperdim += (ncoeff*(ncoeff+1))/2;
  yoffset = nperdim;
  zoffset = 2*nperdim;
  size_peratom_cols = 3*nperdim*atom->ntypes;
  comm_reverse = size_peratom_cols;

  nmax = 0;
  snap = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeSnap::~ComputeSnap()
{
  memory->destroy(snap);
  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(cutsq);
  delete snaptr;
}

/* ---------------------------------------------------------------------- */

void ComputeSnap::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute snap requires a pair style be defined");

  if (cutmax > force->pair->cutforce)
    error->all(FLERR,"Compute snap cutoff is longer than pairwise cutoff");

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"snap") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute snap");
  snaptr->init();
}

/* ---------------------------------------------------------------------- */

void ComputeSnap::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeSnap::compute()
{
  int ntotal = atom->nlocal + atom->nghost;

  invoked_peratom = update->ntimestep;

  // grow snap array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(snap);
    nmax = atom->nmax;
    memory->create(snap,nmax,size_peratom_cols,
                   "snap:snap");
    array = snap;
  }

  // clear local array

  for (int i = 0; i < ntotal; i++)
    for (int icoeff = 0; icoeff < size_peratom_cols; icoeff++) {
      snap[i][icoeff] = 0.0;
    }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  const int inum = list->inum;
  const int* const ilist = list->ilist;
  const int* const numneigh = list->numneigh;
  int** const firstneigh = list->firstneigh;
  int * const type = atom->type;

  // compute sna derivatives for each atom in group
  // use full neighbor list to count atoms less than cutoff

  double** const x = atom->x;
  const int* const mask = atom->mask;

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) {

      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      const int itype = type[i];
      const double radi = radelem[itype];
      const int* const jlist = firstneigh[i];
      const int jnum = numneigh[i];

      // const int typeoffset = threencoeff*(atom->type[i]-1);
      // const int quadraticoffset = threencoeff*atom->ntypes +
      //   threencoeffq*(atom->type[i]-1);
      const int typeoffset = 3*nperdim*(atom->type[i]-1);

      // insure rij, inside, and typej  are of size jnum

      snaptr->grow_rij(jnum);

      // rij[][3] = displacements between atom I and those neighbors
      // inside = indices of neighbors of I within cutoff
      // typej = types of neighbors of I within cutoff
      // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

      int ninside = 0;
      for (int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        j &= NEIGHMASK;

        const double delx = x[j][0] - xtmp;
        const double dely = x[j][1] - ytmp;
        const double delz = x[j][2] - ztmp;
        const double rsq = delx*delx + dely*dely + delz*delz;
        int jtype = type[j];
        if (rsq < cutsq[itype][jtype]&&rsq>1e-20) {
          snaptr->rij[ninside][0] = delx;
          snaptr->rij[ninside][1] = dely;
          snaptr->rij[ninside][2] = delz;
          snaptr->inside[ninside] = j;
          snaptr->wj[ninside] = wjelem[jtype];
          snaptr->rcutij[ninside] = (radi+radelem[jtype])*rcutfac;
          ninside++;
        }
      }

      snaptr->compute_ui(ninside);
      snaptr->compute_zi();
      if (quadraticflag) {
        snaptr->compute_bi();
      }

      for (int jj = 0; jj < ninside; jj++) {
        const int j = snaptr->inside[jj];
        snaptr->compute_duidrj(snaptr->rij[jj],
                                    snaptr->wj[jj],
                                    snaptr->rcutij[jj],jj);
        snaptr->compute_dbidrj();

        // Accumulate -dBi/dRi, -dBi/dRj

        double *snapi = snap[i]+typeoffset;
        double *snapj = snap[j]+typeoffset;

        for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
          snapi[icoeff] += snaptr->dblist[icoeff][0];
          snapi[icoeff+yoffset] += snaptr->dblist[icoeff][1];
          snapi[icoeff+zoffset] += snaptr->dblist[icoeff][2];
          snapj[icoeff] -= snaptr->dblist[icoeff][0];
          snapj[icoeff+yoffset] -= snaptr->dblist[icoeff][1];
          snapj[icoeff+zoffset] -= snaptr->dblist[icoeff][2];
        }

        if (quadraticflag) {
          const int quadraticoffset = ncoeff;
          snapi += quadraticoffset;
          snapj += quadraticoffset;
          int ncount = 0;
          for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
            double bi = snaptr->blist[icoeff];
            double bix = snaptr->dblist[icoeff][0];
            double biy = snaptr->dblist[icoeff][1];
            double biz = snaptr->dblist[icoeff][2];

            // diagonal elements of quadratic matrix

            double dbxtmp = bi*bix;
            double dbytmp = bi*biy;
            double dbztmp = bi*biz;

            snapi[ncount] +=         dbxtmp;
            snapi[ncount+yoffset] += dbytmp;
            snapi[ncount+zoffset] += dbztmp;
            snapj[ncount] -=         dbxtmp;
            snapj[ncount+yoffset] -= dbytmp;
            snapj[ncount+zoffset] -= dbztmp;
            ncount++;

            // upper-triangular elements of quadratic matrix

            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
              double dbxtmp = bi*snaptr->dblist[jcoeff][0]
                + bix*snaptr->blist[jcoeff];
              double dbytmp = bi*snaptr->dblist[jcoeff][1]
                + biy*snaptr->blist[jcoeff];
              double dbztmp = bi*snaptr->dblist[jcoeff][2]
                + biz*snaptr->blist[jcoeff];

              snapi[ncount] +=         dbxtmp;
              snapi[ncount+yoffset] += dbytmp;
              snapi[ncount+zoffset] += dbztmp;
              snapj[ncount] -=         dbxtmp;
              snapj[ncount+yoffset] -= dbytmp;
              snapj[ncount+zoffset] -= dbztmp;
              ncount++;
            }
          }
        }
      }
    }
  }

  // communicate snap contributions between neighbor procs

  comm->reverse_comm_compute(this);

}

/* ---------------------------------------------------------------------- */

int ComputeSnap::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last,icoeff;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    for (icoeff = 0; icoeff < size_peratom_cols; icoeff++)
      buf[m++] = snap[i][icoeff];
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeSnap::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m,icoeff;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (icoeff = 0; icoeff < size_peratom_cols; icoeff++)
      snap[j][icoeff] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double ComputeSnap::memory_usage()
{

  double bytes = nmax*size_peratom_cols * sizeof(double); // snap
  bytes += snaptr->memory_usage();                        // SNA object

  return bytes;
}
