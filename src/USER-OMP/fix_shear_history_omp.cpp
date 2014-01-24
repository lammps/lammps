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

#include "string.h"
#include "stdio.h"
#include "fix_shear_history_omp.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "force.h"
#include "pair.h"
#include "update.h"
#include "modify.h"
#include "error.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;
using namespace FixConst;

/* ----------------------------------------------------------------------
   copy shear partner info from neighbor lists to atom arrays
   so can be exchanged with atoms
------------------------------------------------------------------------- */

void FixShearHistoryOMP::pre_exchange()
{

  const int nlocal = atom->nlocal;
  const int nghost = atom->nghost;
  const int nall = nlocal + nghost;
  const int nthreads = comm->nthreads;
  maxtouch = 0;
  
#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {

#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif

    // each thread works on a fixed chunk of local and ghost atoms.
    const int ldelta = 1 + nlocal/nthreads;
    const int lfrom = tid*ldelta;
    const int lmax = lfrom +ldelta;
    const int lto = (lmax > nlocal) ? nlocal : lmax;

    int i,j,ii,jj,m,n,inum,jnum;
    int *ilist,*jlist,*numneigh,**firstneigh;
    int *touch,**firsttouch;
    double *shear,*allshear,**firstshear;

    // zero npartners for all current atoms and
    // clear page data structures for this thread

    for (i = lfrom; i < lto; i++) npartner[i] = 0;

    MyPage <tagint> &ipg = ipage[tid];
    MyPage <double[3]> &dpg = dpage[tid];
    ipg.reset();
    dpg.reset();

    // 1st loop over neighbor list
    // calculate nparter for each owned atom

    tagint *tag = atom->tag;

    NeighList *list = pair->list;
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    firsttouch = list->listgranhistory->firstneigh;
    firstshear = list->listgranhistory->firstdouble;

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      touch = firsttouch[i];

      for (jj = 0; jj < jnum; jj++) {
        if (touch[jj]) {
          if ((i >= lfrom) && (i < lto))
            npartner[i]++;

          j = jlist[jj];
          j &= NEIGHMASK;
          if ((j >= lfrom) && (j < lto))
            npartner[j]++;
        }
      }
    }

    // get page chunks to store atom IDs and shear history for my atoms

    for (ii = lfrom; ii < lto; ii++) {
      i = ilist[ii];
      n = npartner[i];
      partner[i] = ipg.get(n);
      shearpartner[i] = dpg.get(n);
      if (partner[i] == NULL || shearpartner[i] == NULL)
        error->one(FLERR,"Shear history overflow, boost neigh_modify one");
    }

    // 2nd loop over neighbor list
    // store atom IDs and shear history for my atoms
    // re-zero npartner to use as counter for all my atoms

    for (i = lfrom; i < lto; i++) npartner[i] = 0;

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      jlist = firstneigh[i];
      allshear = firstshear[i];
      jnum = numneigh[i];
      touch = firsttouch[i];

      for (jj = 0; jj < jnum; jj++) {
        if (touch[jj]) {
          shear = &allshear[3*jj];
          j = jlist[jj];
          j &= NEIGHMASK;

          if ((i >= lfrom) && (i < lto)) {
            m = npartner[i];
            partner[i][m] = tag[j];
            shearpartner[i][m][0] = shear[0];
            shearpartner[i][m][1] = shear[1];
            shearpartner[i][m][2] = shear[2];
            npartner[i]++;
          }

          if ((j >= lfrom) && (j < lto)) {
            m = npartner[j];
            partner[j][m] = tag[i];
            shearpartner[j][m][0] = -shear[0];
            shearpartner[j][m][1] = -shear[1];
            shearpartner[j][m][2] = -shear[2];
            npartner[j]++;
          }
        }
      }
    }

    // set maxtouch = max # of partners of any owned atom
    m = 0;
    for (i = lfrom; i < lto; i++)
      m = MAX(m,npartner[i]);

#if defined(_OPENMP)
#pragma omp critical
#endif
    maxtouch = MAX(m,maxtouch);
  }
}
