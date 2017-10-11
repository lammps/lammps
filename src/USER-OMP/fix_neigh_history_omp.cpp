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

#include <string.h>
#include <stdio.h>
#include "fix_neigh_history_omp.h"
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


FixNeighHistoryOMP::FixNeighHistoryOMP(class LAMMPS *lmp,int narg,char **argv)
  : FixNeighHistory(lmp,narg,argv) {

  if (onesided)
    error->all(FLERR,"tri/lj and line/lj are not supported by USER-OMP");
  if (!newton_pair)
    error->all(FLERR,"Newton off for granular is not supported by USER-OMP");
}


/* ----------------------------------------------------------------------
 copy partner info from neighbor data structs (NDS) to atom arrays
 should be called whenever NDS store current history info
   and need to transfer the info to owned atoms
 e.g. when atoms migrate to new procs, new neigh list built, or between runs
   when atoms may be added or deleted (NDS becomes out-of-date)
 the next post_neighbor() will put this info back into new NDS
 called during run before atom exchanges, including for restart files
 called at end of run via post_run()
 do not call during setup of run (setup_pre_exchange)
   b/c there is no guarantee of a current NDS (even on continued run)
 if run command does a 2nd run with pre = no, then no neigh list
   will be built, but old neigh list will still have the info
 
 the USER-OMP version only supports newton on
------------------------------------------------------------------------- */

void FixNeighHistoryOMP::pre_exchange()
{
  const int nthreads = comm->nthreads;
  maxpartner = 0;

#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {

#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif

    int i,j,ii,jj,m,n,inum,jnum;
    int *ilist,*jlist,*numneigh,**firstneigh;
    int *allflags,**firstflag;
    double *allvalues,*onevalues,*jvalues;

    MyPage <tagint> &ipg = ipage_atom[tid];
    MyPage <double> &dpg = dpage_atom[tid];
    ipg.reset();
    dpg.reset();

    // 1st loop over neighbor list
    // calculate npartner for each owned atom

    tagint *tag = atom->tag;

    NeighList *list = pair->list;
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    int nlocal_neigh = 0;
    if (inum) nlocal_neigh = ilist[inum-1] + 1;

    // each thread works on a fixed chunk of local and ghost atoms.
    const int ldelta = 1 + nlocal_neigh/nthreads;
    const int lfrom = tid*ldelta;
    const int lmax = lfrom +ldelta;
    const int lto = (lmax > nlocal_neigh) ? nlocal_neigh : lmax;

    // zero npartners for all current atoms and
    // clear page data structures for this thread

    for (i = lfrom; i < lto; i++) npartner[i] = 0;


    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      allflags = firstflag[i];

      for (jj = 0; jj < jnum; jj++) {
        if (allflags[jj]) {
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

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if ((i >= lfrom) && (i < lto)) {
        n = npartner[i];
        partner[i] = ipg.get(n);
        valuepartner[i] = dpg.get(dnum*n);
        if (partner[i] == NULL || valuepartner[i] == NULL)
          error->one(FLERR,"Neighbor history overflow, boost neigh_modify one");
      }
    }

    // 2nd loop over neighbor list
    // store partner IDs and values for owned+ghost atoms
    // re-zero npartner to use as counter

    for (i = lfrom; i < lto; i++) npartner[i] = 0;

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      allflags = firstflag[i];
      allvalues = firstvalue[i];

      for (jj = 0; jj < jnum; jj++) {
        if (allflags[jj]) {
          onevalues = &allvalues[3*jj];
          j = jlist[jj];
          j &= NEIGHMASK;

          if ((i >= lfrom) && (i < lto)) {
            m = npartner[i]++;
            partner[i][m] = tag[j];
            memcpy(&valuepartner[i][dnum*m],onevalues,dnumbytes);
          }

          if ((j >= lfrom) && (j < lto)) {
            m = npartner[j]++;
            partner[j][m] = tag[i];
            jvalues = &valuepartner[j][dnum*m];
            for (n = 0; n < dnum; n++) jvalues[n] = -onevalues[n];
          }
        }
      }
    }

    // set maxpartner = max # of partners of any owned atom
    maxpartner = m = 0;
    for (i = lfrom; i < lto; i++)
      m = MAX(m,npartner[i]);

#if defined(_OPENMP)
#pragma omp critical
#endif
    maxpartner = MAX(m,maxpartner);
  }
}
