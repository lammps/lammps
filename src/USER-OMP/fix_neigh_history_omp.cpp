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

#include "fix_neigh_history_omp.h"
#include <cstring>
#include "my_page.h"
#include "atom.h"
#include "comm.h"
#include "neigh_list.h"
#include "pair.h"
#include "memory.h"
#include "error.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;
using namespace FixConst;

enum{DEFAULT,NPARTNER,PERPARTNER}; // also set in fix neigh/history


FixNeighHistoryOMP::FixNeighHistoryOMP(class LAMMPS *lmp,int narg,char **argv)
  : FixNeighHistory(lmp,narg,argv) {

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
   onesided and newton on and newton off versions
------------------------------------------------------------------------- */
// below is the pre_exchange() function from the parent class
// void FixNeighHistory::pre_exchange()
// {
//  if (onesided) pre_exchange_onesided();
//  else if (newton_pair) pre_exchange_newton();
//  else pre_exchange_no_newton();
//}

/* ----------------------------------------------------------------------
   onesided version for sphere contact with line/tri particles
   neighbor list has I = sphere, J = line/tri
   only store history info with spheres
------------------------------------------------------------------------- */

void FixNeighHistoryOMP::pre_exchange_onesided()
{
  const int nthreads = comm->nthreads;
  const int nlocal = atom->nlocal;
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
    int *allflags;
    double *allvalues,*onevalues;

    // NOTE: all operations until very end are with:
    //   nlocal_neigh <= current nlocal
    // b/c previous neigh list was built with nlocal_neigh
    // nlocal can be larger if other fixes added atoms at this pre_exchange()

    // clear per-thread paged data structures

    MyPage <tagint> &ipg = ipage_atom[tid];
    MyPage <double> &dpg = dpage_atom[tid];
    ipg.reset();
    dpg.reset();

    // each thread works on a fixed chunk of local and ghost atoms.
    const int ldelta = 1 + nlocal_neigh/nthreads;
    const int lfrom = tid*ldelta;
    const int lmax = lfrom +ldelta;
    const int lto = (lmax > nlocal_neigh) ? nlocal_neigh : lmax;

    // 1st loop over neighbor list, I = sphere, J = tri
    // only calculate npartner for each owned spheres

    for (i = lfrom; i < lto; i++) npartner[i] = 0;

    tagint *tag = atom->tag;
    NeighList *list = pair->list;
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      allflags = firstflag[i];

      for (jj = 0; jj < jnum; jj++)
        if (allflags[jj])
          if ((i >= lfrom) && (i < lto)) npartner[i]++;
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
          onevalues = &allvalues[dnum*jj];
          j = jlist[jj];
          j &= NEIGHMASK;

          if ((i >= lfrom) && (i < lto)) {
            m = npartner[i]++;
            partner[i][m] = tag[j];
            memcpy(&valuepartner[i][dnum*m],onevalues,dnumbytes);
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
    {
      maxpartner = MAX(m,maxpartner);
      comm->maxexchange_fix =MAX(comm->maxexchange_fix,(dnum+1)*maxpartner+1);
    }
  }

  // zero npartner values from previous nlocal_neigh to current nlocal
  for (int i = nlocal_neigh; i < nlocal; ++i) npartner[i] = 0;
}

/* -------------------------------------------------------------------- */

void FixNeighHistoryOMP::pre_exchange_newton()
{
  const int nthreads = comm->nthreads;
  maxpartner = 0;
  for (int i = 0; i < nall_neigh; i++) npartner[i] = 0;

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
    int *allflags;
    double *allvalues,*onevalues,*jvalues;

    MyPage <tagint> &ipg = ipage_atom[tid];
    MyPage <double> &dpg = dpage_atom[tid];
    ipg.reset();
    dpg.reset();

    // 1st loop over neighbor list
    // calculate npartner for each owned+ghost atom

    tagint *tag = atom->tag;

    NeighList *list = pair->list;
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    // each thread works on a fixed chunk of local and ghost atoms.
    const int ldelta = 1 + nlocal_neigh/nthreads;
    const int lfrom = tid*ldelta;
    const int lmax = lfrom +ldelta;
    const int lto = (lmax > nlocal_neigh) ? nlocal_neigh : lmax;

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
#if defined(_OPENMP)
#pragma omp barrier
    {;}

    // perform reverse comm to augment owned npartner counts with ghost counts

#pragma omp master
#endif
    {
      commflag = NPARTNER;
      comm->reverse_comm_fix(this,0);
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

#if defined(_OPENMP)
#pragma omp master
#endif
    {
      for (i = nlocal_neigh; i < nall_neigh; i++) {
        n = npartner[i];
        partner[i] = ipg.get(n);
        valuepartner[i] = dpg.get(dnum*n);
        if (partner[i] == NULL || valuepartner[i] == NULL) {
          error->one(FLERR,"Neighbor history overflow, boost neigh_modify one");
        }
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
          onevalues = &allvalues[dnum*jj];
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
#if defined(_OPENMP)
#pragma omp barrier
    {;}

#pragma omp master
#endif
    {
      // perform reverse comm to augment
      // owned atom partner/valuepartner with ghost info
      // use variable variant b/c size of packed data can be arbitrarily large
      //   if many touching neighbors for large particle

      commflag = PERPARTNER;
      comm->reverse_comm_fix_variable(this);
    }

    // set maxpartner = max # of partners of any owned atom
    m = 0;
    for (i = lfrom; i < lto; i++)
      m = MAX(m,npartner[i]);

#if defined(_OPENMP)
#pragma omp critical
#endif
    {
      maxpartner = MAX(m,maxpartner);
      comm->maxexchange_fix = MAX(comm->maxexchange_fix,(dnum+1)*maxpartner+1);
    }
  }

  // zero npartner values from previous nlocal_neigh to current nlocal

  int nlocal = atom->nlocal;
  for (int i = nlocal_neigh; i < nlocal; i++) npartner[i] = 0;
}

/* -------------------------------------------------------------------- */

void FixNeighHistoryOMP::pre_exchange_no_newton()
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
    int *allflags;
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
          onevalues = &allvalues[dnum*jj];
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
    m = 0;
    for (i = lfrom; i < lto; i++)
      m = MAX(m,npartner[i]);

#if defined(_OPENMP)
#pragma omp critical
#endif
    {
      maxpartner = MAX(m,maxpartner);
      comm->maxexchange_fix = MAX(comm->maxexchange_fix,(dnum+1)*maxpartner+1);
    }
  }
}

/* -------------------------------------------------------------------- */

void FixNeighHistoryOMP::post_neighbor()
{
  const int nthreads = comm->nthreads;
  maxpartner = 0;
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  nlocal_neigh = nlocal;
  nall_neigh = nall;

  // realloc firstflag and firstvalue if needed

  if (maxatom < nlocal) {
    memory->sfree(firstflag);
    memory->sfree(firstvalue);
    maxatom = nall;
    firstflag = (int **)
      memory->smalloc(maxatom*sizeof(int *),"neighbor_history:firstflag");
    firstvalue = (double **)
      memory->smalloc(maxatom*sizeof(double *),"neighbor_history:firstvalue");
  }


#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {

#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif

    int i,j,ii,jj,m,nn,np,inum,jnum,rflag;
    tagint jtag;
    int *ilist,*jlist,*numneigh,**firstneigh;
    int *allflags;
    double *allvalues;

    MyPage <int> &ipg = ipage_neigh[tid];
    MyPage <double> &dpg = dpage_neigh[tid];
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

    // each thread works on a fixed chunk of local and ghost atoms.
    const int ldelta = 1 + inum/nthreads;
    const int lfrom = tid*ldelta;
    const int lmax = lfrom +ldelta;
    const int lto = (lmax > inum) ? inum : lmax;

    for (ii = lfrom; ii < lto; ii++) {
      i = ilist[ii];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      firstflag[i] = allflags = ipg.get(jnum);
      firstvalue[i] = allvalues = dpg.get(jnum*dnum);
      np = npartner[i];
      nn = 0;

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        rflag = sbmask(j);
        j &= NEIGHMASK;
        jlist[jj] = j;

        // rflag = 1 if r < radsum in npair_size() method
        // preserve neigh history info if tag[j] is in old-neigh partner list
        // this test could be more geometrically precise for two sphere/line/tri

        if (rflag) {
          jtag = tag[j];
          for (m = 0; m < np; m++)
            if (partner[i][m] == jtag) break;
          if (m < np) {
            allflags[jj] = 1;
            memcpy(&allvalues[nn],&valuepartner[i][dnum*m],dnumbytes);
          } else {
            allflags[jj] = 0;
            memcpy(&allvalues[nn],zeroes,dnumbytes);
          }
        } else {
          allflags[jj] = 0;
          memcpy(&allvalues[nn],zeroes,dnumbytes);
        }
        nn += dnum;
      }
    }
  }
}
