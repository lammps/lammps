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

#define MAXTOUCH 15

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

  int flag = 0;
#if defined(_OPENMP)
#pragma omp parallel default(none) shared(flag)
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

    const int gdelta = 1 + nghost/nthreads;
    const int gfrom = nlocal + tid*gdelta;
    const int gmax = gfrom + gdelta;
    const int gto = (gmax > nall) ? nall : gmax;


    int i,j,ii,jj,m,inum,jnum;
    int *ilist,*jlist,*numneigh,**firstneigh;
    int *touch,**firsttouch;
    double *shear,*allshear,**firstshear;

    // zero npartners for all current atoms

    for (i = lfrom; i < lto; i++) npartner[i] = 0;

    // copy shear info from neighbor list atoms to atom arrays

    int *tag = atom->tag;

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
      allshear = firstshear[i];
      jnum = numneigh[i];
      touch = firsttouch[i];

      for (jj = 0; jj < jnum; jj++) {
        if (touch[jj]) {
          j = jlist[jj];
          j &= NEIGHMASK;
          shear = &allshear[3*jj];

          if ((i >= lfrom) && (i < lto)) {
            if (npartner[i] < MAXTOUCH) {
              m = npartner[i];
              partner[i][m] = tag[j];
              shearpartner[i][m][0] = shear[0];
              shearpartner[i][m][1] = shear[1];
              shearpartner[i][m][2] = shear[2];
            }
            npartner[i]++;
          }

          if ((j >= lfrom) && (j < lto)) {
            if (npartner[j] < MAXTOUCH) {
              m = npartner[j];
              partner[j][m] = tag[i];
              shearpartner[j][m][0] = -shear[0];
              shearpartner[j][m][1] = -shear[1];
              shearpartner[j][m][2] = -shear[2];
            }
            npartner[j]++;
          }

          if ((j >= gfrom) && (j < gto)) {
            npartner[j]++;
          }
        }
      }
    }

    // test for too many touching neighbors
    int myflag = 0;
    for (i = lfrom; i < lto; i++)
      if (npartner[i] >= MAXTOUCH) myflag = 1;

    if (myflag)
#if defined(_OPENMP)
#pragma omp atomic
#endif
      ++flag;
  }

  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all)
    error->all(FLERR,"Too many touching neighbors - boost MAXTOUCH");
}
