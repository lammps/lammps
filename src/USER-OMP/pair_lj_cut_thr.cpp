/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "math.h"
#include "pair_lj_cut_thr.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"

using namespace LAMMPS_NS;

#if defined(_OPENMP)
#include <omp.h>
#endif

typedef struct { double x,y,z; } dbl3_t;
#if defined(__GNUC__)
#define _noalias __restrict
#else
#define _noalias
#endif

// set loop range thread id, and force array offset for threaded runs.
static inline void loop_setup_thr(int &ifrom, int &ito, int &tid,
                                  int inum, int nthreads)
{
#if defined(_OPENMP)
  tid = omp_get_thread_num();

  // each thread works on a fixed chunk of atoms.
  const int idelta = 1 + inum/nthreads;
  ifrom = tid*idelta;
  ito   = ((ifrom + idelta) > inum) ? inum : ifrom + idelta;
#else
  tid = 0;
  ifrom = 0;
  ito = inum;
#endif
}

/* ---------------------------------------------------------------------- */

PairLJCutThr::PairLJCutThr(LAMMPS *lmp) : PairLJCut(lmp)
{
  respa_enable = 0;
  cut_respa = NULL;
}

/* ---------------------------------------------------------------------- */

void PairLJCutThr::init_style()
{
  if (force->newton_pair)
    error->all(FLERR,"Cannot use newton pair with lj/cut/thr pair style");

  int irequest = neighbor->request(this);

  // request a full neighbor list

  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairLJCutThr::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  if (eflag_either)
    if (vflag_either) eval<1,1>();
    else eval<1,0>();
  else
    if (vflag_either) eval<0,1>();
    else eval<0,0>();

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

// workload distribution summary.
// we have a full neighbor list
// each thread processes a range of the local atoms
// for each atom from the neighbor list we have 2 cases
// 1) the neighbor is in the same chunk of assigned local atoms indices
//    => only process it if i < j 
// 2) the neighbor is outside of the chunk of assigned local atom indices
//    => only compute force on i

template <int EFLAG, int VFLAG>
void PairLJCutThr::eval()
{
  double eg,v0,v1,v2,v3,v4,v5;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  eg=v0=v1=v2=v3=v4=v5=0.0;;

#if defined(_OPENMP)
#pragma omp parallel default(none) reduction(+:eg,v0,v1,v2,v3,v4,v5)
#endif
  {
    int ifrom, ito, tid;
    loop_setup_thr(ifrom, ito, tid, inum, nthreads);

    int i,j,ii,jj,jnum,itype,jtype;
    double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
    double rsq,r2inv,r6inv,forcelj,factor_lj;
    int *ilist,*jlist,*numneigh,**firstneigh;

    evdwl = 0.0;

    const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
    dbl3_t * _noalias const f = (dbl3_t *) atom->f[0];
    const int * _noalias const type = atom->type;
    const double * _noalias const special_lj = force->special_lj;
    double fxtmp,fytmp,fztmp,v[6];

    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    // loop over neighbors of my atoms

    for (ii = ifrom; ii < ito; ++ii) {
      i = ilist[ii];
      xtmp = x[i].x;
      ytmp = x[i].y;
      ztmp = x[i].z;
      itype = type[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      fxtmp=fytmp=fztmp=0.0;

      for (jj = 0; jj < jnum; jj++) {
	j = jlist[jj];
	if ((j < ito) && (i > j)) continue;
	factor_lj = special_lj[sbmask(j)];
	j &= NEIGHMASK;

	delx = xtmp - x[j].x;
	dely = ytmp - x[j].y;
	delz = ztmp - x[j].z;
	rsq = delx*delx + dely*dely + delz*delz;
	jtype = type[j];

	if (rsq < cutsq[itype][jtype]) {
	  r2inv = 1.0/rsq;
	  r6inv = r2inv*r2inv*r2inv;
	  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	  fpair = factor_lj*forcelj*r2inv;

	  fxtmp += delx*fpair;
	  fytmp += dely*fpair;
	  fztmp += delz*fpair;
	  if ((j < ito) && (i < j)) {
	    f[j].x -= delx*fpair;
	    f[j].y -= dely*fpair;
	    f[j].z -= delz*fpair;
	  }

	  if (VFLAG) {
	    fpair *= 0.5;
	    v[0] = delx*delx*fpair;
	    v[1] = dely*dely*fpair;
	    v[2] = delz*delz*fpair;
	    v[3] = delx*dely*fpair;
	    v[4] = delx*delz*fpair;
	    v[5] = dely*delz*fpair;
	  }

	  if (EFLAG) {
	    evdwl = r6inv*(lj3[itype][jtype]*r6inv
			   - lj4[itype][jtype]) - offset[itype][jtype];
	    evdwl *= 0.5*factor_lj;
	  }
	  
	  if (EFLAG || VFLAG) {
	    if (EFLAG) {
	      if (eflag_global) eg += evdwl;
	      if (eflag_atom) eatom[i] += evdwl;
	      
	      if ((j < ito) && (i < j)) {
		if (eflag_global) eg += evdwl;
		if (eflag_atom) eatom[j] += evdwl;
	      }
	    }
	    if (VFLAG) {
	      if (vflag_global) {
		v0 += v[0];
		v1 += v[1];
		v2 += v[2];
		v3 += v[3];
		v4 += v[4];
		v5 += v[5];
	      }
	      
	      if (vflag_atom) {
		vatom[i][0] += v[0];
		vatom[i][1] += v[1];
		vatom[i][2] += v[2];
		vatom[i][3] += v[3];
		vatom[i][4] += v[4];
		vatom[i][5] += v[5];
	      }
	      if ((j < ito) && (i < j)) {
		if (vflag_global) {
		  v0 += v[0];
		  v1 += v[1];
		  v2 += v[2];
		  v3 += v[3];
		  v4 += v[4];
		  v5 += v[5];
		}
		if (vflag_atom) {
		  vatom[j][0] += v[0];
		  vatom[j][1] += v[1];
		  vatom[j][2] += v[2];
		  vatom[j][3] += v[3];
		  vatom[j][4] += v[4];
		  vatom[j][5] += v[5];
		}
	      }
	    }
	  }
	}
      }
      f[i].x += fxtmp;
      f[i].y += fytmp;
      f[i].z += fztmp;
    }
  } // end of omp parallel region
  eng_vdwl += eg;
  virial[0] += v0;
  virial[1] += v1;
  virial[2] += v2;
  virial[3] += v3;
  virial[4] += v4;
  virial[5] += v5;
}

