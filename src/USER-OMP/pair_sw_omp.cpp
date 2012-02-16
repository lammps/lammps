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
#include "pair_sw_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"

#include "suffix.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSWOMP::PairSWOMP(LAMMPS *lmp) :
  PairSW(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairSWOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

    if (evflag) {
      if (eflag) {
	eval<1,1>(ifrom, ito, thr);
      } else {
	eval<1,0>(ifrom, ito, thr);
      }
    } else eval<0,0>(ifrom, ito, thr);

    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG>
void PairSWOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,k,ii,jj,kk,jnum,jnumm1,itag,jtag;
  int itype,jtype,ktype,ijparam,ikparam,ijkparam;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fj[3],fk[3];
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const int * const tag = atom->tag;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double fxtmp,fytmp,fztmp;

  // loop over full neighbor list of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp = fytmp = fztmp = 0.0;

    // two-body interactions, skip half of them

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtag = tag[j];

      if (itag > jtag) {
	if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
	if ((itag+jtag) % 2 == 1) continue;
      } else {
	if (x[j][2] < ztmp) continue;
	if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
	if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      jtype = map[type[j]];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      ijparam = elem2param[itype][jtype][jtype];
      if (rsq > params[ijparam].cutsq) continue;

      twobody(&params[ijparam],rsq,fpair,EFLAG,evdwl);

      fxtmp += delx*fpair;
      fytmp += dely*fpair;
      fztmp += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (EVFLAG) ev_tally_thr(this,i,j,nlocal,/* newton_pair */ 1,
			       evdwl,0.0,fpair,delx,dely,delz,thr);
    }

    jnumm1 = jnum - 1;

    for (jj = 0; jj < jnumm1; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];
      ijparam = elem2param[itype][jtype][jtype];
      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 > params[ijparam].cutsq) continue;

      double fjxtmp,fjytmp,fjztmp;
      fjxtmp = fjytmp = fjztmp = 0.0;

      for (kk = jj+1; kk < jnum; kk++) {
	k = jlist[kk];
	k &= NEIGHMASK;
	ktype = map[type[k]];
	ikparam = elem2param[itype][ktype][ktype];
	ijkparam = elem2param[itype][jtype][ktype];

	delr2[0] = x[k][0] - xtmp;
	delr2[1] = x[k][1] - ytmp;
	delr2[2] = x[k][2] - ztmp;
	rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
	if (rsq2 > params[ikparam].cutsq) continue;

	threebody(&params[ijparam],&params[ikparam],&params[ijkparam],
		  rsq1,rsq2,delr1,delr2,fj,fk,EFLAG,evdwl);

	fxtmp -= fj[0] + fk[0];
	fytmp -= fj[1] + fk[1];
	fztmp -= fj[2] + fk[2];
	fjxtmp += fj[0];
	fjytmp += fj[1];
	fjztmp += fj[2];
	f[k][0] += fk[0];
	f[k][1] += fk[1];
	f[k][2] += fk[2];

	if (EVFLAG) ev_tally3_thr(this,i,j,k,evdwl,0.0,fj,fk,delr1,delr2,thr);
      }
      f[j][0] += fjxtmp;
      f[j][1] += fjytmp;
      f[j][2] += fjztmp;
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairSWOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairSW::memory_usage();

  return bytes;
}
