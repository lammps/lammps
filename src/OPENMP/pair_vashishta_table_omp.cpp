// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_vashishta_table_omp.h"

#include "atom.h"
#include "comm.h"
#include "memory.h"
#include "neigh_list.h"
#include "suffix.h"

#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairVashishtaTableOMP::PairVashishtaTableOMP(LAMMPS *lmp) :
  PairVashishtaTable(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairVashishtaTableOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    if (evflag) {
      if (eflag) {
        eval<1,1>(ifrom, ito, thr);
      } else {
        eval<1,0>(ifrom, ito, thr);
      }
    } else eval<0,0>(ifrom, ito, thr);

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG>
void PairVashishtaTableOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,k,ii,jj,kk,jnum,jnumm1,maxshort_thr;
  tagint itag,jtag;
  int itype,jtype,ktype,ijparam,ikparam,ijkparam;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fj[3],fk[3];
  int *ilist,*jlist,*numneigh,**firstneigh,*neighshort_thr;

  evdwl = 0.0;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const tagint * _noalias const tag = atom->tag;
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double cutshortsq = r0max*r0max;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  maxshort_thr = maxshort;
  memory->create(neighshort_thr,maxshort_thr,"pair_thr:neighshort_thr");

  double fxtmp,fytmp,fztmp;

  // loop over full neighbor list of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    fxtmp = fytmp = fztmp = 0.0;

    // two-body interactions, skip half of them

    jlist = firstneigh[i];
    jnum = numneigh[i];
    int numshort = 0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutshortsq) {
        neighshort_thr[numshort++] = j;
        if (numshort >= maxshort_thr) {
          maxshort_thr += maxshort_thr/2;
          memory->grow(neighshort_thr,maxshort_thr,"pair_thr:neighshort_thr");
        }
      }

      jtag = tag[j];
      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j].z < ztmp) continue;
        if (x[j].z == ztmp && x[j].y < ytmp) continue;
        if (x[j].z == ztmp && x[j].y == ytmp && x[j].x < xtmp) continue;
      }

      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];
      if (rsq >= params[ijparam].cutsq) continue;

      twobody_table(params[ijparam],rsq,fpair,EFLAG,evdwl);

      fxtmp += delx*fpair;
      fytmp += dely*fpair;
      fztmp += delz*fpair;
      f[j].x -= delx*fpair;
      f[j].y -= dely*fpair;
      f[j].z -= delz*fpair;

      if (EVFLAG) ev_tally_thr(this,i,j,nlocal,/* newton_pair */ 1,
                               evdwl,0.0,fpair,delx,dely,delz,thr);
    }

    jnumm1 = numshort - 1;

    for (jj = 0; jj < jnumm1; jj++) {
      j = neighshort_thr[jj];
      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];
      delr1[0] = x[j].x - xtmp;
      delr1[1] = x[j].y - ytmp;
      delr1[2] = x[j].z - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 >= params[ijparam].cutsq2) continue;

      double fjxtmp,fjytmp,fjztmp;
      fjxtmp = fjytmp = fjztmp = 0.0;

      for (kk = jj+1; kk < numshort; kk++) {
        k = neighshort_thr[kk];
        ktype = map[type[k]];
        ikparam = elem3param[itype][ktype][ktype];
        ijkparam = elem3param[itype][jtype][ktype];

        delr2[0] = x[k].x - xtmp;
        delr2[1] = x[k].y - ytmp;
        delr2[2] = x[k].z - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 >= params[ikparam].cutsq2) continue;

        threebody(&params[ijparam],&params[ikparam],&params[ijkparam],
                  rsq1,rsq2,delr1,delr2,fj,fk,EFLAG,evdwl);

        fxtmp -= fj[0] + fk[0];
        fytmp -= fj[1] + fk[1];
        fztmp -= fj[2] + fk[2];
        fjxtmp += fj[0];
        fjytmp += fj[1];
        fjztmp += fj[2];
        f[k].x += fk[0];
        f[k].y += fk[1];
        f[k].z += fk[2];

        if (EVFLAG) ev_tally3_thr(this,i,j,k,evdwl,0.0,fj,fk,delr1,delr2,thr);
      }
      f[j].x += fjxtmp;
      f[j].y += fjytmp;
      f[j].z += fjztmp;
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
  memory->destroy(neighshort_thr);
}

/* ---------------------------------------------------------------------- */

double PairVashishtaTableOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairVashishtaTable::memory_usage();

  return bytes;
}
