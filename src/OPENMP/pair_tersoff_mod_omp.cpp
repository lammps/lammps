// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_tersoff_mod_omp.h"

#include "atom.h"
#include "comm.h"
#include "math_extra.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace MathExtra;

/* ---------------------------------------------------------------------- */

PairTersoffMODOMP::PairTersoffMODOMP(LAMMPS *lmp) :
  PairTersoffMOD(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairTersoffMODOMP::compute(int eflag, int vflag)
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

    if (shift_flag) {
      if (evflag) {
        if (eflag) {
          if (vflag_either) eval<1,1,1,1>(ifrom, ito, thr);
          else eval<1,1,1,0>(ifrom, ito, thr);
        } else {
          if (vflag_either) eval<1,1,0,1>(ifrom, ito, thr);
          else eval<1,1,0,0>(ifrom, ito, thr);
        }
      } else eval<1,0,0,0>(ifrom, ito, thr);
    } else {
      if (evflag) {
        if (eflag) {
          if (vflag_either) eval<0,1,1,1>(ifrom, ito, thr);
          else eval<0,1,1,0>(ifrom, ito, thr);
        } else {
          if (vflag_either) eval<0,1,0,1>(ifrom, ito, thr);
          else eval<0,1,0,0>(ifrom, ito, thr);
        }
      } else eval<0,0,0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int SHIFT_FLAG, int EVFLAG, int EFLAG, int VFLAG_EITHER>
void PairTersoffMODOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,k,ii,jj,kk,jnum;
  tagint itag,jtag;
  int itype,jtype,ktype,iparam_ij,iparam_ijk;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double fforce;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fi[3],fj[3],fk[3];
  double r1_hat[3],r2_hat[3];
  double zeta_ij,prefactor;
  double forceshiftfac;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const tagint * _noalias const tag = atom->tag;
  const int * _noalias const type = atom->type;
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
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
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
        if (x[j].z < ztmp) continue;
        if (x[j].z == ztmp && x[j].y < ytmp) continue;
        if (x[j].z == ztmp && x[j].y == ytmp && x[j].x < xtmp) continue;
      }

      jtype = map[type[j]];

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;

      // shift rsq and store correction for force

      if (SHIFT_FLAG) {
        double rsqtmp = rsq + shift*shift + 2*sqrt(rsq)*shift;
        forceshiftfac = sqrt(rsqtmp/rsq);
        rsq = rsqtmp;
      }

      iparam_ij = elem3param[itype][jtype][jtype];
      if (rsq > params[iparam_ij].cutsq) continue;

      repulsive(&params[iparam_ij],rsq,fpair,EFLAG,evdwl);

      // correct force for shift in rsq

      if (SHIFT_FLAG) fpair *= forceshiftfac;

      fxtmp += delx*fpair;
      fytmp += dely*fpair;
      fztmp += delz*fpair;
      f[j].x -= delx*fpair;
      f[j].y -= dely*fpair;
      f[j].z -= delz*fpair;

      if (EVFLAG) ev_tally_thr(this,i,j,nlocal,/* newton_pair */ 1,
                               evdwl,0.0,fpair,delx,dely,delz,thr);
    }

    // three-body interactions
    // skip immediately if I-J is not within cutoff
    double fjxtmp,fjytmp,fjztmp;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];
      iparam_ij = elem3param[itype][jtype][jtype];

      delr1[0] = x[j].x - xtmp;
      delr1[1] = x[j].y - ytmp;
      delr1[2] = x[j].z - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

      if (SHIFT_FLAG)
        rsq1 += shift*shift + 2*sqrt(rsq1)*shift;

      if (rsq1 > params[iparam_ij].cutsq) continue;

      const double r1inv = 1.0/sqrt(dot3(delr1, delr1));
      scale3(r1inv, delr1, r1_hat);

      // accumulate bondorder zeta for each i-j interaction via loop over k

      fjxtmp = fjytmp = fjztmp = 0.0;
      zeta_ij = 0.0;

      for (kk = 0; kk < jnum; kk++) {
        if (jj == kk) continue;
        k = jlist[kk];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        iparam_ijk = elem3param[itype][jtype][ktype];

        delr2[0] = x[k].x - xtmp;
        delr2[1] = x[k].y - ytmp;
        delr2[2] = x[k].z - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

        if (SHIFT_FLAG)
          rsq2 += shift*shift + 2*sqrt(rsq2)*shift;

        if (rsq2 > params[iparam_ijk].cutsq) continue;

        const double r2inv = 1.0/sqrt(dot3(delr2, delr2));
        scale3(r2inv, delr2, r2_hat);

        zeta_ij += zeta(&params[iparam_ijk],rsq1,rsq2,r1_hat,r2_hat);
      }

      // pairwise force due to zeta

      force_zeta(&params[iparam_ij],rsq1,zeta_ij,fforce,prefactor,EFLAG,evdwl);

      fpair = fforce*r1inv;

      fxtmp += delr1[0]*fpair;
      fytmp += delr1[1]*fpair;
      fztmp += delr1[2]*fpair;
      fjxtmp -= delr1[0]*fpair;
      fjytmp -= delr1[1]*fpair;
      fjztmp -= delr1[2]*fpair;

      if (EVFLAG) ev_tally_thr(this,i,j,nlocal,/* newton_pair */ 1,evdwl,0.0,
                               -fpair,-delr1[0],-delr1[1],-delr1[2],thr);

      // attractive term via loop over k

      for (kk = 0; kk < jnum; kk++) {
        if (jj == kk) continue;
        k = jlist[kk];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        iparam_ijk = elem3param[itype][jtype][ktype];

        delr2[0] = x[k].x - xtmp;
        delr2[1] = x[k].y - ytmp;
        delr2[2] = x[k].z - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

        if (SHIFT_FLAG)
          rsq2 += shift*shift + 2*sqrt(rsq2)*shift;

        if (rsq2 > params[iparam_ijk].cutsq) continue;

        const double r2inv = 1.0/sqrt(dot3(delr2, delr2));
        scale3(r2inv, delr2, r2_hat);

        attractive(&params[iparam_ijk],prefactor,
                   rsq1,rsq2,r1_hat,r2_hat,fi,fj,fk);

        fxtmp += fi[0];
        fytmp += fi[1];
        fztmp += fi[2];
        fjxtmp += fj[0];
        fjytmp += fj[1];
        fjztmp += fj[2];
        f[k].x += fk[0];
        f[k].y += fk[1];
        f[k].z += fk[2];

        if (VFLAG_EITHER) v_tally3_thr(this,i,j,k,fj,fk,delr1,delr2,thr);
      }
      f[j].x += fjxtmp;
      f[j].y += fjytmp;
      f[j].z += fjztmp;
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairTersoffMODOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairTersoffMOD::memory_usage();

  return bytes;
}
