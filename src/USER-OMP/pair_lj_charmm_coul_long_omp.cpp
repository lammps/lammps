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

#include "omp_compat.h"
#include <cmath>
#include "pair_lj_charmm_coul_long_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"

#include "suffix.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCharmmCoulLongOMP::PairLJCharmmCoulLongOMP(LAMMPS *lmp) :
  PairLJCharmmCoulLong(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  cut_respa = NULL;
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulLongOMP::compute(int eflag, int vflag)
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
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, NULL, thr);

    if (evflag) {
      if (eflag) {
        if (force->newton_pair) eval<1,1,1>(ifrom, ito, thr);
        else eval<1,1,0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair) eval<1,0,1>(ifrom, ito, thr);
        else eval<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(ifrom, ito, thr);
      else eval<0,0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLJCharmmCoulLongOMP::eval(int iifrom, int iito, ThrData * const thr)
{

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const double * _noalias const q = atom->q;
  const int * _noalias const type = atom->type;
  const double * _noalias const special_coul = force->special_coul;
  const double * _noalias const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;
  const double inv_denom_lj = 1.0/denom_lj;

  const int * const ilist = list->ilist;
  const int * const numneigh = list->numneigh;
  const int * const * const firstneigh = list->firstneigh;
  const int nlocal = atom->nlocal;

  // loop over neighbors of my atoms

  for (int ii = iifrom; ii < iito; ++ii) {

    const int i = ilist[ii];
    const int itype = type[i];
    const double qtmp = q[i];
    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    double fxtmp,fytmp,fztmp;
    fxtmp=fytmp=fztmp=0.0;

    const int * const jlist = firstneigh[i];
    const int jnum = numneigh[i];
    const double * _noalias const lj1i = lj1[itype];
    const double * _noalias const lj2i = lj2[itype];
    const double * _noalias const lj3i = lj3[itype];
    const double * _noalias const lj4i = lj4[itype];

    for (int jj = 0; jj < jnum; jj++) {
      double forcecoul, forcelj, evdwl, ecoul;
      forcecoul = forcelj = evdwl = ecoul = 0.0;

      const int sbindex = sbmask(jlist[jj]);
      const int j = jlist[jj] & NEIGHMASK;

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq = delx*delx + dely*dely + delz*delz;
      const int jtype = type[j];

      if (rsq < cut_bothsq) {
        const double r2inv = 1.0/rsq;

        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq) {
            const double A1 =  0.254829592;
            const double A2 = -0.284496736;
            const double A3 =  1.421413741;
            const double A4 = -1.453152027;
            const double A5 =  1.061405429;
            const double EWALD_F = 1.12837917;
            const double INV_EWALD_P = 1.0/0.3275911;

            const double r = sqrt(rsq);
            const double grij = g_ewald * r;
            const double expm2 = exp(-grij*grij);
            const double t = INV_EWALD_P / (INV_EWALD_P + grij);
            const double erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
            const double prefactor = qqrd2e * qtmp*q[j]/r;
            forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
            if (EFLAG) ecoul = prefactor*erfc;
            if (sbindex) {
              const double adjust = (1.0-special_coul[sbindex])*prefactor;
              forcecoul -= adjust;
              if (EFLAG) ecoul -= adjust;
            }
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            const int itable = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
            const double fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
            const double table = ftable[itable] + fraction*dftable[itable];
            forcecoul = qtmp*q[j] * table;
            if (EFLAG) ecoul = qtmp*q[j] * (etable[itable] + fraction*detable[itable]);
            if (sbindex) {
              const double table2 = ctable[itable] + fraction*dctable[itable];
              const double prefactor = qtmp*q[j] * table2;
              const double adjust = (1.0-special_coul[sbindex])*prefactor;
              forcecoul -= adjust;
              if (EFLAG) ecoul -= adjust;
            }
          }
        }

        if (rsq < cut_ljsq) {
          const double r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (lj1i[jtype]*r6inv - lj2i[jtype]);
          const double philj = r6inv*(lj3i[jtype]*r6inv-lj4i[jtype]);
          if (EFLAG) evdwl = philj;

          if (rsq > cut_lj_innersq) {
            const double drsq = cut_ljsq - rsq;
            const double cut2 = (rsq - cut_lj_innersq) * drsq;
            const double switch1 = drsq * (drsq*drsq + 3.0*cut2) * inv_denom_lj;
            const double switch2 = 12.0*rsq * cut2 * inv_denom_lj;
            forcelj = forcelj*switch1 + philj*switch2;
            if (EFLAG) evdwl *= switch1;
          }

          if (sbindex) {
            const double factor_lj = special_lj[sbindex];
            forcelj *= factor_lj;
            if (EFLAG) evdwl *= factor_lj;
          }
        }
        const double fpair = (forcecoul + forcelj) * r2inv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if (EVFLAG) ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,
                                 evdwl,ecoul,fpair,delx,dely,delz,thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairLJCharmmCoulLongOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJCharmmCoulLong::memory_usage();

  return bytes;
}
