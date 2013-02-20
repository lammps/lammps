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
#include "pair_coul_long_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"

#include "suffix.h"
using namespace LAMMPS_NS;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairCoulLongOMP::PairCoulLongOMP(LAMMPS *lmp) :
  PairCoulLong(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  cut_respa = NULL;
}

/* ---------------------------------------------------------------------- */

void PairCoulLongOMP::compute(int eflag, int vflag)
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

    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairCoulLongOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itable,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double fraction,table;
  double r,r2inv,rsq,forcecoul,factor_coul;
  double grij,expm2,prefactor,t,erfc;
  int *ilist,*jlist,*numneigh,**firstneigh;

  ecoul = 0.0;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const double * _noalias const q = atom->q;
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_coul = force->special_coul;
  const double qqrd2e = force->qqrd2e;
  double fxtmp,fytmp,fztmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cut_coulsq) {
        r2inv = 1.0/rsq;
        if (!ncoultablebits || rsq <= tabinnersq) {
          r = sqrt(rsq);
          grij = g_ewald * r;
          expm2 = exp(-grij*grij);
          t = 1.0 / (1.0 + EWALD_P*grij);
          erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
          prefactor = qqrd2e * scale[itype][jtype] * qtmp*q[j]/r;
          forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
          if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
        } else {
          union_int_float_t rsq_lookup;
          rsq_lookup.f = rsq;
          itable = rsq_lookup.i & ncoulmask;
          itable >>= ncoulshiftbits;
          fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
          table = ftable[itable] + fraction*dftable[itable];
          forcecoul = scale[itype][jtype] * qtmp*q[j] * table;
          if (factor_coul < 1.0) {
            table = ctable[itable] + fraction*dctable[itable];
            prefactor = scale[itype][jtype] * qtmp*q[j] * table;
            forcecoul -= (1.0-factor_coul)*prefactor;
          }
        }

        fpair = forcecoul * r2inv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if (EFLAG) {
          if (!ncoultablebits || rsq <= tabinnersq)
            ecoul = prefactor*erfc;
          else {
            table = etable[itable] + fraction*detable[itable];
            ecoul = scale[itype][jtype] * qtmp*q[j] * table;
          }
          if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
        }

        if (EVFLAG) ev_tally_thr(this, i,j,nlocal,NEWTON_PAIR,
                                 0.0,ecoul,fpair,delx,dely,delz,thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairCoulLongOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairCoulLong::memory_usage();

  return bytes;
}
