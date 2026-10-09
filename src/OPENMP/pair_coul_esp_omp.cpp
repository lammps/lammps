// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Chun-Yeol You (DGIST)
------------------------------------------------------------------------- */

#include "pair_coul_esp_omp.h"

#include "atom.h"
#include "comm.h"
#include "ewald_const.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"

using namespace LAMMPS_NS;
using namespace EwaldConst;

/* ---------------------------------------------------------------------- */

PairCoulEspOMP::PairCoulEspOMP(LAMMPS *lmp) :
  PairCoulEsp(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  cut_respa = nullptr;
}

/* ---------------------------------------------------------------------- */

void PairCoulEspOMP::compute(int eflag, int vflag)
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
void PairCoulEspOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itable,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double fraction,table;
  double r,r2inv,rsq,forcecoul,factor_coul,prefactor;
  int *ilist,*jlist,*numneigh,**firstneigh;

  ecoul = 0.0;

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
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

          double force_poly_appx = force_poly_coeff[0];
          double force_poly_r = 1.0;
          double r_scal = 2.0*r/cut_coul - 1.0;

          for (int index = 1; index < num_of_force_poly; index++) {
            force_poly_r *= r_scal;
            force_poly_appx += force_poly_coeff[index]*force_poly_r;
          }
          prefactor = qqrd2e * qtmp*q[j]/r;
          forcecoul = prefactor * force_poly_appx;
          if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
        } else {
          union_int_float_t rsq_lookup;
          rsq_lookup.f = rsq;
          itable = rsq_lookup.i & ncoulmask;
          itable >>= ncoulshiftbits;
          fraction = ((double) rsq_lookup.f - rtable[itable]) * drtable[itable];
          table = ftable[itable] + fraction*dftable[itable];
          forcecoul = qtmp*q[j] * table;
          if (factor_coul < 1.0) {
            table = ctable[itable] + fraction*dctable[itable];
            prefactor = qtmp*q[j] * table;
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
          if (!ncoultablebits || rsq <= tabinnersq) {
            double energy_poly_appx = energy_poly_coeff[0];
            double energy_poly_r = 1.0;
            double r_scal = 2.0*r/cut_coul - 1.0;
            for (int index = 1; index < num_of_energy_poly; index++) {
              energy_poly_r *= r_scal;
              energy_poly_appx = energy_poly_appx + energy_poly_coeff[index]*energy_poly_r;
            }
            ecoul = prefactor*energy_poly_appx;
          } else {
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

double PairCoulEspOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairCoulEsp::memory_usage();

  return bytes;
}
