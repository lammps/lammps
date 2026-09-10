// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_coul_long_cs_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"

using namespace LAMMPS_NS;

static constexpr double EWALD_F = 1.12837917;
static constexpr double EWALD_P = 9.95473818e-1;
static constexpr double B0 = -0.1335096380159268;
static constexpr double B1 = -2.57839507e-1;
static constexpr double B2 = -1.37203639e-1;
static constexpr double B3 = -8.88822059e-3;
static constexpr double B4 = -5.80844129e-3;
static constexpr double B5 = 1.14652755e-1;

static constexpr double EPSILON = 1.0e-20;
static constexpr double EPS_EWALD = 1.0e-6;
static constexpr double EPS_EWALD_SQR = 1.0e-12;

/* ---------------------------------------------------------------------- */

PairCoulLongCSOMP::PairCoulLongCSOMP(LAMMPS *lmp) :
  PairCoulLongCS(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  cut_respa = nullptr;
}

/* ---------------------------------------------------------------------- */

void PairCoulLongCSOMP::compute(int eflag, int vflag)
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
void PairCoulLongCSOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itable,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double fraction,table;
  double r,r2inv,rsq,forcecoul,factor_coul;
  double grij,expm2,prefactor,t,erfc,u;
  int *ilist,*jlist,*numneigh,**firstneigh;

  ecoul = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const double * const q = atom->q;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_coul = force->special_coul;
  const double qqrd2e = force->qqrd2e;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cut_coulsq) {
        rsq += EPSILON; // Add Epsilon for case: r = 0; Interaction must be removed by special bond;
        r2inv = 1.0/rsq;
        if (!ncoultablebits || rsq <= tabinnersq) {
          r = sqrt(rsq);
          prefactor = qqrd2e * scale[itype][jtype] * qtmp*q[j];
          if (factor_coul < 1.0) {
            // When bonded parts are being calculated a minimal distance (EPS_EWALD)
            // has to be added to the prefactor and erfc in order to make the
            // used approximation functions for the Ewald correction valid
            grij = g_ewald * (r+EPS_EWALD);
            expm2 = exp(-grij*grij);
            t = 1.0 / (1.0 + EWALD_P*grij);
            u = 1.0 - t;
            erfc = t * (1.+u*(B0+u*(B1+u*(B2+u*(B3+u*(B4+u*B5)))))) * expm2;
            prefactor /= (r+EPS_EWALD);
            forcecoul = prefactor * (erfc + EWALD_F*grij*expm2 - (1.0-factor_coul));
            // Additionally r2inv needs to be accordingly modified since the later
            // scaling of the overall force shall be consistent
            r2inv = 1.0/(rsq + EPS_EWALD_SQR);
          } else {
            grij = g_ewald * r;
            expm2 = exp(-grij*grij);
            t = 1.0 / (1.0 + EWALD_P*grij);
            u = 1.0 - t;
            erfc = t * (1.+u*(B0+u*(B1+u*(B2+u*(B3+u*(B4+u*B5)))))) * expm2;
            prefactor /= r;
            forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
          }
        } else {
          union_int_float_t rsq_lookup;
          rsq_lookup.f = rsq;
          itable = rsq_lookup.i & ncoulmask;
          itable >>= ncoulshiftbits;
          fraction = ((double) rsq_lookup.f - rtable[itable]) * drtable[itable];
          table = ftable[itable] + fraction*dftable[itable];
          forcecoul = scale[itype][jtype] * qtmp*q[j] * table;
          if (factor_coul < 1.0) {
            table = ctable[itable] + fraction*dctable[itable];
            prefactor = scale[itype][jtype] * qtmp*q[j] * table;
            forcecoul -= (1.0-factor_coul)*prefactor;
          }
        }

        fpair = forcecoul * r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
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
  }
}

/* ---------------------------------------------------------------------- */

double PairCoulLongCSOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairCoulLongCS::memory_usage();

  return bytes;
}
