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

#include "pair_peri_lps_omp.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "fix_peri_neigh.h"
#include "force.h"
#include "lattice.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairPeriLPSOMP::PairPeriLPSOMP(LAMMPS *lmp) :
  PairPeriLPS(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairPeriLPSOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  // grow bond forces array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(s0_new);
    memory->destroy(theta);
    nmax = atom->nmax;
    memory->create(s0_new,nmax,"pair:s0_new");
    memory->create(theta,nmax,"pair:theta");
  }

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

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairPeriLPSOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double xtmp0,ytmp0,ztmp0,delx0,dely0,delz0,rsq0;
  double rsq,r,dr,rk,evdwl,fpair,fbond;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double d_ij,delta,stretch;

  evdwl = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  double fxtmp,fytmp,fztmp;

  double *vfrac = atom->vfrac;
  double *s0 = atom->s0;
  double **x0 = atom->x0;
  double **r0   = fix_peri_neigh->r0;
  tagint **partner = fix_peri_neigh->partner;
  int *npartner = fix_peri_neigh->npartner;
  double *wvolume = fix_peri_neigh->wvolume;

  // lc = lattice constant
  // init_style guarantees it's the same in x, y, and z

  double lc = domain->lattice->xlattice;
  double half_lc = 0.5*lc;
  double vfrac_scale = 1.0;

  // short-range forces

  int periodic = (domain->xperiodic || domain->yperiodic || domain->zperiodic);

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // need minimg() for x0 difference since not ghosted

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    xtmp0 = x0[i][0];
    ytmp0 = x0[i][1];
    ztmp0 = x0[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx*delx + dely*dely + delz*delz;
      delx0 = xtmp0 - x0[j][0];
      dely0 = ytmp0 - x0[j][1];
      delz0 = ztmp0 - x0[j][2];
      if (periodic) domain->minimum_image(delx0,dely0,delz0);
      rsq0 = delx0*delx0 + dely0*dely0 + delz0*delz0;
      jtype = type[j];

      r = sqrt(rsq);

      // short-range interaction distance based on initial particle position
      // 0.9 and 1.35 are constants

      d_ij = MIN(0.9*sqrt(rsq0),1.35*lc);

      // short-range contact forces
      // 15 is constant taken from the EMU Theory Manual
      // Silling, 12 May 2005, p 18

      if (r < d_ij) {
        dr = r - d_ij;

        // kshort based upon short-range force constant
        // of the bond-based theory used in PMB model

        double kshort = (15.0 * 18.0 * bulkmodulus[itype][itype]) /
          (MY_PI * cutsq[itype][jtype] * cutsq[itype][jtype]);
        rk = (kshort * vfrac[j]) * (dr / cut[itype][jtype]);

        if (r > 0.0) fpair = -(rk/r);
        else fpair = 0.0;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (EFLAG) evdwl = 0.5*rk*dr;
        if (EVFLAG) ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,evdwl,0.0,
                                 fpair*vfrac[i],delx,dely,delz,thr);
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }

  // wait until all threads are done since we
  // need to distribute the work differently.
  sync_threads();

#if defined(_OPENMP)
  // each thread works on a fixed chunk of atoms.
  const int idelta = 1 + nlocal/comm->nthreads;
  iifrom = thr->get_tid()*idelta;
  iito   = ((iifrom + idelta) > nlocal) ? nlocal : (iifrom + idelta);
#else
  iifrom = 0;
  iito = nlocal;
#endif

  // Compute the dilatation on each particle
  if (iifrom < nlocal)
    compute_dilatation(iifrom, iito);

  // wait until all threads are done before communication
  sync_threads();

#if defined(_OPENMP)
#pragma omp master
#endif
  { // communicate dilatation (theta) of each particle
    comm->forward_comm(this);
    // communicate weighted volume (wvolume) upon every reneighbor
    if (neighbor->ago == 0)
      comm->forward_comm(fix_peri_neigh);
  }

  sync_threads();

  // Volume-dependent part of the energy
  if (EFLAG) {
    for (i = iifrom; i < iito; i++) {
      itype = type[i];
      e_tally_thr(this, i, i, nlocal, NEWTON_PAIR,
                  0.5 * bulkmodulus[itype][itype] * (theta[i] * theta[i]), 0.0, thr);
    }
  }

  // loop over my particles and their partners
  // partner list contains all bond partners, so I-J appears twice
  // if bond already broken, skip this partner
  // first = true if this is first neighbor of particle i

  bool first;
  double omega_minus, omega_plus;

  for (i = iifrom; i < iito; ++i) {
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    xtmp0 = x0[i][0];
    ytmp0 = x0[i][1];
    ztmp0 = x0[i][2];
    itype = type[i];
    jnum = npartner[i];
    first = true;

    for (jj = 0; jj < jnum; jj++) {
      if (partner[i][jj] == 0) continue;
      j = atom->map(partner[i][jj]);

      // check if lost a partner without first breaking bond

      if (j < 0) {
        partner[i][jj] = 0;
        continue;
      }

      // compute force density, add to PD equation of motion

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      if (periodic) domain->minimum_image(delx,dely,delz);
      rsq = delx*delx + dely*dely + delz*delz;
      delx0 = xtmp0 - x0[j][0];
      dely0 = ytmp0 - x0[j][1];
      delz0 = ztmp0 - x0[j][2];
      if (periodic) domain->minimum_image(delx0,dely0,delz0);
      jtype = type[j];
      delta = cut[itype][jtype];
      r = sqrt(rsq);
      dr = r - r0[i][jj];

      // avoid roundoff errors

      if (fabs(dr) < 2.2204e-016) dr = 0.0;

      // scale vfrac[j] if particle j near the horizon

      if ((fabs(r0[i][jj] - delta)) <= half_lc)
        vfrac_scale = (-1.0/(2*half_lc))*(r0[i][jj]) +
          (1.0 + ((delta - half_lc)/(2*half_lc) ) );
      else vfrac_scale = 1.0;

      omega_plus  = influence_function(-1.0*delx0,-1.0*dely0,-1.0*delz0);
      omega_minus = influence_function(delx0,dely0,delz0);
      if ((wvolume[i] > 0.0) && (wvolume[j] > 0.0)) {
        rk = ( (3.0 * bulkmodulus[itype][itype]) -
               (5.0 * shearmodulus[itype][itype]) ) * vfrac[j] * vfrac_scale *
          ( (omega_plus * theta[i] / wvolume[i]) +
            ( omega_minus * theta[j] / wvolume[j] ) ) * r0[i][jj];
        rk +=  15.0 * ( shearmodulus[itype][itype] * vfrac[j] * vfrac_scale ) *
          ( (omega_plus / wvolume[i]) + (omega_minus / wvolume[j]) ) * dr;
      } else rk = 0.0;

      if (r > 0.0) fbond = -(rk/r);
      else fbond = 0.0;

      f[i][0] += delx*fbond;
      f[i][1] += dely*fbond;
      f[i][2] += delz*fbond;

      // since I-J is double counted, set newton off & use 1/2 factor and I,I

      double deviatoric_extension = dr - (theta[i]* r0[i][jj] / 3.0);
      if (EFLAG && (wvolume[i] > 0.0))
        evdwl = 0.5 * 15 * (shearmodulus[itype][itype]/wvolume[i]) *
                   omega_plus*(deviatoric_extension * deviatoric_extension) *
                   vfrac[j] * vfrac_scale;
      else evdwl = 0.0;
      if (EVFLAG) ev_tally_thr(this,i,i,nlocal,0,0.5*evdwl,0.0,
                               0.5*fbond*vfrac[i],delx,dely,delz,thr);

      // find stretch in bond I-J and break if necessary
      // use s0 from previous timestep

      stretch = dr / r0[i][jj];
      if (stretch > MIN(s0[i],s0[j])) partner[i][jj] = 0;

      // update s0 for next timestep

      if (first)
         s0_new[i] = s00[itype][jtype] - (alpha[itype][jtype] * stretch);
      else
         s0_new[i] = MAX(s0_new[i],s00[itype][jtype] - (alpha[itype][jtype] * stretch));

      first = false;
    }
  }

  sync_threads();

  // store new s0 (in parallel)
  if (iifrom < nlocal)
    for (i = iifrom; i < iito; i++) s0[i] = s0_new[i];
}

/* ---------------------------------------------------------------------- */

double PairPeriLPSOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairPeriLPS::memory_usage();

  return bytes;
}
