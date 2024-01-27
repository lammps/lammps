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

#include "pair_brownian_poly_omp.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "fix_wall.h"
#include "force.h"
#include "input.h"
#include "math_const.h"
#include "math_special.h"
#include "neigh_list.h"
#include "random_mars.h"
#include "suffix.h"
#include "update.h"
#include "variable.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

static constexpr double EPSILON = 1.0e-10;

/* ---------------------------------------------------------------------- */

PairBrownianPolyOMP::PairBrownianPolyOMP(LAMMPS *lmp) :
  PairBrownianPoly(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  random_thr = nullptr;
  nthreads = 0;
}

/* ---------------------------------------------------------------------- */

PairBrownianPolyOMP::~PairBrownianPolyOMP()
{
  if (random_thr) {
    for (int i=1; i < nthreads; ++i)
      delete random_thr[i];

    delete[] random_thr;
    random_thr = nullptr;
  }
}

/* ---------------------------------------------------------------------- */

void PairBrownianPolyOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int inum = list->inum;

  // This section of code adjusts R0/RT0/RS0 if necessary due to changes
  // in the volume fraction as a result of fix deform or moving walls

  double dims[3], wallcoord;
  if (flagVF) // Flag for volume fraction corrections
    if (flagdeform || flagwall == 2) { // Possible changes in volume fraction
      if (flagdeform && !flagwall)
        for (int j = 0; j < 3; j++)
          dims[j] = domain->prd[j];
      else if (flagwall == 2 || (flagdeform && flagwall == 1)) {
        double wallhi[3], walllo[3];
        for (int j = 0; j < 3; j++) {
          wallhi[j] = domain->prd[j];
          walllo[j] = 0;
        }
        for (int m = 0; m < wallfix->nwall; m++) {
          int dim = wallfix->wallwhich[m] / 2;
          int side = wallfix->wallwhich[m] % 2;
          if (wallfix->xstyle[m] == FixWall::VARIABLE) {
            wallcoord = input->variable->compute_equal(wallfix->xindex[m]);
          }
          else wallcoord = wallfix->coord0[m];
          if (side == 0) walllo[dim] = wallcoord;
          else wallhi[dim] = wallcoord;
        }
        for (int j = 0; j < 3; j++)
          dims[j] = wallhi[j] - walllo[j];
      }
      double vol_T = dims[0]*dims[1]*dims[2];
      double vol_f = vol_P/vol_T;
      if (flaglog == 0) {
        R0  = 6*MY_PI*mu*rad*(1.0 + 2.16*vol_f);
        RT0 = 8*MY_PI*mu*cube(rad);
        //RS0 = 20.0/3.0*MY_PI*mu*pow(rad,3)*(1.0 + 3.33*vol_f + 2.80*vol_f*vol_f);
      } else {
        R0  = 6*MY_PI*mu*rad*(1.0 + 2.725*vol_f - 6.583*vol_f*vol_f);
        RT0 = 8*MY_PI*mu*cube(rad)*(1.0 + 0.749*vol_f - 2.469*vol_f*vol_f);
        //RS0 = 20.0/3.0*MY_PI*mu*pow(rad,3)*(1.0 + 3.64*vol_f - 6.95*vol_f*vol_f);
      }
    }


  // number of threads has changed. reallocate pool of pRNGs
  if (nthreads != comm->nthreads) {
    if (random_thr) {
      for (int i=1; i < nthreads; ++i)
        delete random_thr[i];

      delete[] random_thr;
    }

    nthreads = comm->nthreads;
    random_thr = new RanMars*[nthreads];
    for (int i=1; i < nthreads; ++i)
      random_thr[i] = nullptr;

    // to ensure full compatibility with the serial BrownianPoly style
    // we use is random number generator instance for thread 0
    random_thr[0] = random;
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

    // generate a random number generator instance for
    // all threads != 0. make sure we use unique seeds.
    if ((tid > 0) && (random_thr[tid] == nullptr))
      random_thr[tid] = new RanMars(Pair::lmp, seed + comm->me
                                    + comm->nprocs*tid);

    if (flaglog) {
      if (evflag)
        eval<1,1>(ifrom, ito, thr);
      else
        eval<1,0>(ifrom, ito, thr);
    } else {
      if (evflag)
        eval<0,1>(ifrom, ito, thr);
      else eval<0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int FLAGLOG, int EVFLAG>
void PairBrownianPolyOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz,tx,ty,tz;
  double rsq,r,h_sep,beta0,beta1,radi,radj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  double * const * const torque = thr->get_torque();
  const double * const radius = atom->radius;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;

  RanMars &rng = *random_thr[thr->get_tid()];

  double vxmu2f = force->vxmu2f;
  double randr;
  double prethermostat;
  double xl[3],a_sq,a_sh,a_pu,Fbmag;
  double p1[3],p2[3],p3[3];

  // scale factor for Brownian moments

  prethermostat = sqrt(24.0*force->boltz*t_target/update->dt);
  prethermostat *= sqrt(force->vxmu2f/force->ftm2v/force->mvv2e);

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // FLD contribution to force and torque due to isotropic terms

    if (flagfld) {
      f[i][0] += prethermostat*sqrt(R0*radi)*(rng.uniform()-0.5);
      f[i][1] += prethermostat*sqrt(R0*radi)*(rng.uniform()-0.5);
      f[i][2] += prethermostat*sqrt(R0*radi)*(rng.uniform()-0.5);
      if (FLAGLOG) {
        const double rad3 = radi*radi*radi;
        torque[i][0] += prethermostat*sqrt(RT0*rad3)*(rng.uniform()-0.5);
        torque[i][1] += prethermostat*sqrt(RT0*rad3)*(rng.uniform()-0.5);
        torque[i][2] += prethermostat*sqrt(RT0*rad3)*(rng.uniform()-0.5);
      }
    }

    if (!flagHI) continue;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      radj = radius[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);

        // scalar resistances a_sq and a_sh

        h_sep = r - radi-radj;

        // if less than minimum gap, use minimum gap instead

        if (r < cut_inner[itype][jtype])
          h_sep = cut_inner[itype][jtype] - radi-radj;

        // scale h_sep by radi

        h_sep = h_sep/radi;
        beta0 = radj/radi;
        beta1 = 1.0 + beta0;

        // scalar resistances

        if (FLAGLOG) {
          a_sq = beta0*beta0/beta1/beta1/h_sep +
            (1.0+7.0*beta0+beta0*beta0)/5.0/cube(beta1)*log(1.0/h_sep);
          a_sq += (1.0+18.0*beta0-29.0*beta0*beta0+18.0*cube(beta0) +
                   powint(beta0,4))/21.0/powint(beta1,4)*h_sep*log(1.0/h_sep);
          a_sq *= 6.0*MY_PI*mu*radi;
          a_sh = 4.0*beta0*(2.0+beta0+2.0*beta0*beta0)/15.0/cube(beta1) *
            log(1.0/h_sep);
          a_sh += 4.0*(16.0-45.0*beta0+58.0*beta0*beta0-45.0*cube(beta0) +
                       16.0*powint(beta0,4))/375.0/powint(beta1,4) *
            h_sep*log(1.0/h_sep);
          a_sh *= 6.0*MY_PI*mu*radi;
          a_pu = beta0*(4.0+beta0)/10.0/beta1/beta1*log(1.0/h_sep);
          a_pu += (32.0-33.0*beta0+83.0*beta0*beta0+43.0 *
                   cube(beta0))/250.0/cube(beta1)*h_sep*log(1.0/h_sep);
          a_pu *= 8.0*MY_PI*mu*cube(radi);

        } else a_sq = 6.0*MY_PI*mu*radi*(beta0*beta0/beta1/beta1/h_sep);

        // generate the Pairwise Brownian Force: a_sq

        Fbmag = prethermostat*sqrt(a_sq);

        // generate a random number

        randr = rng.uniform()-0.5;

        // contribution due to Brownian motion

        fx = Fbmag*randr*delx/r;
        fy = Fbmag*randr*dely/r;
        fz = Fbmag*randr*delz/r;

        // add terms due to a_sh

        if (FLAGLOG) {

          // generate two orthogonal vectors to the line of centers

          p1[0] = delx/r; p1[1] = dely/r; p1[2] = delz/r;
          set_3_orthogonal_vectors(p1,p2,p3);

          // magnitude

          Fbmag = prethermostat*sqrt(a_sh);

          // force in each of the two directions

          randr = rng.uniform()-0.5;
          fx += Fbmag*randr*p2[0];
          fy += Fbmag*randr*p2[1];
          fz += Fbmag*randr*p2[2];

          randr = rng.uniform()-0.5;
          fx += Fbmag*randr*p3[0];
          fy += Fbmag*randr*p3[1];
          fz += Fbmag*randr*p3[2];
        }

        // scale forces to appropriate units

        fx = vxmu2f*fx;
        fy = vxmu2f*fy;
        fz = vxmu2f*fz;

        // sum to total force

        f[i][0] -= fx;
        f[i][1] -= fy;
        f[i][2] -= fz;

        // torque due to the Brownian Force

        if (FLAGLOG) {

          // location of the point of closest approach on I from its center

          xl[0] = -delx/r*radi;
          xl[1] = -dely/r*radi;
          xl[2] = -delz/r*radi;

          // torque = xl_cross_F

          tx = xl[1]*fz - xl[2]*fy;
          ty = xl[2]*fx - xl[0]*fz;
          tz = xl[0]*fy - xl[1]*fx;

          // torque is same on both particles

          torque[i][0] -= tx;
          torque[i][1] -= ty;
          torque[i][2] -= tz;

          // torque due to a_pu

          Fbmag = prethermostat*sqrt(a_pu);

          // force in each direction

          randr = rng.uniform()-0.5;
          tx = Fbmag*randr*p2[0];
          ty = Fbmag*randr*p2[1];
          tz = Fbmag*randr*p2[2];

          randr = rng.uniform()-0.5;
          tx += Fbmag*randr*p3[0];
          ty += Fbmag*randr*p3[1];
          tz += Fbmag*randr*p3[2];

          // torque has opposite sign on two particles

          torque[i][0] -= tx;
          torque[i][1] -= ty;
          torque[i][2] -= tz;

        }

        // set j = nlocal so that only I gets tallied

        if (EVFLAG) ev_tally_xyz(i,nlocal,nlocal,/* newton_pair */ 0,
                                 0.0,0.0,-fx,-fy,-fz,delx,dely,delz);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairBrownianPolyOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairBrownianPoly::memory_usage();
  bytes += (double)nthreads * sizeof(RanMars*);
  bytes += (double)nthreads * sizeof(RanMars);

  return bytes;
}
