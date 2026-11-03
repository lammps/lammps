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

/* ----------------------------------------------------------------------
   the base class settings() re-creates the serial random number generator,
   so the pool of per-thread generators must be re-created on the next call
------------------------------------------------------------------------- */

void PairBrownianPolyOMP::settings(int narg, char **arg)
{
  PairBrownianPoly::settings(narg, arg);

  if (random_thr) {
    for (int i = 1; i < nthreads; ++i) delete random_thr[i];
    delete[] random_thr;
    random_thr = nullptr;
  }
  nthreads = 0;
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

      // R0 and RT0 are multiplied by the radius of each particle (cubed for RT0)

      if (flaglog == 0) {
        R0  = 6*MY_PI*mu*(1.0 + 2.16*vol_f);
        RT0 = 8*MY_PI*mu;
      } else {
        R0  = 6*MY_PI*mu*(1.0 + 2.725*vol_f - 6.583*vol_f*vol_f);
        RT0 = 8*MY_PI*mu*(1.0 + 0.749*vol_f - 2.469*vol_f*vol_f);
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
    if ((tid > 0) && (random_thr[tid] == nullptr)) {
      // wrap around to stay within the valid seed range of RanMars
      const bigint tseed = (bigint) seed + comm->me + (bigint) comm->nprocs * tid;
      random_thr[tid] = new RanMars(Pair::lmp, (int) ((tseed - 1) % 900000000 + 1));
    }

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
  double rsq,r,gap,h_sep,beta0,beta1,radi,radj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  double * const * const torque = thr->get_torque();
  const double * const radius = atom->radius;
  const int * const type = atom->type;
  const tagint * const tag = atom->tag;
  const int nlocal = atom->nlocal;
  const bigint step = update->ntimestep;

  RanMars &rng = *random_thr[thr->get_tid()];

  double vxmu2f = force->vxmu2f;
  double randr;
  double prethermostat;
  double a_sq,a_sh,a_pu,Fbmag;
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

        // the random numbers and the reference frame of a pair do not depend on
        // which of the two atoms is atom i, so that the two atoms get equal and
        // opposite forces.  the atom with the lower ID is the reference atom,
        // csgn = 1 if atom i is the reference atom and -1 otherwise, and the
        // unit vector e points from the other atom to the reference atom.

        const double csgn = (tag[i] < tag[j]) ? 1.0 : -1.0;
        const double rref = (csgn > 0.0) ? radi : radj;
        const double roth = (csgn > 0.0) ? radj : radi;
        const double ex = csgn*delx/r;
        const double ey = csgn*dely/r;
        const double ez = csgn*delz/r;
        const uint64_t key = pair_key(tag[i],tag[j],step,seed);

        // gap between the particle surfaces
        // if less than minimum gap, use minimum gap instead

        if (r < cut_inner[itype][jtype]) gap = cut_inner[itype][jtype] - (radi+radj);
        else gap = r - (radi+radj);

        // scale gap by the radius of the reference atom

        h_sep = gap/rref;
        beta0 = roth/rref;
        beta1 = 1.0 + beta0;

        // scalar resistances as in pair style lubricate/poly
        // log terms use the gap xi = 2*gap/(radi+radj) of Jeffrey and Onishi

        if (FLAGLOG) {
          const double logxi = log(0.5*beta1/h_sep);
          a_sq = beta0*beta0/beta1/beta1/h_sep +
            beta0*(1.0+7.0*beta0+beta0*beta0)/5.0/cube(beta1)*logxi;
          a_sq += (1.0+18.0*beta0-29.0*beta0*beta0+18.0*cube(beta0) +
                   powint(beta0,4))/21.0/powint(beta1,4)*h_sep*logxi;
          a_sq *= 6.0*MY_PI*mu*rref;
          a_sh = 4.0*beta0*(2.0+beta0+2.0*beta0*beta0)/15.0/cube(beta1)*logxi;
          a_sh += 4.0*(16.0-45.0*beta0+58.0*beta0*beta0-45.0*cube(beta0) +
                       16.0*powint(beta0,4))/375.0/powint(beta1,4)*h_sep*logxi;
          a_sh *= 6.0*MY_PI*mu*rref;

          // rotational resistance of atom i

          const double beta0i = radj/radi;
          const double beta1i = 1.0 + beta0i;
          a_pu = 2.0*beta0i/5.0/beta1i*logxi;
          a_pu += 2.0*(8.0+6.0*beta0i+33.0*beta0i*beta0i)/125.0/square(beta1i) *
            gap/radi*logxi;
          a_pu *= 8.0*MY_PI*mu*cube(radi);
        } else a_sq = 6.0*MY_PI*mu*rref*(beta0*beta0/beta1/beta1/h_sep);

        // generate the Pairwise Brownian Force: a_sq

        Fbmag = prethermostat*sqrt(a_sq);

        // generate a random number

        randr = pair_uniform(key,0)-0.5;

        // contribution due to Brownian motion

        fx = Fbmag*randr*ex;
        fy = Fbmag*randr*ey;
        fz = Fbmag*randr*ez;

        // add terms due to a_sh

        if (FLAGLOG) {

          // generate two orthogonal vectors to the line of centers

          p1[0] = ex; p1[1] = ey; p1[2] = ez;
          set_3_orthogonal_vectors(p1,p2,p3);

          // magnitude

          Fbmag = prethermostat*sqrt(a_sh);

          // force in each of the two directions

          randr = pair_uniform(key,1)-0.5;
          fx += Fbmag*randr*p2[0];
          fy += Fbmag*randr*p2[1];
          fz += Fbmag*randr*p2[2];

          randr = pair_uniform(key,2)-0.5;
          fx += Fbmag*randr*p3[0];
          fy += Fbmag*randr*p3[1];
          fz += Fbmag*randr*p3[2];
        }

        // scale forces to appropriate units

        fx = vxmu2f*fx;
        fy = vxmu2f*fy;
        fz = vxmu2f*fz;

        // sum to total force: -F on the reference atom, +F on the other atom

        f[i][0] -= csgn*fx;
        f[i][1] -= csgn*fy;
        f[i][2] -= csgn*fz;

        // torque due to the Brownian Force

        if (FLAGLOG) {

          // torque = xl_cross_F with xl = -csgn*radi*e, the location of the
          // point of closest approach on I from its center, and -csgn*F

          tx = radi*(ey*fz - ez*fy);
          ty = radi*(ez*fx - ex*fz);
          tz = radi*(ex*fy - ey*fx);

          torque[i][0] += tx;
          torque[i][1] += ty;
          torque[i][2] += tz;

          // torque due to a_pu

          Fbmag = prethermostat*sqrt(a_pu);

          // force in each direction

          randr = pair_uniform(key,3)-0.5;
          tx = Fbmag*randr*p2[0];
          ty = Fbmag*randr*p2[1];
          tz = Fbmag*randr*p2[2];

          randr = pair_uniform(key,4)-0.5;
          tx += Fbmag*randr*p3[0];
          ty += Fbmag*randr*p3[1];
          tz += Fbmag*randr*p3[2];

          // torque has opposite sign on two particles

          torque[i][0] -= csgn*tx;
          torque[i][1] -= csgn*ty;
          torque[i][2] -= csgn*tz;
        }

        // set j = nlocal so that only I gets tallied

        if (EVFLAG) ev_tally_xyz_thr(this,i,nlocal,nlocal,/* newton_pair */ 0,
                                     0.0,0.0,-csgn*fx,-csgn*fy,-csgn*fz,delx,dely,delz,thr);
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
