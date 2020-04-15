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
#include "pair_lubricate_poly_omp.h"
#include <cmath>
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "input.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "update.h"
#include "variable.h"
#include "random_mars.h"
#include "fix_wall.h"
#include "fix_deform.h"
#include "math_const.h"

#include "suffix.h"
using namespace LAMMPS_NS;
using namespace MathConst;

// same as fix_wall.cpp

enum{EDGE,CONSTANT,VARIABLE};

/* ---------------------------------------------------------------------- */

PairLubricatePolyOMP::PairLubricatePolyOMP(LAMMPS *lmp) :
  PairLubricatePoly(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairLubricatePolyOMP::~PairLubricatePolyOMP()
{}

/* ---------------------------------------------------------------------- */

void PairLubricatePolyOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;


  // This section of code adjusts R0/RT0/RS0 if necessary due to changes
  // in the volume fraction as a result of fix deform or moving walls

  double dims[3], wallcoord;
  if (flagVF) // Flag for volume fraction corrections
    if (flagdeform || flagwall == 2){ // Possible changes in volume fraction
      if (flagdeform && !flagwall)
        for (int j = 0; j < 3; j++)
          dims[j] = domain->prd[j];
      else if (flagwall == 2 || (flagdeform && flagwall == 1)){
         double wallhi[3], walllo[3];
         for (int j = 0; j < 3; j++){
           wallhi[j] = domain->prd[j];
           walllo[j] = 0;
         }
         for (int m = 0; m < wallfix->nwall; m++){
           int dim = wallfix->wallwhich[m] / 2;
           int side = wallfix->wallwhich[m] % 2;
           if (wallfix->xstyle[m] == VARIABLE){
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
        R0  = 6*MY_PI*mu*(1.0 + 2.16*vol_f);
        RT0 = 8*MY_PI*mu;
        RS0 = 20.0/3.0*MY_PI*mu*(1.0 + 3.33*vol_f + 2.80*vol_f*vol_f);
      } else {
        R0  = 6*MY_PI*mu*(1.0 + 2.725*vol_f - 6.583*vol_f*vol_f);
        RT0 = 8*MY_PI*mu*(1.0 + 0.749*vol_f - 2.469*vol_f*vol_f);
        RS0 = 20.0/3.0*MY_PI*mu*(1.0 + 3.64*vol_f - 6.95*vol_f*vol_f);
      }
    }

  // end of R0 adjustment code


#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, NULL, thr);

    if (flaglog) {
      if (shearing) {
        if (evflag)
          eval<1,1,1>(ifrom, ito, thr);
        else eval<1,1,0>(ifrom, ito, thr);
      } else {
        if (evflag)
          eval<1,0,1>(ifrom, ito, thr);
        else eval<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (shearing) {
        if (evflag)
          eval<0,1,1>(ifrom, ito, thr);
        else eval<0,1,0>(ifrom, ito, thr);
      } else {
        if (evflag)
          eval<0,0,1>(ifrom, ito, thr);
        else eval<0,0,0>(ifrom, ito, thr);
      }
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int FLAGLOG, int SHEARING, int EVFLAG>
void PairLubricatePolyOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz,tx,ty,tz;
  double rsq,r,h_sep,beta0,beta1,radi,radj;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3;
  double vt1,vt2,vt3,wt1,wt2,wt3,wdotn;
  double vRS0;
  double vi[3],vj[3],wi[3],wj[3],xl[3],jl[3];
  double a_sq,a_sh,a_pu;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double lamda[3],vstream[3];

  double vxmu2f = force->vxmu2f;

  double * const * const x = atom->x;
  double * const * const v = atom->v;
  double * const * const f = thr->get_f();
  double * const * const omega = atom->omega;
  double * const * const torque = thr->get_torque();
  const double * const radius = atom->radius;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;

  int overlaps = 0;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // subtract streaming component of velocity, omega, angmom
  // assume fluid streaming velocity = box deformation rate
  // vstream = (ux,uy,uz)
  // ux = h_rate[0]*x + h_rate[5]*y + h_rate[4]*z
  // uy = h_rate[1]*y + h_rate[3]*z
  // uz = h_rate[2]*z
  // omega_new = omega - curl(vstream)/2
  // angmom_new = angmom - I*curl(vstream)/2
  // Ef = (grad(vstream) + (grad(vstream))^T) / 2

  if (shearing) {
    double *h_rate = domain->h_rate;
    double *h_ratelo = domain->h_ratelo;

    for (ii = iifrom; ii < iito; ii++) {
      i = ilist[ii];
      itype = type[i];
      radi = radius[i];
      domain->x2lamda(x[i],lamda);
      vstream[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] +
        h_rate[4]*lamda[2] + h_ratelo[0];
      vstream[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
      vstream[2] = h_rate[2]*lamda[2] + h_ratelo[2];
      v[i][0] -= vstream[0];
      v[i][1] -= vstream[1];
      v[i][2] -= vstream[2];

      omega[i][0] += 0.5*h_rate[3];
      omega[i][1] -= 0.5*h_rate[4];
      omega[i][2] += 0.5*h_rate[5];
    }

    // set Ef from h_rate in strain units

    Ef[0][0] = h_rate[0]/domain->xprd;
    Ef[1][1] = h_rate[1]/domain->yprd;
    Ef[2][2] = h_rate[2]/domain->zprd;
    Ef[0][1] = Ef[1][0] = 0.5 * h_rate[5]/domain->yprd;
    Ef[0][2] = Ef[2][0] = 0.5 * h_rate[4]/domain->zprd;
    Ef[1][2] = Ef[2][1] = 0.5 * h_rate[3]/domain->zprd;

    // copy updated omega to the ghost particles
    // no need to do this if not shearing since comm->ghost_velocity is set

    sync_threads();

    // MPI communication only on master thread
#if defined(_OPENMP)
#pragma omp master
#endif
    { comm->forward_comm_pair(this); }

    sync_threads();
  }

  // R0 adjustment has already been done in this->compute()

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // angular velocity

    wi[0] = omega[i][0];
    wi[1] = omega[i][1];
    wi[2] = omega[i][2];

    // FLD contribution to force and torque due to isotropic terms
    // FLD contribution to stress from isotropic RS0

    if (flagfld) {
      f[i][0] -= vxmu2f*R0*radi*v[i][0];
      f[i][1] -= vxmu2f*R0*radi*v[i][1];
      f[i][2] -= vxmu2f*R0*radi*v[i][2];
      const double rad3 = radi*radi*radi;
      torque[i][0] -= vxmu2f*RT0*rad3*wi[0];
      torque[i][1] -= vxmu2f*RT0*rad3*wi[1];
      torque[i][2] -= vxmu2f*RT0*rad3*wi[2];

      if (SHEARING && vflag_either) {
        vRS0 = -vxmu2f * RS0*rad3;
        v_tally_tensor(i,i,nlocal,/* newton_pair */ 0,
                       vRS0*Ef[0][0],vRS0*Ef[1][1],vRS0*Ef[2][2],
                       vRS0*Ef[0][1],vRS0*Ef[0][2],vRS0*Ef[1][2]);
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
      radj = atom->radius[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);

        // angular momentum = I*omega = 2/5 * M*R^2 * omega

        wj[0] = omega[j][0];
        wj[1] = omega[j][1];
        wj[2] = omega[j][2];

        // xl = point of closest approach on particle i from its center

        xl[0] = -delx/r*radi;
        xl[1] = -dely/r*radi;
        xl[2] = -delz/r*radi;
        jl[0] = -delx/r*radj;
        jl[1] = -dely/r*radj;
        jl[2] = -delz/r*radj;

        // velocity at the point of closest approach on both particles
        // v = v + omega_cross_xl - Ef.xl

        // particle i

        vi[0] = v[i][0] + (wi[1]*xl[2] - wi[2]*xl[1])
                        - (Ef[0][0]*xl[0] + Ef[0][1]*xl[1] + Ef[0][2]*xl[2]);

        vi[1] = v[i][1] + (wi[2]*xl[0] - wi[0]*xl[2])
                        - (Ef[1][0]*xl[0] + Ef[1][1]*xl[1] + Ef[1][2]*xl[2]);

        vi[2] = v[i][2] + (wi[0]*xl[1] - wi[1]*xl[0])
                        - (Ef[2][0]*xl[0] + Ef[2][1]*xl[1] + Ef[2][2]*xl[2]);

        // particle j

        vj[0] = v[j][0] - (wj[1]*jl[2] - wj[2]*jl[1])
                        + (Ef[0][0]*jl[0] + Ef[0][1]*jl[1] + Ef[0][2]*jl[2]);

        vj[1] = v[j][1] - (wj[2]*jl[0] - wj[0]*jl[2])
                        + (Ef[1][0]*jl[0] + Ef[1][1]*jl[1] + Ef[1][2]*jl[2]);

        vj[2] = v[j][2] - (wj[0]*jl[1] - wj[1]*jl[0])
                        + (Ef[2][0]*jl[0] + Ef[2][1]*jl[1] + Ef[2][2]*jl[2]);

        // scalar resistances XA and YA

        h_sep = r - radi-radj;

        // check for overlaps

        if (h_sep < 0.0) overlaps++;

        // if less than the minimum gap use the minimum gap instead

        if (r < cut_inner[itype][jtype])
          h_sep = cut_inner[itype][jtype] - radi-radj;

        // scale h_sep by radi

        h_sep = h_sep/radi;
        beta0 = radj/radi;
        beta1 = 1.0 + beta0;

        // scalar resistances

        if (FLAGLOG) {
          a_sq = beta0*beta0/beta1/beta1/h_sep +
            (1.0+7.0*beta0+beta0*beta0)/5.0/pow(beta1,3.0)*log(1.0/h_sep);
          a_sq += (1.0+18.0*beta0-29.0*beta0*beta0+18.0 *
                   pow(beta0,3.0)+pow(beta0,4.0))/21.0/pow(beta1,4.0) *
            h_sep*log(1.0/h_sep);
          a_sq *= 6.0*MY_PI*mu*radi;
          a_sh = 4.0*beta0*(2.0+beta0+2.0*beta0*beta0)/15.0/pow(beta1,3.0) *
            log(1.0/h_sep);
          a_sh += 4.0*(16.0-45.0*beta0+58.0*beta0*beta0-45.0*pow(beta0,3.0) +
                       16.0*pow(beta0,4.0))/375.0/pow(beta1,4.0) *
            h_sep*log(1.0/h_sep);
          a_sh *= 6.0*MY_PI*mu*radi;
          // old invalid eq for pumping term
          // changed 29Jul16 from eq 9.25 -> 9.27 in Kim and Karilla
//          a_pu = beta0*(4.0+beta0)/10.0/beta1/beta1*log(1.0/h_sep);
//          a_pu += (32.0-33.0*beta0+83.0*beta0*beta0+43.0 *
//                   pow(beta0,3.0))/250.0/pow(beta1,3.0)*h_sep*log(1.0/h_sep);
//          a_pu *= 8.0*MY_PI*mu*pow(radi,3.0);
          a_pu = 2.0*beta0/5.0/beta1*log(1.0/h_sep);
          a_pu += 2.0*(8.0+6.0*beta0+33.0*beta0*beta0)/125.0/beta1/beta1*
                   h_sep*log(1.0/h_sep);
          a_pu *= 8.0*MY_PI*mu*pow(radi,3.0);
        } else a_sq = 6.0*MY_PI*mu*radi*(beta0*beta0/beta1/beta1/h_sep);

        // relative velocity at the point of closest approach
        // includes fluid velocity

        vr1 = vi[0] - vj[0];
        vr2 = vi[1] - vj[1];
        vr3 = vi[2] - vj[2];

        // normal component (vr.n)n

        vnnr = (vr1*delx + vr2*dely + vr3*delz)/r;
        vn1 = vnnr*delx/r;
        vn2 = vnnr*dely/r;
        vn3 = vnnr*delz/r;

        // tangential component vr - (vr.n)n

        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // force due to squeeze type motion

        fx  = a_sq*vn1;
        fy  = a_sq*vn2;
        fz  = a_sq*vn3;

        // force due to all shear kind of motions

        if (FLAGLOG) {
          fx = fx + a_sh*vt1;
          fy = fy + a_sh*vt2;
          fz = fz + a_sh*vt3;
        }

        // scale forces for appropriate units

        fx *= vxmu2f;
        fy *= vxmu2f;
        fz *= vxmu2f;

        // add to total force

        f[i][0] -= fx;
        f[i][1] -= fy;
        f[i][2] -= fz;

        // torque due to this force

        if (FLAGLOG) {
          tx = xl[1]*fz - xl[2]*fy;
          ty = xl[2]*fx - xl[0]*fz;
          tz = xl[0]*fy - xl[1]*fx;

          torque[i][0] -= vxmu2f*tx;
          torque[i][1] -= vxmu2f*ty;
          torque[i][2] -= vxmu2f*tz;

          // torque due to a_pu

          wdotn = ((wi[0]-wj[0])*delx + (wi[1]-wj[1])*dely +
                   (wi[2]-wj[2])*delz)/r;
          wt1 = (wi[0]-wj[0]) - wdotn*delx/r;
          wt2 = (wi[1]-wj[1]) - wdotn*dely/r;
          wt3 = (wi[2]-wj[2]) - wdotn*delz/r;

          tx = a_pu*wt1;
          ty = a_pu*wt2;
          tz = a_pu*wt3;

          torque[i][0] -= vxmu2f*tx;
          torque[i][1] -= vxmu2f*ty;
          torque[i][2] -= vxmu2f*tz;

        }

        if (EVFLAG) ev_tally_xyz(i,nlocal,nlocal, /* newton_pair */ 0,
                                 0.0,0.0,-fx,-fy,-fz,delx,dely,delz);
      }
    }
  }

  // restore streaming component of velocity, omega, angmom

  if (SHEARING) {
    double *h_rate = domain->h_rate;
    double *h_ratelo = domain->h_ratelo;

    for (ii = iifrom; ii < iito; ii++) {
      i = ilist[ii];
      itype = type[i];
      radi = radius[i];

      domain->x2lamda(x[i],lamda);
      vstream[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] +
        h_rate[4]*lamda[2] + h_ratelo[0];
      vstream[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
      vstream[2] = h_rate[2]*lamda[2] + h_ratelo[2];
      v[i][0] += vstream[0];
      v[i][1] += vstream[1];
      v[i][2] += vstream[2];

      omega[i][0] -= 0.5*h_rate[3];
      omega[i][1] += 0.5*h_rate[4];
      omega[i][2] -= 0.5*h_rate[5];
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairLubricatePolyOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLubricatePoly::memory_usage();

  return bytes;
}
