/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Randy Schunk (SNL)
                         Amit Kumar and Michael Bybee (UIUC)
                         Dave Heine (Corning), polydispersity
------------------------------------------------------------------------- */

#include "pair_lubricate_poly.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "fix_deform.h"
#include "fix_wall.h"
#include "input.h"
#include "variable.h"
#include "math_const.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

// same as fix_wall.cpp

enum{EDGE,CONSTANT,VARIABLE};

/* ---------------------------------------------------------------------- */

PairLubricatePoly::PairLubricatePoly(LAMMPS *lmp) : PairLubricate(lmp)
{
  no_virial_fdotr_compute = 1;
}

/* ---------------------------------------------------------------------- */

void PairLubricatePoly::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
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

  ev_init(eflag,vflag);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
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

    for (ii = 0; ii < inum; ii++) {
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

    comm->forward_comm_pair(this);
  }

  // This section of code adjusts R0/RT0/RS0 if necessary due to changes
  // in the volume fraction as a result of fix deform or moving walls

  double dims[3], wallcoord;
  if (flagVF) // Flag for volume fraction corrections
    if (flagdeform || flagwall == 2){ // Possible changes in volume fraction
      if (flagdeform && !flagwall)
        for (j = 0; j < 3; j++)
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

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    radi = radius[i];

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
      const double radi3 = radi*radi*radi;
      torque[i][0] -= vxmu2f*RT0*radi3*wi[0];
      torque[i][1] -= vxmu2f*RT0*radi3*wi[1];
      torque[i][2] -= vxmu2f*RT0*radi3*wi[2];

      if (shearing && vflag_either) {
        vRS0 = -vxmu2f * RS0*radi3;
        v_tally_tensor(i,i,nlocal,newton_pair,
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

        // if less than the minimum gap use the minimum gap instead

        if (r < cut_inner[itype][jtype])
          h_sep = cut_inner[itype][jtype] - radi-radj;

        // scale h_sep by radi

        h_sep = h_sep/radi;
        beta0 = radj/radi;
        beta1 = 1.0 + beta0;

        // scalar resistances

        if (flaglog) {
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

        if (flaglog) {
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

        if (flaglog) {
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

        // set j = nlocal so that only I gets tallied

        if (evflag) ev_tally_xyz(i,nlocal,nlocal,0,
                                 0.0,0.0,-fx,-fy,-fz,delx,dely,delz);
      }
    }
  }

  // restore streaming component of velocity, omega, angmom

  if (shearing) {
    double *h_rate = domain->h_rate;
    double *h_ratelo = domain->h_ratelo;

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      itype = type[i];
      radi = atom->radius[i];

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

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLubricatePoly::init_style()
{
  if (force->newton_pair == 1)
    error->all(FLERR,"Pair lubricate/poly requires newton pair off");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,
               "Pair lubricate/poly requires ghost atoms store velocity");
  if (!atom->sphere_flag)
    error->all(FLERR,"Pair lubricate/poly requires atom style sphere");

  // ensure all particles are finite-size
  // for pair hybrid, should limit test to types using the pair style

  double *radius = atom->radius;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (radius[i] == 0.0)
      error->one(FLERR,"Pair lubricate/poly requires extended particles");

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  // set the isotropic constants that depend on the volume fraction
  // vol_T = total volume

  // check for fix deform, if exists it must use "remap v"
  // If box will change volume, set appropriate flag so that volume
  // and v.f. corrections are re-calculated at every step.

  // if available volume is different from box volume
  // due to walls, set volume appropriately; if walls will
  // move, set appropriate flag so that volume and v.f. corrections
  // are re-calculated at every step.

  shearing = flagdeform = flagwall = 0;
  for (int i = 0; i < modify->nfix; i++){
    if (strcmp(modify->fix[i]->style,"deform") == 0) {
      shearing = flagdeform = 1;
      if (((FixDeform *) modify->fix[i])->remapflag != Domain::V_REMAP)
        error->all(FLERR,"Using pair lubricate with inconsistent "
                   "fix deform remap option");
    }
    if (strstr(modify->fix[i]->style,"wall") != NULL) {
      if (flagwall)
        error->all(FLERR,
                   "Cannot use multiple fix wall commands with "
                   "pair lubricate/poly");
      flagwall = 1; // Walls exist
      wallfix = (FixWall *) modify->fix[i];
      if (wallfix->xflag) flagwall = 2; // Moving walls exist
    }

    if (strstr(modify->fix[i]->style,"wall") != NULL){
      flagwall = 1; // Walls exist
      if (((FixWall *) modify->fix[i])->xflag ) {
        flagwall = 2; // Moving walls exist
        wallfix = (FixWall *) modify->fix[i];
      }
    }
  }

  double vol_T;
  double wallcoord;
  if (!flagwall) vol_T = domain->xprd*domain->yprd*domain->zprd;
  else {
    double wallhi[3], walllo[3];
    for (int j = 0; j < 3; j++){
      wallhi[j] = domain->prd[j];
      walllo[j] = 0;
    }
    for (int m = 0; m < wallfix->nwall; m++){
      int dim = wallfix->wallwhich[m] / 2;
      int side = wallfix->wallwhich[m] % 2;
      if (wallfix->xstyle[m] == VARIABLE){
        wallfix->xindex[m] = input->variable->find(wallfix->xstr[m]);
        //Since fix->wall->init happens after pair->init_style
        wallcoord = input->variable->compute_equal(wallfix->xindex[m]);
      }
      else wallcoord = wallfix->coord0[m];

      if (side == 0) walllo[dim] = wallcoord;
      else wallhi[dim] = wallcoord;
    }
    vol_T = (wallhi[0] - walllo[0]) * (wallhi[1] - walllo[1]) *
      (wallhi[2] - walllo[2]);
  }

  double volP = 0.0;
  for (int i = 0; i < nlocal; i++)
    volP += (4.0/3.0)*MY_PI*pow(atom->radius[i],3.0);
  MPI_Allreduce(&volP,&vol_P,1,MPI_DOUBLE,MPI_SUM,world);

  double vol_f = vol_P/vol_T;

  if (!flagVF) vol_f = 0;

  // set isotropic constants

  if (flaglog == 0) {
    R0  = 6*MY_PI*mu*(1.0 + 2.16*vol_f);
    RT0 = 8*MY_PI*mu;
    RS0 = 20.0/3.0*MY_PI*mu*(1.0 + 3.33*vol_f + 2.80*vol_f*vol_f);
  } else {
    R0  = 6*MY_PI*mu*(1.0 + 2.725*vol_f - 6.583*vol_f*vol_f);
    RT0 = 8*MY_PI*mu*(1.0 + 0.749*vol_f - 2.469*vol_f*vol_f);
    RS0 = 20.0/3.0*MY_PI*mu*(1.0 + 3.64*vol_f - 6.95*vol_f*vol_f);
  }

  // check for fix deform, if exists it must use "remap v"

  shearing = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"deform") == 0) {
      shearing = 1;
      if (((FixDeform *) modify->fix[i])->remapflag != Domain::V_REMAP)
        error->all(FLERR,"Using pair lubricate/poly with inconsistent "
                   "fix deform remap option");
    }

  // set Ef = 0 since used whether shearing or not

  Ef[0][0] = Ef[0][1] = Ef[0][2] = 0.0;
  Ef[1][0] = Ef[1][1] = Ef[1][2] = 0.0;
  Ef[2][0] = Ef[2][1] = Ef[2][2] = 0.0;
}
