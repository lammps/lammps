// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Amit Kumar and Michael Bybee (UIUC)
                         Dave Heine (Corning), polydispersity
------------------------------------------------------------------------- */

#include "pair_brownian_poly.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "fix_wall.h"
#include "force.h"
#include "input.h"
#include "math_const.h"
#include "math_special.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "random_mars.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

// same as fix_wall.cpp

enum{EDGE,CONSTANT,VARIABLE};

/* ---------------------------------------------------------------------- */

PairBrownianPoly::PairBrownianPoly(LAMMPS *lmp) : PairBrownian(lmp)
{
  no_virial_fdotr_compute = 1;
}

/* ---------------------------------------------------------------------- */

void PairBrownianPoly::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz,tx,ty,tz;
  double rsq,r,h_sep,beta0,beta1,radi,radj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double vxmu2f = force->vxmu2f;
  int overlaps = 0;
  double randr;
  double prethermostat;
  double xl[3],a_sq,a_sh,a_pu,Fbmag;
  double p1[3],p2[3],p3[3];

  // this section of code adjusts R0/RT0/RS0 if necessary due to changes
  // in the volume fraction as a result of fix deform or moving walls

  double dims[3], wallcoord;
  if (flagVF) // Flag for volume fraction corrections
    if (flagdeform || flagwall == 2) { // Possible changes in volume fraction
      if (flagdeform && !flagwall)
        for (j = 0; j < 3; j++)
          dims[j] = domain->prd[j];
      else if (flagwall == 2 || (flagdeform && flagwall == 1)) {
        double wallhi[3], walllo[3];
        for (j = 0; j < 3; j++) {
          wallhi[j] = domain->prd[j];
          walllo[j] = 0;
        }
        for (int m = 0; m < wallfix->nwall; m++) {
          int dim = wallfix->wallwhich[m] / 2;
          int side = wallfix->wallwhich[m] % 2;
          if (wallfix->xstyle[m] == VARIABLE) {
            wallcoord = input->variable->compute_equal(wallfix->xindex[m]);
          }
          else wallcoord = wallfix->coord0[m];
          if (side == 0) walllo[dim] = wallcoord;
          else wallhi[dim] = wallcoord;
        }
        for (j = 0; j < 3; j++)
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
  // scale factor for Brownian moments

  prethermostat = sqrt(24.0*force->boltz*t_target/update->dt);
  prethermostat *= sqrt(force->vxmu2f/force->ftm2v/force->mvv2e);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;


  for (ii = 0; ii < inum; ii++) {
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
      f[i][0] += prethermostat*sqrt(R0*radi)*(random->uniform()-0.5);
      f[i][1] += prethermostat*sqrt(R0*radi)*(random->uniform()-0.5);
      f[i][2] += prethermostat*sqrt(R0*radi)*(random->uniform()-0.5);
      if (flaglog) {
        const double radi3 = radi*radi*radi;
        torque[i][0] += prethermostat*sqrt(RT0*radi3) *
          (random->uniform()-0.5);
        torque[i][1] += prethermostat*sqrt(RT0*radi3) *
          (random->uniform()-0.5);
        torque[i][2] += prethermostat*sqrt(RT0*radi3) *
          (random->uniform()-0.5);
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

        // check for overlaps

        if (h_sep < 0.0) overlaps++;

        // if less than minimum gap, use minimum gap instead

        if (r < cut_inner[itype][jtype])
          h_sep = cut_inner[itype][jtype] - radi-radj;

        // scale h_sep by radi

        h_sep = h_sep/radi;
        beta0 = radj/radi;
        beta1 = 1.0 + beta0;

        // scalar resistances

        if (flaglog) {
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

        randr = random->uniform()-0.5;

        // contribution due to Brownian motion

        fx = Fbmag*randr*delx/r;
        fy = Fbmag*randr*dely/r;
        fz = Fbmag*randr*delz/r;

        // add terms due to a_sh

        if (flaglog) {

          // generate two orthogonal vectors to the line of centers

          p1[0] = delx/r; p1[1] = dely/r; p1[2] = delz/r;
          set_3_orthogonal_vectors(p1,p2,p3);

          // magnitude

          Fbmag = prethermostat*sqrt(a_sh);

          // force in each of the two directions

          randr = random->uniform()-0.5;
          fx += Fbmag*randr*p2[0];
          fy += Fbmag*randr*p2[1];
          fz += Fbmag*randr*p2[2];

          randr = random->uniform()-0.5;
          fx += Fbmag*randr*p3[0];
          fy += Fbmag*randr*p3[1];
          fz += Fbmag*randr*p3[2];
        }

        // scale forces to appropriate units

        fx = vxmu2f*fx;
        fy = vxmu2f*fy;
        fz = vxmu2f*fz;

        // sum to total Force

        f[i][0] -= fx;
        f[i][1] -= fy;
        f[i][2] -= fz;

        // torque due to the Brownian Force

        if (flaglog) {

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

          randr = random->uniform()-0.5;
          tx = Fbmag*randr*p2[0];
          ty = Fbmag*randr*p2[1];
          tz = Fbmag*randr*p2[2];

          randr = random->uniform()-0.5;
          tx += Fbmag*randr*p3[0];
          ty += Fbmag*randr*p3[1];
          tz += Fbmag*randr*p3[2];

          // torque has opposite sign on two particles

          torque[i][0] -= tx;
          torque[i][1] -= ty;
          torque[i][2] -= tz;

        }

        // set j = nlocal so that only I gets tallied

        if (evflag) ev_tally_xyz(i,nlocal,nlocal,0,
                                 0.0,0.0,-fx,-fy,-fz,delx,dely,delz);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBrownianPoly::init_style()
{
  if (force->newton_pair == 1)
    error->all(FLERR,"Pair brownian/poly requires newton pair off");
  if (!atom->sphere_flag)
    error->all(FLERR,"Pair brownian/poly requires atom style sphere");

  // insure all particles are finite-size
  // for pair hybrid, should limit test to types using the pair style

  double *radius = atom->radius;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (radius[i] == 0.0)
      error->one(FLERR,"Pair brownian/poly requires extended particles");

  neighbor->add_request(this, NeighConst::REQ_FULL);

  // set the isotropic constants that depend on the volume fraction
  // vol_T = total volume
  // check for fix deform, if exists it must use "remap v"
  // If box will change volume, set appropriate flag so that volume
  // and v.f. corrections are re-calculated at every step.
  //
  // If available volume is different from box volume
  // due to walls, set volume appropriately; if walls will
  // move, set appropriate flag so that volume and v.f. corrections
  // are re-calculated at every step.

  flagdeform = flagwall = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"deform") == 0)
      flagdeform = 1;
    else if (strstr(modify->fix[i]->style,"wall") != nullptr) {
      if (flagwall)
        error->all(FLERR,
                   "Cannot use multiple fix wall commands with pair brownian");
      flagwall = 1; // Walls exist
      wallfix = dynamic_cast<FixWall *>(modify->fix[i]);
      if (wallfix->xflag) flagwall = 2; // Moving walls exist
    }
  }

  // set the isotropic constants that depend on the volume fraction
  // vol_T = total volume

  double vol_T, wallcoord;
  if (!flagwall) vol_T = domain->xprd*domain->yprd*domain->zprd;
  else {
    double wallhi[3], walllo[3];
    for (int j = 0; j < 3; j++) {
      wallhi[j] = domain->prd[j];
      walllo[j] = 0;
    }
    for (int m = 0; m < wallfix->nwall; m++) {
      int dim = wallfix->wallwhich[m] / 2;
      int side = wallfix->wallwhich[m] % 2;
      if (wallfix->xstyle[m] == VARIABLE) {
        wallfix->xindex[m] = input->variable->find(wallfix->xstr[m]);
        // Since fix->wall->init happens after pair->init_style
        wallcoord = input->variable->compute_equal(wallfix->xindex[m]);
      }

      else wallcoord = wallfix->coord0[m];

      if (side == 0) walllo[dim] = wallcoord;
      else wallhi[dim] = wallcoord;
    }
    vol_T = (wallhi[0] - walllo[0]) * (wallhi[1] - walllo[1]) *
      (wallhi[2] - walllo[2]);
  }

  // vol_P = volume of particles, assuming mono-dispersity
  // vol_f = volume fraction

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
  } else {
    R0  = 6*MY_PI*mu*(1.0 + 2.725*vol_f - 6.583*vol_f*vol_f);
    RT0 = 8*MY_PI*mu*(1.0 + 0.749*vol_f - 2.469*vol_f*vol_f);
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBrownianPoly::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    cut_inner[i][j] = mix_distance(cut_inner[i][i],cut_inner[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  cut_inner[j][i] = cut_inner[i][j];
  return cut[i][j];
}
