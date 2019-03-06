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
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_oxdna_hbond.h"
#include "mf_oxdna.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "atom_vec_ellipsoid.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MFOxdna;

// sequence-specific base-pairing strength
// A:0 C:1 G:2 T:3, 5'- (i,j) -3'
static const double alpha[4][4] =
{{1.00000,1.00000,1.00000,0.82915},
 {1.00000,1.00000,1.15413,1.00000},
 {1.00000,1.15413,1.00000,1.00000},
 {0.82915,1.00000,1.00000,1.00000}};

/* ---------------------------------------------------------------------- */

PairOxdnaHbond::PairOxdnaHbond(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairOxdnaHbond::~PairOxdnaHbond()
{
  if (allocated) {

    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(epsilon_hb);
    memory->destroy(a_hb);
    memory->destroy(cut_hb_0);
    memory->destroy(cut_hb_c);
    memory->destroy(cut_hb_lo);
    memory->destroy(cut_hb_hi);
    memory->destroy(cut_hb_lc);
    memory->destroy(cut_hb_hc);
    memory->destroy(b_hb_lo);
    memory->destroy(b_hb_hi);
    memory->destroy(shift_hb);
    memory->destroy(cutsq_hb_hc);

    memory->destroy(a_hb1);
    memory->destroy(theta_hb1_0);
    memory->destroy(dtheta_hb1_ast);
    memory->destroy(b_hb1);
    memory->destroy(dtheta_hb1_c);

    memory->destroy(a_hb2);
    memory->destroy(theta_hb2_0);
    memory->destroy(dtheta_hb2_ast);
    memory->destroy(b_hb2);
    memory->destroy(dtheta_hb2_c);

    memory->destroy(a_hb3);
    memory->destroy(theta_hb3_0);
    memory->destroy(dtheta_hb3_ast);
    memory->destroy(b_hb3);
    memory->destroy(dtheta_hb3_c);

    memory->destroy(a_hb4);
    memory->destroy(theta_hb4_0);
    memory->destroy(dtheta_hb4_ast);
    memory->destroy(b_hb4);
    memory->destroy(dtheta_hb4_c);

    memory->destroy(a_hb7);
    memory->destroy(theta_hb7_0);
    memory->destroy(dtheta_hb7_ast);
    memory->destroy(b_hb7);
    memory->destroy(dtheta_hb7_c);

    memory->destroy(a_hb8);
    memory->destroy(theta_hb8_0);
    memory->destroy(dtheta_hb8_ast);
    memory->destroy(b_hb8);
    memory->destroy(dtheta_hb8_c);

  }
}

/* ----------------------------------------------------------------------
   compute function for oxDNA pair interactions
   hb=hydrogen bonding site
------------------------------------------------------------------------- */

void PairOxdnaHbond::compute(int eflag, int vflag)
{

  double delf[3],delta[3],deltb[3]; // force, torque increment;
  double evdwl,fpair,finc,tpair,factor_lj;
  double delr_hb[3],delr_hb_norm[3],rsq_hb,r_hb,rinv_hb;
  double theta1,t1dir[3],cost1;
  double theta2,t2dir[3],cost2;
  double theta3,t3dir[3],cost3;
  double theta4,t4dir[3],cost4;
  double theta7,t7dir[3],cost7;
  double theta8,t8dir[3],cost8;

  // distance COM-hbonding site
  double d_chb=+0.4;
  // vectors COM-h-bonding site in lab frame
  double ra_chb[3],rb_chb[3];

  // quaternions and Cartesian unit vectors in lab frame
  double *qa,ax[3],ay[3],az[3];
  double *qb,bx[3],by[3],bz[3];

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  int *type = atom->type;

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int *alist,*blist,*numneigh,**firstneigh;
  double *special_lj = force->special_lj;

  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;

  int a,b,ia,ib,anum,bnum,atype,btype;

  double f1,f4t1,f4t4,f4t2,f4t3,f4t7,f4t8;
  double df1,df4t1,df4t4,df4t2,df4t3,df4t7,df4t8;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  anum = list->inum;
  alist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over pair interaction neighbors of my atoms

  for (ia = 0; ia < anum; ia++) {

    a = alist[ia];
    atype = type[a];

    qa=bonus[a].quat;
    MathExtra::q_to_exyz(qa,ax,ay,az);

    ra_chb[0] = d_chb*ax[0];
    ra_chb[1] = d_chb*ax[1];
    ra_chb[2] = d_chb*ax[2];

    blist = firstneigh[a];
    bnum = numneigh[a];

    for (ib = 0; ib < bnum; ib++) {

      b = blist[ib];
      factor_lj = special_lj[sbmask(b)]; // = 0 for nearest neighbors
      b &= NEIGHMASK;

      btype = type[b];

      qb=bonus[b].quat;
      MathExtra::q_to_exyz(qb,bx,by,bz);

      rb_chb[0] = d_chb*bx[0];
      rb_chb[1] = d_chb*bx[1];
      rb_chb[2] = d_chb*bx[2];

      // vector h-bonding site b to a
      delr_hb[0] = x[a][0] + ra_chb[0] - x[b][0] - rb_chb[0];
      delr_hb[1] = x[a][1] + ra_chb[1] - x[b][1] - rb_chb[1];
      delr_hb[2] = x[a][2] + ra_chb[2] - x[b][2] - rb_chb[2];

      rsq_hb = delr_hb[0]*delr_hb[0] + delr_hb[1]*delr_hb[1] + delr_hb[2]*delr_hb[2];
      r_hb = sqrt(rsq_hb);
      rinv_hb = 1.0/r_hb;

      delr_hb_norm[0] = delr_hb[0] * rinv_hb;
      delr_hb_norm[1] = delr_hb[1] * rinv_hb;
      delr_hb_norm[2] = delr_hb[2] * rinv_hb;

      f1 = F1(r_hb, epsilon_hb[atype][btype], a_hb[atype][btype], cut_hb_0[atype][btype],
            cut_hb_lc[atype][btype], cut_hb_hc[atype][btype], cut_hb_lo[atype][btype], cut_hb_hi[atype][btype],
            b_hb_lo[atype][btype], b_hb_hi[atype][btype], shift_hb[atype][btype]);

      // early rejection criterium
      if (f1) {

      cost1 = -1.0*MathExtra::dot3(ax,bx);
      if (cost1 >  1.0) cost1 =  1.0;
      if (cost1 < -1.0) cost1 = -1.0;
      theta1 = acos(cost1);

      f4t1 = F4(theta1, a_hb1[atype][btype], theta_hb1_0[atype][btype], dtheta_hb1_ast[atype][btype],
            b_hb1[atype][btype], dtheta_hb1_c[atype][btype]);

      // early rejection criterium
      if (f4t1) {

      cost2 = -1.0*MathExtra::dot3(ax,delr_hb_norm);
      if (cost2 >  1.0) cost2 =  1.0;
      if (cost2 < -1.0) cost2 = -1.0;
      theta2 = acos(cost2);

      f4t2 = F4(theta2, a_hb2[atype][btype], theta_hb2_0[atype][btype], dtheta_hb2_ast[atype][btype],
            b_hb2[atype][btype], dtheta_hb2_c[atype][btype]);

      // early rejection criterium
      if (f4t2) {

      cost3 = MathExtra::dot3(bx,delr_hb_norm);
      if (cost3 >  1.0) cost3 =  1.0;
      if (cost3 < -1.0) cost3 = -1.0;
      theta3 = acos(cost3);

      f4t3 = F4(theta3, a_hb3[atype][btype], theta_hb3_0[atype][btype], dtheta_hb3_ast[atype][btype],
            b_hb3[atype][btype], dtheta_hb3_c[atype][btype]);

      // early rejection criterium
      if (f4t3) {

      cost4 = MathExtra::dot3(az,bz);
      if (cost4 >  1.0) cost4 =  1.0;
      if (cost4 < -1.0) cost4 = -1.0;
      theta4 = acos(cost4);

      f4t4 = F4(theta4, a_hb4[atype][btype], theta_hb4_0[atype][btype], dtheta_hb4_ast[atype][btype],
            b_hb4[atype][btype], dtheta_hb4_c[atype][btype]);

      // early rejection criterium
      if (f4t4) {

      cost7 = -1.0*MathExtra::dot3(az,delr_hb_norm);
      if (cost7 >  1.0) cost7 =  1.0;
      if (cost7 < -1.0) cost7 = -1.0;
      theta7 = acos(cost7);

      f4t7 = F4(theta7, a_hb7[atype][btype], theta_hb7_0[atype][btype], dtheta_hb7_ast[atype][btype],
            b_hb7[atype][btype], dtheta_hb7_c[atype][btype]);

      // early rejection criterium
      if (f4t7) {

      cost8 = MathExtra::dot3(bz,delr_hb_norm);
      if (cost8 >  1.0) cost8 =  1.0;
      if (cost8 < -1.0) cost8 = -1.0;
      theta8 = acos(cost8);

      f4t8 = F4(theta8, a_hb8[atype][btype], theta_hb8_0[atype][btype], dtheta_hb8_ast[atype][btype],
            b_hb8[atype][btype], dtheta_hb8_c[atype][btype]);

      evdwl = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj;

      // early rejection criterium
      if (evdwl) {

      df1 = DF1(r_hb, epsilon_hb[atype][btype], a_hb[atype][btype], cut_hb_0[atype][btype],
            cut_hb_lc[atype][btype], cut_hb_hc[atype][btype], cut_hb_lo[atype][btype], cut_hb_hi[atype][btype],
            b_hb_lo[atype][btype], b_hb_hi[atype][btype]);

      df4t1 = DF4(theta1, a_hb1[atype][btype], theta_hb1_0[atype][btype], dtheta_hb1_ast[atype][btype],
            b_hb1[atype][btype], dtheta_hb1_c[atype][btype])/sin(theta1);

      df4t2 = DF4(theta2, a_hb2[atype][btype], theta_hb2_0[atype][btype], dtheta_hb2_ast[atype][btype],
            b_hb2[atype][btype], dtheta_hb2_c[atype][btype])/sin(theta2);

      df4t3 = DF4(theta3, a_hb3[atype][btype], theta_hb3_0[atype][btype], dtheta_hb3_ast[atype][btype],
            b_hb3[atype][btype], dtheta_hb3_c[atype][btype])/sin(theta3);

      df4t4 = DF4(theta4, a_hb4[atype][btype], theta_hb4_0[atype][btype], dtheta_hb4_ast[atype][btype],
            b_hb4[atype][btype], dtheta_hb4_c[atype][btype])/sin(theta4);

      df4t7 = DF4(theta7, a_hb7[atype][btype], theta_hb7_0[atype][btype], dtheta_hb7_ast[atype][btype],
            b_hb7[atype][btype], dtheta_hb7_c[atype][btype])/sin(theta7);

      df4t8 = DF4(theta8, a_hb8[atype][btype], theta_hb8_0[atype][btype], dtheta_hb8_ast[atype][btype],
            b_hb8[atype][btype], dtheta_hb8_c[atype][btype])/sin(theta8);

      // force, torque and virial contribution for forces between h-bonding sites

      fpair = 0.0;

      delf[0] = 0.0;
      delf[1] = 0.0;
      delf[2] = 0.0;

      delta[0] = 0.0;
      delta[1] = 0.0;
      delta[2] = 0.0;

      deltb[0] = 0.0;
      deltb[1] = 0.0;
      deltb[2] = 0.0;

      // radial force
      finc  = -df1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj;
      fpair += finc;

      delf[0] += delr_hb[0] * finc;
      delf[1] += delr_hb[1] * finc;
      delf[2] += delr_hb[2] * finc;

      // theta2 force
      if (theta2) {

        finc  = -f1 * f4t1 * df4t2 * f4t3 * f4t4 * f4t7 * f4t8 * rinv_hb * factor_lj;
        fpair += finc;

        delf[0] += (delr_hb_norm[0]*cost2 + ax[0]) * finc;
        delf[1] += (delr_hb_norm[1]*cost2 + ax[1]) * finc;
        delf[2] += (delr_hb_norm[2]*cost2 + ax[2]) * finc;

      }

      // theta3 force
      if (theta3) {

        finc  = -f1 * f4t1 * f4t2 * df4t3 * f4t4 * f4t7 * f4t8 * rinv_hb * factor_lj;
        fpair += finc;

        delf[0] += (delr_hb_norm[0]*cost3 - bx[0]) * finc;
        delf[1] += (delr_hb_norm[1]*cost3 - bx[1]) * finc;
        delf[2] += (delr_hb_norm[2]*cost3 - bx[2]) * finc;

      }

      // theta7 force
      if (theta7) {

        finc  = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * df4t7 * f4t8 * rinv_hb * factor_lj;
        fpair += finc;

        delf[0] += (delr_hb_norm[0]*cost7 + az[0]) * finc;
        delf[1] += (delr_hb_norm[1]*cost7 + az[1]) * finc;
        delf[2] += (delr_hb_norm[2]*cost7 + az[2]) * finc;

      }

      // theta8 force
      if (theta8) {

        finc  = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * df4t8 * rinv_hb * factor_lj;
        fpair += finc;

        delf[0] += (delr_hb_norm[0]*cost8 - bz[0]) * finc;
        delf[1] += (delr_hb_norm[1]*cost8 - bz[1]) * finc;
        delf[2] += (delr_hb_norm[2]*cost8 - bz[2]) * finc;

      }

      // increment forces and torques

      f[a][0] += delf[0];
      f[a][1] += delf[1];
      f[a][2] += delf[2];

      MathExtra::cross3(ra_chb,delf,delta);

      torque[a][0] += delta[0];
      torque[a][1] += delta[1];
      torque[a][2] += delta[2];

      if (newton_pair || b < nlocal) {

        f[b][0] -= delf[0];
        f[b][1] -= delf[1];
        f[b][2] -= delf[2];


        MathExtra::cross3(rb_chb,delf,deltb);

        torque[b][0] -= deltb[0];
        torque[b][1] -= deltb[1];
        torque[b][2] -= deltb[2];

      }

      // increment energy and virial
      if (evflag) ev_tally(a,b,nlocal,newton_pair,evdwl,0.0,fpair,delr_hb[0],delr_hb[1],delr_hb[2]);

      // pure torques not expressible as r x f

      delta[0] = 0.0;
      delta[1] = 0.0;
      delta[2] = 0.0;
      deltb[0] = 0.0;
      deltb[1] = 0.0;
      deltb[2] = 0.0;

      // theta1 torque
      if (theta1) {

        tpair = -f1 * df4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj;
        MathExtra::cross3(ax,bx,t1dir);

        delta[0] += t1dir[0]*tpair;
        delta[1] += t1dir[1]*tpair;
        delta[2] += t1dir[2]*tpair;

        deltb[0] += t1dir[0]*tpair;
        deltb[1] += t1dir[1]*tpair;
        deltb[2] += t1dir[2]*tpair;

      }

      // theta2 torque
      if (theta2) {

        tpair = -f1 * f4t1 * df4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj;
        MathExtra::cross3(ax,delr_hb_norm,t2dir);

        delta[0] += t2dir[0]*tpair;
        delta[1] += t2dir[1]*tpair;
        delta[2] += t2dir[2]*tpair;

      }

      // theta3 torque
      if (theta3) {

        tpair = -f1 * f4t1 * f4t2 * df4t3 * f4t4 * f4t7 * f4t8 * factor_lj;
        MathExtra::cross3(bx,delr_hb_norm,t3dir);

        deltb[0] += t3dir[0]*tpair;
        deltb[1] += t3dir[1]*tpair;
        deltb[2] += t3dir[2]*tpair;

      }

      // theta4 torque
      if (theta4) {

        tpair = -f1 * f4t1 * f4t2 * f4t3 * df4t4 * f4t7 * f4t8 * factor_lj;
        MathExtra::cross3(bz,az,t4dir);

        delta[0] += t4dir[0]*tpair;
        delta[1] += t4dir[1]*tpair;
        delta[2] += t4dir[2]*tpair;

        deltb[0] += t4dir[0]*tpair;
        deltb[1] += t4dir[1]*tpair;
        deltb[2] += t4dir[2]*tpair;

      }

      // theta7 torque
      if (theta7) {

        tpair = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * df4t7 * f4t8 * factor_lj;
        MathExtra::cross3(az,delr_hb_norm,t7dir);

        delta[0] += t7dir[0]*tpair;
        delta[1] += t7dir[1]*tpair;
        delta[2] += t7dir[2]*tpair;

      }

      // theta8 torque
      if (theta8) {

        tpair = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * df4t8 * factor_lj;
        MathExtra::cross3(bz,delr_hb_norm,t8dir);

        deltb[0] += t8dir[0]*tpair;
        deltb[1] += t8dir[1]*tpair;
        deltb[2] += t8dir[2]*tpair;

      }

      // increment torques

      torque[a][0] += delta[0];
      torque[a][1] += delta[1];
      torque[a][2] += delta[2];

      if (newton_pair || b < nlocal) {

        torque[b][0] -= deltb[0];
        torque[b][1] -= deltb[1];
        torque[b][2] -= deltb[2];

      }

      }
      }
      }
      }
      }
      }
      }// end early rejection criteria

    }

  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairOxdnaHbond::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(epsilon_hb,n+1,n+1,"pair:epsilon_hb");
  memory->create(a_hb,n+1,n+1,"pair:a_hb");
  memory->create(cut_hb_0,n+1,n+1,"pair:cut_hb_0");
  memory->create(cut_hb_c,n+1,n+1,"pair:cut_hb_c");
  memory->create(cut_hb_lo,n+1,n+1,"pair:cut_hb_lo");
  memory->create(cut_hb_hi,n+1,n+1,"pair:cut_hb_hi");
  memory->create(cut_hb_lc,n+1,n+1,"pair:cut_hb_lc");
  memory->create(cut_hb_hc,n+1,n+1,"pair:cut_hb_hc");
  memory->create(b_hb_lo,n+1,n+1,"pair:b_hb_lo");
  memory->create(b_hb_hi,n+1,n+1,"pair:b_hb_hi");
  memory->create(shift_hb,n+1,n+1,"pair:shift_hb");
  memory->create(cutsq_hb_hc,n+1,n+1,"pair:cutsq_hb_hc");

  memory->create(a_hb1,n+1,n+1,"pair:a_hb1");
  memory->create(theta_hb1_0,n+1,n+1,"pair:theta_hb1_0");
  memory->create(dtheta_hb1_ast,n+1,n+1,"pair:dtheta_hb1_ast");
  memory->create(b_hb1,n+1,n+1,"pair:b_hb1");
  memory->create(dtheta_hb1_c,n+1,n+1,"pair:dtheta_hb1_c");

  memory->create(a_hb2,n+1,n+1,"pair:a_hb2");
  memory->create(theta_hb2_0,n+1,n+1,"pair:theta_hb2_0");
  memory->create(dtheta_hb2_ast,n+1,n+1,"pair:dtheta_hb2_ast");
  memory->create(b_hb2,n+1,n+1,"pair:b_hb2");
  memory->create(dtheta_hb2_c,n+1,n+1,"pair:dtheta_hb2_c");

  memory->create(a_hb3,n+1,n+1,"pair:a_hb3");
  memory->create(theta_hb3_0,n+1,n+1,"pair:theta_hb3_0");
  memory->create(dtheta_hb3_ast,n+1,n+1,"pair:dtheta_hb3_ast");
  memory->create(b_hb3,n+1,n+1,"pair:b_hb3");
  memory->create(dtheta_hb3_c,n+1,n+1,"pair:dtheta_hb3_c");

  memory->create(a_hb4,n+1,n+1,"pair:a_hb4");
  memory->create(theta_hb4_0,n+1,n+1,"pair:theta_hb4_0");
  memory->create(dtheta_hb4_ast,n+1,n+1,"pair:dtheta_hb4_ast");
  memory->create(b_hb4,n+1,n+1,"pair:b_hb4");
  memory->create(dtheta_hb4_c,n+1,n+1,"pair:dtheta_hb4_c");

  memory->create(a_hb7,n+1,n+1,"pair:a_hb7");
  memory->create(theta_hb7_0,n+1,n+1,"pair:theta_hb7_0");
  memory->create(dtheta_hb7_ast,n+1,n+1,"pair:dtheta_hb7_ast");
  memory->create(b_hb7,n+1,n+1,"pair:b_hb7");
  memory->create(dtheta_hb7_c,n+1,n+1,"pair:dtheta_hb7_c");

  memory->create(a_hb8,n+1,n+1,"pair:a_hb8");
  memory->create(theta_hb8_0,n+1,n+1,"pair:theta_hb8_0");
  memory->create(dtheta_hb8_ast,n+1,n+1,"pair:dtheta_hb8_ast");
  memory->create(b_hb8,n+1,n+1,"pair:b_hb8");
  memory->create(dtheta_hb8_c,n+1,n+1,"pair:dtheta_hb8_c");

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairOxdnaHbond::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairOxdnaHbond::coeff(int narg, char **arg)
{
  int count;

  if (narg != 27) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/hbond");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  // h-bonding interaction
  count = 0;

  double epsilon_hb_one, a_hb_one, cut_hb_0_one, cut_hb_c_one, cut_hb_lo_one, cut_hb_hi_one;
  double b_hb_lo_one, b_hb_hi_one, cut_hb_lc_one, cut_hb_hc_one, tmp, shift_hb_one;

  double a_hb1_one, theta_hb1_0_one, dtheta_hb1_ast_one;
  double b_hb1_one, dtheta_hb1_c_one;

  double a_hb2_one, theta_hb2_0_one, dtheta_hb2_ast_one;
  double b_hb2_one, dtheta_hb2_c_one;

  double a_hb3_one, theta_hb3_0_one, dtheta_hb3_ast_one;
  double b_hb3_one, dtheta_hb3_c_one;

  double a_hb4_one, theta_hb4_0_one, dtheta_hb4_ast_one;
  double b_hb4_one, dtheta_hb4_c_one;

  double a_hb7_one, theta_hb7_0_one, dtheta_hb7_ast_one;
  double b_hb7_one, dtheta_hb7_c_one;

  double a_hb8_one, theta_hb8_0_one, dtheta_hb8_ast_one;
  double b_hb8_one, dtheta_hb8_c_one;

  if (strcmp(arg[2], "seqav") != 0 && strcmp(arg[2], "seqdep") != 0) {
    error->all(FLERR,"Incorrect setting, select seqav or seqdep in oxdna/hbond");
  }
  if (strcmp(arg[2],"seqav")  == 0) seqdepflag = 0;
  if (strcmp(arg[2],"seqdep") == 0) seqdepflag = 1;

  epsilon_hb_one = force->numeric(FLERR,arg[3]);
  a_hb_one = force->numeric(FLERR,arg[4]);
  cut_hb_0_one = force->numeric(FLERR,arg[5]);
  cut_hb_c_one = force->numeric(FLERR,arg[6]);
  cut_hb_lo_one = force->numeric(FLERR,arg[7]);
  cut_hb_hi_one = force->numeric(FLERR,arg[8]);

  a_hb1_one = force->numeric(FLERR,arg[9]);
  theta_hb1_0_one = force->numeric(FLERR,arg[10]);
  dtheta_hb1_ast_one = force->numeric(FLERR,arg[11]);

  a_hb2_one = force->numeric(FLERR,arg[12]);
  theta_hb2_0_one = force->numeric(FLERR,arg[13]);
  dtheta_hb2_ast_one = force->numeric(FLERR,arg[14]);

  a_hb3_one = force->numeric(FLERR,arg[15]);
  theta_hb3_0_one = force->numeric(FLERR,arg[16]);
  dtheta_hb3_ast_one = force->numeric(FLERR,arg[17]);

  a_hb4_one = force->numeric(FLERR,arg[18]);
  theta_hb4_0_one = force->numeric(FLERR,arg[19]);
  dtheta_hb4_ast_one = force->numeric(FLERR,arg[20]);

  a_hb7_one = force->numeric(FLERR,arg[21]);
  theta_hb7_0_one = force->numeric(FLERR,arg[22]);
  dtheta_hb7_ast_one = force->numeric(FLERR,arg[23]);

  a_hb8_one = force->numeric(FLERR,arg[24]);
  theta_hb8_0_one = force->numeric(FLERR,arg[25]);
  dtheta_hb8_ast_one = force->numeric(FLERR,arg[26]);

  b_hb_lo_one = 2*a_hb_one*exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one))*
        2*a_hb_one*exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one))*
        (1-exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one)))*
        (1-exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one)))/
        (4*((1-exp(-a_hb_one*(cut_hb_lo_one -cut_hb_0_one)))*
        (1-exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one)))-
        (1-exp(-a_hb_one*(cut_hb_c_one -cut_hb_0_one)))*
        (1-exp(-a_hb_one*(cut_hb_c_one-cut_hb_0_one)))));

  cut_hb_lc_one = cut_hb_lo_one - a_hb_one*exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one))*
        (1-exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one)))/b_hb_lo_one;

  b_hb_hi_one = 2*a_hb_one*exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one))*
        2*a_hb_one*exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one))*
        (1-exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one)))*
        (1-exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one)))/
        (4*((1-exp(-a_hb_one*(cut_hb_hi_one -cut_hb_0_one)))*
        (1-exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one)))-
        (1-exp(-a_hb_one*(cut_hb_c_one -cut_hb_0_one)))*
        (1-exp(-a_hb_one*(cut_hb_c_one-cut_hb_0_one)))));

  cut_hb_hc_one = cut_hb_hi_one - a_hb_one*exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one))*
        (1-exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one)))/b_hb_hi_one;

  tmp = 1 - exp(-(cut_hb_c_one-cut_hb_0_one) * a_hb_one);
  shift_hb_one = epsilon_hb_one * tmp * tmp;

  b_hb1_one = a_hb1_one*a_hb1_one*dtheta_hb1_ast_one*dtheta_hb1_ast_one/(1-a_hb1_one*dtheta_hb1_ast_one*dtheta_hb1_ast_one);
  dtheta_hb1_c_one = 1/(a_hb1_one*dtheta_hb1_ast_one);

  b_hb2_one = a_hb2_one*a_hb2_one*dtheta_hb2_ast_one*dtheta_hb2_ast_one/(1-a_hb2_one*dtheta_hb2_ast_one*dtheta_hb2_ast_one);
  dtheta_hb2_c_one = 1/(a_hb2_one*dtheta_hb2_ast_one);

  b_hb3_one = a_hb3_one*a_hb3_one*dtheta_hb3_ast_one*dtheta_hb3_ast_one/(1-a_hb3_one*dtheta_hb3_ast_one*dtheta_hb3_ast_one);
  dtheta_hb3_c_one = 1/(a_hb3_one*dtheta_hb3_ast_one);

  b_hb4_one = a_hb4_one*a_hb4_one*dtheta_hb4_ast_one*dtheta_hb4_ast_one/(1-a_hb4_one*dtheta_hb4_ast_one*dtheta_hb4_ast_one);
  dtheta_hb4_c_one = 1/(a_hb4_one*dtheta_hb4_ast_one);

  b_hb7_one = a_hb7_one*a_hb7_one*dtheta_hb7_ast_one*dtheta_hb7_ast_one/(1-a_hb7_one*dtheta_hb7_ast_one*dtheta_hb7_ast_one);
  dtheta_hb7_c_one = 1/(a_hb7_one*dtheta_hb7_ast_one);

  b_hb8_one = a_hb8_one*a_hb8_one*dtheta_hb8_ast_one*dtheta_hb8_ast_one/(1-a_hb8_one*dtheta_hb8_ast_one*dtheta_hb8_ast_one);
  dtheta_hb8_c_one = 1/(a_hb8_one*dtheta_hb8_ast_one);

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {

      epsilon_hb[i][j] = epsilon_hb_one;
      if (seqdepflag) epsilon_hb[i][j] *= alpha[i-1][j-1];
      a_hb[i][j] = a_hb_one;
      cut_hb_0[i][j] = cut_hb_0_one;
      cut_hb_c[i][j] = cut_hb_c_one;
      cut_hb_lo[i][j] = cut_hb_lo_one;
      cut_hb_hi[i][j] = cut_hb_hi_one;
      cut_hb_lc[i][j] = cut_hb_lc_one;
      cut_hb_hc[i][j] = cut_hb_hc_one;
      b_hb_lo[i][j] = b_hb_lo_one;
      b_hb_hi[i][j] = b_hb_hi_one;
      shift_hb[i][j] = shift_hb_one;
      if (seqdepflag) shift_hb[i][j] *= alpha[i-1][j-1];

      a_hb1[i][j] = a_hb1_one;
      theta_hb1_0[i][j] = theta_hb1_0_one;
      dtheta_hb1_ast[i][j] = dtheta_hb1_ast_one;
      b_hb1[i][j] = b_hb1_one;
      dtheta_hb1_c[i][j] = dtheta_hb1_c_one;

      a_hb2[i][j] = a_hb2_one;
      theta_hb2_0[i][j] = theta_hb2_0_one;
      dtheta_hb2_ast[i][j] = dtheta_hb2_ast_one;
      b_hb2[i][j] = b_hb2_one;
      dtheta_hb2_c[i][j] = dtheta_hb2_c_one;

      a_hb3[i][j] = a_hb3_one;
      theta_hb3_0[i][j] = theta_hb3_0_one;
      dtheta_hb3_ast[i][j] = dtheta_hb3_ast_one;
      b_hb3[i][j] = b_hb3_one;
      dtheta_hb3_c[i][j] = dtheta_hb3_c_one;

      a_hb4[i][j] = a_hb4_one;
      theta_hb4_0[i][j] = theta_hb4_0_one;
      dtheta_hb4_ast[i][j] = dtheta_hb4_ast_one;
      b_hb4[i][j] = b_hb4_one;
      dtheta_hb4_c[i][j] = dtheta_hb4_c_one;

      a_hb7[i][j] = a_hb7_one;
      theta_hb7_0[i][j] = theta_hb7_0_one;
      dtheta_hb7_ast[i][j] = dtheta_hb7_ast_one;
      b_hb7[i][j] = b_hb7_one;
      dtheta_hb7_c[i][j] = dtheta_hb7_c_one;

      a_hb8[i][j] = a_hb8_one;
      theta_hb8_0[i][j] = theta_hb8_0_one;
      dtheta_hb8_ast[i][j] = dtheta_hb8_ast_one;
      b_hb8[i][j] = b_hb8_one;
      dtheta_hb8_c[i][j] = dtheta_hb8_c_one;

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/hbond");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairOxdnaHbond::init_style()
{
  int irequest;

  // request regular neighbor lists

  irequest = neighbor->request(this,instance_me);

}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use regular
------------------------------------------------------------------------- */

void PairOxdnaHbond::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  if (id  > 0) error->all(FLERR,"Respa not supported");

}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairOxdnaHbond::init_one(int i, int j)
{

  if (setflag[i][j] == 0) {
    error->all(FLERR,"Coefficient mixing not defined in oxDNA");
  }
  if (offset_flag) {
    error->all(FLERR,"Offset not supported in oxDNA");
  }

  if (seqdepflag) {
    epsilon_hb[j][i] = epsilon_hb[i][j] / alpha[i-1][j-1] * alpha[j-1][i-1];
  }
  else {
    epsilon_hb[j][i] = epsilon_hb[i][j];
  }
  a_hb[j][i] = a_hb[i][j];
  cut_hb_0[j][i] = cut_hb_0[i][j];
  cut_hb_c[j][i] = cut_hb_c[i][j];
  cut_hb_lo[j][i] = cut_hb_lo[i][j];
  cut_hb_hi[j][i] = cut_hb_hi[i][j];
  b_hb_lo[j][i] = b_hb_lo[i][j];
  b_hb_hi[j][i] = b_hb_hi[i][j];
  cut_hb_lc[j][i] = cut_hb_lc[i][j];
  cut_hb_hc[j][i] = cut_hb_hc[i][j];
  if (seqdepflag) {
    shift_hb[j][i] = shift_hb[i][j] / alpha[i-1][j-1] * alpha[j-1][i-1];
  }
  else {
    shift_hb[j][i] = shift_hb[i][j];
  }

  a_hb1[j][i] = a_hb1[i][j];
  theta_hb1_0[j][i] = theta_hb1_0[i][j];
  dtheta_hb1_ast[j][i] = dtheta_hb1_ast[i][j];
  b_hb1[j][i] = b_hb1[i][j];
  dtheta_hb1_c[j][i] = dtheta_hb1_c[i][j];

  a_hb2[j][i] = a_hb2[i][j];
  theta_hb2_0[j][i] = theta_hb2_0[i][j];
  dtheta_hb2_ast[j][i] = dtheta_hb2_ast[i][j];
  b_hb2[j][i] = b_hb2[i][j];
  dtheta_hb2_c[j][i] = dtheta_hb2_c[i][j];

  a_hb3[j][i] = a_hb3[i][j];
  theta_hb3_0[j][i] = theta_hb3_0[i][j];
  dtheta_hb3_ast[j][i] = dtheta_hb3_ast[i][j];
  b_hb3[j][i] = b_hb3[i][j];
  dtheta_hb3_c[j][i] = dtheta_hb3_c[i][j];

  a_hb4[j][i] = a_hb4[i][j];
  theta_hb4_0[j][i] = theta_hb4_0[i][j];
  dtheta_hb4_ast[j][i] = dtheta_hb4_ast[i][j];
  b_hb4[j][i] = b_hb4[i][j];
  dtheta_hb4_c[j][i] = dtheta_hb4_c[i][j];

  a_hb7[j][i] = a_hb7[i][j];
  theta_hb7_0[j][i] = theta_hb7_0[i][j];
  dtheta_hb7_ast[j][i] = dtheta_hb7_ast[i][j];
  b_hb7[j][i] = b_hb7[i][j];
  dtheta_hb7_c[j][i] = dtheta_hb7_c[i][j];

  a_hb8[j][i] = a_hb8[i][j];
  theta_hb8_0[j][i] = theta_hb8_0[i][j];
  dtheta_hb8_ast[j][i] = dtheta_hb8_ast[i][j];
  b_hb8[j][i] = b_hb8[i][j];
  dtheta_hb8_c[j][i] = dtheta_hb8_c[i][j];

  cutsq_hb_hc[i][j] = cut_hb_hc[i][j]*cut_hb_hc[i][j];
  cutsq_hb_hc[j][i] = cutsq_hb_hc[i][j];

  // set the master list distance cutoff
  return cut_hb_hc[i][j];

}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdnaHbond::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {

        fwrite(&epsilon_hb[i][j],sizeof(double),1,fp);
        fwrite(&a_hb[i][j],sizeof(double),1,fp);
        fwrite(&cut_hb_0[i][j],sizeof(double),1,fp);
        fwrite(&cut_hb_c[i][j],sizeof(double),1,fp);
        fwrite(&cut_hb_lo[i][j],sizeof(double),1,fp);
        fwrite(&cut_hb_hi[i][j],sizeof(double),1,fp);
        fwrite(&cut_hb_lc[i][j],sizeof(double),1,fp);
        fwrite(&cut_hb_hc[i][j],sizeof(double),1,fp);
        fwrite(&b_hb_lo[i][j],sizeof(double),1,fp);
        fwrite(&b_hb_hi[i][j],sizeof(double),1,fp);
        fwrite(&shift_hb[i][j],sizeof(double),1,fp);

        fwrite(&a_hb1[i][j],sizeof(double),1,fp);
        fwrite(&theta_hb1_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_hb1_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_hb1[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_hb1_c[i][j],sizeof(double),1,fp);

        fwrite(&a_hb2[i][j],sizeof(double),1,fp);
        fwrite(&theta_hb2_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_hb2_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_hb2[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_hb2_c[i][j],sizeof(double),1,fp);

        fwrite(&a_hb3[i][j],sizeof(double),1,fp);
        fwrite(&theta_hb3_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_hb3_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_hb3[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_hb3_c[i][j],sizeof(double),1,fp);

        fwrite(&a_hb4[i][j],sizeof(double),1,fp);
        fwrite(&theta_hb4_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_hb4_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_hb4[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_hb4_c[i][j],sizeof(double),1,fp);

        fwrite(&a_hb7[i][j],sizeof(double),1,fp);
        fwrite(&theta_hb7_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_hb7_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_hb7[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_hb7_c[i][j],sizeof(double),1,fp);

        fwrite(&a_hb8[i][j],sizeof(double),1,fp);
        fwrite(&theta_hb8_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_hb8_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_hb8[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_hb8_c[i][j],sizeof(double),1,fp);

    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdnaHbond::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {

          fread(&epsilon_hb[i][j],sizeof(double),1,fp);
          fread(&a_hb[i][j],sizeof(double),1,fp);
          fread(&cut_hb_0[i][j],sizeof(double),1,fp);
          fread(&cut_hb_c[i][j],sizeof(double),1,fp);
          fread(&cut_hb_lo[i][j],sizeof(double),1,fp);
          fread(&cut_hb_hi[i][j],sizeof(double),1,fp);
          fread(&cut_hb_lc[i][j],sizeof(double),1,fp);
          fread(&cut_hb_hc[i][j],sizeof(double),1,fp);
          fread(&b_hb_lo[i][j],sizeof(double),1,fp);
          fread(&b_hb_hi[i][j],sizeof(double),1,fp);
          fread(&shift_hb[i][j],sizeof(double),1,fp);

          fread(&a_hb1[i][j],sizeof(double),1,fp);
          fread(&theta_hb1_0[i][j],sizeof(double),1,fp);
          fread(&dtheta_hb1_ast[i][j],sizeof(double),1,fp);
          fread(&b_hb1[i][j],sizeof(double),1,fp);
          fread(&dtheta_hb1_c[i][j],sizeof(double),1,fp);

          fread(&a_hb2[i][j],sizeof(double),1,fp);
          fread(&theta_hb2_0[i][j],sizeof(double),1,fp);
          fread(&dtheta_hb2_ast[i][j],sizeof(double),1,fp);
          fread(&b_hb2[i][j],sizeof(double),1,fp);
          fread(&dtheta_hb2_c[i][j],sizeof(double),1,fp);

          fread(&a_hb3[i][j],sizeof(double),1,fp);
          fread(&theta_hb3_0[i][j],sizeof(double),1,fp);
          fread(&dtheta_hb3_ast[i][j],sizeof(double),1,fp);
          fread(&b_hb3[i][j],sizeof(double),1,fp);
          fread(&dtheta_hb3_c[i][j],sizeof(double),1,fp);

          fread(&a_hb4[i][j],sizeof(double),1,fp);
          fread(&theta_hb4_0[i][j],sizeof(double),1,fp);
          fread(&dtheta_hb4_ast[i][j],sizeof(double),1,fp);
          fread(&b_hb4[i][j],sizeof(double),1,fp);
          fread(&dtheta_hb4_c[i][j],sizeof(double),1,fp);

          fread(&a_hb7[i][j],sizeof(double),1,fp);
          fread(&theta_hb7_0[i][j],sizeof(double),1,fp);
          fread(&dtheta_hb7_ast[i][j],sizeof(double),1,fp);
          fread(&b_hb7[i][j],sizeof(double),1,fp);
          fread(&dtheta_hb7_c[i][j],sizeof(double),1,fp);

          fread(&a_hb8[i][j],sizeof(double),1,fp);
          fread(&theta_hb8_0[i][j],sizeof(double),1,fp);
          fread(&dtheta_hb8_ast[i][j],sizeof(double),1,fp);
          fread(&b_hb8[i][j],sizeof(double),1,fp);
          fread(&dtheta_hb8_c[i][j],sizeof(double),1,fp);

        }

        MPI_Bcast(&epsilon_hb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a_hb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_hb_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_hb_c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_hb_lo[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_hb_hi[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_hb_lc[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_hb_hc[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_hb_lo[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_hb_hi[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&shift_hb[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_hb1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_hb1_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_hb1_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_hb1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_hb1_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_hb2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_hb2_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_hb2_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_hb2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_hb2_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_hb3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_hb3_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_hb3_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_hb3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_hb3_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_hb4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_hb4_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_hb4_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_hb4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_hb4_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_hb7[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_hb7_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_hb7_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_hb7[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_hb7_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_hb8[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_hb8_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_hb8_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_hb8[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_hb8_c[i][j],1,MPI_DOUBLE,0,world);

      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdnaHbond::write_restart_settings(FILE *fp)
{
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdnaHbond::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairOxdnaHbond::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d\
         %g %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         \n",i,
        epsilon_hb[i][i],a_hb[i][i],cut_hb_0[i][i],cut_hb_c[i][i],cut_hb_lo[i][i],cut_hb_hi[i][i],
        cut_hb_lc[i][i],cut_hb_hc[i][i],b_hb_lo[i][i],b_hb_hi[i][i],shift_hb[i][i],
        a_hb1[i][i],theta_hb1_0[i][i],dtheta_hb1_ast[i][i],b_hb1[i][i],dtheta_hb1_c[i][i],
        a_hb2[i][i],theta_hb2_0[i][i],dtheta_hb2_ast[i][i],b_hb2[i][i],dtheta_hb2_c[i][i],
        a_hb3[i][i],theta_hb3_0[i][i],dtheta_hb3_ast[i][i],b_hb3[i][i],dtheta_hb3_c[i][i],
        a_hb4[i][i],theta_hb4_0[i][i],dtheta_hb4_ast[i][i],b_hb4[i][i],dtheta_hb4_c[i][i],
        a_hb7[i][i],theta_hb7_0[i][i],dtheta_hb7_ast[i][i],b_hb7[i][i],dtheta_hb7_c[i][i],
        a_hb8[i][i],theta_hb8_0[i][i],dtheta_hb8_ast[i][i],b_hb8[i][i],dtheta_hb8_c[i][i]);

}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairOxdnaHbond::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d\
         %g %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         \n",i,j,
        epsilon_hb[i][j],a_hb[i][j],cut_hb_0[i][j],cut_hb_c[i][j],cut_hb_lo[i][j],cut_hb_hi[i][j],
        cut_hb_lc[i][j],cut_hb_hc[i][j],b_hb_lo[i][j],b_hb_hi[i][j],shift_hb[i][j],
        a_hb1[i][j],theta_hb1_0[i][j],dtheta_hb1_ast[i][j],b_hb1[i][j],dtheta_hb1_c[i][j],
        a_hb2[i][j],theta_hb2_0[i][j],dtheta_hb2_ast[i][j],b_hb2[i][j],dtheta_hb2_c[i][j],
        a_hb3[i][j],theta_hb3_0[i][j],dtheta_hb3_ast[i][j],b_hb3[i][j],dtheta_hb3_c[i][j],
        a_hb4[i][j],theta_hb4_0[i][j],dtheta_hb4_ast[i][j],b_hb4[i][j],dtheta_hb4_c[i][j],
        a_hb7[i][j],theta_hb7_0[i][j],dtheta_hb7_ast[i][j],b_hb7[i][j],dtheta_hb7_c[i][j],
        a_hb8[i][j],theta_hb8_0[i][j],dtheta_hb8_ast[i][j],b_hb8[i][j],dtheta_hb8_c[i][j]);

}

/* ---------------------------------------------------------------------- */

void *PairOxdnaHbond::extract(const char *str, int &dim)
{
  dim = 2;

  if (strcmp(str,"epsilon_hb") == 0) return (void *) epsilon_hb;
  if (strcmp(str,"a_hb") == 0) return (void *) a_hb;
  if (strcmp(str,"cut_hb_0") == 0) return (void *) cut_hb_0;
  if (strcmp(str,"cut_hb_c") == 0) return (void *) cut_hb_c;
  if (strcmp(str,"cut_hb_lo") == 0) return (void *) cut_hb_lo;
  if (strcmp(str,"cut_hb_hi") == 0) return (void *) cut_hb_hi;
  if (strcmp(str,"cut_hb_lc") == 0) return (void *) cut_hb_lc;
  if (strcmp(str,"cut_hb_hc") == 0) return (void *) cut_hb_hc;
  if (strcmp(str,"b_hb_lo") == 0) return (void *) b_hb_lo;
  if (strcmp(str,"b_hb_hi") == 0) return (void *) b_hb_hi;
  if (strcmp(str,"shift_hb") == 0) return (void *) shift_hb;

  if (strcmp(str,"a_hb1") == 0) return (void *) a_hb1;
  if (strcmp(str,"theta_hb1_0") == 0) return (void *) theta_hb1_0;
  if (strcmp(str,"dtheta_hb1_ast") == 0) return (void *) dtheta_hb1_ast;
  if (strcmp(str,"b_hb1") == 0) return (void *) b_hb1;
  if (strcmp(str,"dtheta_hb1_c") == 0) return (void *) dtheta_hb1_c;

  if (strcmp(str,"a_hb2") == 0) return (void *) a_hb2;
  if (strcmp(str,"theta_hb2_0") == 0) return (void *) theta_hb2_0;
  if (strcmp(str,"dtheta_hb2_ast") == 0) return (void *) dtheta_hb2_ast;
  if (strcmp(str,"b_hb2") == 0) return (void *) b_hb2;
  if (strcmp(str,"dtheta_hb2_c") == 0) return (void *) dtheta_hb2_c;

  if (strcmp(str,"a_hb3") == 0) return (void *) a_hb3;
  if (strcmp(str,"theta_hb3_0") == 0) return (void *) theta_hb3_0;
  if (strcmp(str,"dtheta_hb3_ast") == 0) return (void *) dtheta_hb3_ast;
  if (strcmp(str,"b_hb3") == 0) return (void *) b_hb3;
  if (strcmp(str,"dtheta_hb3_c") == 0) return (void *) dtheta_hb3_c;

  if (strcmp(str,"a_hb4") == 0) return (void *) a_hb4;
  if (strcmp(str,"theta_hb4_0") == 0) return (void *) theta_hb4_0;
  if (strcmp(str,"dtheta_hb4_ast") == 0) return (void *) dtheta_hb4_ast;
  if (strcmp(str,"b_hb4") == 0) return (void *) b_hb4;
  if (strcmp(str,"dtheta_hb4_c") == 0) return (void *) dtheta_hb4_c;

  if (strcmp(str,"a_hb7") == 0) return (void *) a_hb7;
  if (strcmp(str,"theta_hb7_0") == 0) return (void *) theta_hb7_0;
  if (strcmp(str,"dtheta_hb7_ast") == 0) return (void *) dtheta_hb7_ast;
  if (strcmp(str,"b_hb7") == 0) return (void *) b_hb7;
  if (strcmp(str,"dtheta_hb7_c") == 0) return (void *) dtheta_hb7_c;

  if (strcmp(str,"a_hb8") == 0) return (void *) a_hb8;
  if (strcmp(str,"theta_hb8_0") == 0) return (void *) theta_hb8_0;
  if (strcmp(str,"dtheta_hb8_ast") == 0) return (void *) dtheta_hb8_ast;
  if (strcmp(str,"b_hb8") == 0) return (void *) b_hb8;
  if (strcmp(str,"dtheta_hb8_c") == 0) return (void *) dtheta_hb8_c;

  return NULL;
}
