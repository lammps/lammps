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

#include "pair_oxdna_xstk.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "mf_oxdna.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "atom_vec_ellipsoid.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MFOxdna;

/* ---------------------------------------------------------------------- */

PairOxdnaXstk::PairOxdnaXstk(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairOxdnaXstk::~PairOxdnaXstk()
{
  if (allocated) {

    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(k_xst);
    memory->destroy(cut_xst_0);
    memory->destroy(cut_xst_c);
    memory->destroy(cut_xst_lo);
    memory->destroy(cut_xst_hi);
    memory->destroy(cut_xst_lc);
    memory->destroy(cut_xst_hc);
    memory->destroy(cutsq_xst_hc);
    memory->destroy(b_xst_lo);
    memory->destroy(b_xst_hi);

    memory->destroy(a_xst1);
    memory->destroy(theta_xst1_0);
    memory->destroy(dtheta_xst1_ast);
    memory->destroy(b_xst1);
    memory->destroy(dtheta_xst1_c);

    memory->destroy(a_xst2);
    memory->destroy(theta_xst2_0);
    memory->destroy(dtheta_xst2_ast);
    memory->destroy(b_xst2);
    memory->destroy(dtheta_xst2_c);

    memory->destroy(a_xst3);
    memory->destroy(theta_xst3_0);
    memory->destroy(dtheta_xst3_ast);
    memory->destroy(b_xst3);
    memory->destroy(dtheta_xst3_c);

    memory->destroy(a_xst4);
    memory->destroy(theta_xst4_0);
    memory->destroy(dtheta_xst4_ast);
    memory->destroy(b_xst4);
    memory->destroy(dtheta_xst4_c);

    memory->destroy(a_xst7);
    memory->destroy(theta_xst7_0);
    memory->destroy(dtheta_xst7_ast);
    memory->destroy(b_xst7);
    memory->destroy(dtheta_xst7_c);

    memory->destroy(a_xst8);
    memory->destroy(theta_xst8_0);
    memory->destroy(dtheta_xst8_ast);
    memory->destroy(b_xst8);
    memory->destroy(dtheta_xst8_c);

  }
}

/* ----------------------------------------------------------------------
   compute function for oxDNA pair interactions
   hb=hydrogen bonding site

   NOTE: The cross-stacking interaction takes place between hb sites
------------------------------------------------------------------------- */

void PairOxdnaXstk::compute(int eflag, int vflag)
{

  double delf[3],delta[3],deltb[3]; // force, torque increment;
  double evdwl,fpair,finc,tpair,factor_lj;
  double delr_hb[3],delr_hb_norm[3],rsq_hb,r_hb,rinv_hb;
  double theta1,t1dir[3],cost1;
  double theta2,t2dir[3],cost2;
  double theta3,t3dir[3],cost3;
  double theta4,theta4p,t4dir[3],cost4;
  double theta7,theta7p,t7dir[3],cost7;
  double theta8,theta8p,t8dir[3],cost8;

  // distance COM-h-bonding site
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

  double f2,f4t1,f4t4,f4t2,f4t3,f4t7,f4t8;
  double df2,df4t1,df4t4,df4t2,df4t3,df4t7,df4t8,rsint;

  evdwl = 0.0;
  ev_init(eflag,vflag);

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

      f2 = F2(r_hb, k_xst[atype][btype], cut_xst_0[atype][btype],
           cut_xst_lc[atype][btype], cut_xst_hc[atype][btype], cut_xst_lo[atype][btype], cut_xst_hi[atype][btype],
           b_xst_lo[atype][btype], b_xst_hi[atype][btype], cut_xst_c[atype][btype]);

      // early rejection criterium
      if (f2) {

      cost1 = -1.0*MathExtra::dot3(ax,bx);
      if (cost1 >  1.0) cost1 =  1.0;
      if (cost1 < -1.0) cost1 = -1.0;
      theta1 = acos(cost1);

      f4t1 = F4(theta1, a_xst1[atype][btype], theta_xst1_0[atype][btype], dtheta_xst1_ast[atype][btype],
             b_xst1[atype][btype], dtheta_xst1_c[atype][btype]);

      // early rejection criterium
      if (f4t1) {

      cost2 = -1.0*MathExtra::dot3(ax,delr_hb_norm);
      if (cost2 >  1.0) cost2 =  1.0;
      if (cost2 < -1.0) cost2 = -1.0;
      theta2 = acos(cost2);

      f4t2 = F4(theta2, a_xst2[atype][btype], theta_xst2_0[atype][btype], dtheta_xst2_ast[atype][btype],
             b_xst2[atype][btype], dtheta_xst2_c[atype][btype]);

      // early rejection criterium
      if (f4t2) {

      cost3 = MathExtra::dot3(bx,delr_hb_norm);
      if (cost3 >  1.0) cost3 =  1.0;
      if (cost3 < -1.0) cost3 = -1.0;
      theta3 = acos(cost3);

      f4t3 = F4(theta3, a_xst3[atype][btype], theta_xst3_0[atype][btype], dtheta_xst3_ast[atype][btype],
             b_xst3[atype][btype], dtheta_xst3_c[atype][btype]);

      // early rejection criterium
      if (f4t3) {

      cost4 = MathExtra::dot3(az,bz);
      if (cost4 >  1.0) cost4 =  1.0;
      if (cost4 < -1.0) cost4 = -1.0;
      theta4 = acos(cost4);
      theta4p = MY_PI - theta4;

      f4t4 = F4(theta4, a_xst4[atype][btype], theta_xst4_0[atype][btype], dtheta_xst4_ast[atype][btype],
             b_xst4[atype][btype], dtheta_xst4_c[atype][btype]) +
             F4(theta4p, a_xst4[atype][btype], theta_xst4_0[atype][btype], dtheta_xst4_ast[atype][btype],
             b_xst4[atype][btype], dtheta_xst4_c[atype][btype]);

      // early rejection criterium
      if (f4t4) {

      cost7 = -1.0*MathExtra::dot3(az,delr_hb_norm);
      if (cost7 >  1.0) cost7 =  1.0;
      if (cost7 < -1.0) cost7 = -1.0;
      theta7 = acos(cost7);
      theta7p = MY_PI - theta7;

      f4t7 = F4(theta7, a_xst7[atype][btype], theta_xst7_0[atype][btype], dtheta_xst7_ast[atype][btype],
             b_xst7[atype][btype], dtheta_xst7_c[atype][btype]) +
             F4(theta7p, a_xst7[atype][btype], theta_xst7_0[atype][btype], dtheta_xst7_ast[atype][btype],
             b_xst7[atype][btype], dtheta_xst7_c[atype][btype]);

      // early rejection criterium
      if (f4t7) {

      cost8 = MathExtra::dot3(bz,delr_hb_norm);
      if (cost8 >  1.0) cost8 =  1.0;
      if (cost8 < -1.0) cost8 = -1.0;
      theta8 = acos(cost8);
      theta8p = MY_PI -theta8;

      f4t8 = F4(theta8, a_xst8[atype][btype], theta_xst8_0[atype][btype], dtheta_xst8_ast[atype][btype],
             b_xst8[atype][btype], dtheta_xst8_c[atype][btype]) +
             F4(theta8p, a_xst8[atype][btype], theta_xst8_0[atype][btype], dtheta_xst8_ast[atype][btype],
             b_xst8[atype][btype], dtheta_xst8_c[atype][btype]);


      evdwl = f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj;


      // early rejection criterium
      if (evdwl) {

      df2 = DF2(r_hb, k_xst[atype][btype], cut_xst_0[atype][btype],
            cut_xst_lc[atype][btype], cut_xst_hc[atype][btype], cut_xst_lo[atype][btype], cut_xst_hi[atype][btype],
            b_xst_lo[atype][btype], b_xst_hi[atype][btype]);

      df4t1 = DF4(theta1, a_xst1[atype][btype], theta_xst1_0[atype][btype], dtheta_xst1_ast[atype][btype],
              b_xst1[atype][btype], dtheta_xst1_c[atype][btype])/sin(theta1);

      df4t2 = DF4(theta2, a_xst2[atype][btype], theta_xst2_0[atype][btype], dtheta_xst2_ast[atype][btype],
              b_xst2[atype][btype], dtheta_xst2_c[atype][btype])/sin(theta2);

      df4t3 = DF4(theta3, a_xst3[atype][btype], theta_xst3_0[atype][btype], dtheta_xst3_ast[atype][btype],
              b_xst3[atype][btype], dtheta_xst3_c[atype][btype])/sin(theta3);

      rsint = 1.0/sin(theta4);
      df4t4 = DF4(theta4, a_xst4[atype][btype], theta_xst4_0[atype][btype], dtheta_xst4_ast[atype][btype],
              b_xst4[atype][btype], dtheta_xst4_c[atype][btype])*rsint -
              DF4(theta4p, a_xst4[atype][btype], theta_xst4_0[atype][btype], dtheta_xst4_ast[atype][btype],
              b_xst4[atype][btype], dtheta_xst4_c[atype][btype])*rsint;

      rsint = 1.0/sin(theta7);
      df4t7 = DF4(theta7, a_xst7[atype][btype], theta_xst7_0[atype][btype], dtheta_xst7_ast[atype][btype],
              b_xst7[atype][btype], dtheta_xst7_c[atype][btype])*rsint -
              DF4(theta7p, a_xst7[atype][btype], theta_xst7_0[atype][btype], dtheta_xst7_ast[atype][btype],
              b_xst7[atype][btype], dtheta_xst7_c[atype][btype])*rsint;

      rsint = 1.0/sin(theta8);
      df4t8 = DF4(theta8, a_xst8[atype][btype], theta_xst8_0[atype][btype], dtheta_xst8_ast[atype][btype],
              b_xst8[atype][btype], dtheta_xst8_c[atype][btype])*rsint -
              DF4(theta8p, a_xst8[atype][btype], theta_xst8_0[atype][btype], dtheta_xst8_ast[atype][btype],
              b_xst8[atype][btype], dtheta_xst8_c[atype][btype])*rsint;

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
      finc  = -df2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 * rinv_hb *factor_lj;
      fpair += finc;

      delf[0] += delr_hb[0] * finc;
      delf[1] += delr_hb[1] * finc;
      delf[2] += delr_hb[2] * finc;

      // theta2 force
      if (theta2) {

        finc  = -f2 * f4t1 * df4t2 * f4t3 * f4t4 * f4t7 * f4t8 * rinv_hb * factor_lj;
        fpair += finc;

        delf[0] += (delr_hb_norm[0]*cost2 + ax[0]) * finc;
        delf[1] += (delr_hb_norm[1]*cost2 + ax[1]) * finc;
        delf[2] += (delr_hb_norm[2]*cost2 + ax[2]) * finc;

      }

      // theta3 force
      if (theta3) {

        finc  = -f2 * f4t1 * f4t2 * df4t3 * f4t4 * f4t7 * f4t8 * rinv_hb * factor_lj;
        fpair += finc;

        delf[0] += (delr_hb_norm[0]*cost3 - bx[0]) * finc;
        delf[1] += (delr_hb_norm[1]*cost3 - bx[1]) * finc;
        delf[2] += (delr_hb_norm[2]*cost3 - bx[2]) * finc;

      }

      // theta7 force
      if (theta7) {

        finc  = -f2 * f4t1 * f4t2 * f4t3 * f4t4 * df4t7 * f4t8 * rinv_hb * factor_lj;
        fpair += finc;

        delf[0] += (delr_hb_norm[0]*cost7 + az[0]) * finc;
        delf[1] += (delr_hb_norm[1]*cost7 + az[1]) * finc;
        delf[2] += (delr_hb_norm[2]*cost7 + az[2]) * finc;

      }

      // theta8 force
      if (theta8) {

        finc  = -f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * df4t8 * rinv_hb * factor_lj;
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

        tpair = -f2 * df4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj;
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

        tpair = -f2 * f4t1 * df4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj;
        MathExtra::cross3(ax,delr_hb_norm,t2dir);

        delta[0] += t2dir[0]*tpair;
        delta[1] += t2dir[1]*tpair;
        delta[2] += t2dir[2]*tpair;

      }

      // theta3 torque
      if (theta3) {

        tpair = -f2 * f4t1 * f4t2 * df4t3 * f4t4 * f4t7 * f4t8 * factor_lj;
        MathExtra::cross3(bx,delr_hb_norm,t3dir);

        deltb[0] += t3dir[0]*tpair;
        deltb[1] += t3dir[1]*tpair;
        deltb[2] += t3dir[2]*tpair;

      }

      // theta4 torque
      if (theta4 && theta4p) {

        tpair = -f2 * f4t1 * f4t2 * f4t3 * df4t4 * f4t7 * f4t8 * factor_lj;
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

        tpair = -f2 * f4t1 * f4t2 * f4t3 * f4t4 * df4t7 * f4t8 * factor_lj;
        MathExtra::cross3(az,delr_hb_norm,t7dir);

        delta[0] += t7dir[0]*tpair;
        delta[1] += t7dir[1]*tpair;
        delta[2] += t7dir[2]*tpair;

      }

      // theta8 torque
      if (theta8) {

        tpair = -f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * df4t8 * factor_lj;
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

void PairOxdnaXstk::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(k_xst,n+1,n+1,"pair:k_xst");
  memory->create(cut_xst_0,n+1,n+1,"pair:cut_xst_0");
  memory->create(cut_xst_c,n+1,n+1,"pair:cut_xst_c");
  memory->create(cut_xst_lo,n+1,n+1,"pair:cut_xst_lo");
  memory->create(cut_xst_hi,n+1,n+1,"pair:cut_xst_hi");
  memory->create(cut_xst_lc,n+1,n+1,"pair:cut_xst_lc");
  memory->create(cut_xst_hc,n+1,n+1,"pair:cut_xst_hc");
  memory->create(b_xst_lo,n+1,n+1,"pair:b_xst_lo");
  memory->create(b_xst_hi,n+1,n+1,"pair:b_xst_hi");
  memory->create(cutsq_xst_hc,n+1,n+1,"pair:cutsq_xst_hc");

  memory->create(a_xst1,n+1,n+1,"pair:a_xst1");
  memory->create(theta_xst1_0,n+1,n+1,"pair:theta_xst1_0");
  memory->create(dtheta_xst1_ast,n+1,n+1,"pair:dtheta_xst1_ast");
  memory->create(b_xst1,n+1,n+1,"pair:b_xst1");
  memory->create(dtheta_xst1_c,n+1,n+1,"pair:dtheta_xst1_c");

  memory->create(a_xst2,n+1,n+1,"pair:a_xst2");
  memory->create(theta_xst2_0,n+1,n+1,"pair:theta_xst2_0");
  memory->create(dtheta_xst2_ast,n+1,n+1,"pair:dtheta_xst2_ast");
  memory->create(b_xst2,n+1,n+1,"pair:b_xst2");
  memory->create(dtheta_xst2_c,n+1,n+1,"pair:dtheta_xst2_c");

  memory->create(a_xst3,n+1,n+1,"pair:a_xst3");
  memory->create(theta_xst3_0,n+1,n+1,"pair:theta_xst3_0");
  memory->create(dtheta_xst3_ast,n+1,n+1,"pair:dtheta_xst3_ast");
  memory->create(b_xst3,n+1,n+1,"pair:b_xst3");
  memory->create(dtheta_xst3_c,n+1,n+1,"pair:dtheta_xst3_c");

  memory->create(a_xst4,n+1,n+1,"pair:a_xst4");
  memory->create(theta_xst4_0,n+1,n+1,"pair:theta_xst4_0");
  memory->create(dtheta_xst4_ast,n+1,n+1,"pair:dtheta_xst4_ast");
  memory->create(b_xst4,n+1,n+1,"pair:b_xst4");
  memory->create(dtheta_xst4_c,n+1,n+1,"pair:dtheta_xst4_c");

  memory->create(a_xst7,n+1,n+1,"pair:a_xst7");
  memory->create(theta_xst7_0,n+1,n+1,"pair:theta_xst7_0");
  memory->create(dtheta_xst7_ast,n+1,n+1,"pair:dtheta_xst7_ast");
  memory->create(b_xst7,n+1,n+1,"pair:b_xst7");
  memory->create(dtheta_xst7_c,n+1,n+1,"pair:dtheta_xst7_c");

  memory->create(a_xst8,n+1,n+1,"pair:a_xst8");
  memory->create(theta_xst8_0,n+1,n+1,"pair:theta_xst8_0");
  memory->create(dtheta_xst8_ast,n+1,n+1,"pair:dtheta_xst8_ast");
  memory->create(b_xst8,n+1,n+1,"pair:b_xst8");
  memory->create(dtheta_xst8_c,n+1,n+1,"pair:dtheta_xst8_c");

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairOxdnaXstk::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairOxdnaXstk::coeff(int narg, char **arg)
{
  int count;

  if (narg != 25) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/xstk");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  // cross-stacking interaction
  count = 0;

  double k_xst_one, cut_xst_0_one, cut_xst_c_one, cut_xst_lo_one, cut_xst_hi_one;
  double b_xst_lo_one, b_xst_hi_one, cut_xst_lc_one, cut_xst_hc_one;

  double a_xst1_one, theta_xst1_0_one, dtheta_xst1_ast_one;
  double b_xst1_one, dtheta_xst1_c_one;

  double a_xst2_one, theta_xst2_0_one, dtheta_xst2_ast_one;
  double b_xst2_one, dtheta_xst2_c_one;

  double a_xst3_one, theta_xst3_0_one, dtheta_xst3_ast_one;
  double b_xst3_one, dtheta_xst3_c_one;

  double a_xst4_one, theta_xst4_0_one, dtheta_xst4_ast_one;
  double b_xst4_one, dtheta_xst4_c_one;

  double a_xst7_one, theta_xst7_0_one, dtheta_xst7_ast_one;
  double b_xst7_one, dtheta_xst7_c_one;

  double a_xst8_one, theta_xst8_0_one, dtheta_xst8_ast_one;
  double b_xst8_one, dtheta_xst8_c_one;

  k_xst_one = force->numeric(FLERR,arg[2]);
  cut_xst_0_one = force->numeric(FLERR,arg[3]);
  cut_xst_c_one = force->numeric(FLERR,arg[4]);
  cut_xst_lo_one = force->numeric(FLERR,arg[5]);
  cut_xst_hi_one = force->numeric(FLERR,arg[6]);

  a_xst1_one = force->numeric(FLERR,arg[7]);
  theta_xst1_0_one = force->numeric(FLERR,arg[8]);
  dtheta_xst1_ast_one = force->numeric(FLERR,arg[9]);

  a_xst2_one = force->numeric(FLERR,arg[10]);
  theta_xst2_0_one = force->numeric(FLERR,arg[11]);
  dtheta_xst2_ast_one = force->numeric(FLERR,arg[12]);

  a_xst3_one = force->numeric(FLERR,arg[13]);
  theta_xst3_0_one = force->numeric(FLERR,arg[14]);
  dtheta_xst3_ast_one = force->numeric(FLERR,arg[15]);

  a_xst4_one = force->numeric(FLERR,arg[16]);
  theta_xst4_0_one = force->numeric(FLERR,arg[17]);
  dtheta_xst4_ast_one = force->numeric(FLERR,arg[18]);

  a_xst7_one = force->numeric(FLERR,arg[19]);
  theta_xst7_0_one = force->numeric(FLERR,arg[20]);
  dtheta_xst7_ast_one = force->numeric(FLERR,arg[21]);

  a_xst8_one = force->numeric(FLERR,arg[22]);
  theta_xst8_0_one = force->numeric(FLERR,arg[23]);
  dtheta_xst8_ast_one = force->numeric(FLERR,arg[24]);


  b_xst_lo_one = 0.25 * (cut_xst_lo_one - cut_xst_0_one) * (cut_xst_lo_one - cut_xst_0_one)/
        (0.5 * (cut_xst_lo_one - cut_xst_0_one) * (cut_xst_lo_one - cut_xst_0_one) -
        k_xst_one * 0.5 * (cut_xst_0_one -cut_xst_c_one) * (cut_xst_0_one - cut_xst_c_one)/k_xst_one);

  cut_xst_lc_one = cut_xst_lo_one - 0.5 * (cut_xst_lo_one - cut_xst_0_one)/b_xst_lo_one;;

  b_xst_hi_one = 0.25 * (cut_xst_hi_one - cut_xst_0_one) * (cut_xst_hi_one - cut_xst_0_one)/
        (0.5 * (cut_xst_hi_one - cut_xst_0_one) * (cut_xst_hi_one - cut_xst_0_one) -
        k_xst_one * 0.5 * (cut_xst_0_one -cut_xst_c_one) * (cut_xst_0_one - cut_xst_c_one)/k_xst_one);

  cut_xst_hc_one = cut_xst_hi_one - 0.5* (cut_xst_hi_one - cut_xst_0_one)/b_xst_hi_one;


  b_xst1_one = a_xst1_one*a_xst1_one*dtheta_xst1_ast_one*dtheta_xst1_ast_one/(1-a_xst1_one*dtheta_xst1_ast_one*dtheta_xst1_ast_one);
  dtheta_xst1_c_one = 1/(a_xst1_one*dtheta_xst1_ast_one);

  b_xst2_one = a_xst2_one*a_xst2_one*dtheta_xst2_ast_one*dtheta_xst2_ast_one/(1-a_xst2_one*dtheta_xst2_ast_one*dtheta_xst2_ast_one);
  dtheta_xst2_c_one = 1/(a_xst2_one*dtheta_xst2_ast_one);

  b_xst3_one = a_xst3_one*a_xst3_one*dtheta_xst3_ast_one*dtheta_xst3_ast_one/(1-a_xst3_one*dtheta_xst3_ast_one*dtheta_xst3_ast_one);
  dtheta_xst3_c_one = 1/(a_xst3_one*dtheta_xst3_ast_one);

  b_xst4_one = a_xst4_one*a_xst4_one*dtheta_xst4_ast_one*dtheta_xst4_ast_one/(1-a_xst4_one*dtheta_xst4_ast_one*dtheta_xst4_ast_one);
  dtheta_xst4_c_one = 1/(a_xst4_one*dtheta_xst4_ast_one);

  b_xst7_one = a_xst7_one*a_xst7_one*dtheta_xst7_ast_one*dtheta_xst7_ast_one/(1-a_xst7_one*dtheta_xst7_ast_one*dtheta_xst7_ast_one);
  dtheta_xst7_c_one = 1/(a_xst7_one*dtheta_xst7_ast_one);

  b_xst8_one = a_xst8_one*a_xst8_one*dtheta_xst8_ast_one*dtheta_xst8_ast_one/(1-a_xst8_one*dtheta_xst8_ast_one*dtheta_xst8_ast_one);
  dtheta_xst8_c_one = 1/(a_xst8_one*dtheta_xst8_ast_one);

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {

      k_xst[i][j] = k_xst_one;
      cut_xst_0[i][j] = cut_xst_0_one;
      cut_xst_c[i][j] = cut_xst_c_one;
      cut_xst_lo[i][j] = cut_xst_lo_one;
      cut_xst_hi[i][j] = cut_xst_hi_one;
      cut_xst_lc[i][j] = cut_xst_lc_one;
      cut_xst_hc[i][j] = cut_xst_hc_one;
      b_xst_lo[i][j] = b_xst_lo_one;
      b_xst_hi[i][j] = b_xst_hi_one;

      a_xst1[i][j] = a_xst1_one;
      theta_xst1_0[i][j] = theta_xst1_0_one;
      dtheta_xst1_ast[i][j] = dtheta_xst1_ast_one;
      b_xst1[i][j] = b_xst1_one;
      dtheta_xst1_c[i][j] = dtheta_xst1_c_one;

      a_xst2[i][j] = a_xst2_one;
      theta_xst2_0[i][j] = theta_xst2_0_one;
      dtheta_xst2_ast[i][j] = dtheta_xst2_ast_one;
      b_xst2[i][j] = b_xst2_one;
      dtheta_xst2_c[i][j] = dtheta_xst2_c_one;

      a_xst3[i][j] = a_xst3_one;
      theta_xst3_0[i][j] = theta_xst3_0_one;
      dtheta_xst3_ast[i][j] = dtheta_xst3_ast_one;
      b_xst3[i][j] = b_xst3_one;
      dtheta_xst3_c[i][j] = dtheta_xst3_c_one;

      a_xst4[i][j] = a_xst4_one;
      theta_xst4_0[i][j] = theta_xst4_0_one;
      dtheta_xst4_ast[i][j] = dtheta_xst4_ast_one;
      b_xst4[i][j] = b_xst4_one;
      dtheta_xst4_c[i][j] = dtheta_xst4_c_one;

      a_xst7[i][j] = a_xst7_one;
      theta_xst7_0[i][j] = theta_xst7_0_one;
      dtheta_xst7_ast[i][j] = dtheta_xst7_ast_one;
      b_xst7[i][j] = b_xst7_one;
      dtheta_xst7_c[i][j] = dtheta_xst7_c_one;

      a_xst8[i][j] = a_xst8_one;
      theta_xst8_0[i][j] = theta_xst8_0_one;
      dtheta_xst8_ast[i][j] = dtheta_xst8_ast_one;
      b_xst8[i][j] = b_xst8_one;
      dtheta_xst8_c[i][j] = dtheta_xst8_c_one;

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/xstk");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairOxdnaXstk::init_style()
{
  int irequest;

  // request regular neighbor lists

  irequest = neighbor->request(this,instance_me);

}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use regular
------------------------------------------------------------------------- */

void PairOxdnaXstk::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  if (id  > 0) error->all(FLERR,"Respa not supported");

}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairOxdnaXstk::init_one(int i, int j)
{

  if (setflag[i][j] == 0) {
    error->all(FLERR,"Coefficient mixing not defined in oxDNA");
  }
  if (offset_flag) {
    error->all(FLERR,"Offset not supported in oxDNA");
  }

  k_xst[j][i] = k_xst[i][j];
  cut_xst_0[j][i] = cut_xst_0[i][j];
  cut_xst_c[j][i] = cut_xst_c[i][j];
  cut_xst_lo[j][i] = cut_xst_lo[i][j];
  cut_xst_hi[j][i] = cut_xst_hi[i][j];
  b_xst_lo[j][i] = b_xst_lo[i][j];
  b_xst_hi[j][i] = b_xst_hi[i][j];
  cut_xst_lc[j][i] = cut_xst_lc[i][j];
  cut_xst_hc[j][i] = cut_xst_hc[i][j];

  a_xst1[j][i] = a_xst1[i][j];
  theta_xst1_0[j][i] = theta_xst1_0[i][j];
  dtheta_xst1_ast[j][i] = dtheta_xst1_ast[i][j];
  b_xst1[j][i] = b_xst1[i][j];
  dtheta_xst1_c[j][i] = dtheta_xst1_c[i][j];

  a_xst2[j][i] = a_xst2[i][j];
  theta_xst2_0[j][i] = theta_xst2_0[i][j];
  dtheta_xst2_ast[j][i] = dtheta_xst2_ast[i][j];
  b_xst2[j][i] = b_xst2[i][j];
  dtheta_xst2_c[j][i] = dtheta_xst2_c[i][j];

  a_xst3[j][i] = a_xst3[i][j];
  theta_xst3_0[j][i] = theta_xst3_0[i][j];
  dtheta_xst3_ast[j][i] = dtheta_xst3_ast[i][j];
  b_xst3[j][i] = b_xst3[i][j];
  dtheta_xst3_c[j][i] = dtheta_xst3_c[i][j];

  a_xst4[j][i] = a_xst4[i][j];
  theta_xst4_0[j][i] = theta_xst4_0[i][j];
  dtheta_xst4_ast[j][i] = dtheta_xst4_ast[i][j];
  b_xst4[j][i] = b_xst4[i][j];
  dtheta_xst4_c[j][i] = dtheta_xst4_c[i][j];

  a_xst7[j][i] = a_xst7[i][j];
  theta_xst7_0[j][i] = theta_xst7_0[i][j];
  dtheta_xst7_ast[j][i] = dtheta_xst7_ast[i][j];
  b_xst7[j][i] = b_xst7[i][j];
  dtheta_xst7_c[j][i] = dtheta_xst7_c[i][j];

  a_xst8[j][i] = a_xst8[i][j];
  theta_xst8_0[j][i] = theta_xst8_0[i][j];
  dtheta_xst8_ast[j][i] = dtheta_xst8_ast[i][j];
  b_xst8[j][i] = b_xst8[i][j];
  dtheta_xst8_c[j][i] = dtheta_xst8_c[i][j];

  cutsq_xst_hc[i][j] = cut_xst_hc[i][j]*cut_xst_hc[i][j];
  cutsq_xst_hc[j][i] = cutsq_xst_hc[i][j];

  // set the master list distance cutoff
  return cut_xst_hc[i][j];

}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdnaXstk::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {

        fwrite(&k_xst[i][j],sizeof(double),1,fp);
        fwrite(&cut_xst_0[i][j],sizeof(double),1,fp);
        fwrite(&cut_xst_c[i][j],sizeof(double),1,fp);
        fwrite(&cut_xst_lo[i][j],sizeof(double),1,fp);
        fwrite(&cut_xst_hi[i][j],sizeof(double),1,fp);
        fwrite(&cut_xst_lc[i][j],sizeof(double),1,fp);
        fwrite(&cut_xst_hc[i][j],sizeof(double),1,fp);
        fwrite(&b_xst_lo[i][j],sizeof(double),1,fp);
        fwrite(&b_xst_hi[i][j],sizeof(double),1,fp);

        fwrite(&a_xst1[i][j],sizeof(double),1,fp);
        fwrite(&theta_xst1_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_xst1_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_xst1[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_xst1_c[i][j],sizeof(double),1,fp);

        fwrite(&a_xst2[i][j],sizeof(double),1,fp);
        fwrite(&theta_xst2_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_xst2_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_xst2[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_xst2_c[i][j],sizeof(double),1,fp);

        fwrite(&a_xst3[i][j],sizeof(double),1,fp);
        fwrite(&theta_xst3_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_xst3_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_xst3[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_xst3_c[i][j],sizeof(double),1,fp);

        fwrite(&a_xst4[i][j],sizeof(double),1,fp);
        fwrite(&theta_xst4_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_xst4_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_xst4[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_xst4_c[i][j],sizeof(double),1,fp);

        fwrite(&a_xst7[i][j],sizeof(double),1,fp);
        fwrite(&theta_xst7_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_xst7_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_xst7[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_xst7_c[i][j],sizeof(double),1,fp);

        fwrite(&a_xst8[i][j],sizeof(double),1,fp);
        fwrite(&theta_xst8_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_xst8_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_xst8[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_xst8_c[i][j],sizeof(double),1,fp);

    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdnaXstk::read_restart(FILE *fp)
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

          fread(&k_xst[i][j],sizeof(double),1,fp);
          fread(&cut_xst_0[i][j],sizeof(double),1,fp);
          fread(&cut_xst_c[i][j],sizeof(double),1,fp);
          fread(&cut_xst_lo[i][j],sizeof(double),1,fp);
          fread(&cut_xst_hi[i][j],sizeof(double),1,fp);
          fread(&cut_xst_lc[i][j],sizeof(double),1,fp);
          fread(&cut_xst_hc[i][j],sizeof(double),1,fp);
          fread(&b_xst_lo[i][j],sizeof(double),1,fp);
          fread(&b_xst_hi[i][j],sizeof(double),1,fp);

          fread(&a_xst1[i][j],sizeof(double),1,fp);
          fread(&theta_xst1_0[i][j],sizeof(double),1,fp);
          fread(&dtheta_xst1_ast[i][j],sizeof(double),1,fp);
          fread(&b_xst1[i][j],sizeof(double),1,fp);
          fread(&dtheta_xst1_c[i][j],sizeof(double),1,fp);

          fread(&a_xst2[i][j],sizeof(double),1,fp);
          fread(&theta_xst2_0[i][j],sizeof(double),1,fp);
          fread(&dtheta_xst2_ast[i][j],sizeof(double),1,fp);
          fread(&b_xst2[i][j],sizeof(double),1,fp);
          fread(&dtheta_xst2_c[i][j],sizeof(double),1,fp);

          fread(&a_xst3[i][j],sizeof(double),1,fp);
          fread(&theta_xst3_0[i][j],sizeof(double),1,fp);
          fread(&dtheta_xst3_ast[i][j],sizeof(double),1,fp);
          fread(&b_xst3[i][j],sizeof(double),1,fp);
          fread(&dtheta_xst3_c[i][j],sizeof(double),1,fp);

          fread(&a_xst4[i][j],sizeof(double),1,fp);
          fread(&theta_xst4_0[i][j],sizeof(double),1,fp);
          fread(&dtheta_xst4_ast[i][j],sizeof(double),1,fp);
          fread(&b_xst4[i][j],sizeof(double),1,fp);
          fread(&dtheta_xst4_c[i][j],sizeof(double),1,fp);

          fread(&a_xst7[i][j],sizeof(double),1,fp);
          fread(&theta_xst7_0[i][j],sizeof(double),1,fp);
          fread(&dtheta_xst7_ast[i][j],sizeof(double),1,fp);
          fread(&b_xst7[i][j],sizeof(double),1,fp);
          fread(&dtheta_xst7_c[i][j],sizeof(double),1,fp);

          fread(&a_xst8[i][j],sizeof(double),1,fp);
          fread(&theta_xst8_0[i][j],sizeof(double),1,fp);
          fread(&dtheta_xst8_ast[i][j],sizeof(double),1,fp);
          fread(&b_xst8[i][j],sizeof(double),1,fp);
          fread(&dtheta_xst8_c[i][j],sizeof(double),1,fp);

        }

        MPI_Bcast(&k_xst[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_xst_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_xst_c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_xst_lo[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_xst_hi[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_xst_lc[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_xst_hc[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_xst_lo[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_xst_hi[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_xst1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_xst1_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_xst1_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_xst1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_xst1_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_xst2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_xst2_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_xst2_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_xst2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_xst2_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_xst3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_xst3_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_xst3_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_xst3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_xst3_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_xst4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_xst4_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_xst4_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_xst4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_xst4_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_xst7[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_xst7_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_xst7_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_xst7[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_xst7_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_xst8[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_xst8_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_xst8_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_xst8[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_xst8_c[i][j],1,MPI_DOUBLE,0,world);

      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdnaXstk::write_restart_settings(FILE *fp)
{
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdnaXstk::read_restart_settings(FILE *fp)
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

void PairOxdnaXstk::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d\
         %g %g %g %g %g\
         %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         \n",i,
        k_xst[i][i],cut_xst_0[i][i],cut_xst_c[i][i],cut_xst_lo[i][i],cut_xst_hi[i][i],
        cut_xst_lc[i][i],cut_xst_hc[i][i],b_xst_lo[i][i],b_xst_hi[i][i],
        a_xst1[i][i],theta_xst1_0[i][i],dtheta_xst1_ast[i][i],b_xst1[i][i],dtheta_xst1_c[i][i],
        a_xst2[i][i],theta_xst2_0[i][i],dtheta_xst2_ast[i][i],b_xst2[i][i],dtheta_xst2_c[i][i],
        a_xst3[i][i],theta_xst3_0[i][i],dtheta_xst3_ast[i][i],b_xst3[i][i],dtheta_xst3_c[i][i],
        a_xst4[i][i],theta_xst4_0[i][i],dtheta_xst4_ast[i][i],b_xst4[i][i],dtheta_xst4_c[i][i],
        a_xst7[i][i],theta_xst7_0[i][i],dtheta_xst7_ast[i][i],b_xst7[i][i],dtheta_xst7_c[i][i],
        a_xst8[i][i],theta_xst8_0[i][i],dtheta_xst8_ast[i][i],b_xst8[i][i],dtheta_xst8_c[i][i]);

}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairOxdnaXstk::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d\
         %g %g %g %g %g\
         %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         \n",i,j,
        k_xst[i][j],cut_xst_0[i][j],cut_xst_c[i][j],cut_xst_lo[i][j],cut_xst_hi[i][j],
        cut_xst_lc[i][j],cut_xst_hc[i][j],b_xst_lo[i][j],b_xst_hi[i][j],
        a_xst1[i][j],theta_xst1_0[i][j],dtheta_xst1_ast[i][j],b_xst1[i][j],dtheta_xst1_c[i][j],
        a_xst2[i][j],theta_xst2_0[i][j],dtheta_xst2_ast[i][j],b_xst2[i][j],dtheta_xst2_c[i][j],
        a_xst3[i][j],theta_xst3_0[i][j],dtheta_xst3_ast[i][j],b_xst3[i][j],dtheta_xst3_c[i][j],
        a_xst4[i][j],theta_xst4_0[i][j],dtheta_xst4_ast[i][j],b_xst4[i][j],dtheta_xst4_c[i][j],
        a_xst7[i][j],theta_xst7_0[i][j],dtheta_xst7_ast[i][j],b_xst7[i][j],dtheta_xst7_c[i][j],
        a_xst8[i][j],theta_xst8_0[i][j],dtheta_xst8_ast[i][j],b_xst8[i][j],dtheta_xst8_c[i][j]);

}

/* ---------------------------------------------------------------------- */

void *PairOxdnaXstk::extract(const char *str, int &dim)
{
  dim = 2;

  if (strcmp(str,"k_xst") == 0) return (void *) k_xst;
  if (strcmp(str,"cut_xst_0") == 0) return (void *) cut_xst_0;
  if (strcmp(str,"cut_xst_c") == 0) return (void *) cut_xst_c;
  if (strcmp(str,"cut_xst_lo") == 0) return (void *) cut_xst_lo;
  if (strcmp(str,"cut_xst_hi") == 0) return (void *) cut_xst_hi;
  if (strcmp(str,"cut_xst_lc") == 0) return (void *) cut_xst_lc;
  if (strcmp(str,"cut_xst_hc") == 0) return (void *) cut_xst_hc;
  if (strcmp(str,"b_xst_lo") == 0) return (void *) b_xst_lo;
  if (strcmp(str,"b_xst_hi") == 0) return (void *) b_xst_hi;

  if (strcmp(str,"a_xst1") == 0) return (void *) a_xst1;
  if (strcmp(str,"theta_xst1_0") == 0) return (void *) theta_xst1_0;
  if (strcmp(str,"dtheta_xst1_ast") == 0) return (void *) dtheta_xst1_ast;
  if (strcmp(str,"b_xst1") == 0) return (void *) b_xst1;
  if (strcmp(str,"dtheta_xst1_c") == 0) return (void *) dtheta_xst1_c;

  if (strcmp(str,"a_xst2") == 0) return (void *) a_xst2;
  if (strcmp(str,"theta_xst2_0") == 0) return (void *) theta_xst2_0;
  if (strcmp(str,"dtheta_xst2_ast") == 0) return (void *) dtheta_xst2_ast;
  if (strcmp(str,"b_xst2") == 0) return (void *) b_xst2;
  if (strcmp(str,"dtheta_xst2_c") == 0) return (void *) dtheta_xst2_c;

  if (strcmp(str,"a_xst3") == 0) return (void *) a_xst3;
  if (strcmp(str,"theta_xst3_0") == 0) return (void *) theta_xst3_0;
  if (strcmp(str,"dtheta_xst3_ast") == 0) return (void *) dtheta_xst3_ast;
  if (strcmp(str,"b_xst3") == 0) return (void *) b_xst3;
  if (strcmp(str,"dtheta_xst3_c") == 0) return (void *) dtheta_xst3_c;

  if (strcmp(str,"a_xst4") == 0) return (void *) a_xst4;
  if (strcmp(str,"theta_xst4_0") == 0) return (void *) theta_xst4_0;
  if (strcmp(str,"dtheta_xst4_ast") == 0) return (void *) dtheta_xst4_ast;
  if (strcmp(str,"b_xst4") == 0) return (void *) b_xst4;
  if (strcmp(str,"dtheta_xst4_c") == 0) return (void *) dtheta_xst4_c;

  if (strcmp(str,"a_xst7") == 0) return (void *) a_xst7;
  if (strcmp(str,"theta_xst7_0") == 0) return (void *) theta_xst7_0;
  if (strcmp(str,"dtheta_xst7_ast") == 0) return (void *) dtheta_xst7_ast;
  if (strcmp(str,"b_xst7") == 0) return (void *) b_xst7;
  if (strcmp(str,"dtheta_xst7_c") == 0) return (void *) dtheta_xst7_c;

  if (strcmp(str,"a_xst8") == 0) return (void *) a_xst8;
  if (strcmp(str,"theta_xst8_0") == 0) return (void *) theta_xst8_0;
  if (strcmp(str,"dtheta_xst8_ast") == 0) return (void *) dtheta_xst8_ast;
  if (strcmp(str,"b_xst8") == 0) return (void *) b_xst8;
  if (strcmp(str,"dtheta_xst8_c") == 0) return (void *) dtheta_xst8_c;

  return NULL;
}
