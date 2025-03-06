// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include "pair_oxdna_coaxstk.h"

#include "atom.h"
#include "comm.h"
#include "constants_oxdna.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "mf_oxdna.h"
#include "neigh_list.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MFOxdna;

/* ---------------------------------------------------------------------- */

PairOxdnaCoaxstk::PairOxdnaCoaxstk(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  writedata = 1;
  trim_flag = 0;
}

/* ---------------------------------------------------------------------- */

PairOxdnaCoaxstk::~PairOxdnaCoaxstk()
{
  if (allocated) {

    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(k_cxst);
    memory->destroy(cut_cxst_0);
    memory->destroy(cut_cxst_c);
    memory->destroy(cut_cxst_lo);
    memory->destroy(cut_cxst_hi);
    memory->destroy(cut_cxst_lc);
    memory->destroy(cut_cxst_hc);
    memory->destroy(cutsq_cxst_hc);
    memory->destroy(b_cxst_lo);
    memory->destroy(b_cxst_hi);

    memory->destroy(a_cxst1);
    memory->destroy(theta_cxst1_0);
    memory->destroy(dtheta_cxst1_ast);
    memory->destroy(b_cxst1);
    memory->destroy(dtheta_cxst1_c);

    memory->destroy(a_cxst4);
    memory->destroy(theta_cxst4_0);
    memory->destroy(dtheta_cxst4_ast);
    memory->destroy(b_cxst4);
    memory->destroy(dtheta_cxst4_c);

    memory->destroy(a_cxst5);
    memory->destroy(theta_cxst5_0);
    memory->destroy(dtheta_cxst5_ast);
    memory->destroy(b_cxst5);
    memory->destroy(dtheta_cxst5_c);

    memory->destroy(a_cxst6);
    memory->destroy(theta_cxst6_0);
    memory->destroy(dtheta_cxst6_ast);
    memory->destroy(b_cxst6);
    memory->destroy(dtheta_cxst6_c);

    memory->destroy(a_cxst3p);
    memory->destroy(cosphi_cxst3p_ast);
    memory->destroy(b_cxst3p);
    memory->destroy(cosphi_cxst3p_c);
    memory->destroy(a_cxst4p);
    memory->destroy(cosphi_cxst4p_ast);
    memory->destroy(b_cxst4p);
    memory->destroy(cosphi_cxst4p_c);

  }
}

/* ----------------------------------------------------------------------
   compute function for oxDNA pair interactions
   st=stacking site
------------------------------------------------------------------------- */

void PairOxdnaCoaxstk::compute(int eflag, int vflag)
{
  double delf[3],delt[3],delta[3],deltb[3]; // force, torque increment;
  double evdwl,finc,tpair,factor_lj;
  double v1tmp[3],v2tmp[3],v3tmp[3];
  double delr_ss[3],delr_ss_norm[3],rsq_ss,r_ss,rinv_ss;
  double delr_st[3],delr_st_norm[3],rsq_st,r_st,rinv_st;
  double theta1,theta1p,t1dir[3],cost1;
  double theta4,t4dir[3],cost4;
  double theta5,theta5p,t5dir[3],cost5;
  double theta6,theta6p,t6dir[3],cost6;
  double cosphi3;

  double gamma,gammacub,rinv_ss_cub,fac;
  double aybx,azbx,rax,ray,raz,rbx;
  double dcdr,dcdrbx;
  double dcdaxbx,dcdaybx,dcdazbx;
  double dcdrax,dcdray,dcdraz;

  // distances COM-backbone site, COM-stacking site
  double d_cs = ConstantsOxdna::get_d_cs();
  double d_cst = ConstantsOxdna::get_d_cst();
  // vectors COM-backbone site, COM-stacking site in lab frame
  double ra_cs[3],ra_cst[3];
  double rb_cs[3],rb_cst[3];
  // Cartesian unit vectors in lab frame
  double ax[3],ay[3],az[3];
  double bx[3],bz[3];

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  int *type = atom->type;

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int *alist,*blist,*numneigh,**firstneigh;
  double *special_lj = force->special_lj;

  int a,b,ia,ib,anum,bnum,atype,btype;

  double f2,f4t1,f4t4,f4t5,f4t6,f5c3;
  double df2,df4t1,df4t4,df4t5,df4t6,df5c3,rsint;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  anum = list->inum;
  alist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // n(x/y/z)_xtrct = extracted local unit vectors in lab frame from oxdna_excv
  int dim;
  nx_xtrct = (double **) force->pair->extract("nx",dim);
  ny_xtrct = (double **) force->pair->extract("ny",dim);
  nz_xtrct = (double **) force->pair->extract("nz",dim);

  // loop over pair interaction neighbors of my atoms

  for (ia = 0; ia < anum; ia++) {

    a = alist[ia];
    atype = type[a];

    ax[0] = nx_xtrct[a][0];
    ax[1] = nx_xtrct[a][1];
    ax[2] = nx_xtrct[a][2];
    // a(y/z) not needed here as oxDNA(1) co-linear

    // vector COM a - stacking site a
    ra_cst[0] = d_cst*ax[0];
    ra_cst[1] = d_cst*ax[1];
    ra_cst[2] = d_cst*ax[2];

    // vector COM a - backbone site a
    ra_cs[0] = d_cs*ax[0];
    ra_cs[1] = d_cs*ax[1];
    ra_cs[2] = d_cs*ax[2];

    blist = firstneigh[a];
    bnum = numneigh[a];

    for (ib = 0; ib < bnum; ib++) {

      b = blist[ib];
      factor_lj = special_lj[sbmask(b)]; // = 0 for nearest neighbors
      b &= NEIGHMASK;

      btype = type[b];

      bx[0] = nx_xtrct[b][0];
      bx[1] = nx_xtrct[b][1];
      bx[2] = nx_xtrct[b][2];
      // b(y/z) not needed here as oxDNA(1) co-linear

      // vector COM b - stacking site b
      rb_cst[0] = d_cst*bx[0];
      rb_cst[1] = d_cst*bx[1];
      rb_cst[2] = d_cst*bx[2];

      // vector stacking site b to a
      delr_st[0] = x[a][0] + ra_cst[0] - x[b][0] - rb_cst[0];
      delr_st[1] = x[a][1] + ra_cst[1] - x[b][1] - rb_cst[1];
      delr_st[2] = x[a][2] + ra_cst[2] - x[b][2] - rb_cst[2];

      rsq_st = delr_st[0]*delr_st[0] + delr_st[1]*delr_st[1] + delr_st[2]*delr_st[2];
      r_st = sqrt(rsq_st);
      rinv_st = 1.0/r_st;

      delr_st_norm[0] = delr_st[0] * rinv_st;
      delr_st_norm[1] = delr_st[1] * rinv_st;
      delr_st_norm[2] = delr_st[2] * rinv_st;

      // vector COM b - backbone site b
      rb_cs[0] = d_cs*bx[0];
      rb_cs[1] = d_cs*bx[1];
      rb_cs[2] = d_cs*bx[2];

      // vector backbone site b to a
      delr_ss[0] = (x[a][0] + ra_cs[0] - x[b][0] - rb_cs[0]);
      delr_ss[1] = (x[a][1] + ra_cs[1] - x[b][1] - rb_cs[1]);
      delr_ss[2] = (x[a][2] + ra_cs[2] - x[b][2] - rb_cs[2]);

      rsq_ss = delr_ss[0]*delr_ss[0] + delr_ss[1]*delr_ss[1] + delr_ss[2]*delr_ss[2];
      r_ss = sqrt(rsq_ss);
      rinv_ss = 1.0/r_ss;

      delr_ss_norm[0] = delr_ss[0] * rinv_ss;
      delr_ss_norm[1] = delr_ss[1] * rinv_ss;
      delr_ss_norm[2] = delr_ss[2] * rinv_ss;

      cost1 = -1.0*MathExtra::dot3(ax,bx);
      if (cost1 >  1.0) cost1 =  1.0;
      if (cost1 < -1.0) cost1 = -1.0;
      theta1 = acos(cost1);
      theta1p = 2 * MY_PI - theta1;

      f4t1 = F4(theta1, a_cxst1[atype][btype], theta_cxst1_0[atype][btype], dtheta_cxst1_ast[atype][btype],
             b_cxst1[atype][btype], dtheta_cxst1_c[atype][btype]) +
             F4(theta1p, a_cxst1[atype][btype], theta_cxst1_0[atype][btype], dtheta_cxst1_ast[atype][btype],
             b_cxst1[atype][btype], dtheta_cxst1_c[atype][btype]);

      // early rejection criterium
      if (f4t1) {

      az[0] = nz_xtrct[a][0];
      az[1] = nz_xtrct[a][1];
      az[2] = nz_xtrct[a][2];
      bz[0] = nz_xtrct[b][0];
      bz[1] = nz_xtrct[b][1];
      bz[2] = nz_xtrct[b][2];

      cost4 = MathExtra::dot3(az,bz);
      if (cost4 >  1.0) cost4 =  1.0;
      if (cost4 < -1.0) cost4 = -1.0;
      theta4 = acos(cost4);

      f4t4 = F4(theta4, a_cxst4[atype][btype], theta_cxst4_0[atype][btype], dtheta_cxst4_ast[atype][btype],
             b_cxst4[atype][btype], dtheta_cxst4_c[atype][btype]);

      // early rejection criterium
      if (f4t4) {

      cost5 = MathExtra::dot3(delr_st_norm,az);
      if (cost5 >  1.0) cost5 =  1.0;
      if (cost5 < -1.0) cost5 = -1.0;
      theta5 = acos(cost5);
      theta5p = MY_PI - theta5;

      f4t5 = F4(theta5, a_cxst5[atype][btype], theta_cxst5_0[atype][btype], dtheta_cxst5_ast[atype][btype],
             b_cxst5[atype][btype], dtheta_cxst5_c[atype][btype]) +
             F4(theta5p, a_cxst5[atype][btype], theta_cxst5_0[atype][btype], dtheta_cxst5_ast[atype][btype],
             b_cxst5[atype][btype], dtheta_cxst5_c[atype][btype]);

      // early rejection criterium
      if (f4t5) {

      cost6 = MathExtra::dot3(delr_st_norm,bz);
      if (cost6 >  1.0) cost6 =  1.0;
      if (cost6 < -1.0) cost6 = -1.0;
      theta6 = acos(cost6);
      theta6p = MY_PI - theta6;

      f4t6 = F4(theta6, a_cxst6[atype][btype], theta_cxst6_0[atype][btype], dtheta_cxst6_ast[atype][btype],
             b_cxst6[atype][btype], dtheta_cxst6_c[atype][btype]) +
             F4(theta6p, a_cxst6[atype][btype], theta_cxst6_0[atype][btype], dtheta_cxst6_ast[atype][btype],
             b_cxst6[atype][btype], dtheta_cxst6_c[atype][btype]);

      MathExtra::cross3(delr_ss_norm,ax,v1tmp);
      cosphi3 = MathExtra::dot3(delr_st_norm,v1tmp);
      if (cosphi3 >  1.0) cosphi3 =  1.0;
      if (cosphi3 < -1.0) cosphi3 = -1.0;

      f2 = F2(r_st, k_cxst[atype][btype], cut_cxst_0[atype][btype],
           cut_cxst_lc[atype][btype], cut_cxst_hc[atype][btype], cut_cxst_lo[atype][btype], cut_cxst_hi[atype][btype],
           b_cxst_lo[atype][btype], b_cxst_hi[atype][btype], cut_cxst_c[atype][btype]);


      f5c3 = F5(cosphi3, a_cxst3p[atype][btype], cosphi_cxst3p_ast[atype][btype], b_cxst3p[atype][btype],
             cosphi_cxst3p_c[atype][btype]);

      evdwl = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5c3 * f5c3 * factor_lj;

      // early rejection criterium
      if (evdwl) {

      df2 = DF2(r_st, k_cxst[atype][btype], cut_cxst_0[atype][btype],
            cut_cxst_lc[atype][btype], cut_cxst_hc[atype][btype], cut_cxst_lo[atype][btype], cut_cxst_hi[atype][btype],
            b_cxst_lo[atype][btype], b_cxst_hi[atype][btype]);

      rsint = 1.0/sin(theta1);
      df4t1 = DF4(theta1, a_cxst1[atype][btype], theta_cxst1_0[atype][btype], dtheta_cxst1_ast[atype][btype],
              b_cxst1[atype][btype], dtheta_cxst1_c[atype][btype])*rsint -
              DF4(theta1p, a_cxst1[atype][btype], theta_cxst1_0[atype][btype], dtheta_cxst1_ast[atype][btype],
              b_cxst1[atype][btype], dtheta_cxst1_c[atype][btype])*rsint;

      df4t4 = DF4(theta4, a_cxst4[atype][btype], theta_cxst4_0[atype][btype], dtheta_cxst4_ast[atype][btype],
              b_cxst4[atype][btype], dtheta_cxst4_c[atype][btype])/sin(theta4);

      rsint = 1.0/sin(theta5);
      df4t5 = DF4(theta5, a_cxst5[atype][btype], theta_cxst5_0[atype][btype], dtheta_cxst5_ast[atype][btype],
              b_cxst5[atype][btype], dtheta_cxst5_c[atype][btype])*rsint -
              DF4(theta5p, a_cxst5[atype][btype], theta_cxst5_0[atype][btype], dtheta_cxst5_ast[atype][btype],
              b_cxst5[atype][btype], dtheta_cxst5_c[atype][btype])*rsint;

      rsint = 1.0/sin(theta6);
      df4t6 = DF4(theta6, a_cxst6[atype][btype], theta_cxst6_0[atype][btype], dtheta_cxst6_ast[atype][btype],
              b_cxst6[atype][btype], dtheta_cxst6_c[atype][btype])*rsint -
              DF4(theta6p, a_cxst6[atype][btype], theta_cxst6_0[atype][btype], dtheta_cxst6_ast[atype][btype],
              b_cxst6[atype][btype], dtheta_cxst6_c[atype][btype])*rsint;

      df5c3 = DF5(cosphi3, a_cxst3p[atype][btype], cosphi_cxst3p_ast[atype][btype], b_cxst3p[atype][btype],
              cosphi_cxst3p_c[atype][btype]);


     // force, torque and virial contribution for forces between stacking sites

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
      finc  = -df2 * f4t1 * f4t4 * f4t5 * f4t6 * f5c3 * f5c3 * rinv_st * factor_lj;

      delf[0] += delr_st[0] * finc;
      delf[1] += delr_st[1] * finc;
      delf[2] += delr_st[2] * finc;

      // theta5 force
      if (theta5 && theta5p) {

        finc   = -f2 * f4t1 * f4t4 * df4t5 * f4t6 * f5c3 * f5c3 * rinv_st * factor_lj;

        delf[0] += (delr_st_norm[0]*cost5 - az[0]) * finc;
        delf[1] += (delr_st_norm[1]*cost5 - az[1]) * finc;
        delf[2] += (delr_st_norm[2]*cost5 - az[2]) * finc;

      }

      // theta6 force
      if (theta6 && theta6p) {

        finc   = -f2 * f4t1* f4t4 * f4t5 * df4t6 * f5c3 * f5c3 * rinv_st * factor_lj;

        delf[0] += (delr_st_norm[0]*cost6 - bz[0]) * finc;
        delf[1] += (delr_st_norm[1]*cost6 - bz[1]) * finc;
        delf[2] += (delr_st_norm[2]*cost6 - bz[2]) * finc;

      }


      // cosphi3 and cosphi4 (=cosphi3) force and virial
      if (cosphi3) {

        ay[0] = ny_xtrct[a][0];
        ay[1] = ny_xtrct[a][1];
        ay[2] = ny_xtrct[a][2];

        finc   = -f2 * f4t1* f4t4 * f4t5 * f4t6 * 2.0 * f5c3 * df5c3 * factor_lj;

        gamma = d_cs - d_cst;
        gammacub = gamma * gamma * gamma;
        rinv_ss_cub = rinv_ss * rinv_ss * rinv_ss;
        aybx = MathExtra::dot3(ay,bx);
        azbx = MathExtra::dot3(az,bx);
        rax = MathExtra::dot3(delr_st_norm,ax);
        ray = MathExtra::dot3(delr_st_norm,ay);
        raz = MathExtra::dot3(delr_st_norm,az);
        rbx = MathExtra::dot3(delr_st_norm,bx);

        fac = (raz * aybx - ray * azbx);

        dcdr    = -gamma * fac * (gamma * (rax - rbx) + r_st) * rinv_ss_cub;
        dcdaxbx =  gammacub * fac * rinv_ss_cub;
        dcdaybx =  gamma * raz * rinv_ss;
        dcdazbx = -gamma * ray * rinv_ss;
        dcdrax  = -gamma*gamma * fac * r_st * rinv_ss_cub;
        dcdray  = -gamma * azbx * rinv_ss;
        dcdraz  =  gamma * aybx * rinv_ss;
        dcdrbx  =  gamma*gamma * fac * r_st * rinv_ss_cub;

        delf[0] += (delr_st_norm[0] * dcdr + ((ax[0] - delr_st_norm[0] * rax) * dcdrax +
                                              (ay[0] - delr_st_norm[0] * ray) * dcdray +
                                              (az[0] - delr_st_norm[0] * raz) * dcdraz +
                                              (bx[0] - delr_st_norm[0] * rbx) * dcdrbx) * rinv_st) * finc * factor_lj;

        delf[1] += (delr_st_norm[1] * dcdr + ((ax[1] - delr_st_norm[1] * rax) * dcdrax +
                                              (ay[1] - delr_st_norm[1] * ray) * dcdray +
                                              (az[1] - delr_st_norm[1] * raz) * dcdraz +
                                              (bx[1] - delr_st_norm[1] * rbx) * dcdrbx) * rinv_st) * finc * factor_lj;

        delf[2] += (delr_st_norm[2] * dcdr + ((ax[2] - delr_st_norm[2] * rax) * dcdrax +
                                              (ay[2] - delr_st_norm[2] * ray) * dcdray +
                                              (az[2] - delr_st_norm[2] * raz) * dcdraz +
                                              (bx[2] - delr_st_norm[2] * rbx) * dcdrbx) * rinv_st) * finc * factor_lj;

      }

      // increment forces and torques

      f[a][0] += delf[0];
      f[a][1] += delf[1];
      f[a][2] += delf[2];

      MathExtra::cross3(ra_cst,delf,delta);

      torque[a][0] += delta[0];
      torque[a][1] += delta[1];
      torque[a][2] += delta[2];

      if (newton_pair || b < nlocal) {

        f[b][0] -= delf[0];
        f[b][1] -= delf[1];
        f[b][2] -= delf[2];

        MathExtra::cross3(rb_cst,delf,deltb);

        torque[b][0] -= deltb[0];
        torque[b][1] -= deltb[1];
        torque[b][2] -= deltb[2];

      }

      // increment energy and virial
      // NOTE: The virial is calculated on the 'molecular' basis.
      // (see G. Ciccotti and J.P. Ryckaert, Comp. Phys. Rep. 4, 345-392 (1986))

      if (evflag) ev_tally_xyz(a,b,nlocal,newton_pair,evdwl,0.0,
          delf[0],delf[1],delf[2],x[a][0]-x[b][0],x[a][1]-x[b][1],x[a][2]-x[b][2]);

      // pure torques not expressible as r x f

      delta[0] = 0.0;
      delta[1] = 0.0;
      delta[2] = 0.0;
      deltb[0] = 0.0;
      deltb[1] = 0.0;
      deltb[2] = 0.0;

      // theta1 torque
      if (theta1 && theta1p) {

        tpair = -f2 * df4t1 * f4t4 * f4t5 * f4t6 * f5c3 * f5c3 * factor_lj;
        MathExtra::cross3(ax,bx,t1dir);

        delta[0] += t1dir[0]*tpair;
        delta[1] += t1dir[1]*tpair;
        delta[2] += t1dir[2]*tpair;

        deltb[0] += t1dir[0]*tpair;
        deltb[1] += t1dir[1]*tpair;
        deltb[2] += t1dir[2]*tpair;

      }

      // theta4 torque
      if (theta4) {

        tpair = -f2 * f4t1 * df4t4 * f4t5 * f4t6 * f5c3 * f5c3 * factor_lj;
        MathExtra::cross3(bz,az,t4dir);

        delta[0] += t4dir[0]*tpair;
        delta[1] += t4dir[1]*tpair;
        delta[2] += t4dir[2]*tpair;

        deltb[0] += t4dir[0]*tpair;
        deltb[1] += t4dir[1]*tpair;
        deltb[2] += t4dir[2]*tpair;

      }

      // theta5 torque
      if (theta5 && theta5p) {

        tpair = -f2 * f4t1 * f4t4 * df4t5 * f4t6 * f5c3 * f5c3 * factor_lj;
        MathExtra::cross3(delr_st_norm,az,t5dir);

        delta[0] += t5dir[0] * tpair;
        delta[1] += t5dir[1] * tpair;
        delta[2] += t5dir[2] * tpair;

      }

      // theta6 torque
      if (theta6 && theta6p) {

        tpair = -f2 * f4t1 * f4t4 * f4t5 * df4t6 * f5c3 * f5c3 * factor_lj;
        MathExtra::cross3(delr_st_norm,bz,t6dir);

        deltb[0] -= t6dir[0] * tpair;
        deltb[1] -= t6dir[1] * tpair;
        deltb[2] -= t6dir[2] * tpair;

      }

      // Full cosphi3 and cosphi4 (=cosphi3) contribution to the torque
      if (cosphi3) {

        gamma = d_cs - d_cst;
        gammacub = gamma * gamma * gamma;
        rinv_ss_cub = rinv_ss * rinv_ss * rinv_ss;
        aybx = MathExtra::dot3(ay,bx);
        azbx = MathExtra::dot3(az,bx);
        rax = MathExtra::dot3(delr_st_norm,ax);
        ray = MathExtra::dot3(delr_st_norm,ay);
        raz = MathExtra::dot3(delr_st_norm,az);
        rbx = MathExtra::dot3(delr_st_norm,bx);

        fac = (raz * aybx - ray * azbx);

        dcdr    = -gamma * fac * (gamma * (rax - rbx) + r_st) * rinv_ss_cub;
        dcdaxbx =  gammacub * fac * rinv_ss_cub;
        dcdaybx =  gamma * raz * rinv_ss;
        dcdazbx = -gamma * ray * rinv_ss;
        dcdrax  = -gamma*gamma * fac * r_st * rinv_ss_cub;
        dcdray  = -gamma * azbx * rinv_ss;
        dcdraz  =  gamma * aybx * rinv_ss;
        dcdrbx  =  gamma*gamma * fac * r_st * rinv_ss_cub;

        tpair   = -f2 * f4t1 * f4t4 * f4t5 * f4t6 * 2.0 * f5c3 * df5c3 * factor_lj;

        MathExtra::cross3(ax,bx,v1tmp);
        MathExtra::cross3(ay,bx,v2tmp);
        MathExtra::cross3(az,bx,v3tmp);

        delt[0] = (v1tmp[0] * dcdaxbx + v2tmp[0] * dcdaybx + v3tmp[0] * dcdazbx) * tpair;
        delt[1] = (v1tmp[1] * dcdaxbx + v2tmp[1] * dcdaybx + v3tmp[1] * dcdazbx) * tpair;
        delt[2] = (v1tmp[2] * dcdaxbx + v2tmp[2] * dcdaybx + v3tmp[2] * dcdazbx) * tpair;

        delta[0] += delt[0];
        delta[1] += delt[1];
        delta[2] += delt[2];
        deltb[0] += delt[0];
        deltb[1] += delt[1];
        deltb[2] += delt[2];

        MathExtra::cross3(ax,delr_st_norm,v1tmp);
        MathExtra::cross3(ay,delr_st_norm,v2tmp);
        MathExtra::cross3(az,delr_st_norm,v3tmp);

        delta[0] += (v1tmp[0] * dcdrax + v2tmp[0] * dcdray + v3tmp[0] * dcdraz) * tpair;
        delta[1] += (v1tmp[1] * dcdrax + v2tmp[1] * dcdray + v3tmp[1] * dcdraz) * tpair;
        delta[2] += (v1tmp[2] * dcdrax + v2tmp[2] * dcdray + v3tmp[2] * dcdraz) * tpair;

        MathExtra::cross3(bx,delr_st_norm,v1tmp);

        deltb[0] -= v1tmp[0] * dcdrbx * tpair;
        deltb[1] -= v1tmp[1] * dcdrbx * tpair;
        deltb[2] -= v1tmp[2] * dcdrbx * tpair;

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
      }// end early rejection criteria


    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairOxdnaCoaxstk::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(k_cxst,n+1,n+1,"pair:k_cxst");
  memory->create(cut_cxst_0,n+1,n+1,"pair:cut_cxst_0");
  memory->create(cut_cxst_c,n+1,n+1,"pair:cut_cxst_c");
  memory->create(cut_cxst_lo,n+1,n+1,"pair:cut_cxst_lo");
  memory->create(cut_cxst_hi,n+1,n+1,"pair:cut_cxst_hi");
  memory->create(cut_cxst_lc,n+1,n+1,"pair:cut_cxst_lc");
  memory->create(cut_cxst_hc,n+1,n+1,"pair:cut_cxst_hc");
  memory->create(b_cxst_lo,n+1,n+1,"pair:b_cxst_lo");
  memory->create(b_cxst_hi,n+1,n+1,"pair:b_cxst_hi");
  memory->create(cutsq_cxst_hc,n+1,n+1,"pair:cutsq_cxst_hc");

  memory->create(a_cxst1,n+1,n+1,"pair:a_cxst1");
  memory->create(theta_cxst1_0,n+1,n+1,"pair:theta_cxst1_0");
  memory->create(dtheta_cxst1_ast,n+1,n+1,"pair:dtheta_cxst1_ast");
  memory->create(b_cxst1,n+1,n+1,"pair:b_cxst1");
  memory->create(dtheta_cxst1_c,n+1,n+1,"pair:dtheta_cxst1_c");

  memory->create(a_cxst4,n+1,n+1,"pair:a_cxst4");
  memory->create(theta_cxst4_0,n+1,n+1,"pair:theta_cxst4_0");
  memory->create(dtheta_cxst4_ast,n+1,n+1,"pair:dtheta_cxst4_ast");
  memory->create(b_cxst4,n+1,n+1,"pair:b_cxst4");
  memory->create(dtheta_cxst4_c,n+1,n+1,"pair:dtheta_cxst4_c");

  memory->create(a_cxst5,n+1,n+1,"pair:a_cxst5");
  memory->create(theta_cxst5_0,n+1,n+1,"pair:theta_cxst5_0");
  memory->create(dtheta_cxst5_ast,n+1,n+1,"pair:dtheta_cxst5_ast");
  memory->create(b_cxst5,n+1,n+1,"pair:b_cxst5");
  memory->create(dtheta_cxst5_c,n+1,n+1,"pair:dtheta_cxst5_c");

  memory->create(a_cxst6,n+1,n+1,"pair:a_cxst6");
  memory->create(theta_cxst6_0,n+1,n+1,"pair:theta_cxst6_0");
  memory->create(dtheta_cxst6_ast,n+1,n+1,"pair:dtheta_cxst6_ast");
  memory->create(b_cxst6,n+1,n+1,"pair:b_cxst6");
  memory->create(dtheta_cxst6_c,n+1,n+1,"pair:dtheta_cxst6_c");

  memory->create(a_cxst3p,n+1,n+1,"pair:a_cxst3p");
  memory->create(cosphi_cxst3p_ast,n+1,n+1,"pair:cosphi_cxst3p_ast");
  memory->create(b_cxst3p,n+1,n+1,"pair:b_cxst3p");
  memory->create(cosphi_cxst3p_c,n+1,n+1,"pair:cosphi_cxst3p_c");
  memory->create(a_cxst4p,n+1,n+1,"pair:a_cxst4p");
  memory->create(cosphi_cxst4p_ast,n+1,n+1,"pair:cosphi_cxst4p_ast");
  memory->create(b_cxst4p,n+1,n+1,"pair:b_cxst4p");
  memory->create(cosphi_cxst4p_c,n+1,n+1,"pair:cosphi_cxst4p_c");

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairOxdnaCoaxstk::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairOxdnaCoaxstk::coeff(int narg, char **arg)
{
  int count;

  if (narg != 3 && narg != 23) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/coaxstk");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  // cross-stacking interaction
  count = 0;

  double k_cxst_one, cut_cxst_0_one, cut_cxst_c_one, cut_cxst_lo_one, cut_cxst_hi_one;
  double b_cxst_lo_one, b_cxst_hi_one, cut_cxst_lc_one, cut_cxst_hc_one;

  double a_cxst1_one, theta_cxst1_0_one, dtheta_cxst1_ast_one;
  double b_cxst1_one, dtheta_cxst1_c_one;

  double a_cxst4_one, theta_cxst4_0_one, dtheta_cxst4_ast_one;
  double b_cxst4_one, dtheta_cxst4_c_one;

  double a_cxst5_one, theta_cxst5_0_one, dtheta_cxst5_ast_one;
  double b_cxst5_one, dtheta_cxst5_c_one;

  double a_cxst6_one, theta_cxst6_0_one, dtheta_cxst6_ast_one;
  double b_cxst6_one, dtheta_cxst6_c_one;

  double a_cxst3p_one, cosphi_cxst3p_ast_one, b_cxst3p_one, cosphi_cxst3p_c_one;
  double a_cxst4p_one, cosphi_cxst4p_ast_one, b_cxst4p_one, cosphi_cxst4p_c_one;

  if (narg == 23) {
    k_cxst_one = utils::numeric(FLERR,arg[2],false,lmp);
    cut_cxst_0_one = utils::numeric(FLERR,arg[3],false,lmp);
    cut_cxst_c_one = utils::numeric(FLERR,arg[4],false,lmp);
    cut_cxst_lo_one = utils::numeric(FLERR,arg[5],false,lmp);
    cut_cxst_hi_one = utils::numeric(FLERR,arg[6],false,lmp);

    a_cxst1_one = utils::numeric(FLERR,arg[7],false,lmp);
    theta_cxst1_0_one = utils::numeric(FLERR,arg[8],false,lmp);
    dtheta_cxst1_ast_one = utils::numeric(FLERR,arg[9],false,lmp);

    a_cxst4_one = utils::numeric(FLERR,arg[10],false,lmp);
    theta_cxst4_0_one = utils::numeric(FLERR,arg[11],false,lmp);
    dtheta_cxst4_ast_one = utils::numeric(FLERR,arg[12],false,lmp);

    a_cxst5_one = utils::numeric(FLERR,arg[13],false,lmp);
    theta_cxst5_0_one = utils::numeric(FLERR,arg[14],false,lmp);
    dtheta_cxst5_ast_one = utils::numeric(FLERR,arg[15],false,lmp);

    a_cxst6_one = utils::numeric(FLERR,arg[16],false,lmp);
    theta_cxst6_0_one = utils::numeric(FLERR,arg[17],false,lmp);
    dtheta_cxst6_ast_one = utils::numeric(FLERR,arg[18],false,lmp);

    a_cxst3p_one = utils::numeric(FLERR,arg[19],false,lmp);
    cosphi_cxst3p_ast_one = utils::numeric(FLERR,arg[20],false,lmp);
    a_cxst4p_one = utils::numeric(FLERR,arg[21],false,lmp);
    cosphi_cxst4p_ast_one = utils::numeric(FLERR,arg[22],false,lmp);
  } else {
    if (comm->me == 0) { // read values from potential file
      PotentialFileReader reader(lmp, arg[2], "oxdna potential", " (coaxstk)");
      char * line;
      std::string iloc, jloc, potential_name;

      while ((line = reader.next_line())) {
        try {
          ValueTokenizer values(line);
          iloc = values.next_string();
          jloc = values.next_string();
          potential_name = values.next_string();
          if (iloc == arg[0] && jloc == arg[1] && potential_name == "coaxstk") {
            k_cxst_one = values.next_double();
            cut_cxst_0_one = values.next_double();
            cut_cxst_c_one = values.next_double();
            cut_cxst_lo_one = values.next_double();
            cut_cxst_hi_one = values.next_double();

            a_cxst1_one = values.next_double();
            theta_cxst1_0_one = values.next_double();
            dtheta_cxst1_ast_one = values.next_double();

            a_cxst4_one = values.next_double();
            theta_cxst4_0_one = values.next_double();
            dtheta_cxst4_ast_one = values.next_double();

            a_cxst5_one = values.next_double();
            theta_cxst5_0_one = values.next_double();
            dtheta_cxst5_ast_one = values.next_double();

            a_cxst6_one = values.next_double();
            theta_cxst6_0_one = values.next_double();
            dtheta_cxst6_ast_one = values.next_double();

            a_cxst3p_one = values.next_double();
            cosphi_cxst3p_ast_one = values.next_double();
            a_cxst4p_one = values.next_double();
            cosphi_cxst4p_ast_one = values.next_double();

            break;
          } else continue;
        } catch (std::exception &e) {
          error->one(FLERR, "Problem parsing oxDNA potential file: {}", e.what());
        }
      }
      if ((iloc != arg[0]) || (jloc != arg[1]) || (potential_name != "coaxstk"))
        error->one(FLERR, "No corresponding coaxstk potential found in file {} for pair type {} {}",
                   arg[2], arg[0], arg[1]);
    }

    MPI_Bcast(&k_cxst_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cut_cxst_0_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cut_cxst_c_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cut_cxst_lo_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cut_cxst_hi_one, 1, MPI_DOUBLE, 0, world);

    MPI_Bcast(&a_cxst1_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&theta_cxst1_0_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&dtheta_cxst1_ast_one, 1, MPI_DOUBLE, 0, world);

    MPI_Bcast(&a_cxst4_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&theta_cxst4_0_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&dtheta_cxst4_ast_one, 1, MPI_DOUBLE, 0, world);

    MPI_Bcast(&a_cxst5_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&theta_cxst5_0_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&dtheta_cxst5_ast_one, 1, MPI_DOUBLE, 0, world);

    MPI_Bcast(&a_cxst6_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&theta_cxst6_0_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&dtheta_cxst6_ast_one, 1, MPI_DOUBLE, 0, world);

    MPI_Bcast(&a_cxst3p_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cosphi_cxst3p_ast_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&a_cxst4p_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cosphi_cxst4p_ast_one, 1, MPI_DOUBLE, 0, world);
  }

  b_cxst_lo_one = 0.25 * (cut_cxst_lo_one - cut_cxst_0_one) * (cut_cxst_lo_one - cut_cxst_0_one)/
        (0.5 * (cut_cxst_lo_one - cut_cxst_0_one) * (cut_cxst_lo_one - cut_cxst_0_one) -
        k_cxst_one * 0.5 * (cut_cxst_0_one -cut_cxst_c_one) * (cut_cxst_0_one - cut_cxst_c_one)/k_cxst_one);

  cut_cxst_lc_one = cut_cxst_lo_one - 0.5 * (cut_cxst_lo_one - cut_cxst_0_one)/b_cxst_lo_one;

  b_cxst_hi_one = 0.25 * (cut_cxst_hi_one - cut_cxst_0_one) * (cut_cxst_hi_one - cut_cxst_0_one)/
        (0.5 * (cut_cxst_hi_one - cut_cxst_0_one) * (cut_cxst_hi_one - cut_cxst_0_one) -
        k_cxst_one * 0.5 * (cut_cxst_0_one -cut_cxst_c_one) * (cut_cxst_0_one - cut_cxst_c_one)/k_cxst_one);

  cut_cxst_hc_one = cut_cxst_hi_one - 0.5* (cut_cxst_hi_one - cut_cxst_0_one)/b_cxst_hi_one;


  b_cxst1_one = a_cxst1_one*a_cxst1_one*dtheta_cxst1_ast_one*dtheta_cxst1_ast_one/(1-a_cxst1_one*dtheta_cxst1_ast_one*dtheta_cxst1_ast_one);
  dtheta_cxst1_c_one = 1/(a_cxst1_one*dtheta_cxst1_ast_one);

  b_cxst4_one = a_cxst4_one*a_cxst4_one*dtheta_cxst4_ast_one*dtheta_cxst4_ast_one/(1-a_cxst4_one*dtheta_cxst4_ast_one*dtheta_cxst4_ast_one);
  dtheta_cxst4_c_one = 1/(a_cxst4_one*dtheta_cxst4_ast_one);

  b_cxst5_one = a_cxst5_one*a_cxst5_one*dtheta_cxst5_ast_one*dtheta_cxst5_ast_one/(1-a_cxst5_one*dtheta_cxst5_ast_one*dtheta_cxst5_ast_one);
  dtheta_cxst5_c_one = 1/(a_cxst5_one*dtheta_cxst5_ast_one);

  b_cxst6_one = a_cxst6_one*a_cxst6_one*dtheta_cxst6_ast_one*dtheta_cxst6_ast_one/(1-a_cxst6_one*dtheta_cxst6_ast_one*dtheta_cxst6_ast_one);
  dtheta_cxst6_c_one = 1/(a_cxst6_one*dtheta_cxst6_ast_one);

  b_cxst3p_one = a_cxst3p_one*a_cxst3p_one*cosphi_cxst3p_ast_one*cosphi_cxst3p_ast_one/(1-a_cxst3p_one*cosphi_cxst3p_ast_one*cosphi_cxst3p_ast_one);
  cosphi_cxst3p_c_one=1/(a_cxst3p_one*cosphi_cxst3p_ast_one);

  b_cxst4p_one = a_cxst4p_one*a_cxst4p_one*cosphi_cxst4p_ast_one*cosphi_cxst4p_ast_one/(1-a_cxst4p_one*cosphi_cxst4p_ast_one*cosphi_cxst4p_ast_one);
  cosphi_cxst4p_c_one=1/(a_cxst4p_one*cosphi_cxst4p_ast_one);

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {

      k_cxst[i][j] = k_cxst_one;
      cut_cxst_0[i][j] = cut_cxst_0_one;
      cut_cxst_c[i][j] = cut_cxst_c_one;
      cut_cxst_lo[i][j] = cut_cxst_lo_one;
      cut_cxst_hi[i][j] = cut_cxst_hi_one;
      cut_cxst_lc[i][j] = cut_cxst_lc_one;
      cut_cxst_hc[i][j] = cut_cxst_hc_one;
      b_cxst_lo[i][j] = b_cxst_lo_one;
      b_cxst_hi[i][j] = b_cxst_hi_one;

      a_cxst1[i][j] = a_cxst1_one;
      theta_cxst1_0[i][j] = theta_cxst1_0_one;
      dtheta_cxst1_ast[i][j] = dtheta_cxst1_ast_one;
      b_cxst1[i][j] = b_cxst1_one;
      dtheta_cxst1_c[i][j] = dtheta_cxst1_c_one;

      a_cxst4[i][j] = a_cxst4_one;
      theta_cxst4_0[i][j] = theta_cxst4_0_one;
      dtheta_cxst4_ast[i][j] = dtheta_cxst4_ast_one;
      b_cxst4[i][j] = b_cxst4_one;
      dtheta_cxst4_c[i][j] = dtheta_cxst4_c_one;

      a_cxst5[i][j] = a_cxst5_one;
      theta_cxst5_0[i][j] = theta_cxst5_0_one;
      dtheta_cxst5_ast[i][j] = dtheta_cxst5_ast_one;
      b_cxst5[i][j] = b_cxst5_one;
      dtheta_cxst5_c[i][j] = dtheta_cxst5_c_one;

      a_cxst6[i][j] = a_cxst6_one;
      theta_cxst6_0[i][j] = theta_cxst6_0_one;
      dtheta_cxst6_ast[i][j] = dtheta_cxst6_ast_one;
      b_cxst6[i][j] = b_cxst6_one;
      dtheta_cxst6_c[i][j] = dtheta_cxst6_c_one;

      a_cxst3p[i][j] = a_cxst3p_one;
      cosphi_cxst3p_ast[i][j] = cosphi_cxst3p_ast_one;
      b_cxst3p[i][j] = b_cxst3p_one;
      cosphi_cxst3p_c[i][j] = cosphi_cxst3p_c_one;

      a_cxst4p[i][j] = a_cxst4p_one;
      cosphi_cxst4p_ast[i][j] = cosphi_cxst4p_ast_one;
      b_cxst4p[i][j] = b_cxst4p_one;
      cosphi_cxst4p_c[i][j] = cosphi_cxst4p_c_one;

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/coaxstk");

}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use regular
------------------------------------------------------------------------- */

void PairOxdnaCoaxstk::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  if (id  > 0) error->all(FLERR,"Respa not supported");

}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairOxdnaCoaxstk::init_one(int i, int j)
{

  if (setflag[i][j] == 0) {
    error->all(FLERR,"Coefficient mixing not defined in oxDNA");
  }
  if (offset_flag) {
    error->all(FLERR,"Offset not supported in oxDNA");
  }

  k_cxst[j][i] = k_cxst[i][j];
  cut_cxst_0[j][i] = cut_cxst_0[i][j];
  cut_cxst_c[j][i] = cut_cxst_c[i][j];
  cut_cxst_lo[j][i] = cut_cxst_lo[i][j];
  cut_cxst_hi[j][i] = cut_cxst_hi[i][j];
  b_cxst_lo[j][i] = b_cxst_lo[i][j];
  b_cxst_hi[j][i] = b_cxst_hi[i][j];
  cut_cxst_lc[j][i] = cut_cxst_lc[i][j];
  cut_cxst_hc[j][i] = cut_cxst_hc[i][j];

  a_cxst1[j][i] = a_cxst1[i][j];
  theta_cxst1_0[j][i] = theta_cxst1_0[i][j];
  dtheta_cxst1_ast[j][i] = dtheta_cxst1_ast[i][j];
  b_cxst1[j][i] = b_cxst1[i][j];
  dtheta_cxst1_c[j][i] = dtheta_cxst1_c[i][j];

  a_cxst4[j][i] = a_cxst4[i][j];
  theta_cxst4_0[j][i] = theta_cxst4_0[i][j];
  dtheta_cxst4_ast[j][i] = dtheta_cxst4_ast[i][j];
  b_cxst4[j][i] = b_cxst4[i][j];
  dtheta_cxst4_c[j][i] = dtheta_cxst4_c[i][j];

  a_cxst5[j][i] = a_cxst5[i][j];
  theta_cxst5_0[j][i] = theta_cxst5_0[i][j];
  dtheta_cxst5_ast[j][i] = dtheta_cxst5_ast[i][j];
  b_cxst5[j][i] = b_cxst5[i][j];
  dtheta_cxst5_c[j][i] = dtheta_cxst5_c[i][j];

  a_cxst6[j][i] = a_cxst6[i][j];
  theta_cxst6_0[j][i] = theta_cxst6_0[i][j];
  dtheta_cxst6_ast[j][i] = dtheta_cxst6_ast[i][j];
  b_cxst6[j][i] = b_cxst6[i][j];
  dtheta_cxst6_c[j][i] = dtheta_cxst6_c[i][j];

  a_cxst3p[j][i] = a_cxst3p[i][j];
  cosphi_cxst3p_ast[j][i] = cosphi_cxst3p_ast[i][j];
  b_cxst3p[j][i] = b_cxst3p[i][j];
  cosphi_cxst3p_c[j][i] = cosphi_cxst3p_c[i][j];

  a_cxst4p[j][i] = a_cxst4p[i][j];
  cosphi_cxst4p_ast[j][i] = cosphi_cxst4p_ast[i][j];
  b_cxst4p[j][i] = b_cxst4p[i][j];
  cosphi_cxst4p_c[j][i] = cosphi_cxst4p_c[i][j];

  cutsq_cxst_hc[i][j] = cut_cxst_hc[i][j]*cut_cxst_hc[i][j];
  cutsq_cxst_hc[j][i] = cutsq_cxst_hc[i][j];

  // set the master list distance cutoff
  return cut_cxst_hc[i][j];

}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdnaCoaxstk::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {

        fwrite(&k_cxst[i][j],sizeof(double),1,fp);
        fwrite(&cut_cxst_0[i][j],sizeof(double),1,fp);
        fwrite(&cut_cxst_c[i][j],sizeof(double),1,fp);
        fwrite(&cut_cxst_lo[i][j],sizeof(double),1,fp);
        fwrite(&cut_cxst_hi[i][j],sizeof(double),1,fp);
        fwrite(&cut_cxst_lc[i][j],sizeof(double),1,fp);
        fwrite(&cut_cxst_hc[i][j],sizeof(double),1,fp);
        fwrite(&b_cxst_lo[i][j],sizeof(double),1,fp);
        fwrite(&b_cxst_hi[i][j],sizeof(double),1,fp);

        fwrite(&a_cxst1[i][j],sizeof(double),1,fp);
        fwrite(&theta_cxst1_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_cxst1_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_cxst1[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_cxst1_c[i][j],sizeof(double),1,fp);

        fwrite(&a_cxst4[i][j],sizeof(double),1,fp);
        fwrite(&theta_cxst4_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_cxst4_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_cxst4[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_cxst4_c[i][j],sizeof(double),1,fp);

        fwrite(&a_cxst5[i][j],sizeof(double),1,fp);
        fwrite(&theta_cxst5_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_cxst5_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_cxst5[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_cxst5_c[i][j],sizeof(double),1,fp);

        fwrite(&a_cxst6[i][j],sizeof(double),1,fp);
        fwrite(&theta_cxst6_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_cxst6_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_cxst6[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_cxst6_c[i][j],sizeof(double),1,fp);

        fwrite(&a_cxst3p[i][j],sizeof(double),1,fp);
        fwrite(&cosphi_cxst3p_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_cxst3p[i][j],sizeof(double),1,fp);
        fwrite(&cosphi_cxst3p_c[i][j],sizeof(double),1,fp);
        fwrite(&a_cxst4p[i][j],sizeof(double),1,fp);
        fwrite(&cosphi_cxst4p_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_cxst4p[i][j],sizeof(double),1,fp);
        fwrite(&cosphi_cxst4p_c[i][j],sizeof(double),1,fp);

    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdnaCoaxstk::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {

          utils::sfread(FLERR,&k_cxst[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_cxst_0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_cxst_c[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_cxst_lo[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_cxst_hi[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_cxst_lc[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_cxst_hc[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_cxst_lo[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_cxst_hi[i][j],sizeof(double),1,fp,nullptr,error);

          utils::sfread(FLERR,&a_cxst1[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&theta_cxst1_0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_cxst1_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_cxst1[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_cxst1_c[i][j],sizeof(double),1,fp,nullptr,error);

          utils::sfread(FLERR,&a_cxst4[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&theta_cxst4_0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_cxst4_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_cxst4[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_cxst4_c[i][j],sizeof(double),1,fp,nullptr,error);

          utils::sfread(FLERR,&a_cxst5[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&theta_cxst5_0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_cxst5_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_cxst5[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_cxst5_c[i][j],sizeof(double),1,fp,nullptr,error);

          utils::sfread(FLERR,&a_cxst6[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&theta_cxst6_0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_cxst6_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_cxst6[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_cxst6_c[i][j],sizeof(double),1,fp,nullptr,error);

          utils::sfread(FLERR,&a_cxst3p[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cosphi_cxst3p_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_cxst3p[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cosphi_cxst3p_c[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&a_cxst4p[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cosphi_cxst4p_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_cxst4p[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cosphi_cxst4p_c[i][j],sizeof(double),1,fp,nullptr,error);

        }

        MPI_Bcast(&k_cxst[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_cxst_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_cxst_c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_cxst_lo[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_cxst_hi[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_cxst_lc[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_cxst_hc[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_cxst_lo[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_cxst_hi[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_cxst1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_cxst1_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_cxst1_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_cxst1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_cxst1_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_cxst4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_cxst4_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_cxst4_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_cxst4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_cxst4_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_cxst5[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_cxst5_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_cxst5_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_cxst5[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_cxst5_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_cxst6[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_cxst6_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_cxst6_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_cxst6[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_cxst6_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_cxst3p[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cosphi_cxst3p_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_cxst3p[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cosphi_cxst3p_c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a_cxst4p[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cosphi_cxst4p_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_cxst4p[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cosphi_cxst4p_c[i][j],1,MPI_DOUBLE,0,world);

      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdnaCoaxstk::write_restart_settings(FILE *fp)
{
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdnaCoaxstk::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&tail_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairOxdnaCoaxstk::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d\
         %g %g %g %g %g\
         %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g\
         %g %g %g %g\
         \n",i,
        k_cxst[i][i],cut_cxst_0[i][i],cut_cxst_c[i][i],cut_cxst_lo[i][i],cut_cxst_hi[i][i],
        cut_cxst_lc[i][i],cut_cxst_hc[i][i],b_cxst_lo[i][i],b_cxst_hi[i][i],
        a_cxst1[i][i],theta_cxst1_0[i][i],dtheta_cxst1_ast[i][i],b_cxst1[i][i],dtheta_cxst1_c[i][i],
        a_cxst4[i][i],theta_cxst4_0[i][i],dtheta_cxst4_ast[i][i],b_cxst4[i][i],dtheta_cxst4_c[i][i],
        a_cxst5[i][i],theta_cxst5_0[i][i],dtheta_cxst5_ast[i][i],b_cxst5[i][i],dtheta_cxst5_c[i][i],
        a_cxst6[i][i],theta_cxst6_0[i][i],dtheta_cxst6_ast[i][i],b_cxst6[i][i],dtheta_cxst6_c[i][i],
        a_cxst3p[i][i],cosphi_cxst3p_ast[i][i],b_cxst3p[i][i], cosphi_cxst3p_c[i][i],
        a_cxst4p[i][i],cosphi_cxst4p_ast[i][i],b_cxst4p[i][i], cosphi_cxst4p_c[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairOxdnaCoaxstk::write_data_all(FILE *fp)
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
         %g %g %g %g\
         %g %g %g %g\
         \n",i,j,
        k_cxst[i][j],cut_cxst_0[i][j],cut_cxst_c[i][j],cut_cxst_lo[i][j],cut_cxst_hi[i][j],
        cut_cxst_lc[i][j],cut_cxst_hc[i][j],b_cxst_lo[i][j],b_cxst_hi[i][j],
        a_cxst1[i][j],theta_cxst1_0[i][j],dtheta_cxst1_ast[i][j],b_cxst1[i][j],dtheta_cxst1_c[i][j],
        a_cxst4[i][j],theta_cxst4_0[i][j],dtheta_cxst4_ast[i][j],b_cxst4[i][j],dtheta_cxst4_c[i][j],
        a_cxst5[i][j],theta_cxst5_0[i][j],dtheta_cxst5_ast[i][j],b_cxst5[i][j],dtheta_cxst5_c[i][j],
        a_cxst6[i][j],theta_cxst6_0[i][j],dtheta_cxst6_ast[i][j],b_cxst6[i][j],dtheta_cxst6_c[i][j],
        a_cxst3p[i][j],cosphi_cxst3p_ast[i][j],b_cxst3p[i][j],cosphi_cxst3p_c[i][j],
        a_cxst4p[i][j],cosphi_cxst4p_ast[i][j],b_cxst4p[i][j],cosphi_cxst4p_c[i][j]);

}

/* ---------------------------------------------------------------------- */

void *PairOxdnaCoaxstk::extract(const char *str, int &dim)
{
  dim = 2;

  if (strcmp(str,"k_cxst") == 0) return (void *) k_cxst;
  if (strcmp(str,"cut_cxst_0") == 0) return (void *) cut_cxst_0;
  if (strcmp(str,"cut_cxst_c") == 0) return (void *) cut_cxst_c;
  if (strcmp(str,"cut_cxst_lo") == 0) return (void *) cut_cxst_lo;
  if (strcmp(str,"cut_cxst_hi") == 0) return (void *) cut_cxst_hi;
  if (strcmp(str,"cut_cxst_lc") == 0) return (void *) cut_cxst_lc;
  if (strcmp(str,"cut_cxst_hc") == 0) return (void *) cut_cxst_hc;
  if (strcmp(str,"b_cxst_lo") == 0) return (void *) b_cxst_lo;
  if (strcmp(str,"b_cxst_hi") == 0) return (void *) b_cxst_hi;

  if (strcmp(str,"a_cxst1") == 0) return (void *) a_cxst1;
  if (strcmp(str,"theta_cxst1_0") == 0) return (void *) theta_cxst1_0;
  if (strcmp(str,"dtheta_cxst1_ast") == 0) return (void *) dtheta_cxst1_ast;
  if (strcmp(str,"b_cxst1") == 0) return (void *) b_cxst1;
  if (strcmp(str,"dtheta_cxst1_c") == 0) return (void *) dtheta_cxst1_c;

  if (strcmp(str,"a_cxst4") == 0) return (void *) a_cxst4;
  if (strcmp(str,"theta_cxst4_0") == 0) return (void *) theta_cxst4_0;
  if (strcmp(str,"dtheta_cxst4_ast") == 0) return (void *) dtheta_cxst4_ast;
  if (strcmp(str,"b_cxst4") == 0) return (void *) b_cxst4;
  if (strcmp(str,"dtheta_cxst4_c") == 0) return (void *) dtheta_cxst4_c;

  if (strcmp(str,"a_cxst5") == 0) return (void *) a_cxst5;
  if (strcmp(str,"theta_cxst5_0") == 0) return (void *) theta_cxst5_0;
  if (strcmp(str,"dtheta_cxst5_ast") == 0) return (void *) dtheta_cxst5_ast;
  if (strcmp(str,"b_cxst5") == 0) return (void *) b_cxst5;
  if (strcmp(str,"dtheta_cxst5_c") == 0) return (void *) dtheta_cxst5_c;

  if (strcmp(str,"a_cxst6") == 0) return (void *) a_cxst6;
  if (strcmp(str,"theta_cxst6_0") == 0) return (void *) theta_cxst6_0;
  if (strcmp(str,"dtheta_cxst6_ast") == 0) return (void *) dtheta_cxst6_ast;
  if (strcmp(str,"b_cxst6") == 0) return (void *) b_cxst6;
  if (strcmp(str,"dtheta_cxst6_c") == 0) return (void *) dtheta_cxst6_c;

  if (strcmp(str,"a_cxst3p") == 0) return (void *) a_cxst3p;
  if (strcmp(str,"cosphi_cxst3p_ast") == 0) return (void *) cosphi_cxst3p_ast;
  if (strcmp(str,"b_cxst3p") == 0) return (void *) b_cxst3p;
  if (strcmp(str,"cosphi_cxst3p_c") == 0) return (void *) cosphi_cxst3p_c;

  if (strcmp(str,"a_cxst4p") == 0) return (void *) a_cxst4p;
  if (strcmp(str,"cosphi_cxst4p_ast") == 0) return (void *) cosphi_cxst4p_ast;
  if (strcmp(str,"b_cxst4p") == 0) return (void *) b_cxst4p;
  if (strcmp(str,"cosphi_cxst4p_c") == 0) return (void *) cosphi_cxst4p_c;

  return nullptr;
}
