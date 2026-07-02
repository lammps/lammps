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

#include "pair_oxdna3_xstk.h"
#include "constants_oxdna.h"
#include "mf_oxdna.h"
#include "nucleotide_oxdna.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix_oxdna_lrf.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"
#include "math_special.h"

#include <cmath>
#include <cstring>
#include <cassert>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;
using namespace MFOxdna;

/* ---------------------------------------------------------------------- */

PairOxdna3Xstk::PairOxdna3Xstk(LAMMPS *lmp) :
    Pair(lmp), k_xst(nullptr), cut_xst_0_33(nullptr), cut_xst_c_33(nullptr),
    cut_xst_lo_33(nullptr), cut_xst_hi_33(nullptr), cut_xst_lc_33(nullptr), cut_xst_hc_33(nullptr),
    cutsq_xst_hc_33(nullptr), cut_xst_0_55(nullptr), cut_xst_c_55(nullptr), cut_xst_lo_55(nullptr),
    cut_xst_hi_55(nullptr), cut_xst_lc_55(nullptr), cut_xst_hc_55(nullptr),
    cutsq_xst_hc_55(nullptr), b_xst_lo(nullptr), b_xst_hi(nullptr), a_xst1(nullptr),
    theta_xst1_0(nullptr), dtheta_xst1_ast(nullptr), b_xst1(nullptr), dtheta_xst1_c(nullptr),
    a_xst2(nullptr), theta_xst2_0(nullptr), dtheta_xst2_ast(nullptr), b_xst2(nullptr),
    dtheta_xst2_c(nullptr), a_xst3(nullptr), theta_xst3_0(nullptr), dtheta_xst3_ast(nullptr),
    b_xst3(nullptr), dtheta_xst3_c(nullptr), a_xst4_33(nullptr), theta_xst4_0_33(nullptr),
    dtheta_xst4_ast_33(nullptr), b_xst4_33(nullptr), dtheta_xst4_c_33(nullptr), a_xst4_55(nullptr),
    theta_xst4_0_55(nullptr), dtheta_xst4_ast_55(nullptr), b_xst4_55(nullptr),
    dtheta_xst4_c_55(nullptr), a_xst7(nullptr), theta_xst7_0_33(nullptr), theta_xst7_0_55(nullptr),
    dtheta_xst7_ast(nullptr), b_xst7(nullptr), dtheta_xst7_c(nullptr), a_xst8(nullptr),
    theta_xst8_0_33(nullptr), theta_xst8_0_55(nullptr), dtheta_xst8_ast(nullptr), b_xst8(nullptr),
    dtheta_xst8_c(nullptr), nxyz_xtrct(nullptr), fix_lrf(nullptr)
{
  single_enable = 0;
  writedata = 0;
  trim_flag = 0;
}

/* ---------------------------------------------------------------------- */

PairOxdna3Xstk::~PairOxdna3Xstk()
{
  if (allocated) {

    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(k_xst);

    memory->destroy(cut_xst_0_33);
    memory->destroy(cut_xst_c_33);
    memory->destroy(cut_xst_lo_33);
    memory->destroy(cut_xst_hi_33);
    memory->destroy(cut_xst_lc_33);
    memory->destroy(cut_xst_hc_33);
    memory->destroy(cutsq_xst_hc_33);

    memory->destroy(cut_xst_0_55);
    memory->destroy(cut_xst_c_55);
    memory->destroy(cut_xst_lo_55);
    memory->destroy(cut_xst_hi_55);
    memory->destroy(cut_xst_lc_55);
    memory->destroy(cut_xst_hc_55);
    memory->destroy(cutsq_xst_hc_55);

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

    memory->destroy(a_xst4_33);
    memory->destroy(theta_xst4_0_33);
    memory->destroy(dtheta_xst4_ast_33);
    memory->destroy(b_xst4_33);
    memory->destroy(dtheta_xst4_c_33);

    memory->destroy(a_xst4_55);
    memory->destroy(theta_xst4_0_55);
    memory->destroy(dtheta_xst4_ast_55);
    memory->destroy(b_xst4_55);
    memory->destroy(dtheta_xst4_c_55);

    memory->destroy(a_xst7);
    memory->destroy(theta_xst7_0_33);
    memory->destroy(theta_xst7_0_55);
    memory->destroy(dtheta_xst7_ast);
    memory->destroy(b_xst7);
    memory->destroy(dtheta_xst7_c);

    memory->destroy(a_xst8);
    memory->destroy(theta_xst8_0_33);
    memory->destroy(theta_xst8_0_55);
    memory->destroy(dtheta_xst8_ast);
    memory->destroy(b_xst8);
    memory->destroy(dtheta_xst8_c);

  }
}

/* --------------------------------------------------------------
   compute vector COM-hydrogen bonding interaction site in oxDNA3
   A=1, C=2, G=3, T=0
----------------------------------------------------------------- */
inline void PairOxdna3Xstk::compute_base_site(int type, double e1[3],
  double /*e2*/[3], double /*e3*/[3], double rbs[3]) const
{
  NucleotideOxdna3 oxdna3;
  switch (type) {
    case 0:
      oxdna3.base_site<0>(e1, nullptr, nullptr, rbs);
      break;
    case 1:
      oxdna3.base_site<1>(e1, nullptr, nullptr, rbs);
      break;
    case 2:
      oxdna3.base_site<2>(e1, nullptr, nullptr, rbs);
      break;
    case 3:
      oxdna3.base_site<3>(e1, nullptr, nullptr, rbs);
      break;
  }
}

/* ----------------------------------------------------------------------
   compute function for oxDNA pair interactions
   hb=hydrogen bonding site

   NOTE: The cross-stacking interaction takes place between hb sites
------------------------------------------------------------------------- */

void PairOxdna3Xstk::compute(int eflag, int vflag)
{
  double delf[3],delta[3],deltb[3]; // force, torque increment;
  double evdwl,finc,tpair,factor_lj;
  double delr_bsbs[3],delr_bsbs_norm[3],rsq_bsbs,r_bsbs,rinv_bsbs;
  double theta1,t1dir[3],cost1;
  double theta2,t2dir[3],cost2;
  double theta3,t3dir[3],cost3;
  double theta4,t4dir[3],cost4;
  double theta7,t7dir[3],cost7;
  double theta8,t8dir[3],cost8;

  // vectors COM-h-bonding site in lab frame
  double ra_cbs[3],rb_cbs[3];
  // Cartesian unit vectors in lab frame
  double ax[3],ay[3],az[3];
  double bx[3],by[3],bz[3];

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  int *type = atom->type;

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int *alist,*blist,*numneigh,**firstneigh;
  double *special_lj = force->special_lj;

  tagint *id3p = atom->id3p;
  tagint *id5p = atom->id5p;

  int a,b,ia,ib,anum,bnum;
  int a3ptype,atype,a5ptype,b3ptype,btype,b5ptype;

  double f2_33,f2_55,f4t1,f4t2,f4t3,f4t4_33,f4t4_55,f4t7_33,f4t7_55,f4t8_33,f4t8_55;
  double df2_33,df2_55,df4t1,df4t2,df4t3,df4t4_33,df4t4_55,df4t7_33,df4t7_55,df4t8_33,df4t8_55;
  double rsint;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  anum = list->inum;
  alist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // nxyz_xtrct = extracted local unit vectors in lab frame from fix OXDNA/LRF
  nxyz_xtrct = fix_lrf->array_atom;

  // loop over pair interaction neighbors of my atoms

  for (ia = 0; ia < anum; ia++) {

    a = alist[ia];
    atype = type[a];

    if (id3p[a] != -1) {
      a3ptype = type[atom->map(id3p[a])];
    }
    else a3ptype = 0;

    if (id5p[a] != -1) {
      a5ptype = type[atom->map(id5p[a])];
    }
    else a5ptype = 0;

    ax[0] = nxyz_xtrct[a][0];
    ax[1] = nxyz_xtrct[a][1];
    ax[2] = nxyz_xtrct[a][2];

    // vector COM - base site a
    compute_base_site(atype%4,ax,ay,az,ra_cbs);

    blist = firstneigh[a];
    bnum = numneigh[a];

    for (ib = 0; ib < bnum; ib++) {

      b = blist[ib];
      factor_lj = special_lj[sbmask(b)]; // = 0 for nearest neighbors
      b &= NEIGHMASK;

      btype = type[b];

      if (id3p[b] != -1) {
        b3ptype = type[atom->map(id3p[b])];
      }
      else b3ptype = 0;

      if (id5p[b] != -1) {
        b5ptype = type[atom->map(id5p[b])];
      }
      else b5ptype = 0;

      bx[0] = nxyz_xtrct[b][0];
      bx[1] = nxyz_xtrct[b][1];
      bx[2] = nxyz_xtrct[b][2];

      // vector COM - base site b
      compute_base_site(btype%4,bx,by,bz,rb_cbs);

      // vector h-bonding site b to a
      delr_bsbs[0] = x[a][0] + ra_cbs[0] - x[b][0] - rb_cbs[0];
      delr_bsbs[1] = x[a][1] + ra_cbs[1] - x[b][1] - rb_cbs[1];
      delr_bsbs[2] = x[a][2] + ra_cbs[2] - x[b][2] - rb_cbs[2];

      rsq_bsbs = delr_bsbs[0]*delr_bsbs[0] + delr_bsbs[1]*delr_bsbs[1] + delr_bsbs[2]*delr_bsbs[2];
      r_bsbs = sqrt(rsq_bsbs);
      rinv_bsbs = 1.0/r_bsbs;

      delr_bsbs_norm[0] = delr_bsbs[0] * rinv_bsbs;
      delr_bsbs_norm[1] = delr_bsbs[1] * rinv_bsbs;
      delr_bsbs_norm[2] = delr_bsbs[2] * rinv_bsbs;

      f2_33 = F2(r_bsbs, k_xst[atype][btype], cut_xst_0_33[a3ptype][atype][btype][b3ptype],
              cut_xst_lc_33[a3ptype][atype][btype][b3ptype], cut_xst_hc_33[a3ptype][atype][btype][b3ptype],
              cut_xst_lo_33[a3ptype][atype][btype][b3ptype], cut_xst_hi_33[a3ptype][atype][btype][b3ptype],
              b_xst_lo[atype][btype], b_xst_hi[atype][btype], cut_xst_c_33[a3ptype][atype][btype][b3ptype]);

      f2_55 = F2(r_bsbs, k_xst[atype][btype], cut_xst_0_55[a5ptype][atype][btype][b5ptype],
              cut_xst_lc_55[a5ptype][atype][btype][b5ptype], cut_xst_hc_55[a5ptype][atype][btype][b5ptype],
              cut_xst_lo_55[a5ptype][atype][btype][b5ptype], cut_xst_hi_55[a5ptype][atype][btype][b5ptype],
              b_xst_lo[atype][btype], b_xst_hi[atype][btype], cut_xst_c_55[a5ptype][atype][btype][b5ptype]);

      // early rejection criterium
      if ((f2_33 != 0.0) || (f2_55 != 0.0)) {

      cost1 = -1.0*MathExtra::dot3(ax,bx);
      if (cost1 >  1.0) cost1 =  1.0;
      if (cost1 < -1.0) cost1 = -1.0;
      theta1 = acos(cost1);

      f4t1 = F4(theta1, a_xst1[atype][btype], theta_xst1_0[atype][btype], dtheta_xst1_ast[atype][btype],
             b_xst1[atype][btype], dtheta_xst1_c[atype][btype]);

      // early rejection criterium
      if (f4t1 != 0.0) {

      cost2 = -1.0*MathExtra::dot3(ax,delr_bsbs_norm);
      if (cost2 >  1.0) cost2 =  1.0;
      if (cost2 < -1.0) cost2 = -1.0;
      theta2 = acos(cost2);

      f4t2 = F4(theta2, a_xst2[atype][btype], theta_xst2_0[atype][btype], dtheta_xst2_ast[atype][btype],
             b_xst2[atype][btype], dtheta_xst2_c[atype][btype]);

      // early rejection criterium
      if (f4t2 != 0.0) {

      cost3 = MathExtra::dot3(bx,delr_bsbs_norm);
      if (cost3 >  1.0) cost3 =  1.0;
      if (cost3 < -1.0) cost3 = -1.0;
      theta3 = acos(cost3);

      f4t3 = F4(theta3, a_xst3[atype][btype], theta_xst3_0[atype][btype], dtheta_xst3_ast[atype][btype],
             b_xst3[atype][btype], dtheta_xst3_c[atype][btype]);

      // early rejection criterium
      if (f4t3 != 0.0) {

      az[0] = nxyz_xtrct[a][6];
      az[1] = nxyz_xtrct[a][7];
      az[2] = nxyz_xtrct[a][8];
      bz[0] = nxyz_xtrct[b][6];
      bz[1] = nxyz_xtrct[b][7];
      bz[2] = nxyz_xtrct[b][8];

      cost4 = MathExtra::dot3(az,bz);
      if (cost4 >  1.0) cost4 =  1.0;
      if (cost4 < -1.0) cost4 = -1.0;
      theta4 = acos(cost4);

      f4t4_33 = F4(theta4, a_xst4_33[a3ptype][atype][btype][b3ptype], theta_xst4_0_33[a3ptype][atype][btype][b3ptype],
                  dtheta_xst4_ast_33[a3ptype][atype][btype][b3ptype], b_xst4_33[a3ptype][atype][btype][b3ptype],
                  dtheta_xst4_c_33[a3ptype][atype][btype][b3ptype]);

      f4t4_55 = F4(theta4, a_xst4_55[a5ptype][atype][btype][b5ptype], theta_xst4_0_55[a5ptype][atype][btype][b5ptype],
                  dtheta_xst4_ast_55[a5ptype][atype][btype][b5ptype], b_xst4_55[a5ptype][atype][btype][b5ptype],
                  dtheta_xst4_c_55[a5ptype][atype][btype][b5ptype]);

      // early rejection criterium
      if ((f4t4_33 != 0.0) || (f4t4_55 != 0.0)) {

      cost7 = -1.0*MathExtra::dot3(az,delr_bsbs_norm);
      if (cost7 >  1.0) cost7 =  1.0;
      if (cost7 < -1.0) cost7 = -1.0;
      theta7 = acos(cost7);

      f4t7_33 = F4(theta7, a_xst7[atype][btype], theta_xst7_0_33[atype][btype], dtheta_xst7_ast[atype][btype],
                 b_xst7[atype][btype], dtheta_xst7_c[atype][btype]);

      f4t7_55 = F4(theta7, a_xst7[atype][btype], theta_xst7_0_55[atype][btype], dtheta_xst7_ast[atype][btype],
                 b_xst7[atype][btype], dtheta_xst7_c[atype][btype]);

      // early rejection criterium
      if ((f4t7_33 != 0.0) || (f4t7_55 != 0.0)) {

      cost8 = MathExtra::dot3(bz,delr_bsbs_norm);
      if (cost8 >  1.0) cost8 =  1.0;
      if (cost8 < -1.0) cost8 = -1.0;
      theta8 = acos(cost8);

      f4t8_33 = F4(theta8, a_xst8[atype][btype], theta_xst8_0_33[atype][btype], dtheta_xst8_ast[atype][btype],
                 b_xst8[atype][btype], dtheta_xst8_c[atype][btype]);

      f4t8_55 = F4(theta8, a_xst8[atype][btype], theta_xst8_0_55[atype][btype], dtheta_xst8_ast[atype][btype],
                 b_xst8[atype][btype], dtheta_xst8_c[atype][btype]);

      evdwl = f4t1 * f4t2 * f4t3 * (f2_33 * f4t4_33 * f4t7_33 * f4t8_33 + f2_55 * f4t4_55 * f4t7_55 * f4t8_55) * factor_lj;

      // early rejection criterium
      if (evdwl != 0.0) {

      df2_33 = DF2(r_bsbs, k_xst[atype][btype], cut_xst_0_33[a3ptype][atype][btype][b3ptype],
                 cut_xst_lc_33[a3ptype][atype][btype][b3ptype], cut_xst_hc_33[a3ptype][atype][btype][b3ptype],
                 cut_xst_lo_33[a3ptype][atype][btype][b3ptype], cut_xst_hi_33[a3ptype][atype][btype][b3ptype],
                 b_xst_lo[atype][btype], b_xst_hi[atype][btype]);

      df2_55 = DF2(r_bsbs, k_xst[atype][btype], cut_xst_0_55[a5ptype][atype][btype][b5ptype],
                 cut_xst_lc_55[a5ptype][atype][btype][b5ptype], cut_xst_hc_55[a5ptype][atype][btype][b5ptype],
                 cut_xst_lo_55[a5ptype][atype][btype][b5ptype], cut_xst_hi_55[a5ptype][atype][btype][b5ptype],
                 b_xst_lo[atype][btype], b_xst_hi[atype][btype]);

      df4t1 = DF4(theta1, a_xst1[atype][btype], theta_xst1_0[atype][btype], dtheta_xst1_ast[atype][btype],
                b_xst1[atype][btype], dtheta_xst1_c[atype][btype])/sin(theta1);

      df4t2 = DF4(theta2, a_xst2[atype][btype], theta_xst2_0[atype][btype], dtheta_xst2_ast[atype][btype],
                b_xst2[atype][btype], dtheta_xst2_c[atype][btype])/sin(theta2);

      df4t3 = DF4(theta3, a_xst3[atype][btype], theta_xst3_0[atype][btype], dtheta_xst3_ast[atype][btype],
                b_xst3[atype][btype], dtheta_xst3_c[atype][btype])/sin(theta3);

      rsint = 1.0/sin(theta4);
      df4t4_33 = DF4(theta4, a_xst4_33[a3ptype][atype][btype][b3ptype], theta_xst4_0_33[a3ptype][atype][btype][b3ptype],
                  dtheta_xst4_ast_33[a3ptype][atype][btype][b3ptype], b_xst4_33[a3ptype][atype][btype][b3ptype],
                  dtheta_xst4_c_33[a3ptype][atype][btype][b3ptype])*rsint;

      df4t4_55 = DF4(theta4, a_xst4_55[a5ptype][atype][btype][b5ptype], theta_xst4_0_55[a5ptype][atype][btype][b5ptype],
                  dtheta_xst4_ast_55[a5ptype][atype][btype][b5ptype], b_xst4_55[a5ptype][atype][btype][b5ptype],
                  dtheta_xst4_c_55[a5ptype][atype][btype][b5ptype])*rsint;

      rsint = 1.0/sin(theta7);
      df4t7_33 = DF4(theta7, a_xst7[atype][btype], theta_xst7_0_33[atype][btype], dtheta_xst7_ast[atype][btype],
                  b_xst7[atype][btype], dtheta_xst7_c[atype][btype])*rsint;

      df4t7_55 = DF4(theta7, a_xst7[atype][btype], theta_xst7_0_55[atype][btype], dtheta_xst7_ast[atype][btype],
                  b_xst7[atype][btype], dtheta_xst7_c[atype][btype])*rsint;

      rsint = 1.0/sin(theta8);
      df4t8_33 = DF4(theta8, a_xst8[atype][btype], theta_xst8_0_33[atype][btype], dtheta_xst8_ast[atype][btype],
                  b_xst8[atype][btype], dtheta_xst8_c[atype][btype])*rsint;

      df4t8_55 = DF4(theta8, a_xst8[atype][btype], theta_xst8_0_55[atype][btype], dtheta_xst8_ast[atype][btype],
                  b_xst8[atype][btype], dtheta_xst8_c[atype][btype])*rsint;

      // force, torque and virial contribution for forces between h-bonding sites

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
      finc  = -f4t1 * f4t2 * f4t3 * (df2_33 * f4t4_33 * f4t7_33 * f4t8_33 + df2_55 * f4t4_55 * f4t7_55 * f4t8_55) * rinv_bsbs *factor_lj;

      delf[0] += delr_bsbs[0] * finc;
      delf[1] += delr_bsbs[1] * finc;
      delf[2] += delr_bsbs[2] * finc;

      // theta2 force
      if (theta2 != 0.0) {

        finc  = -f4t1 * df4t2 * f4t3 * (f2_33 * f4t4_33 * f4t7_33 * f4t8_33 + f2_55 * f4t4_55 * f4t7_55 * f4t8_55) * rinv_bsbs * factor_lj;

        delf[0] += (delr_bsbs_norm[0]*cost2 + ax[0]) * finc;
        delf[1] += (delr_bsbs_norm[1]*cost2 + ax[1]) * finc;
        delf[2] += (delr_bsbs_norm[2]*cost2 + ax[2]) * finc;

      }

      // theta3 force
      if (theta3 != 0.0) {

        finc  = -f4t1 * f4t2 * df4t3 * (f2_33 * f4t4_33 * f4t7_33 * f4t8_33 + f2_55 * f4t4_55 * f4t7_55 * f4t8_55) * rinv_bsbs * factor_lj;

        delf[0] += (delr_bsbs_norm[0]*cost3 - bx[0]) * finc;
        delf[1] += (delr_bsbs_norm[1]*cost3 - bx[1]) * finc;
        delf[2] += (delr_bsbs_norm[2]*cost3 - bx[2]) * finc;

      }

      // theta7 force
      if (theta7 != 0.0) {

        finc  = -f4t1 * f4t2 * f4t3 * (f2_33 * f4t4_33 * df4t7_33 * f4t8_33 + f2_55 * f4t4_55 * df4t7_55 * f4t8_55) * rinv_bsbs * factor_lj;

        delf[0] += (delr_bsbs_norm[0]*cost7 + az[0]) * finc;
        delf[1] += (delr_bsbs_norm[1]*cost7 + az[1]) * finc;
        delf[2] += (delr_bsbs_norm[2]*cost7 + az[2]) * finc;

      }

      // theta8 force
      if (theta8 != 0.0) {

        finc  = -f4t1 * f4t2 * f4t3 * (f2_33 * f4t4_33 * f4t7_33 * df4t8_33 + f2_55 * f4t4_55 * f4t7_55 * df4t8_55) * rinv_bsbs * factor_lj;

        delf[0] += (delr_bsbs_norm[0]*cost8 - bz[0]) * finc;
        delf[1] += (delr_bsbs_norm[1]*cost8 - bz[1]) * finc;
        delf[2] += (delr_bsbs_norm[2]*cost8 - bz[2]) * finc;

      }

      // increment forces and torques

      f[a][0] += delf[0];
      f[a][1] += delf[1];
      f[a][2] += delf[2];

      MathExtra::cross3(ra_cbs,delf,delta);

      torque[a][0] += delta[0];
      torque[a][1] += delta[1];
      torque[a][2] += delta[2];

      if (newton_pair || b < nlocal) {

        f[b][0] -= delf[0];
        f[b][1] -= delf[1];
        f[b][2] -= delf[2];


        MathExtra::cross3(rb_cbs,delf,deltb);

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
      if (theta1 != 0.0) {

        tpair = -df4t1 * f4t2 * f4t3 * (f2_33 * f4t4_33 * f4t7_33 * f4t8_33 + f2_55 * f4t4_55 * f4t7_55 * f4t8_55) * factor_lj;
        MathExtra::cross3(ax,bx,t1dir);

        delta[0] += t1dir[0]*tpair;
        delta[1] += t1dir[1]*tpair;
        delta[2] += t1dir[2]*tpair;

        deltb[0] += t1dir[0]*tpair;
        deltb[1] += t1dir[1]*tpair;
        deltb[2] += t1dir[2]*tpair;

      }

      // theta2 torque
      if (theta2 != 0.0) {

        tpair = -f4t1 * df4t2 * f4t3 * (f2_33 * f4t4_33 * f4t7_33 * f4t8_33 + f2_55 * f4t4_55 * f4t7_55 * f4t8_55) * factor_lj;
        MathExtra::cross3(ax,delr_bsbs_norm,t2dir);

        delta[0] += t2dir[0]*tpair;
        delta[1] += t2dir[1]*tpair;
        delta[2] += t2dir[2]*tpair;

      }

      // theta3 torque
      if (theta3 != 0.0) {

        tpair = -f4t1 * f4t2 * df4t3 * (f2_33 * f4t4_33 * f4t7_33 * f4t8_33 + f2_55 * f4t4_55 * f4t7_55 * f4t8_55) * factor_lj;
        MathExtra::cross3(bx,delr_bsbs_norm,t3dir);

        deltb[0] += t3dir[0]*tpair;
        deltb[1] += t3dir[1]*tpair;
        deltb[2] += t3dir[2]*tpair;

      }

      // theta4 torque
      if (theta4 != 0.0) {

        tpair = -f4t1 * f4t2 * f4t3 * (f2_33 * df4t4_33 * f4t7_33 * f4t8_33 + f2_55 * df4t4_55 * f4t7_55 * f4t8_55) * factor_lj;
        MathExtra::cross3(bz,az,t4dir);

        delta[0] += t4dir[0]*tpair;
        delta[1] += t4dir[1]*tpair;
        delta[2] += t4dir[2]*tpair;

        deltb[0] += t4dir[0]*tpair;
        deltb[1] += t4dir[1]*tpair;
        deltb[2] += t4dir[2]*tpair;

      }

      // theta7 torque
      if (theta7 != 0.0) {

        tpair = -f4t1 * f4t2 * f4t3 * (f2_33 * f4t4_33 * df4t7_33 * f4t8_33 + f2_55 * f4t4_55 * df4t7_55 * f4t8_55) * factor_lj;
        MathExtra::cross3(az,delr_bsbs_norm,t7dir);

        delta[0] += t7dir[0]*tpair;
        delta[1] += t7dir[1]*tpair;
        delta[2] += t7dir[2]*tpair;

      }

      // theta8 torque
      if (theta8 != 0.0) {

        tpair = -f4t1 * f4t2 * f4t3 * (f2_33 * f4t4_33 * f4t7_33 * df4t8_33 + f2_55 * f4t4_55 * f4t7_55 * df4t8_55) * factor_lj;
        MathExtra::cross3(bz,delr_bsbs_norm,t8dir);

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

void PairOxdna3Xstk::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(k_xst,n+1,n+1,"pair:k_xst");

  memory->create(cut_xst_0_33,n+1,n+1,n+1,n+1,"pair:cut_xst_0");
  memory->create(cut_xst_c_33,n+1,n+1,n+1,n+1,"pair:cut_xst_c");
  memory->create(cut_xst_lo_33,n+1,n+1,n+1,n+1,"pair:cut_xst_lo");
  memory->create(cut_xst_hi_33,n+1,n+1,n+1,n+1,"pair:cut_xst_hi");
  memory->create(cut_xst_lc_33,n+1,n+1,n+1,n+1,"pair:cut_xst_lc");
  memory->create(cut_xst_hc_33,n+1,n+1,n+1,n+1,"pair:cut_xst_hc");
  memory->create(cutsq_xst_hc_33,n+1,n+1,n+1,n+1,"pair:cutsq_xst_hc");

  memory->create(cut_xst_0_55,n+1,n+1,n+1,n+1,"pair:cut_xst_0");
  memory->create(cut_xst_c_55,n+1,n+1,n+1,n+1,"pair:cut_xst_c");
  memory->create(cut_xst_lo_55,n+1,n+1,n+1,n+1,"pair:cut_xst_lo");
  memory->create(cut_xst_hi_55,n+1,n+1,n+1,n+1,"pair:cut_xst_hi");
  memory->create(cut_xst_lc_55,n+1,n+1,n+1,n+1,"pair:cut_xst_lc");
  memory->create(cut_xst_hc_55,n+1,n+1,n+1,n+1,"pair:cut_xst_hc");
  memory->create(cutsq_xst_hc_55,n+1,n+1,n+1,n+1,"pair:cutsq_xst_hc");

  memory->create(b_xst_lo,n+1,n+1,"pair:b_xst_lo");
  memory->create(b_xst_hi,n+1,n+1,"pair:b_xst_hi");

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

  memory->create(a_xst4_33,n+1,n+1,n+1,n+1,"pair:a_xst4");
  memory->create(theta_xst4_0_33,n+1,n+1,n+1,n+1,"pair:theta_xst4_0");
  memory->create(dtheta_xst4_ast_33,n+1,n+1,n+1,n+1,"pair:dtheta_xst4_ast");
  memory->create(b_xst4_33,n+1,n+1,n+1,n+1,"pair:b_xst4");
  memory->create(dtheta_xst4_c_33,n+1,n+1,n+1,n+1,"pair:dtheta_xst4_c");

  memory->create(a_xst4_55,n+1,n+1,n+1,n+1,"pair:a_xst4");
  memory->create(theta_xst4_0_55,n+1,n+1,n+1,n+1,"pair:theta_xst4_0");
  memory->create(dtheta_xst4_ast_55,n+1,n+1,n+1,n+1,"pair:dtheta_xst4_ast");
  memory->create(b_xst4_55,n+1,n+1,n+1,n+1,"pair:b_xst4");
  memory->create(dtheta_xst4_c_55,n+1,n+1,n+1,n+1,"pair:dtheta_xst4_c");

  memory->create(a_xst7,n+1,n+1,"pair:a_xst7");
  memory->create(theta_xst7_0_33,n+1,n+1,"pair:theta_xst7_0_33");
  memory->create(theta_xst7_0_55,n+1,n+1,"pair:theta_xst7_0_55");
  memory->create(dtheta_xst7_ast,n+1,n+1,"pair:dtheta_xst7_ast");
  memory->create(b_xst7,n+1,n+1,"pair:b_xst7");
  memory->create(dtheta_xst7_c,n+1,n+1,"pair:dtheta_xst7_c");

  memory->create(a_xst8,n+1,n+1,"pair:a_xst8");
  memory->create(theta_xst8_0_33,n+1,n+1,"pair:theta_xst8_0_33");
  memory->create(theta_xst8_0_55,n+1,n+1,"pair:theta_xst8_0_55");
  memory->create(dtheta_xst8_ast,n+1,n+1,"pair:dtheta_xst8_ast");
  memory->create(b_xst8,n+1,n+1,"pair:b_xst8");
  memory->create(dtheta_xst8_c,n+1,n+1,"pair:dtheta_xst8_c");

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairOxdna3Xstk::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairOxdna3Xstk::coeff(int narg, char **arg)
{
  int count;

  if (narg != 3) error->all(FLERR,"Incorrect args for pair coefficients in oxdna3/xstk" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi,nlo,nhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  assert((ilo == jlo) & (ihi == jhi));
  nlo = ilo;
  nhi = ihi;

  if (nhi > 4) error->all(FLERR, "pair oxdna3/xstk does not support more than 4 atom types for A, C, G and T");

  // cross-stacking interaction
  count = 0;

  double k_xst_one;
  double b_xst_lo_one, b_xst_hi_one;

  double a_xst1_one, theta_xst1_0_one, dtheta_xst1_ast_one;
  double b_xst1_one, dtheta_xst1_c_one;

  double a_xst2_one, theta_xst2_0_one, dtheta_xst2_ast_one;
  double b_xst2_one, dtheta_xst2_c_one;

  double a_xst3_one, theta_xst3_0_one, dtheta_xst3_ast_one;
  double b_xst3_one, dtheta_xst3_c_one;

  double a_xst7_one, theta_xst7_0_33_one, theta_xst7_0_55_one, dtheta_xst7_ast_one;
  double b_xst7_one, dtheta_xst7_c_one;

  double a_xst8_one, theta_xst8_0_33_one, theta_xst8_0_55_one, dtheta_xst8_ast_one;
  double b_xst8_one, dtheta_xst8_c_one;


  for (int i = 0; i <= nhi; i++) { // type 0 for terminal j
    for (int j = 0; j <= nhi; j++) {
      for (int k = 0; k <= nhi; k++) {
        for (int l = 0; l <= nhi; l++) { // type 0 for terminal k
          cut_xst_0_33[i][j][k][l] = 0.0;
          cut_xst_0_55[i][j][k][l] = 0.0;
          cut_xst_c_33[i][j][k][l] = 0.0;
          cut_xst_c_55[i][j][k][l] = 0.0;
          cut_xst_lo_33[i][j][k][l] = 0.0;
          cut_xst_lo_55[i][j][k][l] = 0.0;
          cut_xst_hi_33[i][j][k][l] = 0.0;
          cut_xst_hi_55[i][j][k][l] = 0.0;

          a_xst4_33[i][j][k][l] = 0.0;
          theta_xst4_0_33[i][j][k][l] = 0.0;
          dtheta_xst4_ast_33[i][j][k][l] = 0.0;

          a_xst4_55[i][j][k][l] = 0.0;
          theta_xst4_0_55[i][j][k][l] = 0.0;
          dtheta_xst4_ast_55[i][j][k][l] = 0.0;
        }
      }
    }
  }

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, arg[2], "oxdna3 potential", " (xstk)");
    reader.set_bufsize(131072);
    char * line;
    std::string iloc, jloc, potential_name;

    while ((line = reader.next_line())) {
      try {
        ValueTokenizer values(line);
        iloc = values.next_string();
        jloc = values.next_string();
        potential_name = values.next_string();
        if (iloc == arg[0] && jloc == arg[1] && potential_name == "xstk") {
          k_xst_one = values.next_double();

          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  cut_xst_0_33[i][j][k][l] = values.next_double();
                  cut_xst_0_33[i][j][k][0] += cut_xst_0_33[i][j][k][l];
                  cut_xst_0_33[0][j][k][l] += cut_xst_0_33[i][j][k][l];
                  cut_xst_0_33[0][j][k][0] += cut_xst_0_33[i][j][k][l];
                }
              }
            }
          }

          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  cut_xst_c_33[i][j][k][l] = values.next_double();
                  cut_xst_c_33[i][j][k][0] += cut_xst_c_33[i][j][k][l];
                  cut_xst_c_33[0][j][k][l] += cut_xst_c_33[i][j][k][l];
                  cut_xst_c_33[0][j][k][0] += cut_xst_c_33[i][j][k][l];
                }
              }
            }
          }
          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  cut_xst_lo_33[i][j][k][l] = values.next_double();
                  cut_xst_lo_33[i][j][k][0] += cut_xst_lo_33[i][j][k][l];
                  cut_xst_lo_33[0][j][k][l] += cut_xst_lo_33[i][j][k][l];
                  cut_xst_lo_33[0][j][k][0] += cut_xst_lo_33[i][j][k][l];
                }
              }
            }
          }
          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  cut_xst_hi_33[i][j][k][l] = values.next_double();
                  cut_xst_hi_33[i][j][k][0] += cut_xst_hi_33[i][j][k][l];
                  cut_xst_hi_33[0][j][k][l] += cut_xst_hi_33[i][j][k][l];
                  cut_xst_hi_33[0][j][k][0] += cut_xst_hi_33[i][j][k][l];
                }
              }
            }
          }

          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  cut_xst_0_55[i][j][k][l] = values.next_double();
                  cut_xst_0_55[i][j][k][0] += cut_xst_0_55[i][j][k][l];
                  cut_xst_0_55[0][j][k][l] += cut_xst_0_55[i][j][k][l];
                  cut_xst_0_55[0][j][k][0] += cut_xst_0_55[i][j][k][l];
                }
              }
            }
          }
          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  cut_xst_c_55[i][j][k][l] = values.next_double();
                  cut_xst_c_55[i][j][k][0] += cut_xst_c_55[i][j][k][l];
                  cut_xst_c_55[0][j][k][l] += cut_xst_c_55[i][j][k][l];
                  cut_xst_c_55[0][j][k][0] += cut_xst_c_55[i][j][k][l];
                }
              }
            }
          }
          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  cut_xst_lo_55[i][j][k][l] = values.next_double();
                  cut_xst_lo_55[i][j][k][0] += cut_xst_lo_55[i][j][k][l];
                  cut_xst_lo_55[0][j][k][l] += cut_xst_lo_55[i][j][k][l];
                  cut_xst_lo_55[0][j][k][0] += cut_xst_lo_55[i][j][k][l];
                }
              }
            }
          }
          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  cut_xst_hi_55[i][j][k][l] = values.next_double();
                  cut_xst_hi_55[i][j][k][0] += cut_xst_hi_55[i][j][k][l];
                  cut_xst_hi_55[0][j][k][l] += cut_xst_hi_55[i][j][k][l];
                  cut_xst_hi_55[0][j][k][0] += cut_xst_hi_55[i][j][k][l];
                }
              }
            }
          }

          a_xst1_one = values.next_double();
          theta_xst1_0_one = values.next_double();
          dtheta_xst1_ast_one = values.next_double();

          a_xst2_one = values.next_double();
          theta_xst2_0_one = values.next_double();
          dtheta_xst2_ast_one = values.next_double();

          a_xst3_one = values.next_double();
          theta_xst3_0_one = values.next_double();
          dtheta_xst3_ast_one = values.next_double();

          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  a_xst4_33[i][j][k][l] = values.next_double();
                  a_xst4_33[i][j][k][0] += a_xst4_33[i][j][k][l];
                  a_xst4_33[0][j][k][l] += a_xst4_33[i][j][k][l];
                  a_xst4_33[0][j][k][0] += a_xst4_33[i][j][k][l];
                }
              }
            }
          }
          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  theta_xst4_0_33[i][j][k][l] = values.next_double();
                  theta_xst4_0_33[i][j][k][0] += theta_xst4_0_33[i][j][k][l];
                  theta_xst4_0_33[0][j][k][l] += theta_xst4_0_33[i][j][k][l];
                  theta_xst4_0_33[0][j][k][0] += theta_xst4_0_33[i][j][k][l];
                }
              }
            }
          }
          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  dtheta_xst4_ast_33[i][j][k][l] = values.next_double();
                  dtheta_xst4_ast_33[i][j][k][0] += dtheta_xst4_ast_33[i][j][k][l];
                  dtheta_xst4_ast_33[0][j][k][l] += dtheta_xst4_ast_33[i][j][k][l];
                  dtheta_xst4_ast_33[0][j][k][0] += dtheta_xst4_ast_33[i][j][k][l];
                }
              }
            }
          }
          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  a_xst4_55[i][j][k][l] = values.next_double();
                  a_xst4_55[i][j][k][0] += a_xst4_55[i][j][k][l];
                  a_xst4_55[0][j][k][l] += a_xst4_55[i][j][k][l];
                  a_xst4_55[0][j][k][0] += a_xst4_55[i][j][k][l];
                }
              }
            }
          }
          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  theta_xst4_0_55[i][j][k][l] = values.next_double();
                  theta_xst4_0_55[i][j][k][0] += theta_xst4_0_55[i][j][k][l];
                  theta_xst4_0_55[0][j][k][l] += theta_xst4_0_55[i][j][k][l];
                  theta_xst4_0_55[0][j][k][0] += theta_xst4_0_55[i][j][k][l];
                }
              }
            }
          }
          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  dtheta_xst4_ast_55[i][j][k][l] = values.next_double();
                  dtheta_xst4_ast_55[i][j][k][0] += dtheta_xst4_ast_55[i][j][k][l];
                  dtheta_xst4_ast_55[0][j][k][l] += dtheta_xst4_ast_55[i][j][k][l];
                  dtheta_xst4_ast_55[0][j][k][0] += dtheta_xst4_ast_55[i][j][k][l];
                }
              }
            }
          }

          a_xst7_one = values.next_double();
          theta_xst7_0_33_one = values.next_double();
          theta_xst7_0_55_one = values.next_double();
          dtheta_xst7_ast_one = values.next_double();

          a_xst8_one = values.next_double();
          theta_xst8_0_33_one = values.next_double();
          theta_xst8_0_55_one = values.next_double();
          dtheta_xst8_ast_one = values.next_double();

          break;
        } else continue;
      } catch (std::exception &e) {
        error->one(FLERR, "Problem parsing oxDNA3 potential file: {}", e.what());
      }
    }
    if ((iloc != arg[0]) || (jloc != arg[1]) || (potential_name != "xstk"))
      error->one(FLERR, "No corresponding xstk potential found in file {} for pair type {} {}",
                 arg[2], arg[0], arg[1]);
  }

  // calculate sequence-averaged parameters for terminal base step j-k
  for (int i = nlo; i <= nhi; i++) {
    for (int j = nlo; j <= nhi; j++) {
      for (int k = nlo; k <= nhi; k++) {
        cut_xst_0_33[i][j][k][0] /= nhi;
        cut_xst_c_33[i][j][k][0] /= nhi;
        cut_xst_lo_33[i][j][k][0] /= nhi;
        cut_xst_hi_33[i][j][k][0] /= nhi;

        cut_xst_0_55[i][j][k][0] /= nhi;
        cut_xst_c_55[i][j][k][0] /= nhi;
        cut_xst_lo_55[i][j][k][0] /= nhi;
        cut_xst_hi_55[i][j][k][0] /= nhi;

        a_xst4_33[i][j][k][0] /= nhi;
        theta_xst4_0_33[i][j][k][0] /= nhi;
        dtheta_xst4_ast_33[i][j][k][0] /= nhi;

        a_xst4_55[i][j][k][0] /= nhi;
        theta_xst4_0_55[i][j][k][0] /= nhi;
        dtheta_xst4_ast_55[i][j][k][0] /= nhi;
      }
    }
  }
  for (int j = nlo; j <= nhi; j++) {
    for (int k = nlo; k <= nhi; k++) {
      for (int l = nlo; l <= nhi; l++) {
        cut_xst_0_33[0][j][k][l] /= nhi;
        cut_xst_c_33[0][j][k][l] /= nhi;
        cut_xst_lo_33[0][j][k][l] /= nhi;
        cut_xst_hi_33[0][j][k][l] /= nhi;

        cut_xst_0_55[0][j][k][l] /= nhi;
        cut_xst_c_55[0][j][k][l] /= nhi;
        cut_xst_lo_55[0][j][k][l] /= nhi;
        cut_xst_hi_55[0][j][k][l] /= nhi;

        a_xst4_33[0][j][k][l] /= nhi;
        theta_xst4_0_33[0][j][k][l] /= nhi;
        dtheta_xst4_ast_33[0][j][k][l] /= nhi;

        a_xst4_55[0][j][k][l] /= nhi;
        theta_xst4_0_55[0][j][k][l] /= nhi;
        dtheta_xst4_ast_55[0][j][k][l] /= nhi;
      }
    }
  }
  for (int j = nlo; j <= nhi; j++) {
    for (int k = nlo; k <= nhi; k++) {
      cut_xst_0_33[0][j][k][0] /= powint(nhi,2);
      cut_xst_c_33[0][j][k][0] /= powint(nhi,2);
      cut_xst_lo_33[0][j][k][0] /= powint(nhi,2);
      cut_xst_hi_33[0][j][k][0] /= powint(nhi,2);

      cut_xst_0_55[0][j][k][0] /= powint(nhi,2);
      cut_xst_c_55[0][j][k][0] /= powint(nhi,2);
      cut_xst_lo_55[0][j][k][0] /= powint(nhi,2);
      cut_xst_hi_55[0][j][k][0] /= powint(nhi,2);

      a_xst4_33[0][j][k][0] /= powint(nhi,2);
      theta_xst4_0_33[0][j][k][0] /= powint(nhi,2);
      dtheta_xst4_ast_33[0][j][k][0] /= powint(nhi,2);

      a_xst4_55[0][j][k][0] /= powint(nhi,2);
      theta_xst4_0_55[0][j][k][0] /= powint(nhi,2);
      dtheta_xst4_ast_55[0][j][k][0] /= powint(nhi,2);
    }
  }

  MPI_Bcast(&k_xst_one, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(&cut_xst_0_33[0][0][0][0], 625, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_xst_c_33[0][0][0][0], 625, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_xst_lo_33[0][0][0][0], 625, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_xst_hi_33[0][0][0][0], 625, MPI_DOUBLE, 0, world);

  MPI_Bcast(&cut_xst_0_55[0][0][0][0], 625, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_xst_c_55[0][0][0][0], 625, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_xst_lo_55[0][0][0][0], 625, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_xst_hi_55[0][0][0][0], 625, MPI_DOUBLE, 0, world);

  MPI_Bcast(&a_xst1_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_xst1_0_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_xst1_ast_one, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(&a_xst2_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_xst2_0_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_xst2_ast_one, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(&a_xst3_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_xst3_0_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_xst3_ast_one, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(&a_xst4_33[0][0][0][0], 625, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_xst4_0_33[0][0][0][0], 625, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_xst4_ast_33[0][0][0][0], 625, MPI_DOUBLE, 0, world);

  MPI_Bcast(&a_xst4_55[0][0][0][0], 625, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_xst4_0_55[0][0][0][0], 625, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_xst4_ast_55[0][0][0][0], 625, MPI_DOUBLE, 0, world);

  MPI_Bcast(&a_xst7_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_xst7_0_33_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_xst7_0_55_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_xst7_ast_one, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(&a_xst8_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_xst8_0_33_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_xst8_0_55_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_xst8_ast_one, 1, MPI_DOUBLE, 0, world);

  // smoothing - determined through continuity and differentiability

  // smoothing strength coincidentally identical for all pairs, hence use AAAA tetramer value below
  b_xst_lo_one = 0.25 * (cut_xst_lo_33[1][1][1][1] - cut_xst_0_33[1][1][1][1]) * (cut_xst_lo_33[1][1][1][1] - cut_xst_0_33[1][1][1][1])/
        (0.5 * (cut_xst_lo_33[1][1][1][1] - cut_xst_0_33[1][1][1][1]) * (cut_xst_lo_33[1][1][1][1] - cut_xst_0_33[1][1][1][1]) -
        k_xst_one * 0.5 * (cut_xst_0_33[1][1][1][1] -cut_xst_c_33[1][1][1][1]) * (cut_xst_0_33[1][1][1][1] - cut_xst_c_33[1][1][1][1])/k_xst_one);

  // smoothing strength coincidentally identical for all pairs, hence use AAAA tetramer value below
  b_xst_hi_one = 0.25 * (cut_xst_hi_33[1][1][1][1] - cut_xst_0_33[1][1][1][1]) * (cut_xst_hi_33[1][1][1][1] - cut_xst_0_33[1][1][1][1])/
        (0.5 * (cut_xst_hi_33[1][1][1][1] - cut_xst_0_33[1][1][1][1]) * (cut_xst_hi_33[1][1][1][1] - cut_xst_0_33[1][1][1][1]) -
        k_xst_one * 0.5 * (cut_xst_0_33[1][1][1][1] -cut_xst_c_33[1][1][1][1]) * (cut_xst_0_33[1][1][1][1] - cut_xst_c_33[1][1][1][1])/k_xst_one);

  b_xst1_one = a_xst1_one*a_xst1_one*dtheta_xst1_ast_one*dtheta_xst1_ast_one/(1-a_xst1_one*dtheta_xst1_ast_one*dtheta_xst1_ast_one);
  dtheta_xst1_c_one = 1/(a_xst1_one*dtheta_xst1_ast_one);

  b_xst2_one = a_xst2_one*a_xst2_one*dtheta_xst2_ast_one*dtheta_xst2_ast_one/(1-a_xst2_one*dtheta_xst2_ast_one*dtheta_xst2_ast_one);
  dtheta_xst2_c_one = 1/(a_xst2_one*dtheta_xst2_ast_one);

  b_xst3_one = a_xst3_one*a_xst3_one*dtheta_xst3_ast_one*dtheta_xst3_ast_one/(1-a_xst3_one*dtheta_xst3_ast_one*dtheta_xst3_ast_one);
  dtheta_xst3_c_one = 1/(a_xst3_one*dtheta_xst3_ast_one);

  b_xst7_one = a_xst7_one*a_xst7_one*dtheta_xst7_ast_one*dtheta_xst7_ast_one/(1-a_xst7_one*dtheta_xst7_ast_one*dtheta_xst7_ast_one);
  dtheta_xst7_c_one = 1/(a_xst7_one*dtheta_xst7_ast_one);

  b_xst8_one = a_xst8_one*a_xst8_one*dtheta_xst8_ast_one*dtheta_xst8_ast_one/(1-a_xst8_one*dtheta_xst8_ast_one*dtheta_xst8_ast_one);
  dtheta_xst8_c_one = 1/(a_xst8_one*dtheta_xst8_ast_one);

  // parameters depending on base step
  for (int i = nlo; i <= nhi; i++) {
    for (int j = nlo; j <= nhi; j++) {

      k_xst[i][j] = k_xst_one;

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

      a_xst7[i][j] = a_xst7_one;
      theta_xst7_0_33[i][j] = theta_xst7_0_33_one;
      theta_xst7_0_55[i][j] = theta_xst7_0_55_one;
      dtheta_xst7_ast[i][j] = dtheta_xst7_ast_one;
      b_xst7[i][j] = b_xst7_one;
      dtheta_xst7_c[i][j] = dtheta_xst7_c_one;

      a_xst8[i][j] = a_xst8_one;
      theta_xst8_0_33[i][j] = theta_xst8_0_33_one;
      theta_xst8_0_55[i][j] = theta_xst8_0_55_one;
      dtheta_xst8_ast[i][j] = dtheta_xst8_ast_one;
      b_xst8[i][j] = b_xst8_one;
      dtheta_xst8_c[i][j] = dtheta_xst8_c_one;

    }
  }

  // parameters depending on tetramer
  for (int i = 0; i <= nhi; i++) { // type 0 for terminal j
    for (int j = nlo; j <= nhi; j++) {
      for (int k = nlo; k <= nhi; k++) {
        for (int l = 0; l <= nhi; l++) { // type 0 for terminal k

          cut_xst_lc_33[i][j][k][l] = cut_xst_lo_33[i][j][k][l] -
               0.5*(cut_xst_lo_33[i][j][k][l] - cut_xst_0_33[i][j][k][l])/b_xst_lo_one;
          cut_xst_hc_33[i][j][k][l] = cut_xst_hi_33[i][j][k][l] -
               0.5*(cut_xst_hi_33[i][j][k][l] - cut_xst_0_33[i][j][k][l])/b_xst_hi_one;

          cut_xst_lc_55[i][j][k][l] = cut_xst_lo_55[i][j][k][l] -
               0.5*(cut_xst_lo_55[i][j][k][l] - cut_xst_0_55[i][j][k][l])/b_xst_lo_one;
          cut_xst_hc_55[i][j][k][l] = cut_xst_hi_55[i][j][k][l] -
               0.5*(cut_xst_hi_55[i][j][k][l] - cut_xst_0_55[i][j][k][l])/b_xst_hi_one;

          b_xst4_33[i][j][k][l] = a_xst4_33[i][j][k][l]*a_xst4_33[i][j][k][l]*dtheta_xst4_ast_33[i][j][k][l]
              *dtheta_xst4_ast_33[i][j][k][l]/(1-a_xst4_33[i][j][k][l]
              *dtheta_xst4_ast_33[i][j][k][l]*dtheta_xst4_ast_33[i][j][k][l]);
          dtheta_xst4_c_33[i][j][k][l] = 1/(a_xst4_33[i][j][k][l]*dtheta_xst4_ast_33[i][j][k][l]);

          b_xst4_55[i][j][k][l] = a_xst4_55[i][j][k][l]*a_xst4_55[i][j][k][l]*dtheta_xst4_ast_55[i][j][k][l]
              *dtheta_xst4_ast_55[i][j][k][l]/(1-a_xst4_55[i][j][k][l]
              *dtheta_xst4_ast_55[i][j][k][l]*dtheta_xst4_ast_55[i][j][k][l]);
          dtheta_xst4_c_55[i][j][k][l] = 1/(a_xst4_55[i][j][k][l]*dtheta_xst4_ast_55[i][j][k][l]);

        }

        setflag[j][k] = 1;
        count++;

       }
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/xstk" + utils::errorurl(21));

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */
void PairOxdna3Xstk::init_style()
{
  fix_lrf = nullptr;
  auto fixes = modify->get_fix_by_style("^OXDNA/LRF");
  if (fixes.size() == 0) error->all(FLERR, "Fix OXDNA/LRF not found. Ensure pair oxdna/excv is present");
  else fix_lrf = dynamic_cast<FixOxdnaLRF *>(fixes[0]);

  neighbor->add_request(this, NeighConst::REQ_DEFAULT);
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use regular
------------------------------------------------------------------------- */

void PairOxdna3Xstk::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  if (id  > 0) error->all(FLERR,"Respa not supported");

}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairOxdna3Xstk::init_one(int i, int j)
{

  if (setflag[i][j] == 0) {
    error->all(FLERR,"Coefficient mixing not defined in oxDNA");
  }
  if (offset_flag) {
    error->all(FLERR,"Offset not supported in oxDNA");
  }

  k_xst[j][i] = k_xst[i][j];

  b_xst_lo[j][i] = b_xst_lo[i][j];
  b_xst_hi[j][i] = b_xst_hi[i][j];

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

  a_xst7[j][i] = a_xst7[i][j];
  theta_xst7_0_33[j][i] = theta_xst7_0_33[i][j];
  theta_xst7_0_55[j][i] = theta_xst7_0_55[i][j];
  dtheta_xst7_ast[j][i] = dtheta_xst7_ast[i][j];
  b_xst7[j][i] = b_xst7[i][j];
  dtheta_xst7_c[j][i] = dtheta_xst7_c[i][j];

  a_xst8[j][i] = a_xst8[i][j];
  theta_xst8_0_33[j][i] = theta_xst8_0_33[i][j];
  theta_xst8_0_55[j][i] = theta_xst8_0_55[i][j];
  dtheta_xst8_ast[j][i] = dtheta_xst8_ast[i][j];
  b_xst8[j][i] = b_xst8[i][j];
  dtheta_xst8_c[j][i] = dtheta_xst8_c[i][j];

  // set the master list distance cutoff
  double cut_max=0.0;

  for (int a=0; a<atom->ntypes; a++) {
    for (int b=0; b<atom->ntypes; b++) {
      cut_max = MAX(cut_xst_hc_33[a][i][j][b],cut_max);
      cut_max = MAX(cut_xst_hc_55[a][i][j][b],cut_max);
    }
  }

  return cut_max;

}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdna3Xstk::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {

        fwrite(&k_xst[i][j],sizeof(double),1,fp);

        fwrite(&cut_xst_0_33[i][j],sizeof(double),1,fp);
        fwrite(&cut_xst_c_33[i][j],sizeof(double),1,fp);
        fwrite(&cut_xst_lo_33[i][j],sizeof(double),1,fp);
        fwrite(&cut_xst_hi_33[i][j],sizeof(double),1,fp);
        fwrite(&cut_xst_lc_33[i][j],sizeof(double),1,fp);
        fwrite(&cut_xst_hc_33[i][j],sizeof(double),1,fp);

        fwrite(&cut_xst_0_55[i][j],sizeof(double),1,fp);
        fwrite(&cut_xst_c_55[i][j],sizeof(double),1,fp);
        fwrite(&cut_xst_lo_55[i][j],sizeof(double),1,fp);
        fwrite(&cut_xst_hi_55[i][j],sizeof(double),1,fp);
        fwrite(&cut_xst_lc_55[i][j],sizeof(double),1,fp);
        fwrite(&cut_xst_hc_55[i][j],sizeof(double),1,fp);

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

        fwrite(&a_xst7[i][j],sizeof(double),1,fp);
        fwrite(&theta_xst7_0_33[i][j],sizeof(double),1,fp);
        fwrite(&theta_xst7_0_55[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_xst7_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_xst7[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_xst7_c[i][j],sizeof(double),1,fp);

        fwrite(&a_xst8[i][j],sizeof(double),1,fp);
        fwrite(&theta_xst8_0_33[i][j],sizeof(double),1,fp);
        fwrite(&theta_xst8_0_55[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_xst8_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_xst8[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_xst8_c[i][j],sizeof(double),1,fp);

    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdna3Xstk::read_restart(FILE *fp)
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

          utils::sfread(FLERR,&k_xst[i][j],sizeof(double),1,fp,nullptr,error);

          utils::sfread(FLERR,&cut_xst_0_33[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_xst_c_33[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_xst_lo_33[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_xst_hi_33[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_xst_lc_33[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_xst_hc_33[i][j],sizeof(double),1,fp,nullptr,error);

          utils::sfread(FLERR,&cut_xst_0_55[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_xst_c_55[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_xst_lo_55[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_xst_hi_55[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_xst_lc_55[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_xst_hc_55[i][j],sizeof(double),1,fp,nullptr,error);

          utils::sfread(FLERR,&b_xst_lo[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_xst_hi[i][j],sizeof(double),1,fp,nullptr,error);

          utils::sfread(FLERR,&a_xst1[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&theta_xst1_0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_xst1_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_xst1[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_xst1_c[i][j],sizeof(double),1,fp,nullptr,error);

          utils::sfread(FLERR,&a_xst2[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&theta_xst2_0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_xst2_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_xst2[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_xst2_c[i][j],sizeof(double),1,fp,nullptr,error);

          utils::sfread(FLERR,&a_xst3[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&theta_xst3_0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_xst3_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_xst3[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_xst3_c[i][j],sizeof(double),1,fp,nullptr,error);

          utils::sfread(FLERR,&a_xst7[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&theta_xst7_0_33[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&theta_xst7_0_55[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_xst7_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_xst7[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_xst7_c[i][j],sizeof(double),1,fp,nullptr,error);

          utils::sfread(FLERR,&a_xst8[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&theta_xst8_0_33[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&theta_xst8_0_55[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_xst8_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_xst8[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_xst8_c[i][j],sizeof(double),1,fp,nullptr,error);

        }

        MPI_Bcast(&k_xst[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&cut_xst_0_33[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_xst_c_33[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_xst_lo_33[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_xst_hi_33[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_xst_lc_33[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_xst_hc_33[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&cut_xst_0_55[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_xst_c_55[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_xst_lo_55[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_xst_hi_55[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_xst_lc_55[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_xst_hc_55[i][j],1,MPI_DOUBLE,0,world);

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

        MPI_Bcast(&a_xst7[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_xst7_0_33[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_xst7_0_55[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_xst7_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_xst7[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_xst7_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_xst8[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_xst8_0_33[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_xst8_0_55[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_xst8_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_xst8[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_xst8_c[i][j],1,MPI_DOUBLE,0,world);

      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdna3Xstk::write_restart_settings(FILE *fp)
{
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdna3Xstk::read_restart_settings(FILE *fp)
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

/* ---------------------------------------------------------------------- */

void *PairOxdna3Xstk::extract(const char *str, int &dim)
{
  dim = 2;

  if (strcmp(str,"k_xst") == 0) return (void *) k_xst;

  if (strcmp(str,"cut_xst_0_33") == 0) return (void *) cut_xst_0_33;
  if (strcmp(str,"cut_xst_c_33") == 0) return (void *) cut_xst_c_33;
  if (strcmp(str,"cut_xst_lo_33") == 0) return (void *) cut_xst_lo_33;
  if (strcmp(str,"cut_xst_hi_33") == 0) return (void *) cut_xst_hi_33;
  if (strcmp(str,"cut_xst_lc_33") == 0) return (void *) cut_xst_lc_33;
  if (strcmp(str,"cut_xst_hc_33") == 0) return (void *) cut_xst_hc_33;

  if (strcmp(str,"cut_xst_0_55") == 0) return (void *) cut_xst_0_55;
  if (strcmp(str,"cut_xst_c_55") == 0) return (void *) cut_xst_c_55;
  if (strcmp(str,"cut_xst_lo_55") == 0) return (void *) cut_xst_lo_55;
  if (strcmp(str,"cut_xst_hi_55") == 0) return (void *) cut_xst_hi_55;
  if (strcmp(str,"cut_xst_lc_55") == 0) return (void *) cut_xst_lc_55;
  if (strcmp(str,"cut_xst_hc_55") == 0) return (void *) cut_xst_hc_55;

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

  if (strcmp(str,"a_xst7") == 0) return (void *) a_xst7;
  if (strcmp(str,"theta_xst7_0_33") == 0) return (void *) theta_xst7_0_33;
  if (strcmp(str,"theta_xst7_0_55") == 0) return (void *) theta_xst7_0_55;
  if (strcmp(str,"dtheta_xst7_ast") == 0) return (void *) dtheta_xst7_ast;
  if (strcmp(str,"b_xst7") == 0) return (void *) b_xst7;
  if (strcmp(str,"dtheta_xst7_c") == 0) return (void *) dtheta_xst7_c;

  if (strcmp(str,"a_xst8") == 0) return (void *) a_xst8;
  if (strcmp(str,"theta_xst8_0_33") == 0) return (void *) theta_xst8_0_33;
  if (strcmp(str,"theta_xst8_0_55") == 0) return (void *) theta_xst8_0_55;
  if (strcmp(str,"dtheta_xst8_ast") == 0) return (void *) dtheta_xst8_ast;
  if (strcmp(str,"b_xst8") == 0) return (void *) b_xst8;
  if (strcmp(str,"dtheta_xst8_c") == 0) return (void *) dtheta_xst8_c;

  return nullptr;
}
