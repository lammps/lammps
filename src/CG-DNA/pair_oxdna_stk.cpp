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

#include "pair_oxdna_stk.h"
#include "nucleotide_oxdna.h"

#include "atom.h"
#include "comm.h"
#include "constants_oxdna.h"
#include "error.h"
#include "fix_oxdna_lrf.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "mf_oxdna.h"
#include "modify.h"
#include "neighbor.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>
#include <cassert>
#include <exception>

using namespace LAMMPS_NS;
using namespace MFOxdna;

/* ---------------------------------------------------------------------- */

PairOxdnaStk::PairOxdnaStk(LAMMPS *lmp) :
    Pair(lmp), epsilon_st(nullptr), a_st(nullptr), cut_st_0(nullptr), cut_st_c(nullptr),
    cut_st_lo(nullptr), cut_st_hi(nullptr), cut_st_lc(nullptr), cut_st_hc(nullptr),
    b_st_lo(nullptr), b_st_hi(nullptr), shift_st(nullptr), cutsq_st_hc(nullptr), a_st4(nullptr),
    theta_st4_0(nullptr), dtheta_st4_ast(nullptr), b_st4(nullptr), dtheta_st4_c(nullptr),
    a_st5(nullptr), theta_st5_0(nullptr), dtheta_st5_ast(nullptr), b_st5(nullptr),
    dtheta_st5_c(nullptr), a_st6(nullptr), theta_st6_0(nullptr), dtheta_st6_ast(nullptr),
    b_st6(nullptr), dtheta_st6_c(nullptr), a_st1(nullptr), cosphi_st1_ast(nullptr), b_st1(nullptr),
    cosphi_st1_c(nullptr), a_st2(nullptr), cosphi_st2_ast(nullptr), b_st2(nullptr),
    cosphi_st2_c(nullptr), nxyz_xtrct(nullptr), fix_lrf(nullptr)
{
  // sequence-specific stacking strength
  // A:0 C:1 G:2 T:3, 3'- [i][j] -5'

  eta_st[0][0] = 1.11960;
  eta_st[1][0] = 1.00852;
  eta_st[2][0] = 0.96950;
  eta_st[3][0] = 0.99632;

  eta_st[0][1] = 1.01889;
  eta_st[1][1] = 0.97804;
  eta_st[2][1] = 1.02681;
  eta_st[3][1] = 0.96950;

  eta_st[0][2] = 0.98169;
  eta_st[1][2] = 1.05913;
  eta_st[2][2] = 0.97804;
  eta_st[3][2] = 1.00852;

  eta_st[0][3] = 0.94694;
  eta_st[1][3] = 0.98169;
  eta_st[2][3] = 1.01889;
  eta_st[3][3] = 0.96383;

  single_enable = 0;
  writedata = 0;
  trim_flag = 0;
}

/* ---------------------------------------------------------------------- */

PairOxdnaStk::~PairOxdnaStk()
{
  if (allocated) {

    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(epsilon_st);
    memory->destroy(a_st);
    memory->destroy(cut_st_0);
    memory->destroy(cut_st_c);
    memory->destroy(cut_st_lo);
    memory->destroy(cut_st_hi);
    memory->destroy(cut_st_lc);
    memory->destroy(cut_st_hc);
    memory->destroy(b_st_lo);
    memory->destroy(b_st_hi);
    memory->destroy(shift_st);
    memory->destroy(cutsq_st_hc);

    memory->destroy(a_st4);
    memory->destroy(theta_st4_0);
    memory->destroy(dtheta_st4_ast);
    memory->destroy(b_st4);
    memory->destroy(dtheta_st4_c);

    memory->destroy(a_st5);
    memory->destroy(theta_st5_0);
    memory->destroy(dtheta_st5_ast);
    memory->destroy(b_st5);
    memory->destroy(dtheta_st5_c);

    memory->destroy(a_st6);
    memory->destroy(theta_st6_0);
    memory->destroy(dtheta_st6_ast);
    memory->destroy(b_st6);
    memory->destroy(dtheta_st6_c);

    memory->destroy(a_st1);
    memory->destroy(cosphi_st1_ast);
    memory->destroy(b_st1);
    memory->destroy(cosphi_st1_c);
    memory->destroy(a_st2);
    memory->destroy(cosphi_st2_ast);
    memory->destroy(b_st2);
    memory->destroy(cosphi_st2_c);

  }
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators

   NOTE: Although this is a pair style interaction, the algorithm below
   follows the virial incrementation of the bond style. This is because
   the bond topology is used in the main compute loop.
------------------------------------------------------------------------- */

void PairOxdnaStk::ev_tally_xyz(int i, int j, int nlocal, int newton_bond,
                    double evdwl,
                    double fx, double fy, double fz,
                    double delx, double dely, double delz)
{
  double evdwlhalf,v[6];

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) eng_vdwl += evdwl;
      else {
        evdwlhalf = 0.5*evdwl;
        if (i < nlocal) eng_vdwl += evdwlhalf;
        if (j < nlocal) eng_vdwl += evdwlhalf;
      }
    }
    if (eflag_atom) {
      evdwlhalf = 0.5*evdwl;
      if (newton_bond || i < nlocal) eatom[i] += evdwlhalf;
      if (newton_bond || j < nlocal) eatom[j] += evdwlhalf;
    }
  }

  if (vflag_either) {
    v[0] = delx*fx;
    v[1] = dely*fy;
    v[2] = delz*fz;
    v[3] = delx*fy;
    v[4] = delx*fz;
    v[5] = dely*fz;

    if (vflag_global) {
      if (newton_bond) {
        virial[0] += v[0];
        virial[1] += v[1];
        virial[2] += v[2];
        virial[3] += v[3];
        virial[4] += v[4];
        virial[5] += v[5];
      } else {
        if (i < nlocal) {
          virial[0] += 0.5*v[0];
          virial[1] += 0.5*v[1];
          virial[2] += 0.5*v[2];
          virial[3] += 0.5*v[3];
          virial[4] += 0.5*v[4];
          virial[5] += 0.5*v[5];
        }
        if (j < nlocal) {
          virial[0] += 0.5*v[0];
          virial[1] += 0.5*v[1];
          virial[2] += 0.5*v[2];
          virial[3] += 0.5*v[3];
          virial[4] += 0.5*v[4];
          virial[5] += 0.5*v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        vatom[i][0] += 0.5*v[0];
        vatom[i][1] += 0.5*v[1];
        vatom[i][2] += 0.5*v[2];
        vatom[i][3] += 0.5*v[3];
        vatom[i][4] += 0.5*v[4];
        vatom[i][5] += 0.5*v[5];
      }
      if (newton_bond || j < nlocal) {
        vatom[j][0] += 0.5*v[0];
        vatom[j][1] += 0.5*v[1];
        vatom[j][2] += 0.5*v[2];
        vatom[j][3] += 0.5*v[3];
        vatom[j][4] += 0.5*v[4];
        vatom[j][5] += 0.5*v[5];
      }
    }
  }
}

/* ----------------------------------------------------------------------
    compute vector COM-backbone interaction site in oxDNA
------------------------------------------------------------------------- */
inline void PairOxdnaStk::compute_backbone_site(double e1[3], double /*e2*/[3],
    double /*e3*/[3], double rbk[3]) const
{
  NucleotideOxdna1 oxdna1;
  oxdna1.backbone_site(e1, nullptr, nullptr, rbk);
}

/* ----------------------------------------------------------------------
    compute vector COM-stacking interaction site in oxDNA
------------------------------------------------------------------------- */
inline void PairOxdnaStk::compute_stacking_site(double e1[3], double /*e2*/[3],
    double /*e3*/[3], double rstk[3]) const
{
  NucleotideOxdna1 oxdna1;
  oxdna1.stacking_site(e1, nullptr, nullptr, rstk);
}

/* ----------------------------------------------------------------------
   compute function for oxDNA pair interactions
   s=sugar-phosphate backbone site, b=base site, st=stacking site
------------------------------------------------------------------------- */

void PairOxdnaStk::compute(int eflag, int vflag)
{
  double delf[3],delta[3],deltb[3]; // force, torque increment;
  double evdwl,finc,tpair;
  double delr_bkbk[3],delr_bkbk_norm[3],rsq_bkbk,r_bkbk,rinv_bkbk;
  double delr_stkstk[3],delr_stkstk_norm[3],rsq_stkstk,r_stkstk,rinv_stkstk;
  double theta4,t4dir[3],cost4;
  double theta5p,t5pdir[3],cost5p;
  double theta6p,t6pdir[3],cost6p;
  double cosphi1,cosphi2,cosphi1dir[3],cosphi2dir[3];

  // vectors COM-backbone site, COM-stacking site in lab frame
  double ra_cbk[3],ra_cstk[3];
  double rb_cbk[3],rb_cstk[3];
  // Cartesian unit vectors in lab frame
  double ax[3],ay[3],az[3];
  double bx[3],by[3],bz[3];

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  int *type = atom->type;

  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;

  tagint *id3p = atom->id3p;
  tagint *id5p = atom->id5p;

  int a,b,btemp,in,a3ptype,atype,btype,b5ptype;

  double f1,f4t4,f4t5,f4t6,f5c1,f5c2;
  double df1,df4t4,df4t5,df4t6,df5c1,df5c2;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  // nxyz_xtrct = extracted local unit vectors in lab frame from fix OXDNA/LRF
  nxyz_xtrct = fix_lrf->array_atom;

  // loop over stacking interaction neighbors using bond topology

  for (in = 0; in < nbondlist; in++) {

    a = bondlist[in][0];
    b = bondlist[in][1];

    // directionality test: a -> b is 3' -> 5'
    if(atom->tag[b] != id5p[a]) {

      btemp = b;
      b = a;
      a = btemp;

    }

    // a now in 3' direction, b in 5' direction

    ax[0] = nxyz_xtrct[a][0];
    ax[1] = nxyz_xtrct[a][1];
    ax[2] = nxyz_xtrct[a][2];
    bx[0] = nxyz_xtrct[b][0];
    bx[1] = nxyz_xtrct[b][1];
    bx[2] = nxyz_xtrct[b][2];
    // (a/b)y/z not needed here as oxDNA(1) co-linear

    // vector COM a - stacking site a
    compute_stacking_site(ax,ay,az,ra_cstk);

    // vector COM b - stacking site b
    compute_stacking_site(bx,by,bz,rb_cstk);

    // vector stacking site a to b
    delr_stkstk[0] = x[b][0] + rb_cstk[0] - x[a][0] - ra_cstk[0];
    delr_stkstk[1] = x[b][1] + rb_cstk[1] - x[a][1] - ra_cstk[1];
    delr_stkstk[2] = x[b][2] + rb_cstk[2] - x[a][2] - ra_cstk[2];

    // determine tetramer types
    // 3'neighbor a - a - b - 5'neighbor b

    if (id3p[a] != -1) {
      a3ptype = type[atom->map(id3p[a])];
    }
    else a3ptype = 0;

    atype = type[a];
    btype = type[b];

    if (id5p[b] != -1) {
      b5ptype = type[atom->map(id5p[b])];
    }
    else b5ptype = 0;

    rsq_stkstk = delr_stkstk[0]*delr_stkstk[0] + delr_stkstk[1]*delr_stkstk[1] + delr_stkstk[2]*delr_stkstk[2];
    r_stkstk = sqrt(rsq_stkstk);
    rinv_stkstk = 1.0/r_stkstk;

    delr_stkstk_norm[0] = delr_stkstk[0] * rinv_stkstk;
    delr_stkstk_norm[1] = delr_stkstk[1] * rinv_stkstk;
    delr_stkstk_norm[2] = delr_stkstk[2] * rinv_stkstk;

    // vector COM a - backbone site a
    compute_backbone_site(ax,ay,az,ra_cbk);

    // vector COM b - backbone site b
    compute_backbone_site(bx,by,bz,rb_cbk);

    // vector backbone site b to a
    delr_bkbk[0] = (x[b][0] + rb_cbk[0] - x[a][0] - ra_cbk[0]);
    delr_bkbk[1] = (x[b][1] + rb_cbk[1] - x[a][1] - ra_cbk[1]);
    delr_bkbk[2] = (x[b][2] + rb_cbk[2] - x[a][2] - ra_cbk[2]);

    rsq_bkbk = delr_bkbk[0]*delr_bkbk[0] + delr_bkbk[1]*delr_bkbk[1] + delr_bkbk[2]*delr_bkbk[2];
    r_bkbk = sqrt(rsq_bkbk);
    rinv_bkbk = 1.0/r_bkbk;

    delr_bkbk_norm[0] = delr_bkbk[0] * rinv_bkbk;
    delr_bkbk_norm[1] = delr_bkbk[1] * rinv_bkbk;
    delr_bkbk_norm[2] = delr_bkbk[2] * rinv_bkbk;

    f1 = F1(r_stkstk, epsilon_st[atype][btype], a_st[atype][btype], cut_st_0[a3ptype][atype][btype][b5ptype],
        cut_st_lc[a3ptype][atype][btype][b5ptype], cut_st_hc[a3ptype][atype][btype][b5ptype],
        cut_st_lo[a3ptype][atype][btype][b5ptype], cut_st_hi[a3ptype][atype][btype][b5ptype],
        b_st_lo[atype][btype], b_st_hi[atype][btype], shift_st[a3ptype][atype][btype][b5ptype]);

    // early rejection criterium
    if (f1 != 0.0) {

    az[0] = nxyz_xtrct[a][6];
    az[1] = nxyz_xtrct[a][7];
    az[2] = nxyz_xtrct[a][8];
    bz[0] = nxyz_xtrct[b][6];
    bz[1] = nxyz_xtrct[b][7];
    bz[2] = nxyz_xtrct[b][8];

    // theta4 angle and correction
    cost4 = MathExtra::dot3(bz,az);
    if (cost4 >  1.0) cost4 =  1.0;
    if (cost4 < -1.0) cost4 = -1.0;
    theta4 = acos(cost4);

    f4t4 = F4(theta4, a_st4[a3ptype][atype][btype][b5ptype], theta_st4_0[atype][btype],
        dtheta_st4_ast[a3ptype][atype][btype][b5ptype], b_st4[a3ptype][atype][btype][b5ptype],
        dtheta_st4_c[a3ptype][atype][btype][b5ptype]);

    // early rejection criterium
    if (f4t4 != 0.0) {

    // theta5 angle and correction
    cost5p  = MathExtra::dot3(delr_stkstk_norm,bz);
    if (cost5p >  1.0) cost5p =  1.0;
    if (cost5p < -1.0) cost5p = -1.0;
    theta5p = acos(cost5p);

    f4t5 = F4(theta5p, a_st5[atype][btype], theta_st5_0[atype][btype], dtheta_st5_ast[atype][btype],
        b_st5[atype][btype], dtheta_st5_c[atype][btype]);

    // early rejection criterium
    if (f4t5 != 0.0) {

    ay[0] = nxyz_xtrct[a][3];
    ay[1] = nxyz_xtrct[a][4];
    ay[2] = nxyz_xtrct[a][5];
    by[0] = nxyz_xtrct[b][3];
    by[1] = nxyz_xtrct[b][4];
    by[2] = nxyz_xtrct[b][5];

    cost6p = MathExtra::dot3(delr_stkstk_norm,az);
    if (cost6p >  1.0) cost6p =  1.0;
    if (cost6p < -1.0) cost6p = -1.0;
    theta6p = acos(cost6p);

    cosphi1 = MathExtra::dot3(delr_bkbk_norm,by);
    if (cosphi1 >  1.0) cosphi1 =  1.0;
    if (cosphi1 < -1.0) cosphi1 = -1.0;

    cosphi2 = MathExtra::dot3(delr_bkbk_norm,ay);
    if (cosphi2 >  1.0) cosphi2 =  1.0;
    if (cosphi2 < -1.0) cosphi2 = -1.0;

    f4t6 = F4(theta6p, a_st6[atype][btype], theta_st6_0[atype][btype], dtheta_st6_ast[atype][btype],
        b_st6[atype][btype], dtheta_st6_c[atype][btype]);

    f5c1 = F5(-cosphi1, a_st1[atype][btype], -cosphi_st1_ast[atype][btype], b_st1[atype][btype],
        -cosphi_st1_c[atype][btype]);

    f5c2 = F5(-cosphi2, a_st2[atype][btype], -cosphi_st2_ast[atype][btype], b_st2[atype][btype],
        -cosphi_st2_c[atype][btype]);

    evdwl = f1 * f4t4 * f4t5 * f4t6 * f5c1 * f5c2;

    // early rejection criterium
    if (evdwl != 0.0) {

    df1 = DF1(r_stkstk, epsilon_st[atype][btype], a_st[atype][btype], cut_st_0[a3ptype][atype][btype][b5ptype],
        cut_st_lc[a3ptype][atype][btype][b5ptype], cut_st_hc[a3ptype][atype][btype][b5ptype], cut_st_lo[a3ptype][atype][btype][b5ptype],
        cut_st_hi[a3ptype][atype][btype][b5ptype], b_st_lo[atype][btype], b_st_hi[atype][btype]);

    df4t4 = DF4(theta4, a_st4[a3ptype][atype][btype][b5ptype], theta_st4_0[atype][btype],
        dtheta_st4_ast[a3ptype][atype][btype][b5ptype], b_st4[a3ptype][atype][btype][b5ptype],
        dtheta_st4_c[a3ptype][atype][btype][b5ptype])/sin(theta4);

    df4t5 = DF4(theta5p, a_st5[atype][btype], theta_st5_0[atype][btype], dtheta_st5_ast[atype][btype],
        b_st5[atype][btype], dtheta_st5_c[atype][btype])/sin(theta5p);

    df4t6 = DF4(theta6p, a_st6[atype][btype], theta_st6_0[atype][btype], dtheta_st6_ast[atype][btype],
        b_st6[atype][btype], dtheta_st6_c[atype][btype])/sin(theta6p);

    df5c1 = DF5(-cosphi1, a_st1[atype][btype], -cosphi_st1_ast[atype][btype], b_st1[atype][btype],
        -cosphi_st1_c[atype][btype]);

    df5c2 = DF5(-cosphi2, a_st2[atype][btype], -cosphi_st2_ast[atype][btype], b_st2[atype][btype],
        -cosphi_st2_c[atype][btype]);


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
    finc  = -df1 * f4t4 * f4t5 * f4t6 * f5c1 * f5c2;

    delf[0] += delr_stkstk[0] * finc;
    delf[1] += delr_stkstk[1] * finc;
    delf[2] += delr_stkstk[2] * finc;

    // theta5p force
    if (theta5p != 0.0) {

      finc   = -f1 * f4t4 * df4t5 * f4t6 * f5c1 * f5c2 * rinv_stkstk;

      delf[0] += (delr_stkstk_norm[0]*cost5p - bz[0]) * finc;
      delf[1] += (delr_stkstk_norm[1]*cost5p - bz[1]) * finc;
      delf[2] += (delr_stkstk_norm[2]*cost5p - bz[2]) * finc;

    }

    // theta6p force
    if (theta6p != 0.0) {

      finc   = -f1 * f4t4 * f4t5 * df4t6 * f5c1 * f5c2 * rinv_stkstk;

      delf[0] += (delr_stkstk_norm[0]*cost6p - az[0]) * finc;
      delf[1] += (delr_stkstk_norm[1]*cost6p - az[1]) * finc;
      delf[2] += (delr_stkstk_norm[2]*cost6p - az[2]) * finc;

    }

    // increment forces and torques
    if (newton_bond || a < nlocal) {

      f[a][0] -= delf[0];
      f[a][1] -= delf[1];
      f[a][2] -= delf[2];

      MathExtra::cross3(ra_cstk,delf,delta);

    }
    if (newton_bond || b < nlocal) {

      f[b][0] += delf[0];
      f[b][1] += delf[1];
      f[b][2] += delf[2];

      MathExtra::cross3(rb_cstk,delf,deltb);

    }

    if (newton_bond || a < nlocal) {

      torque[a][0] -= delta[0];
      torque[a][1] -= delta[1];
      torque[a][2] -= delta[2];

    }
    if (newton_bond || b < nlocal) {

      torque[b][0] += deltb[0];
      torque[b][1] += deltb[1];
      torque[b][2] += deltb[2];

    }

    // increment energy and virial
    // NOTE: The virial is calculated on the 'molecular' basis.
    // (see G. Ciccotti and J.P. Ryckaert, Comp. Phys. Rep. 4, 345-392 (1986))

    if (evflag) ev_tally_xyz(a,b,nlocal,newton_bond,evdwl,
        delf[0],delf[1],delf[2],x[b][0]-x[a][0],x[b][1]-x[a][1],x[b][2]-x[a][2]);

    // force, torque and virial contribution for forces between backbone sites

    delf[0] = 0.0;
    delf[1] = 0.0;
    delf[2] = 0.0;

    delta[0] = 0.0;
    delta[1] = 0.0;
    delta[2] = 0.0;

    deltb[0] = 0.0;
    deltb[1] = 0.0;
    deltb[2] = 0.0;

    // cosphi1 force
    if (cosphi1 != 0.0) {

      finc   = -f1 * f4t4 * f4t5 * f4t6 * df5c1 * f5c2 * rinv_bkbk;

      delf[0] += (delr_bkbk_norm[0]*cosphi1 - by[0]) * finc;
      delf[1] += (delr_bkbk_norm[1]*cosphi1 - by[1]) * finc;
      delf[2] += (delr_bkbk_norm[2]*cosphi1 - by[2]) * finc;

    }

    // cosphi2 force
    if (cosphi2 != 0.0) {

      finc   = -f1 * f4t4 * f4t5 * f4t6 * f5c1 * df5c2 * rinv_bkbk;

      delf[0] += (delr_bkbk_norm[0]*cosphi2 - ay[0]) * finc;
      delf[1] += (delr_bkbk_norm[1]*cosphi2 - ay[1]) * finc;
      delf[2] += (delr_bkbk_norm[2]*cosphi2 - ay[2]) * finc;

    }

    // increment forces and torques
    if (newton_bond || a < nlocal) {

      f[a][0] -= delf[0];
      f[a][1] -= delf[1];
      f[a][2] -= delf[2];

      MathExtra::cross3(ra_cbk,delf,delta);

    }
    if (newton_bond || b < nlocal) {

      f[b][0] += delf[0];
      f[b][1] += delf[1];
      f[b][2] += delf[2];

      MathExtra::cross3(rb_cbk,delf,deltb);

    }

    if (newton_bond || a < nlocal) {

      torque[a][0] -= delta[0];
      torque[a][1] -= delta[1];
      torque[a][2] -= delta[2];

    }
    if (newton_bond || b < nlocal) {

      torque[b][0] += deltb[0];
      torque[b][1] += deltb[1];
      torque[b][2] += deltb[2];

    }

    // increment virial only
    if (evflag) ev_tally_xyz(a,b,nlocal,newton_bond,0.0,
        delf[0],delf[1],delf[2],x[b][0]-x[a][0],x[b][1]-x[a][1],x[b][2]-x[a][2]);

    // pure torques not expressible as r x f

    delta[0] = 0.0;
    delta[1] = 0.0;
    delta[2] = 0.0;
    deltb[0] = 0.0;
    deltb[1] = 0.0;
    deltb[2] = 0.0;

    // theta4 torque
    if (theta4 != 0.0) {

      tpair = -f1 * df4t4 * f4t5 * f4t6 * f5c1 * f5c2;
      MathExtra::cross3(az,bz,t4dir);

      delta[0] += t4dir[0]*tpair;
      delta[1] += t4dir[1]*tpair;
      delta[2] += t4dir[2]*tpair;

      deltb[0] += t4dir[0]*tpair;
      deltb[1] += t4dir[1]*tpair;
      deltb[2] += t4dir[2]*tpair;

    }

    // theta5p torque
    if (theta5p != 0.0) {

      tpair = -f1 * f4t4 * df4t5 * f4t6 * f5c1 * f5c2;
      MathExtra::cross3(delr_stkstk_norm,bz,t5pdir);

      deltb[0] += t5pdir[0] * tpair;
      deltb[1] += t5pdir[1] * tpair;
      deltb[2] += t5pdir[2] * tpair;

    }

    // theta6p torque
    if (theta6p != 0.0) {

      tpair = -f1 * f4t4 * f4t5 * df4t6 * f5c1 * f5c2;
      MathExtra::cross3(delr_stkstk_norm,az,t6pdir);

      delta[0] -= t6pdir[0] * tpair;
      delta[1] -= t6pdir[1] * tpair;
      delta[2] -= t6pdir[2] * tpair;

    }

    // cosphi1 torque
    if (cosphi1 != 0.0) {

      tpair   = -f1 * f4t4 * f4t5 * f4t6 * df5c1 * f5c2;
      MathExtra::cross3(delr_bkbk_norm,by,cosphi1dir);

      deltb[0] += cosphi1dir[0] * tpair;
      deltb[1] += cosphi1dir[1] * tpair;
      deltb[2] += cosphi1dir[2] * tpair;

    }

    // cosphi2 torque
    if (cosphi2 != 0.0) {

      tpair   = -f1 * f4t4 * f4t5 * f4t6 * f5c1 * df5c2;
      MathExtra::cross3(delr_bkbk_norm,ay,cosphi2dir);

      delta[0] -= cosphi2dir[0] * tpair;
      delta[1] -= cosphi2dir[1] * tpair;
      delta[2] -= cosphi2dir[2] * tpair;

    }

    // increment torques
    if (newton_bond || a < nlocal) {

      torque[a][0] -= delta[0];
      torque[a][1] -= delta[1];
      torque[a][2] -= delta[2];

    }
    if (newton_bond || b < nlocal) {

      torque[b][0] += deltb[0];
      torque[b][1] += deltb[1];
      torque[b][2] += deltb[2];

    }

    }
    }
    }
    }
    // end early rejection criteria

  }
  // end stacking interaction

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairOxdnaStk::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(epsilon_st,n+1,n+1,"pair:epsilon_st");
  memory->create(a_st,n+1,n+1,"pair:a_st");
  memory->create(cut_st_0,n+1,n+1,n+1,n+1,"pair:cut_st_0");
  memory->create(cut_st_c,n+1,n+1,n+1,n+1,"pair:cut_st_c");
  memory->create(cut_st_lo,n+1,n+1,n+1,n+1,"pair:cut_st_lo");
  memory->create(cut_st_hi,n+1,n+1,n+1,n+1,"pair:cut_st_hi");
  memory->create(cut_st_lc,n+1,n+1,n+1,n+1,"pair:cut_st_lc");
  memory->create(cut_st_hc,n+1,n+1,n+1,n+1,"pair:cut_st_hc");
  memory->create(b_st_lo,n+1,n+1,"pair:b_st_lo");
  memory->create(b_st_hi,n+1,n+1,"pair:b_st_hi");
  memory->create(shift_st,n+1,n+1,n+1,n+1,"pair:shift_st");
  memory->create(cutsq_st_hc,n+1,n+1,n+1,n+1,"pair:cutsq_st_hc");

  memory->create(a_st4,n+1,n+1,n+1,n+1,"pair:a_st4");
  memory->create(theta_st4_0,n+1,n+1,"pair:theta_st4_0");
  memory->create(dtheta_st4_ast,n+1,n+1,n+1,n+1,"pair:dtheta_st4_ast");
  memory->create(b_st4,n+1,n+1,n+1,n+1,"pair:b_st4");
  memory->create(dtheta_st4_c,n+1,n+1,n+1,n+1,"pair:dtheta_st4_c");

  memory->create(a_st5,n+1,n+1,"pair:a_st5");
  memory->create(theta_st5_0,n+1,n+1,"pair:theta_st5_0");
  memory->create(dtheta_st5_ast,n+1,n+1,"pair:dtheta_st5_ast");
  memory->create(b_st5,n+1,n+1,"pair:b_st5");
  memory->create(dtheta_st5_c,n+1,n+1,"pair:dtheta_st5_c");

  memory->create(a_st6,n+1,n+1,"pair:a_st6");
  memory->create(theta_st6_0,n+1,n+1,"pair:theta_st6_0");
  memory->create(dtheta_st6_ast,n+1,n+1,"pair:dtheta_st6_ast");
  memory->create(b_st6,n+1,n+1,"pair:b_st6");
  memory->create(dtheta_st6_c,n+1,n+1,"pair:dtheta_st6_c");

  memory->create(a_st1,n+1,n+1,"pair:a_st1");
  memory->create(cosphi_st1_ast,n+1,n+1,"pair:cosphi_st1_ast");
  memory->create(b_st1,n+1,n+1,"pair:b_st1");
  memory->create(cosphi_st1_c,n+1,n+1,"pair:cosphi_st1_c");
  memory->create(a_st2,n+1,n+1,"pair:a_st2");
  memory->create(cosphi_st2_ast,n+1,n+1,"pair:cosphi_st2_ast");
  memory->create(b_st2,n+1,n+1,"pair:b_st2");
  memory->create(cosphi_st2_c,n+1,n+1,"pair:cosphi_st2_c");

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairOxdnaStk::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

}

/* ----------------------------------------------------------------------
   return temperature dependent oxDNA stacking strength
------------------------------------------------------------------------- */

double PairOxdnaStk::stacking_strength(double xi_st, double kappa_st, double T)
{
  double eps;

  eps = xi_st + kappa_st * T;

  return eps;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairOxdnaStk::coeff(int narg, char **arg)
{
  int count;

  if (narg != 5 && narg != 24) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/stk" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi,nlo,nhi,imod4,jmod4,kmod4;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  assert((ilo == jlo) & (ihi == jhi));
  nlo = ilo;
  nhi = ihi;

  // stacking interaction
  count = 0;

  double T, epsilon_st_one, xi_st_one, kappa_st_one, a_st_one, b_st_lo_one, b_st_hi_one;
  double cut_st_0_one, cut_st_c_one, cut_st_lo_one, cut_st_hi_one;
  double cut_st_lc_one, cut_st_hc_one, tmp, shift_st_one;

  double a_st4_one, theta_st4_0_one, dtheta_st4_ast_one;
  double b_st4_one, dtheta_st4_c_one;

  double a_st5_one, theta_st5_0_one, dtheta_st5_ast_one;
  double b_st5_one, dtheta_st5_c_one;

  double a_st6_one, theta_st6_0_one, dtheta_st6_ast_one;
  double b_st6_one, dtheta_st6_c_one;

  double a_st1_one, cosphi_st1_ast_one, b_st1_one, cosphi_st1_c_one;
  double a_st2_one, cosphi_st2_ast_one, b_st2_one, cosphi_st2_c_one;

  if (strcmp(arg[2], "seqav") != 0 && strcmp(arg[2], "seqdep") != 0) {
    error->all(FLERR,"Incorrect setting, select seqav or seqdep in oxdna/stk");
  }
  if (strcmp(arg[2],"seqav")  == 0) seqdepflag = 0;
  if (strcmp(arg[2],"seqdep") == 0) seqdepflag = 1;

  T = utils::numeric(FLERR,arg[3],false,lmp);

  if (narg == 24) { // values are listed in input
    xi_st_one = utils::numeric(FLERR,arg[4],false,lmp);
    kappa_st_one = utils::numeric(FLERR,arg[5],false,lmp);

    epsilon_st_one = stacking_strength(xi_st_one, kappa_st_one, T);

    a_st_one = utils::numeric(FLERR,arg[6],false,lmp);
    cut_st_0_one = utils::numeric(FLERR,arg[7],false,lmp);
    cut_st_c_one = utils::numeric(FLERR,arg[8],false,lmp);
    cut_st_lo_one = utils::numeric(FLERR,arg[9],false,lmp);
    cut_st_hi_one = utils::numeric(FLERR,arg[10],false,lmp);

    a_st4_one = utils::numeric(FLERR,arg[11],false,lmp);
    theta_st4_0_one = utils::numeric(FLERR,arg[12],false,lmp);
    dtheta_st4_ast_one = utils::numeric(FLERR,arg[13],false,lmp);
    a_st5_one = utils::numeric(FLERR,arg[14],false,lmp);
    theta_st5_0_one = utils::numeric(FLERR,arg[15],false,lmp);
    dtheta_st5_ast_one = utils::numeric(FLERR,arg[16],false,lmp);
    a_st6_one = utils::numeric(FLERR,arg[17],false,lmp);
    theta_st6_0_one = utils::numeric(FLERR,arg[18],false,lmp);
    dtheta_st6_ast_one = utils::numeric(FLERR,arg[19],false,lmp);
    a_st1_one = utils::numeric(FLERR,arg[20],false,lmp);
    cosphi_st1_ast_one = utils::numeric(FLERR,arg[21],false,lmp);
    a_st2_one = utils::numeric(FLERR,arg[22],false,lmp);
    cosphi_st2_ast_one = utils::numeric(FLERR,arg[23],false,lmp);
  } else { // read values from potential file
    if (comm->me == 0) {
      PotentialFileReader reader(lmp, arg[4], "oxdna potential", " (stk)");
      char * line;
      std::string iloc, jloc, potential_name;

      while ((line = reader.next_line())) {
        try {
          ValueTokenizer values(line);
          iloc = values.next_string();
          jloc = values.next_string();
          potential_name = values.next_string();
          if (iloc == arg[0] && jloc == arg[1] && potential_name == "stk") {

            xi_st_one = values.next_double();
            kappa_st_one = values.next_double();

            epsilon_st_one = stacking_strength(xi_st_one, kappa_st_one, T);

            a_st_one = values.next_double();
            cut_st_0_one = values.next_double();
            cut_st_c_one = values.next_double();
            cut_st_lo_one = values.next_double();
            cut_st_hi_one = values.next_double();

            a_st4_one = values.next_double();
            theta_st4_0_one = values.next_double();
            dtheta_st4_ast_one = values.next_double();
            a_st5_one = values.next_double();
            theta_st5_0_one = values.next_double();
            dtheta_st5_ast_one = values.next_double();
            a_st6_one = values.next_double();
            theta_st6_0_one = values.next_double();
            dtheta_st6_ast_one = values.next_double();
            a_st1_one = values.next_double();
            cosphi_st1_ast_one = values.next_double();
            a_st2_one = values.next_double();
            cosphi_st2_ast_one = values.next_double();

            break;
          } else continue;
        } catch (std::exception &e) {
          error->one(FLERR, "Problem parsing oxDNA potential file: {}", e.what());
        }
      }
      if ((iloc != arg[0]) || (jloc != arg[1]) || (potential_name != "stk"))
        error->one(FLERR, "No corresponding stk potential found in file {} for pair type {} {}",
                   arg[4], arg[0], arg[1]);
    }

    MPI_Bcast(&epsilon_st_one, 1, MPI_DOUBLE, 0, world);

    MPI_Bcast(&a_st_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cut_st_0_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cut_st_c_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cut_st_lo_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cut_st_hi_one, 1, MPI_DOUBLE, 0, world);

    MPI_Bcast(&a_st4_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&theta_st4_0_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&dtheta_st4_ast_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&a_st5_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&theta_st5_0_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&dtheta_st5_ast_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&a_st6_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&theta_st6_0_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&dtheta_st6_ast_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&a_st1_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cosphi_st1_ast_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&a_st2_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cosphi_st2_ast_one, 1, MPI_DOUBLE, 0, world);
  }

  b_st_lo_one = 2*a_st_one*exp(-a_st_one*(cut_st_lo_one-cut_st_0_one))*
        2*a_st_one*exp(-a_st_one*(cut_st_lo_one-cut_st_0_one))*
        (1-exp(-a_st_one*(cut_st_lo_one-cut_st_0_one)))*
        (1-exp(-a_st_one*(cut_st_lo_one-cut_st_0_one)))/
        (4*((1-exp(-a_st_one*(cut_st_lo_one -cut_st_0_one)))*
        (1-exp(-a_st_one*(cut_st_lo_one-cut_st_0_one)))-
        (1-exp(-a_st_one*(cut_st_c_one -cut_st_0_one)))*
        (1-exp(-a_st_one*(cut_st_c_one-cut_st_0_one)))));

  cut_st_lc_one = cut_st_lo_one - a_st_one*exp(-a_st_one*(cut_st_lo_one-cut_st_0_one))*
        (1-exp(-a_st_one*(cut_st_lo_one-cut_st_0_one)))/b_st_lo_one;

  b_st_hi_one = 2*a_st_one*exp(-a_st_one*(cut_st_hi_one-cut_st_0_one))*
        2*a_st_one*exp(-a_st_one*(cut_st_hi_one-cut_st_0_one))*
        (1-exp(-a_st_one*(cut_st_hi_one-cut_st_0_one)))*
        (1-exp(-a_st_one*(cut_st_hi_one-cut_st_0_one)))/
        (4*((1-exp(-a_st_one*(cut_st_hi_one -cut_st_0_one)))*
        (1-exp(-a_st_one*(cut_st_hi_one-cut_st_0_one)))-
        (1-exp(-a_st_one*(cut_st_c_one -cut_st_0_one)))*
        (1-exp(-a_st_one*(cut_st_c_one-cut_st_0_one)))));

  cut_st_hc_one = cut_st_hi_one - a_st_one*exp(-a_st_one*(cut_st_hi_one-cut_st_0_one))*
        (1-exp(-a_st_one*(cut_st_hi_one-cut_st_0_one)))/b_st_hi_one;

  tmp = 1 - exp(-(cut_st_c_one-cut_st_0_one) * a_st_one);
  shift_st_one = epsilon_st_one * tmp * tmp;

  b_st4_one = a_st4_one*a_st4_one*dtheta_st4_ast_one*dtheta_st4_ast_one/(1-a_st4_one*dtheta_st4_ast_one*dtheta_st4_ast_one);
  dtheta_st4_c_one = 1/(a_st4_one*dtheta_st4_ast_one);

  b_st5_one = a_st5_one*a_st5_one*dtheta_st5_ast_one*dtheta_st5_ast_one/(1-a_st5_one*dtheta_st5_ast_one*dtheta_st5_ast_one);
  dtheta_st5_c_one = 1/(a_st5_one*dtheta_st5_ast_one);

  b_st6_one = a_st6_one*a_st6_one*dtheta_st6_ast_one*dtheta_st6_ast_one/(1-a_st6_one*dtheta_st6_ast_one*dtheta_st6_ast_one);
  dtheta_st6_c_one = 1/(a_st6_one*dtheta_st6_ast_one);

  b_st1_one = a_st1_one*a_st1_one*cosphi_st1_ast_one*cosphi_st1_ast_one/(1-a_st1_one*cosphi_st1_ast_one*cosphi_st1_ast_one);
  cosphi_st1_c_one = 1/(a_st1_one*cosphi_st1_ast_one);

  b_st2_one = a_st2_one*a_st2_one*cosphi_st2_ast_one*cosphi_st2_ast_one/(1-a_st2_one*cosphi_st2_ast_one*cosphi_st2_ast_one);
  cosphi_st2_c_one = 1/(a_st2_one*cosphi_st2_ast_one);

  // parameters, uniform or depending on base step
  for (int i = nlo; i <= nhi; i++) {
    imod4 = i%4;
    if (imod4 == 0) imod4 = 4;

    for (int j = nlo; j <= nhi; j++) {
      jmod4 = j%4;
      if (jmod4 == 0) jmod4 = 4;

      epsilon_st[i][j] = epsilon_st_one;
      if (seqdepflag) {
        epsilon_st[i][j] *= eta_st[imod4-1][jmod4-1];
      }

      a_st[i][j] = a_st_one;
      b_st_lo[i][j] = b_st_lo_one;
      b_st_hi[i][j] = b_st_hi_one;

      theta_st4_0[i][j] = theta_st4_0_one;

      a_st5[i][j] = a_st5_one;
      theta_st5_0[i][j] = theta_st5_0_one;
      dtheta_st5_ast[i][j] = dtheta_st5_ast_one;
      b_st5[i][j] = b_st5_one;
      dtheta_st5_c[i][j] = dtheta_st5_c_one;

      a_st6[i][j] = a_st6_one;
      theta_st6_0[i][j] = theta_st6_0_one;
      dtheta_st6_ast[i][j] = dtheta_st6_ast_one;
      b_st6[i][j] = b_st6_one;
      dtheta_st6_c[i][j] = dtheta_st6_c_one;

      a_st1[i][j] = a_st1_one;
      cosphi_st1_ast[i][j] = cosphi_st1_ast_one;
      b_st1[i][j] = b_st1_one;
      cosphi_st1_c[i][j] = cosphi_st1_c_one;

      a_st2[i][j] = a_st2_one;
      cosphi_st2_ast[i][j] = cosphi_st2_ast_one;
      b_st2[i][j] = b_st2_one;
      cosphi_st2_c[i][j] = cosphi_st2_c_one;

    }
  }

  // parameters depending on dummy tetramer
  for (int i = 0; i <= nhi; i++) { // type 0 for terminal j

    for (int j = nlo; j <= nhi; j++) {
      jmod4 = j%4;
      if (jmod4 == 0) jmod4 = 4;

      for (int k = nlo; k <= nhi; k++) {
        kmod4 = k%4;
        if (kmod4 == 0) kmod4 = 4;

        for (int l = 0; l <= nhi; l++) { // type 0 for terminal k
          cut_st_0[i][j][k][l] = cut_st_0_one;
          cut_st_c[i][j][k][l] = cut_st_c_one;
          cut_st_lo[i][j][k][l] = cut_st_lo_one;
          cut_st_hi[i][j][k][l] = cut_st_hi_one;
          cut_st_lc[i][j][k][l] = cut_st_lc_one;
          cut_st_hc[i][j][k][l] = cut_st_hc_one;
          cutsq_st_hc[i][j][k][l] = cut_st_hc[i][j][k][l]*cut_st_hc[i][j][k][l];
          shift_st[i][j][k][l] = shift_st_one;
          if (seqdepflag) {
            shift_st[i][j][k][l] *= eta_st[jmod4-1][kmod4-1];
          }
          a_st4[i][j][k][l] = a_st4_one;
          dtheta_st4_ast[i][j][k][l] = dtheta_st4_ast_one;
          b_st4[i][j][k][l] = b_st4_one;
          dtheta_st4_c[i][j][k][l] = dtheta_st4_c_one;
        }
      }

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/stk" + utils::errorurl(21));

}

/* ----------------------------------------------------------------------
   atom_style hybrid bond ellipsoid oxdna required
------------------------------------------------------------------------- */

void PairOxdnaStk::init_style()
{
  if (!atom->style_match("oxdna")) {
    error->all(FLERR,"Must use 'atom_style hybrid bond ellipsoid oxdna' with pair style oxdna/stk, oxdna2/stk or oxrna2/stk");
  }

  fix_lrf = nullptr;
  auto fixes = modify->get_fix_by_style("^OXDNA/LRF");
  if (fixes.size() == 0) error->all(FLERR, "Fix OXDNA/LRF not found. Ensure pair oxdna/excv is present");
  else fix_lrf = dynamic_cast<FixOxdnaLRF *>(fixes[0]);
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use regular
------------------------------------------------------------------------- */

void PairOxdnaStk::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  if (id  > 0) error->all(FLERR,"Respa not supported");

}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairOxdnaStk::init_one(int i, int j)
{

  if (setflag[i][j] == 0) {
    error->all(FLERR,"Coefficient mixing not defined in oxDNA");
  }
  if (offset_flag) {
    error->all(FLERR,"Offset not supported in oxDNA");
  }

  // set the master list distance cutoff
  double cut_max=0.0;

  for (int a=0; a<atom->ntypes; a++) {
    for (int b=0; b<atom->ntypes; b++) {
      cut_max = MAX(cut_st_hc[a][i][j][b],cut_max);
    }
  }

  return cut_max;

}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdnaStk::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {

        fwrite(&epsilon_st[i][j],sizeof(double),1,fp);
        fwrite(&a_st[i][j],sizeof(double),1,fp);
        fwrite(&cut_st_0[i][j],sizeof(double),1,fp);
        fwrite(&cut_st_c[i][j],sizeof(double),1,fp);
        fwrite(&cut_st_lo[i][j],sizeof(double),1,fp);
        fwrite(&cut_st_hi[i][j],sizeof(double),1,fp);
        fwrite(&cut_st_lc[i][j],sizeof(double),1,fp);
        fwrite(&cut_st_hc[i][j],sizeof(double),1,fp);
        fwrite(&b_st_lo[i][j],sizeof(double),1,fp);
        fwrite(&b_st_hi[i][j],sizeof(double),1,fp);
        fwrite(&shift_st[i][j],sizeof(double),1,fp);

        fwrite(&a_st4[i][j],sizeof(double),1,fp);
        fwrite(&theta_st4_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_st4_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_st4[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_st4_c[i][j],sizeof(double),1,fp);

        fwrite(&a_st5[i][j],sizeof(double),1,fp);
        fwrite(&theta_st5_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_st5_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_st5[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_st5_c[i][j],sizeof(double),1,fp);

        fwrite(&a_st6[i][j],sizeof(double),1,fp);
        fwrite(&theta_st6_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_st6_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_st6[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_st6_c[i][j],sizeof(double),1,fp);

        fwrite(&a_st1[i][j],sizeof(double),1,fp);
        fwrite(&cosphi_st1_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_st1[i][j],sizeof(double),1,fp);
        fwrite(&cosphi_st1_c[i][j],sizeof(double),1,fp);
        fwrite(&a_st2[i][j],sizeof(double),1,fp);
        fwrite(&cosphi_st2_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_st2[i][j],sizeof(double),1,fp);
        fwrite(&cosphi_st2_c[i][j],sizeof(double),1,fp);

    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdnaStk::read_restart(FILE *fp)
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

          utils::sfread(FLERR,&epsilon_st[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&a_st[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_st_0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_st_c[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_st_lo[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_st_hi[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_st_lc[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_st_hc[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_st_lo[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_st_hi[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&shift_st[i][j],sizeof(double),1,fp,nullptr,error);

          utils::sfread(FLERR,&a_st4[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&theta_st4_0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_st4_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_st4[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_st4_c[i][j],sizeof(double),1,fp,nullptr,error);

          utils::sfread(FLERR,&a_st5[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&theta_st5_0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_st5_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_st5[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_st5_c[i][j],sizeof(double),1,fp,nullptr,error);

          utils::sfread(FLERR,&a_st6[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&theta_st6_0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_st6_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_st6[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&dtheta_st6_c[i][j],sizeof(double),1,fp,nullptr,error);

          utils::sfread(FLERR,&a_st1[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cosphi_st1_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_st1[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cosphi_st1_c[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&a_st2[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cosphi_st2_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_st2[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cosphi_st2_c[i][j],sizeof(double),1,fp,nullptr,error);

        }

        MPI_Bcast(&epsilon_st[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a_st[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_st_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_st_c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_st_lo[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_st_hi[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_st_lc[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_st_hc[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st_lo[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st_hi[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&shift_st[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_st4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_st4_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_st4_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_st4_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_st5[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_st5_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_st5_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st5[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_st5_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_st6[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_st6_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_st6_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st6[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_st6_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_st1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cosphi_st1_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cosphi_st1_c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a_st2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cosphi_st2_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cosphi_st2_c[i][j],1,MPI_DOUBLE,0,world);

      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdnaStk::write_restart_settings(FILE *fp)
{
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdnaStk::read_restart_settings(FILE *fp)
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

void *PairOxdnaStk::extract(const char *str, int &dim)
{
  dim = 2;

  if (strcmp(str,"epsilon_st") == 0) return (void *) epsilon_st;
  if (strcmp(str,"a_st") == 0) return (void *) a_st;
  if (strcmp(str,"cut_st_0") == 0) return (void *) cut_st_0;
  if (strcmp(str,"cut_st_c") == 0) return (void *) cut_st_c;
  if (strcmp(str,"cut_st_lo") == 0) return (void *) cut_st_lo;
  if (strcmp(str,"cut_st_hi") == 0) return (void *) cut_st_hi;
  if (strcmp(str,"cut_st_lc") == 0) return (void *) cut_st_lc;
  if (strcmp(str,"cut_st_hc") == 0) return (void *) cut_st_hc;
  if (strcmp(str,"b_st_lo") == 0) return (void *) b_st_lo;
  if (strcmp(str,"b_st_hi") == 0) return (void *) b_st_hi;
  if (strcmp(str,"shift_st") == 0) return (void *) shift_st;

  if (strcmp(str,"a_st4") == 0) return (void *) a_st4;
  if (strcmp(str,"theta_st4_0") == 0) return (void *) theta_st4_0;
  if (strcmp(str,"dtheta_st4_ast") == 0) return (void *) dtheta_st4_ast;
  if (strcmp(str,"b_st4") == 0) return (void *) b_st4;
  if (strcmp(str,"dtheta_st4_c") == 0) return (void *) dtheta_st4_c;

  if (strcmp(str,"a_st5") == 0) return (void *) a_st5;
  if (strcmp(str,"theta_st5_0") == 0) return (void *) theta_st5_0;
  if (strcmp(str,"dtheta_st5_ast") == 0) return (void *) dtheta_st5_ast;
  if (strcmp(str,"b_st5") == 0) return (void *) b_st5;
  if (strcmp(str,"dtheta_st5_c") == 0) return (void *) dtheta_st5_c;

  if (strcmp(str,"a_st6") == 0) return (void *) a_st6;
  if (strcmp(str,"theta_st6_0") == 0) return (void *) theta_st6_0;
  if (strcmp(str,"dtheta_st6_ast") == 0) return (void *) dtheta_st6_ast;
  if (strcmp(str,"b_st6") == 0) return (void *) b_st6;
  if (strcmp(str,"dtheta_st6_c") == 0) return (void *) dtheta_st6_c;

  if (strcmp(str,"a_st1") == 0) return (void *) a_st1;
  if (strcmp(str,"cosphi_st1_ast") == 0) return (void *) cosphi_st1_ast;
  if (strcmp(str,"b_st1") == 0) return (void *) b_st1;
  if (strcmp(str,"cosphi_st1_c") == 0) return (void *) cosphi_st1_c;

  if (strcmp(str,"a_st2") == 0) return (void *) a_st2;
  if (strcmp(str,"cosphi_st2_ast") == 0) return (void *) cosphi_st2_ast;
  if (strcmp(str,"b_st2") == 0) return (void *) b_st2;
  if (strcmp(str,"cosphi_st2_c") == 0) return (void *) cosphi_st2_c;

  return nullptr;
}
