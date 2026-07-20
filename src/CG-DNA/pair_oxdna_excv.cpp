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

#include "pair_oxdna_excv.h"
#include "constants_oxdna.h"
#include "nucleotide_oxdna.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "error.h"
#include "fix_oxdna_lrf.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "mf_oxdna.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>
#include <cassert>
#include <exception>

using namespace LAMMPS_NS;
using namespace MFOxdna;

/* ---------------------------------------------------------------------- */

PairOxdnaExcv::PairOxdnaExcv(LAMMPS *lmp) :
    Pair(lmp), epsilon_bkbk(nullptr), sigma_bkbk(nullptr), cut_bkbk_ast(nullptr),
    cutsq_bkbk_ast(nullptr), lj1_bkbk(nullptr), lj2_bkbk(nullptr), b_bkbk(nullptr),
    cut_bkbk_c(nullptr), cutsq_bkbk_c(nullptr), epsilon_bkbs(nullptr), sigma_bkbs(nullptr),
    cut_bkbs_ast(nullptr), cutsq_bkbs_ast(nullptr), lj1_bkbs(nullptr), lj2_bkbs(nullptr),
    b_bkbs(nullptr), cut_bkbs_c(nullptr), cutsq_bkbs_c(nullptr), epsilon_bsbs(nullptr),
    sigma_bsbs(nullptr), cut_bsbs_ast(nullptr), cutsq_bsbs_ast(nullptr), lj1_bsbs(nullptr),
    lj2_bsbs(nullptr), b_bsbs(nullptr), cut_bsbs_c(nullptr), cutsq_bsbs_c(nullptr),
    sigma4_bsbs(nullptr), cut4_bsbs_ast(nullptr), cut4sq_bsbs_ast(nullptr), lj14_bsbs(nullptr),
    lj24_bsbs(nullptr), b4_bsbs(nullptr), cut4_bsbs_c(nullptr), cut4sq_bsbs_c(nullptr),
    nxyz_xtrct(nullptr), fix_lrf(nullptr)
{
  single_enable = 0;
  writedata = 0;

  trim_flag = 0;
}

/* ---------------------------------------------------------------------- */

PairOxdnaExcv::~PairOxdnaExcv()
{

  if (fix_lrf) modify->delete_fix(fix_lrf->id);

  if (allocated) {

    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(epsilon_bkbk);
    memory->destroy(sigma_bkbk);
    memory->destroy(cut_bkbk_ast);
    memory->destroy(b_bkbk);
    memory->destroy(cut_bkbk_c);
    memory->destroy(lj1_bkbk);
    memory->destroy(lj2_bkbk);
    memory->destroy(cutsq_bkbk_ast);
    memory->destroy(cutsq_bkbk_c);

    memory->destroy(epsilon_bkbs);
    memory->destroy(sigma_bkbs);
    memory->destroy(cut_bkbs_ast);
    memory->destroy(b_bkbs);
    memory->destroy(cut_bkbs_c);
    memory->destroy(lj1_bkbs);
    memory->destroy(lj2_bkbs);
    memory->destroy(cutsq_bkbs_ast);
    memory->destroy(cutsq_bkbs_c);

    memory->destroy(epsilon_bsbs);
    memory->destroy(sigma_bsbs);
    memory->destroy(cut_bsbs_ast);
    memory->destroy(b_bsbs);
    memory->destroy(cut_bsbs_c);
    memory->destroy(lj1_bsbs);
    memory->destroy(lj2_bsbs);
    memory->destroy(cutsq_bsbs_ast);
    memory->destroy(cutsq_bsbs_c);

    memory->destroy(sigma4_bsbs);
    memory->destroy(cut4_bsbs_ast);
    memory->destroy(b4_bsbs);
    memory->destroy(cut4_bsbs_c);
    memory->destroy(lj14_bsbs);
    memory->destroy(lj24_bsbs);
    memory->destroy(cut4sq_bsbs_ast);
    memory->destroy(cut4sq_bsbs_c);
  }
}

/* ----------------------------------------------------------------------
    compute vector COM-sugar-phosphate backbone interaction site in oxDNA
------------------------------------------------------------------------- */
inline void PairOxdnaExcv::compute_backbone_site(double e1[3], double /*e2*/[3],
    double /*e3*/[3], double rbk[3]) const
{
  NucleotideOxdna1 oxdna1;
  oxdna1.backbone_site(e1, nullptr, nullptr, rbk);
}

/* ---------------------------------------------------------------------
    compute vector COM-base site in oxDNA/oxDNA2
    identical templates for A=1, C=2, G=3, T=0
------------------------------------------------------------------------ */
inline void PairOxdnaExcv::compute_base_site(int /*type*/, double e1[3],
  double /*e2*/[3], double /*e3*/[3], double rbs[3]) const
{
  NucleotideOxdna1 oxdna1;
  oxdna1.base_site<0>(e1, nullptr, nullptr, rbs);
}

/* ----------------------------------------------------------------------
   compute function for oxDNA pair interactions
   s=sugar-phosphate backbone site, b=base site, st=stacking site
------------------------------------------------------------------------- */

void PairOxdnaExcv::compute(int eflag, int vflag)
{
  double delf[3],delta[3],deltb[3]; // force, torque increment;
  double evdwl,fpair,factor_lj;
  double rtmp_bk[3],rtmp_bs[3];
  double delr_bkbk[3],rsq_bkbk,delr_bkbs[3],rsq_bkbs;
  double delr_bsbk[3],rsq_bsbk,delr_bsbs[3],rsq_bsbs;

  // vectors COM-backbone site, COM-base site in lab frame
  double ra_cbk[3],ra_cbs[3];
  double rb_cbk[3],rb_cbs[3];
  // Cartesian unit vectors in lab frame
  double ax[3],ay[3],az[3];
  double bx[3],by[3],bz[3];

  double *special_lj = force->special_lj;

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  int *type = atom->type;

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int *alist,*blist,*numneigh,**firstneigh;

  int a,b,ia,ib,anum,bnum,atype,btype;
  tagint *id3p = atom->id3p;
  tagint *id5p = atom->id5p;
  int _3ptype,_5ptype;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  anum = list->inum;
  alist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // nxyz_xtrct = extracted local unit vectors in lab frame from fix oxdna/lrf
  nxyz_xtrct = fix_lrf->array_atom;

  // loop over pair interaction neighbors of my atoms

  for (ia = 0; ia < anum; ia++) {

    a = alist[ia];
    atype = type[a];

    ax[0] = nxyz_xtrct[a][0];
    ax[1] = nxyz_xtrct[a][1];
    ax[2] = nxyz_xtrct[a][2];
    ay[0] = nxyz_xtrct[a][3];
    ay[1] = nxyz_xtrct[a][4];
    ay[2] = nxyz_xtrct[a][5];
    az[0] = nxyz_xtrct[a][6];
    az[1] = nxyz_xtrct[a][7];
    az[2] = nxyz_xtrct[a][8];

    // vector COM - backbone site a
    compute_backbone_site(ax,ay,az,ra_cbk);

    // vector COM - base site a
    compute_base_site(atype%4, ax,ay,az,ra_cbs);

    rtmp_bk[0] = x[a][0] + ra_cbk[0];
    rtmp_bk[1] = x[a][1] + ra_cbk[1];
    rtmp_bk[2] = x[a][2] + ra_cbk[2];

    rtmp_bs[0] = x[a][0] + ra_cbs[0];
    rtmp_bs[1] = x[a][1] + ra_cbs[1];
    rtmp_bs[2] = x[a][2] + ra_cbs[2];

    blist = firstneigh[a];
    bnum = numneigh[a];

    for (ib = 0; ib < bnum; ib++) {

      b = blist[ib];
      factor_lj = special_lj[sbmask(b)]; // = 0 for nearest neighbors
      b &= NEIGHMASK;

      btype = type[b];

      bx[0] = nxyz_xtrct[b][0];
      bx[1] = nxyz_xtrct[b][1];
      bx[2] = nxyz_xtrct[b][2];
      by[0] = nxyz_xtrct[b][3];
      by[1] = nxyz_xtrct[b][4];
      by[2] = nxyz_xtrct[b][5];
      bz[0] = nxyz_xtrct[b][6];
      bz[1] = nxyz_xtrct[b][7];
      bz[2] = nxyz_xtrct[b][8];

      // vector COM - backbone site b
      compute_backbone_site(bx,by,bz,rb_cbk);

      // vector COM - base site b
      compute_base_site(btype%4, bx,by,bz,rb_cbs);

      // vector backbone site b to a
      delr_bkbk[0] = rtmp_bk[0] - (x[b][0] + rb_cbk[0]);
      delr_bkbk[1] = rtmp_bk[1] - (x[b][1] + rb_cbk[1]);
      delr_bkbk[2] = rtmp_bk[2] - (x[b][2] + rb_cbk[2]);
      rsq_bkbk = delr_bkbk[0]*delr_bkbk[0] + delr_bkbk[1]*delr_bkbk[1] + delr_bkbk[2]*delr_bkbk[2];

      // vector base site b to backbone site a
      delr_bkbs[0] =  rtmp_bk[0] - (x[b][0] + rb_cbs[0]);
      delr_bkbs[1] =  rtmp_bk[1] - (x[b][1] + rb_cbs[1]);
      delr_bkbs[2] =  rtmp_bk[2] - (x[b][2] + rb_cbs[2]);
      rsq_bkbs = delr_bkbs[0]*delr_bkbs[0] + delr_bkbs[1]*delr_bkbs[1] + delr_bkbs[2]*delr_bkbs[2];

      // vector backbone site b to base site a
      delr_bsbk[0] = rtmp_bs[0] - (x[b][0] + rb_cbk[0]);
      delr_bsbk[1] = rtmp_bs[1] - (x[b][1] + rb_cbk[1]);
      delr_bsbk[2] = rtmp_bs[2] - (x[b][2] + rb_cbk[2]);
      rsq_bsbk = delr_bsbk[0]*delr_bsbk[0] + delr_bsbk[1]*delr_bsbk[1] + delr_bsbk[2]*delr_bsbk[2];

      // vector base site b to a
      delr_bsbs[0] = rtmp_bs[0] - (x[b][0] + rb_cbs[0]);
      delr_bsbs[1] = rtmp_bs[1] - (x[b][1] + rb_cbs[1]);
      delr_bsbs[2] = rtmp_bs[2] - (x[b][2] + rb_cbs[2]);
      rsq_bsbs = delr_bsbs[0]*delr_bsbs[0] + delr_bsbs[1]*delr_bsbs[1] + delr_bsbs[2]*delr_bsbs[2];

      // excluded volume interaction

      // backbone-backbone

      // interaction coefficients depend on base step
      if (rsq_bkbk < cutsq_bkbk_c[atype][btype]) {

        evdwl = F3(rsq_bkbk,cutsq_bkbk_ast[atype][btype],cut_bkbk_c[atype][btype],lj1_bkbk[atype][btype],
                        lj2_bkbk[atype][btype],epsilon_bkbk[atype][btype],b_bkbk[atype][btype],fpair);

        // knock out nearest-neighbor interaction between backbone sites on same strand
        fpair *= factor_lj;
        evdwl *= factor_lj;

        delf[0] = delr_bkbk[0]*fpair;
        delf[1] = delr_bkbk[1]*fpair;
        delf[2] = delr_bkbk[2]*fpair;

        // increment energy and virial
        // NOTE: The virial is calculated on the 'molecular' basis.
        // (see G. Ciccotti and J.P. Ryckaert, Comp. Phys. Rep. 4, 345-392 (1986))

        if (evflag) ev_tally_xyz(a,b,nlocal,newton_pair,evdwl,0.0,
            delf[0],delf[1],delf[2],x[a][0]-x[b][0],x[a][1]-x[b][1],x[a][2]-x[b][2]);

        f[a][0] += delf[0];
        f[a][1] += delf[1];
        f[a][2] += delf[2];

        MathExtra::cross3(ra_cbk,delf,delta);

        torque[a][0] += delta[0];
        torque[a][1] += delta[1];
        torque[a][2] += delta[2];

        if (newton_pair || b < nlocal) {

          f[b][0] -= delf[0];
          f[b][1] -= delf[1];
          f[b][2] -= delf[2];

          MathExtra::cross3(rb_cbk,delf,deltb);

          torque[b][0] -= deltb[0];
          torque[b][1] -= deltb[1];
          torque[b][2] -= deltb[2];

        }
      }

      // backbone-base

      if (rsq_bkbs < cutsq_bkbs_c[atype][btype]) {
        evdwl = F3(rsq_bkbs,cutsq_bkbs_ast[atype][btype],cut_bkbs_c[atype][btype],lj1_bkbs[atype][btype],
                        lj2_bkbs[atype][btype],epsilon_bkbs[atype][btype],b_bkbs[atype][btype],fpair);

        delf[0] = delr_bkbs[0]*fpair;
        delf[1] = delr_bkbs[1]*fpair;
        delf[2] = delr_bkbs[2]*fpair;

        // increment energy and virial
        if (evflag) ev_tally_xyz(a,b,nlocal,newton_pair,evdwl,0.0,
            delf[0],delf[1],delf[2],x[a][0]-x[b][0],x[a][1]-x[b][1],x[a][2]-x[b][2]);

        f[a][0] += delf[0];
        f[a][1] += delf[1];
        f[a][2] += delf[2];

        MathExtra::cross3(ra_cbk,delf,delta);

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
      }

      // base-backbone

      if (rsq_bsbk < cutsq_bkbs_c[atype][btype]) {
        evdwl = F3(rsq_bsbk,cutsq_bkbs_ast[atype][btype],cut_bkbs_c[atype][btype],lj1_bkbs[atype][btype],
                       lj2_bkbs[atype][btype],epsilon_bkbs[atype][btype],b_bkbs[atype][btype],fpair);

        delf[0] = delr_bsbk[0]*fpair;
        delf[1] = delr_bsbk[1]*fpair;
        delf[2] = delr_bsbk[2]*fpair;

        // increment energy and virial
        if (evflag) ev_tally_xyz(a,b,nlocal,newton_pair,evdwl,0.0,
            delf[0],delf[1],delf[2],x[a][0]-x[b][0],x[a][1]-x[b][1],x[a][2]-x[b][2]);

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

          MathExtra::cross3(rb_cbk,delf,deltb);

          torque[b][0] -= deltb[0];
          torque[b][1] -= deltb[1];
          torque[b][2] -= deltb[2];

        }
      }

      // base-base

      evdwl = 0.0;

      // bond-topology: a-b nearest-neighbors on same strand, arguments depends on tetramer
      if ((atom->tag[a] == id3p[b]) && (atom->tag[b] == id5p[a])) { // a -> b is 3' -> 5'

        // determine type of 3'-partner of a and 5'-partner of b
        if (id3p[a] != -1) {
          _3ptype = type[atom->map(id3p[a])];
        }
        else _3ptype = 0;

        if (id5p[b] != -1) {
          _5ptype = type[atom->map(id5p[b])];
        }
        else _5ptype = 0;

        if (rsq_bsbs < cut4sq_bsbs_c[_3ptype][atype][btype][_5ptype]) {
          evdwl = F3(rsq_bsbs,cut4sq_bsbs_ast[_3ptype][atype][btype][_5ptype],cut4_bsbs_c[_3ptype][atype][btype][_5ptype],
                          lj14_bsbs[_3ptype][atype][btype][_5ptype],lj24_bsbs[_3ptype][atype][btype][_5ptype],
                          epsilon_bsbs[atype][btype],b4_bsbs[_3ptype][atype][btype][_5ptype],fpair);
        }

      }
      // bond-topology: a-b nearest-neighbors on same strand, arguments depends on tetramer
      else if ((atom->tag[a] == id5p[b]) && (atom->tag[b] == id3p[a])) { // b -> a is 3' -> 5'

        // determine type of 3'-partner of b and 5'-partner of a
        if (id3p[b] != -1) {
          _3ptype = type[atom->map(id3p[b])];
        }
        else _3ptype = 0;

        if (id5p[a] != -1) {
          _5ptype = type[atom->map(id5p[a])];
        }
        else _5ptype = 0;

        if (rsq_bsbs < cut4sq_bsbs_c[_3ptype][btype][atype][_5ptype]) {
          evdwl = F3(rsq_bsbs,cut4sq_bsbs_ast[_3ptype][btype][atype][_5ptype],cut4_bsbs_c[_3ptype][btype][atype][_5ptype],
                          lj14_bsbs[_3ptype][btype][atype][_5ptype],lj24_bsbs[_3ptype][btype][atype][_5ptype],
                          epsilon_bsbs[btype][atype],b4_bsbs[_3ptype][btype][atype][_5ptype],fpair);
        }
      }
      else {
        if (rsq_bsbs < cutsq_bsbs_c[atype][btype]) {
          evdwl = F3(rsq_bsbs,cutsq_bsbs_ast[atype][btype],cut_bsbs_c[atype][btype],lj1_bsbs[atype][btype],
                          lj2_bsbs[atype][btype],epsilon_bsbs[atype][btype],b_bsbs[atype][btype],fpair);
        }
      }

      if (evdwl != 0.0) {

        delf[0] = delr_bsbs[0]*fpair;
        delf[1] = delr_bsbs[1]*fpair;
        delf[2] = delr_bsbs[2]*fpair;

        // increment energy and virial
        if (evflag) ev_tally_xyz(a,b,nlocal,newton_pair,evdwl,0.0,
            delf[0],delf[1],delf[2],x[a][0]-x[b][0],x[a][1]-x[b][1],x[a][2]-x[b][2]);

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

      }
      // end excluded volume interaction

    }
  }

  if (vflag_fdotr) virial_fdotr_compute();

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairOxdnaExcv::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(epsilon_bkbk,n+1,n+1,"pair:epsilon_bkbk");
  memory->create(sigma_bkbk,n+1,n+1,"pair:sigma_bkbk");
  memory->create(cut_bkbk_ast,n+1,n+1,"pair:cut_bkbk_ast");
  memory->create(b_bkbk,n+1,n+1,"pair:b_bkbk");
  memory->create(cut_bkbk_c,n+1,n+1,"pair:cut_bkbk_c");
  memory->create(lj1_bkbk,n+1,n+1,"pair:lj1_bkbk");
  memory->create(lj2_bkbk,n+1,n+1,"pair:lj2_bkbk");
  memory->create(cutsq_bkbk_ast,n+1,n+1,"pair:cutsq_bkbk_ast");
  memory->create(cutsq_bkbk_c,n+1,n+1,"pair:cutsq_bkbk_c");

  memory->create(epsilon_bkbs,n+1,n+1,"pair:epsilon_bkbs");
  memory->create(sigma_bkbs,n+1,n+1,"pair:sigma_bkbs");
  memory->create(cut_bkbs_ast,n+1,n+1,"pair:cut_bkbs_ast");
  memory->create(b_bkbs,n+1,n+1,"pair:b_bkbs");
  memory->create(cut_bkbs_c,n+1,n+1,"pair:cut_bkbs_c");
  memory->create(lj1_bkbs,n+1,n+1,"pair:lj1_bkbs");
  memory->create(lj2_bkbs,n+1,n+1,"pair:lj2_bkbs");
  memory->create(cutsq_bkbs_ast,n+1,n+1,"pair:cutsq_bkbs_ast");
  memory->create(cutsq_bkbs_c,n+1,n+1,"pair:cutsq_bkbs_c");

  memory->create(epsilon_bsbs,n+1,n+1,"pair:epsilon_bsbs");
  memory->create(sigma_bsbs,n+1,n+1,"pair:sigma_bsbs");
  memory->create(cut_bsbs_ast,n+1,n+1,"pair:cut_bsbs_ast");
  memory->create(b_bsbs,n+1,n+1,"pair:b_bsbs");
  memory->create(cut_bsbs_c,n+1,n+1,"pair:cut_bsbs_c");
  memory->create(lj1_bsbs,n+1,n+1,"pair:lj1_bsbs");
  memory->create(lj2_bsbs,n+1,n+1,"pair:lj2_bsbs");
  memory->create(cutsq_bsbs_ast,n+1,n+1,"pair:cutsq_bsbs_ast");
  memory->create(cutsq_bsbs_c,n+1,n+1,"pair:cutsq_bsbs_c");

  memory->create(sigma4_bsbs,n+1,n+1,n+1,n+1,"pair:sigma4_bsbs");
  memory->create(cut4_bsbs_ast,n+1,n+1,n+1,n+1,"pair:cut4_bsbs_ast");
  memory->create(b4_bsbs,n+1,n+1,n+1,n+1,"pair:b4_bsbs");
  memory->create(cut4_bsbs_c,n+1,n+1,n+1,n+1,"pair:cut4_bsbs_c");
  memory->create(lj14_bsbs,n+1,n+1,n+1,n+1,"pair:lj14_bsbs");
  memory->create(lj24_bsbs,n+1,n+1,n+1,n+1,"pair:lj24_bsbs");
  memory->create(cut4sq_bsbs_ast,n+1,n+1,n+1,n+1,"pair:cut4sq_bsbs_ast");
  memory->create(cut4sq_bsbs_c,n+1,n+1,n+1,n+1,"pair:cut4sq_bsbs_c");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairOxdnaExcv::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairOxdnaExcv::coeff(int narg, char **arg)
{
  int count;

  if (narg != 3 && narg != 11) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/excv" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi,nlo,nhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  assert((ilo == jlo) & (ihi == jhi));
  nlo = ilo;
  nhi = ihi;

  double epsilon_bkbk_one, sigma_bkbk_one;
  double cut_bkbk_ast_one, cut_bkbk_c_one, b_bkbk_one;

  double epsilon_bkbs_one, sigma_bkbs_one;
  double cut_bkbs_ast_one, cut_bkbs_c_one, b_bkbs_one;

  double epsilon_bsbs_one, sigma_bsbs_one;
  double cut_bsbs_ast_one, cut_bsbs_c_one, b_bsbs_one;

  if (narg == 11) {
    // Excluded volume interaction
    // LJ parameters
    epsilon_bkbk_one = utils::numeric(FLERR,arg[2],false,lmp);
    sigma_bkbk_one = utils::numeric(FLERR,arg[3],false,lmp);
    cut_bkbk_ast_one = utils::numeric(FLERR,arg[4],false,lmp);

    // LJ parameters
    epsilon_bkbs_one = utils::numeric(FLERR,arg[5],false,lmp);
    sigma_bkbs_one = utils::numeric(FLERR,arg[6],false,lmp);
    cut_bkbs_ast_one = utils::numeric(FLERR,arg[7],false,lmp);

    // LJ parameters
    epsilon_bsbs_one = utils::numeric(FLERR,arg[8],false,lmp);
    sigma_bsbs_one = utils::numeric(FLERR,arg[9],false,lmp);
    cut_bsbs_ast_one = utils::numeric(FLERR,arg[10],false,lmp);
  }
  else {
    if (comm->me == 0) {
      PotentialFileReader reader(lmp, arg[2], "oxdna potential", " (excv)");
      char * line;
      std::string iloc, jloc, potential_name;

      while ((line = reader.next_line())) {
        try {
          ValueTokenizer values(line);
          iloc = values.next_string();
          jloc = values.next_string();
          potential_name = values.next_string();
          if (iloc == arg[0] && jloc == arg[1] && potential_name == "excv") {
            // Excluded volume interaction
            // LJ backbone-backbone parameters
            epsilon_bkbk_one = values.next_double();
            sigma_bkbk_one = values.next_double();
            cut_bkbk_ast_one = values.next_double();

            // LJ backbone-base parameters
            epsilon_bkbs_one = values.next_double();
            sigma_bkbs_one = values.next_double();
            cut_bkbs_ast_one = values.next_double();

            // LJ base-base parameters
            epsilon_bsbs_one = values.next_double();
            sigma_bsbs_one = values.next_double();
            cut_bsbs_ast_one = values.next_double();

            break;
          } else continue;
        } catch (std::exception &e) {
          error->one(FLERR, "Problem parsing oxdna potential file: {}", e.what());
        }
      }
      if ((iloc != arg[0]) || (jloc != arg[1]) || (potential_name != "excv"))
        error->one(FLERR, "No corresponding excv potential found in file {} for pair type {} {}",
                   arg[2], arg[0], arg[1]);
    }

    MPI_Bcast(&epsilon_bkbk_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&sigma_bkbk_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cut_bkbk_ast_one, 1, MPI_DOUBLE, 0, world);

    MPI_Bcast(&epsilon_bkbs_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&sigma_bkbs_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cut_bkbs_ast_one, 1, MPI_DOUBLE, 0, world);

    MPI_Bcast(&epsilon_bsbs_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&sigma_bsbs_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cut_bsbs_ast_one, 1, MPI_DOUBLE, 0, world);
  }

  // backbone-backbone
  count = 0;

  // smoothing - determined through continuity and differentiability
  b_bkbk_one = 4.0/sigma_bkbk_one
      *(6.0*pow(sigma_bkbk_one/cut_bkbk_ast_one,7)-12.0*pow(sigma_bkbk_one/cut_bkbk_ast_one,13))
      *4.0/sigma_bkbk_one*(6.0*pow(sigma_bkbk_one/cut_bkbk_ast_one,7)-12.0*pow(sigma_bkbk_one/cut_bkbk_ast_one,13))
      /4.0/(4.0*(pow(sigma_bkbk_one/cut_bkbk_ast_one,12)-pow(sigma_bkbk_one/cut_bkbk_ast_one,6)));

  cut_bkbk_c_one = cut_bkbk_ast_one
      - 2.0*4.0*(pow(sigma_bkbk_one/cut_bkbk_ast_one,12)-pow(sigma_bkbk_one/cut_bkbk_ast_one,6))
      /(4.0/sigma_bkbk_one*(6.0*pow(sigma_bkbk_one/cut_bkbk_ast_one,7)-12.0*pow(sigma_bkbk_one/cut_bkbk_ast_one,13)));

  // backbone-backbone parameters depending on base step
  for (int i = nlo; i <= nhi; i++) {
    for (int j = nlo; j <= nhi; j++) {
      epsilon_bkbk[i][j] = epsilon_bkbk_one;
      sigma_bkbk[i][j] = sigma_bkbk_one;
      cut_bkbk_ast[i][j] = cut_bkbk_ast_one;
      b_bkbk[i][j] = b_bkbk_one;
      cut_bkbk_c[i][j] = cut_bkbk_c_one;
      lj1_bkbk[i][j] = 4.0 * epsilon_bkbk[i][j] * pow(sigma_bkbk[i][j],12.0);
      lj2_bkbk[i][j] = 4.0 * epsilon_bkbk[i][j] * pow(sigma_bkbk[i][j],6.0);
      cutsq_bkbk_ast[i][j] = cut_bkbk_ast[i][j]*cut_bkbk_ast[i][j];
      cutsq_bkbk_c[i][j]  = cut_bkbk_c[i][j]*cut_bkbk_c[i][j];
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/excv" + utils::errorurl(21));

  // backbone-base
  count = 0;

  // smoothing - determined through continuity and differentiability
  b_bkbs_one = 4.0/sigma_bkbs_one
      *(6.0*pow(sigma_bkbs_one/cut_bkbs_ast_one,7)-12.0*pow(sigma_bkbs_one/cut_bkbs_ast_one,13))
      *4.0/sigma_bkbs_one*(6.0*pow(sigma_bkbs_one/cut_bkbs_ast_one,7)-12.0*pow(sigma_bkbs_one/cut_bkbs_ast_one,13))
      /4.0/(4.0*(pow(sigma_bkbs_one/cut_bkbs_ast_one,12)-pow(sigma_bkbs_one/cut_bkbs_ast_one,6)));

  cut_bkbs_c_one = cut_bkbs_ast_one
      - 2.0*4.0*(pow(sigma_bkbs_one/cut_bkbs_ast_one,12)-pow(sigma_bkbs_one/cut_bkbs_ast_one,6))
      /(4.0/sigma_bkbs_one*(6.0*pow(sigma_bkbs_one/cut_bkbs_ast_one,7)-12.0*pow(sigma_bkbs_one/cut_bkbs_ast_one,13)));

  // backbone-base parameters depending on base step
  for (int i = nlo; i <= nhi; i++) {
    for (int j = nlo; j <= nhi; j++) {
      epsilon_bkbs[i][j] = epsilon_bkbs_one;
      sigma_bkbs[i][j] = sigma_bkbs_one;
      cut_bkbs_ast[i][j] = cut_bkbs_ast_one;
      b_bkbs[i][j] = b_bkbs_one;
      cut_bkbs_c[i][j] = cut_bkbs_c_one;
      lj1_bkbs[i][j] = 4.0 * epsilon_bkbs[i][j] * pow(sigma_bkbs[i][j],12.0);
      lj2_bkbs[i][j] = 4.0 * epsilon_bkbs[i][j] * pow(sigma_bkbs[i][j],6.0);
      cutsq_bkbs_ast[i][j] = cut_bkbs_ast[i][j]*cut_bkbs_ast[i][j];
      cutsq_bkbs_c[i][j]  = cut_bkbs_c[i][j]*cut_bkbs_c[i][j];
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/excv" + utils::errorurl(21));

  // base-base
  count = 0;

  // smoothing - determined through continuity and differentiability
  b_bsbs_one = 4.0/sigma_bsbs_one
      *(6.0*pow(sigma_bsbs_one/cut_bsbs_ast_one,7)-12.0*pow(sigma_bsbs_one/cut_bsbs_ast_one,13))
      *4.0/sigma_bsbs_one*(6.0*pow(sigma_bsbs_one/cut_bsbs_ast_one,7)-12.0*pow(sigma_bsbs_one/cut_bsbs_ast_one,13))
      /4.0/(4.0*(pow(sigma_bsbs_one/cut_bsbs_ast_one,12)-pow(sigma_bsbs_one/cut_bsbs_ast_one,6)));

  cut_bsbs_c_one = cut_bsbs_ast_one
      - 2.0*4.0*(pow(sigma_bsbs_one/cut_bsbs_ast_one,12)-pow(sigma_bsbs_one/cut_bsbs_ast_one,6))
      /(4.0/sigma_bsbs_one*(6.0*pow(sigma_bsbs_one/cut_bsbs_ast_one,7)-12.0*pow(sigma_bsbs_one/cut_bsbs_ast_one,13)));

  // base-base parameters depending on base step
  for (int i = nlo; i <= nhi; i++) {
    for (int j = nlo; j <= nhi; j++) {
      epsilon_bsbs[i][j] = epsilon_bsbs_one;
      sigma_bsbs[i][j] = sigma_bsbs_one;
      cut_bsbs_ast[i][j] = cut_bsbs_ast_one;
      b_bsbs[i][j] = b_bsbs_one;
      cut_bsbs_c[i][j] = cut_bsbs_c_one;
      lj1_bsbs[i][j] = 4.0 * epsilon_bsbs[i][j] * pow(sigma_bsbs[i][j],12.0);
      lj2_bsbs[i][j] = 4.0 * epsilon_bsbs[i][j] * pow(sigma_bsbs[i][j],6.0);
      cutsq_bsbs_ast[i][j] = cut_bsbs_ast[i][j]*cut_bsbs_ast[i][j];
      cutsq_bsbs_c[i][j]  = cut_bsbs_c[i][j]*cut_bsbs_c[i][j];
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/excv");

  // base-base parameters depending on tetramer
  count = 0;

  for (int i = 0; i <= nhi; i++) { // type 0 for terminal j
    for (int j = nlo; j <= nhi; j++) {
      for (int k = nlo; k <= nhi; k++) {
        for (int l = 0; l <= nhi; l++) { // type 0 for terminal k
          sigma4_bsbs[i][j][k][l] = sigma_bsbs_one;
          cut4_bsbs_ast[i][j][k][l] = cut_bsbs_ast_one;
          b4_bsbs[i][j][k][l] = b_bsbs_one;
          cut4_bsbs_c[i][j][k][l] = cut_bsbs_c_one;
          cut4sq_bsbs_ast[i][j][k][l] = cut4_bsbs_ast[i][j][k][l]*cut4_bsbs_ast[i][j][k][l];
          cut4sq_bsbs_c[i][j][k][l]  = cut4_bsbs_c[i][j][k][l]*cut4_bsbs_c[i][j][k][l];
          lj14_bsbs[i][j][k][l] = 4.0 * epsilon_bsbs[j][k] * pow(sigma4_bsbs[i][j][k][l],12.0);
          lj24_bsbs[i][j][k][l] = 4.0 * epsilon_bsbs[j][k] * pow(sigma4_bsbs[i][j][k][l],6.0);
          count++;
       }
      }
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/excv");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */
void PairOxdnaExcv::init_style()
{
  // ensure fix oxdna/lrf is added for backward-compatability
  if (!fix_lrf)
    fix_lrf = dynamic_cast<FixOxdnaLRF *>(modify->add_fix("lrf all OXDNA/LRF"));

  neighbor->add_request(this, NeighConst::REQ_DEFAULT);
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use regular
------------------------------------------------------------------------- */

void PairOxdnaExcv::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  if (id  > 0) error->all(FLERR,"Respa not supported");

}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairOxdnaExcv::init_one(int i, int j)
{

  if (setflag[i][j] == 0) {
    error->all(FLERR,"Coefficient mixing not defined in oxdna");
  }
  if (offset_flag) {
    error->all(FLERR,"Offset not supported in oxdna");
  }

  // set the master list distance cutoff
  return cut_bkbk_c[i][j];

}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdnaExcv::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {

        fwrite(&epsilon_bkbk[i][j],sizeof(double),1,fp);
        fwrite(&sigma_bkbk[i][j],sizeof(double),1,fp);
        fwrite(&cut_bkbk_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_bkbk[i][j],sizeof(double),1,fp);
        fwrite(&cut_bkbk_c[i][j],sizeof(double),1,fp);
        fwrite(&epsilon_bkbs[i][j],sizeof(double),1,fp);
        fwrite(&sigma_bkbs[i][j],sizeof(double),1,fp);
        fwrite(&cut_bkbs_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_bkbs[i][j],sizeof(double),1,fp);
        fwrite(&cut_bkbs_c[i][j],sizeof(double),1,fp);
        fwrite(&epsilon_bsbs[i][j],sizeof(double),1,fp);
        fwrite(&sigma_bsbs[i][j],sizeof(double),1,fp);
        fwrite(&cut_bsbs_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_bsbs[i][j],sizeof(double),1,fp);
        fwrite(&cut_bsbs_c[i][j],sizeof(double),1,fp);

    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdnaExcv::read_restart(FILE *fp)
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

          utils::sfread(FLERR,&epsilon_bkbk[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&sigma_bkbk[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_bkbk_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_bkbk[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_bkbk_c[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&epsilon_bkbs[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&sigma_bkbs[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_bkbs_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_bkbs[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_bkbs_c[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&epsilon_bsbs[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&sigma_bsbs[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_bsbs_ast[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b_bsbs[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_bsbs_c[i][j],sizeof(double),1,fp,nullptr,error);

         }

        MPI_Bcast(&epsilon_bkbk[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma_bkbk[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_bkbk_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_bkbk[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_bkbk_c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&epsilon_bkbs[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma_bkbs[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_bkbs_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_bkbs[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_bkbs_c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&epsilon_bsbs[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma_bsbs[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_bsbs_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_bsbs[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_bsbs_c[i][j],1,MPI_DOUBLE,0,world);

      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdnaExcv::write_restart_settings(FILE *fp)
{
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdnaExcv::read_restart_settings(FILE *fp)
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

void *PairOxdnaExcv::extract(const char *str, int &dim)
{
  dim = 2;

  if (strcmp(str,"epsilon_bkbk") == 0) return (void *) epsilon_bkbk;
  if (strcmp(str,"sigma_bkbk") == 0) return (void *) sigma_bkbk;
  if (strcmp(str,"cut_bkbk_ast") == 0) return (void *) cut_bkbk_ast;
  if (strcmp(str,"b_bkbk") == 0) return (void *) b_bkbk;
  if (strcmp(str,"cut_bkbk_c") == 0) return (void *) cut_bkbk_c;

  if (strcmp(str,"epsilon_bkbs") == 0) return (void *) epsilon_bkbs;
  if (strcmp(str,"sigma_bkbs") == 0) return (void *) sigma_bkbs;
  if (strcmp(str,"cut_bkbs_ast") == 0) return (void *) cut_bkbs_ast;
  if (strcmp(str,"b_bkbs") == 0) return (void *) b_bkbs;
  if (strcmp(str,"cut_bkbs_c") == 0) return (void *) cut_bkbs_c;
  if (strcmp(str,"sigma4_bkbs") == 0) return (void *) sigma_bkbs;

  if (strcmp(str,"epsilon_bsbs") == 0) return (void *) epsilon_bsbs;
  if (strcmp(str,"sigma_bsbs") == 0) return (void *) sigma_bsbs;
  if (strcmp(str,"cut_bsbs_ast") == 0) return (void *) cut_bsbs_ast;
  if (strcmp(str,"b_bsbs") == 0) return (void *) b_bsbs;
  if (strcmp(str,"cut_bsbs_c") == 0) return (void *) cut_bsbs_c;

  return nullptr;
}
