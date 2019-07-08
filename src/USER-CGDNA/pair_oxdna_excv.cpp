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

#include "pair_oxdna_excv.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "mf_oxdna.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "atom_vec_ellipsoid.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace MFOxdna;

/* ---------------------------------------------------------------------- */

PairOxdnaExcv::PairOxdnaExcv(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairOxdnaExcv::~PairOxdnaExcv()
{
  if (allocated) {

    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(epsilon_ss);
    memory->destroy(sigma_ss);
    memory->destroy(cut_ss_ast);
    memory->destroy(b_ss);
    memory->destroy(cut_ss_c);
    memory->destroy(lj1_ss);
    memory->destroy(lj2_ss);
    memory->destroy(cutsq_ss_ast);
    memory->destroy(cutsq_ss_c);

    memory->destroy(epsilon_sb);
    memory->destroy(sigma_sb);
    memory->destroy(cut_sb_ast);
    memory->destroy(b_sb);
    memory->destroy(cut_sb_c);
    memory->destroy(lj1_sb);
    memory->destroy(lj2_sb);
    memory->destroy(cutsq_sb_ast);
    memory->destroy(cutsq_sb_c);

    memory->destroy(epsilon_bb);
    memory->destroy(sigma_bb);
    memory->destroy(cut_bb_ast);
    memory->destroy(b_bb);
    memory->destroy(cut_bb_c);
    memory->destroy(lj1_bb);
    memory->destroy(lj2_bb);
    memory->destroy(cutsq_bb_ast);
    memory->destroy(cutsq_bb_c);

  }
}

/* ----------------------------------------------------------------------
    compute vector COM-excluded volume interaction sites in oxDNA
------------------------------------------------------------------------- */
void PairOxdnaExcv::compute_interaction_sites(double e1[3],
  double /*e2*/[3], double rs[3], double rb[3])
{
  double d_cs=-0.4, d_cb=+0.4;

  rs[0] = d_cs*e1[0];
  rs[1] = d_cs*e1[1];
  rs[2] = d_cs*e1[2];

  rb[0] = d_cb*e1[0];
  rb[1] = d_cb*e1[1];
  rb[2] = d_cb*e1[2];

}

/* ----------------------------------------------------------------------
   compute function for oxDNA pair interactions
   s=sugar-phosphate backbone site, b=base site, st=stacking site
------------------------------------------------------------------------- */

void PairOxdnaExcv::compute(int eflag, int vflag)
{

  double delf[3],delta[3],deltb[3]; // force, torque increment;
  double evdwl,fpair,factor_lj;
  double rtmp_s[3],rtmp_b[3];
  double delr_ss[3],rsq_ss,delr_sb[3],rsq_sb;
  double delr_bs[3],rsq_bs,delr_bb[3],rsq_bb;

  // vectors COM-backbone site, COM-base site in lab frame
  double ra_cs[3],ra_cb[3];
  double rb_cs[3],rb_cb[3];

  // quaternions and Cartesian unit vectors in lab frame
  double *qa,ax[3],ay[3],az[3];
  double *qb,bx[3],by[3],bz[3];
  double *special_lj = force->special_lj;

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  int *type = atom->type;

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int *alist,*blist,*numneigh,**firstneigh;

  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;

  int a,b,ia,ib,anum,bnum,atype,btype;

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

    // vector COM - backbone and base site a
    compute_interaction_sites(ax,ay,ra_cs,ra_cb);

    rtmp_s[0] = x[a][0] + ra_cs[0];
    rtmp_s[1] = x[a][1] + ra_cs[1];
    rtmp_s[2] = x[a][2] + ra_cs[2];

    rtmp_b[0] = x[a][0] + ra_cb[0];
    rtmp_b[1] = x[a][1] + ra_cb[1];
    rtmp_b[2] = x[a][2] + ra_cb[2];

    blist = firstneigh[a];
    bnum = numneigh[a];

    for (ib = 0; ib < bnum; ib++) {

      b = blist[ib];
      factor_lj = special_lj[sbmask(b)]; // = 0 for nearest neighbors
      b &= NEIGHMASK;

      btype = type[b];

      qb=bonus[b].quat;
      MathExtra::q_to_exyz(qb,bx,by,bz);

      // vector COM - backbone and base site b
      compute_interaction_sites(bx,by,rb_cs,rb_cb);

      // vector backbone site b to a
      delr_ss[0] = rtmp_s[0] - (x[b][0] + rb_cs[0]);
      delr_ss[1] = rtmp_s[1] - (x[b][1] + rb_cs[1]);
      delr_ss[2] = rtmp_s[2] - (x[b][2] + rb_cs[2]);
      rsq_ss = delr_ss[0]*delr_ss[0] + delr_ss[1]*delr_ss[1] + delr_ss[2]*delr_ss[2];

      // vector base site b to backbone site a
      delr_sb[0] =  rtmp_s[0] - (x[b][0] + rb_cb[0]);
      delr_sb[1] =  rtmp_s[1] - (x[b][1] + rb_cb[1]);
      delr_sb[2] =  rtmp_s[2] - (x[b][2] + rb_cb[2]);
      rsq_sb = delr_sb[0]*delr_sb[0] + delr_sb[1]*delr_sb[1] + delr_sb[2]*delr_sb[2];

      // vector backbone site b to base site a
      delr_bs[0] = rtmp_b[0] - (x[b][0] + rb_cs[0]);
      delr_bs[1] = rtmp_b[1] - (x[b][1] + rb_cs[1]);
      delr_bs[2] = rtmp_b[2] - (x[b][2] + rb_cs[2]);
      rsq_bs = delr_bs[0]*delr_bs[0] + delr_bs[1]*delr_bs[1] + delr_bs[2]*delr_bs[2];

      // vector base site b to a
      delr_bb[0] = rtmp_b[0] - (x[b][0] + rb_cb[0]);
      delr_bb[1] = rtmp_b[1] - (x[b][1] + rb_cb[1]);
      delr_bb[2] = rtmp_b[2] - (x[b][2] + rb_cb[2]);
      rsq_bb = delr_bb[0]*delr_bb[0] + delr_bb[1]*delr_bb[1] + delr_bb[2]*delr_bb[2];

      // excluded volume interaction

      // backbone-backbone
      if (rsq_ss < cutsq_ss_c[atype][btype]) {

        evdwl = F3(rsq_ss,cutsq_ss_ast[atype][btype],cut_ss_c[atype][btype],lj1_ss[atype][btype],
                        lj2_ss[atype][btype],epsilon_ss[atype][btype],b_ss[atype][btype],fpair);

        // knock out nearest-neighbor interaction between ss
        fpair *= factor_lj;
        evdwl *= factor_lj;

        // increment energy and virial
        if (evflag) ev_tally(a,b,nlocal,newton_pair,
                evdwl,0.0,fpair,delr_ss[0],delr_ss[1],delr_ss[2]);

        delf[0] = delr_ss[0]*fpair;
        delf[1] = delr_ss[1]*fpair;
        delf[2] = delr_ss[2]*fpair;

        f[a][0] += delf[0];
        f[a][1] += delf[1];
        f[a][2] += delf[2];

        MathExtra::cross3(ra_cs,delf,delta);

        torque[a][0] += delta[0];
        torque[a][1] += delta[1];
        torque[a][2] += delta[2];

        if (newton_pair || b < nlocal) {

          f[b][0] -= delf[0];
          f[b][1] -= delf[1];
          f[b][2] -= delf[2];

          MathExtra::cross3(rb_cs,delf,deltb);

          torque[b][0] -= deltb[0];
          torque[b][1] -= deltb[1];
          torque[b][2] -= deltb[2];

        }

      }


      // backbone-base
      if (rsq_sb < cutsq_sb_c[atype][btype]) {

        evdwl = F3(rsq_sb,cutsq_sb_ast[atype][btype],cut_sb_c[atype][btype],lj1_sb[atype][btype],
                        lj2_sb[atype][btype],epsilon_sb[atype][btype],b_sb[atype][btype],fpair);

        // increment energy and virial
        if (evflag) ev_tally(a,b,nlocal,newton_pair,
                evdwl,0.0,fpair,delr_sb[0],delr_sb[1],delr_sb[2]);

        delf[0] = delr_sb[0]*fpair;
        delf[1] = delr_sb[1]*fpair;
        delf[2] = delr_sb[2]*fpair;

        f[a][0] += delf[0];
        f[a][1] += delf[1];
        f[a][2] += delf[2];

        MathExtra::cross3(ra_cs,delf,delta);

        torque[a][0] += delta[0];
        torque[a][1] += delta[1];
        torque[a][2] += delta[2];

        if (newton_pair || b < nlocal) {

          f[b][0] -= delf[0];
          f[b][1] -= delf[1];
          f[b][2] -= delf[2];

          MathExtra::cross3(rb_cb,delf,deltb);

          torque[b][0] -= deltb[0];
          torque[b][1] -= deltb[1];
          torque[b][2] -= deltb[2];

        }

      }

      // base-backbone
      if (rsq_bs < cutsq_sb_c[atype][btype]) {

        evdwl = F3(rsq_bs,cutsq_sb_ast[atype][btype],cut_sb_c[atype][btype],lj1_sb[atype][btype],
                        lj2_sb[atype][btype],epsilon_sb[atype][btype],b_sb[atype][btype],fpair);

        // increment energy and virial
        if (evflag) ev_tally(a,b,nlocal,newton_pair,
                evdwl,0.0,fpair,delr_bs[0],delr_bs[1],delr_bs[2]);

        delf[0] = delr_bs[0]*fpair;
        delf[1] = delr_bs[1]*fpair;
        delf[2] = delr_bs[2]*fpair;

        f[a][0] += delf[0];
        f[a][1] += delf[1];
        f[a][2] += delf[2];

        MathExtra::cross3(ra_cb,delf,delta);

        torque[a][0] += delta[0];
        torque[a][1] += delta[1];
        torque[a][2] += delta[2];

        if (newton_pair || b < nlocal) {

          f[b][0] -= delf[0];
          f[b][1] -= delf[1];
          f[b][2] -= delf[2];

          MathExtra::cross3(rb_cs,delf,deltb);

          torque[b][0] -= deltb[0];
          torque[b][1] -= deltb[1];
          torque[b][2] -= deltb[2];

        }

      }

      // base-base
      if (rsq_bb < cutsq_bb_c[atype][btype]) {

        evdwl = F3(rsq_bb,cutsq_bb_ast[atype][btype],cut_bb_c[atype][btype],lj1_bb[atype][btype],
                        lj2_bb[atype][btype],epsilon_bb[atype][btype],b_bb[atype][btype],fpair);

        // increment energy and virial
        if (evflag) ev_tally(a,b,nlocal,newton_pair,
                evdwl,0.0,fpair,delr_bb[0],delr_bb[1],delr_bb[2]);

        delf[0] = delr_bb[0]*fpair;
        delf[1] = delr_bb[1]*fpair;
        delf[2] = delr_bb[2]*fpair;

        f[a][0] += delf[0];
        f[a][1] += delf[1];
        f[a][2] += delf[2];

        MathExtra::cross3(ra_cb,delf,delta);

        torque[a][0] += delta[0];
        torque[a][1] += delta[1];
        torque[a][2] += delta[2];

        if (newton_pair || b < nlocal) {

          f[b][0] -= delf[0];
          f[b][1] -= delf[1];
          f[b][2] -= delf[2];

          MathExtra::cross3(rb_cb,delf,deltb);

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

  memory->create(epsilon_ss,n+1,n+1,"pair:epsilon_ss");
  memory->create(sigma_ss,n+1,n+1,"pair:sigma_ss");
  memory->create(cut_ss_ast,n+1,n+1,"pair:cut_ss_ast");
  memory->create(b_ss,n+1,n+1,"pair:b_ss");
  memory->create(cut_ss_c,n+1,n+1,"pair:cut_ss_c");
  memory->create(lj1_ss,n+1,n+1,"pair:lj1_ss");
  memory->create(lj2_ss,n+1,n+1,"pair:lj2_ss");
  memory->create(cutsq_ss_ast,n+1,n+1,"pair:cutsq_ss_ast");
  memory->create(cutsq_ss_c,n+1,n+1,"pair:cutsq_ss_c");

  memory->create(epsilon_sb,n+1,n+1,"pair:epsilon_sb");
  memory->create(sigma_sb,n+1,n+1,"pair:sigma_sb");
  memory->create(cut_sb_ast,n+1,n+1,"pair:cut_sb_ast");
  memory->create(b_sb,n+1,n+1,"pair:b_sb");
  memory->create(cut_sb_c,n+1,n+1,"pair:cut_sb_c");
  memory->create(lj1_sb,n+1,n+1,"pair:lj1_sb");
  memory->create(lj2_sb,n+1,n+1,"pair:lj2_sb");
  memory->create(cutsq_sb_ast,n+1,n+1,"pair:cutsq_sb_ast");
  memory->create(cutsq_sb_c,n+1,n+1,"pair:cutsq_sb_c");

  memory->create(epsilon_bb,n+1,n+1,"pair:epsilon_bb");
  memory->create(sigma_bb,n+1,n+1,"pair:sigma_bb");
  memory->create(cut_bb_ast,n+1,n+1,"pair:cut_bb_ast");
  memory->create(b_bb,n+1,n+1,"pair:b_bb");
  memory->create(cut_bb_c,n+1,n+1,"pair:cut_bb_c");
  memory->create(lj1_bb,n+1,n+1,"pair:lj1_bb");
  memory->create(lj2_bb,n+1,n+1,"pair:lj2_bb");
  memory->create(cutsq_bb_ast,n+1,n+1,"pair:cutsq_bb_ast");
  memory->create(cutsq_bb_c,n+1,n+1,"pair:cutsq_bb_c");

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

  if (narg != 11) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/excv");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  count = 0;

  double epsilon_ss_one, sigma_ss_one;
  double cut_ss_ast_one, cut_ss_c_one, b_ss_one;

  double epsilon_sb_one, sigma_sb_one;
  double cut_sb_ast_one, cut_sb_c_one, b_sb_one;

  double epsilon_bb_one, sigma_bb_one;
  double cut_bb_ast_one, cut_bb_c_one, b_bb_one;

  // Excluded volume interaction
  // LJ parameters
  epsilon_ss_one = force->numeric(FLERR,arg[2]);
  sigma_ss_one = force->numeric(FLERR,arg[3]);
  cut_ss_ast_one = force->numeric(FLERR,arg[4]);

  // smoothing - determined through continuity and differentiability
  b_ss_one = 4.0/sigma_ss_one
      *(6.0*pow(sigma_ss_one/cut_ss_ast_one,7)-12.0*pow(sigma_ss_one/cut_ss_ast_one,13))
      *4.0/sigma_ss_one*(6.0*pow(sigma_ss_one/cut_ss_ast_one,7)-12.0*pow(sigma_ss_one/cut_ss_ast_one,13))
      /4.0/(4.0*(pow(sigma_ss_one/cut_ss_ast_one,12)-pow(sigma_ss_one/cut_ss_ast_one,6)));

  cut_ss_c_one = cut_ss_ast_one
      - 2.0*4.0*(pow(sigma_ss_one/cut_ss_ast_one,12)-pow(sigma_ss_one/cut_ss_ast_one,6))
      /(4.0/sigma_ss_one*(6.0*pow(sigma_ss_one/cut_ss_ast_one,7)-12.0*pow(sigma_ss_one/cut_ss_ast_one,13)));

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon_ss[i][j] = epsilon_ss_one;
      sigma_ss[i][j] = sigma_ss_one;
      cut_ss_ast[i][j] = cut_ss_ast_one;
      b_ss[i][j] = b_ss_one;
      cut_ss_c[i][j] = cut_ss_c_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/excv");

  count = 0;

  // LJ parameters
  epsilon_sb_one = force->numeric(FLERR,arg[5]);
  sigma_sb_one = force->numeric(FLERR,arg[6]);
  cut_sb_ast_one = force->numeric(FLERR,arg[7]);

  // smoothing - determined through continuity and differentiability
  b_sb_one = 4.0/sigma_sb_one
      *(6.0*pow(sigma_sb_one/cut_sb_ast_one,7)-12.0*pow(sigma_sb_one/cut_sb_ast_one,13))
      *4.0/sigma_sb_one*(6.0*pow(sigma_sb_one/cut_sb_ast_one,7)-12.0*pow(sigma_sb_one/cut_sb_ast_one,13))
      /4.0/(4.0*(pow(sigma_sb_one/cut_sb_ast_one,12)-pow(sigma_sb_one/cut_sb_ast_one,6)));

  cut_sb_c_one = cut_sb_ast_one
      - 2.0*4.0*(pow(sigma_sb_one/cut_sb_ast_one,12)-pow(sigma_sb_one/cut_sb_ast_one,6))
      /(4.0/sigma_sb_one*(6.0*pow(sigma_sb_one/cut_sb_ast_one,7)-12.0*pow(sigma_sb_one/cut_sb_ast_one,13)));

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon_sb[i][j] = epsilon_sb_one;
      sigma_sb[i][j] = sigma_sb_one;
      cut_sb_ast[i][j] = cut_sb_ast_one;
      b_sb[i][j] = b_sb_one;
      cut_sb_c[i][j] = cut_sb_c_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/excv");

  count = 0;

  // LJ parameters
  epsilon_bb_one = force->numeric(FLERR,arg[8]);
  sigma_bb_one = force->numeric(FLERR,arg[9]);
  cut_bb_ast_one = force->numeric(FLERR,arg[10]);

  // smoothing - determined through continuity and differentiability
  b_bb_one = 4.0/sigma_bb_one
      *(6.0*pow(sigma_bb_one/cut_bb_ast_one,7)-12.0*pow(sigma_bb_one/cut_bb_ast_one,13))
      *4.0/sigma_bb_one*(6.0*pow(sigma_bb_one/cut_bb_ast_one,7)-12.0*pow(sigma_bb_one/cut_bb_ast_one,13))
      /4.0/(4.0*(pow(sigma_bb_one/cut_bb_ast_one,12)-pow(sigma_bb_one/cut_bb_ast_one,6)));

  cut_bb_c_one = cut_bb_ast_one
      - 2.0*4.0*(pow(sigma_bb_one/cut_bb_ast_one,12)-pow(sigma_bb_one/cut_bb_ast_one,6))
      /(4.0/sigma_bb_one*(6.0*pow(sigma_bb_one/cut_bb_ast_one,7)-12.0*pow(sigma_bb_one/cut_bb_ast_one,13)));

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon_bb[i][j] = epsilon_bb_one;
      sigma_bb[i][j] = sigma_bb_one;
      cut_bb_ast[i][j] = cut_bb_ast_one;
      b_bb[i][j] = b_bb_one;
      cut_bb_c[i][j] = cut_bb_c_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/excv");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairOxdnaExcv::init_style()
{
  int irequest;

  // request regular neighbor lists

  irequest = neighbor->request(this,instance_me);

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
    error->all(FLERR,"Coefficient mixing not defined in oxDNA");
  }
  if (offset_flag) {
    error->all(FLERR,"Offset not supported in oxDNA");
  }

  epsilon_ss[j][i] = epsilon_ss[i][j];
  sigma_ss[j][i] = sigma_ss[i][j];
  cut_ss_ast[j][i] = cut_ss_ast[i][j];
  cut_ss_c[j][i] = cut_ss_c[i][j];
  b_ss[j][i] = b_ss[i][j];

  epsilon_sb[j][i] = epsilon_sb[i][j];
  sigma_sb[j][i] = sigma_sb[i][j];
  cut_sb_ast[j][i] = cut_sb_ast[i][j];
  cut_sb_c[j][i] = cut_sb_c[i][j];
  b_sb[j][i] = b_sb[i][j];

  epsilon_bb[j][i] = epsilon_bb[i][j];
  sigma_bb[j][i] = sigma_bb[i][j];
  cut_bb_ast[j][i] = cut_bb_ast[i][j];
  cut_bb_c[j][i] = cut_bb_c[i][j];
  b_bb[j][i] = b_bb[i][j];

  // excluded volume auxiliary parameters

  lj1_ss[i][j] = 4.0 * epsilon_ss[i][j] * pow(sigma_ss[i][j],12.0);
  lj2_ss[i][j] = 4.0 * epsilon_ss[i][j] * pow(sigma_ss[i][j],6.0);

  lj1_sb[i][j] = 4.0 * epsilon_sb[i][j] * pow(sigma_sb[i][j],12.0);
  lj2_sb[i][j] = 4.0 * epsilon_sb[i][j] * pow(sigma_sb[i][j],6.0);

  lj1_bb[i][j] = 4.0 * epsilon_bb[i][j] * pow(sigma_bb[i][j],12.0);
  lj2_bb[i][j] = 4.0 * epsilon_bb[i][j] * pow(sigma_bb[i][j],6.0);

  lj1_ss[j][i] = lj1_ss[i][j];
  lj2_ss[j][i] = lj2_ss[i][j];

  lj1_sb[j][i] = lj1_sb[i][j];
  lj2_sb[j][i] = lj2_sb[i][j];

  lj1_bb[j][i] = lj1_bb[i][j];
  lj2_bb[j][i] = lj2_bb[i][j];

  cutsq_ss_ast[i][j] = cut_ss_ast[i][j]*cut_ss_ast[i][j];
  cutsq_ss_c[i][j]  = cut_ss_c[i][j]*cut_ss_c[i][j];

  cutsq_sb_ast[i][j] = cut_sb_ast[i][j]*cut_sb_ast[i][j];
  cutsq_sb_c[i][j]  = cut_sb_c[i][j]*cut_sb_c[i][j];

  cutsq_bb_ast[i][j] = cut_bb_ast[i][j]*cut_bb_ast[i][j];
  cutsq_bb_c[i][j]  = cut_bb_c[i][j]*cut_bb_c[i][j];

  cutsq_ss_ast[j][i] = cutsq_ss_ast[i][j];
  cutsq_ss_c[j][i]  = cutsq_ss_c[i][j];

  cutsq_sb_ast[j][i] = cutsq_sb_ast[i][j];
  cutsq_sb_c[j][i]  = cutsq_sb_c[i][j];

  cutsq_bb_ast[j][i] = cutsq_bb_ast[i][j];
  cutsq_bb_c[j][i]  = cutsq_bb_c[i][j];

  // set the master list distance cutoff
  return cut_ss_c[i][j];

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

        fwrite(&epsilon_ss[i][j],sizeof(double),1,fp);
        fwrite(&sigma_ss[i][j],sizeof(double),1,fp);
        fwrite(&cut_ss_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_ss[i][j],sizeof(double),1,fp);
        fwrite(&cut_ss_c[i][j],sizeof(double),1,fp);
        fwrite(&epsilon_sb[i][j],sizeof(double),1,fp);
        fwrite(&sigma_sb[i][j],sizeof(double),1,fp);
        fwrite(&cut_sb_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_sb[i][j],sizeof(double),1,fp);
        fwrite(&cut_sb_c[i][j],sizeof(double),1,fp);
        fwrite(&epsilon_bb[i][j],sizeof(double),1,fp);
        fwrite(&sigma_bb[i][j],sizeof(double),1,fp);
        fwrite(&cut_bb_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_bb[i][j],sizeof(double),1,fp);
        fwrite(&cut_bb_c[i][j],sizeof(double),1,fp);

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
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {

          fread(&epsilon_ss[i][j],sizeof(double),1,fp);
          fread(&sigma_ss[i][j],sizeof(double),1,fp);
          fread(&cut_ss_ast[i][j],sizeof(double),1,fp);
          fread(&b_ss[i][j],sizeof(double),1,fp);
          fread(&cut_ss_c[i][j],sizeof(double),1,fp);
          fread(&epsilon_sb[i][j],sizeof(double),1,fp);
          fread(&sigma_sb[i][j],sizeof(double),1,fp);
          fread(&cut_sb_ast[i][j],sizeof(double),1,fp);
          fread(&b_sb[i][j],sizeof(double),1,fp);
          fread(&cut_sb_c[i][j],sizeof(double),1,fp);
          fread(&epsilon_bb[i][j],sizeof(double),1,fp);
          fread(&sigma_bb[i][j],sizeof(double),1,fp);
          fread(&cut_bb_ast[i][j],sizeof(double),1,fp);
          fread(&b_bb[i][j],sizeof(double),1,fp);
          fread(&cut_bb_c[i][j],sizeof(double),1,fp);

         }

        MPI_Bcast(&epsilon_ss[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma_ss[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_ss_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_ss[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_ss_c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&epsilon_sb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma_sb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_sb_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_sb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_sb_c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&epsilon_bb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma_bb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_bb_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_bb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_bb_c[i][j],1,MPI_DOUBLE,0,world);

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

void PairOxdnaExcv::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         \n",i,
        epsilon_ss[i][i],sigma_ss[i][i],cut_ss_ast[i][i],b_ss[i][i],cut_ss_c[i][i],
        epsilon_sb[i][i],sigma_sb[i][i],cut_sb_ast[i][i],b_sb[i][i],cut_sb_c[i][i],
        epsilon_bb[i][i],sigma_bb[i][i],cut_bb_ast[i][i],b_bb[i][i],cut_bb_c[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairOxdnaExcv::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         \n",i,j,
        epsilon_ss[i][j],sigma_ss[i][j],cut_ss_ast[i][j],b_ss[i][j],cut_ss_c[i][j],
        epsilon_sb[i][j],sigma_sb[i][j],cut_sb_ast[i][j],b_sb[i][j],cut_sb_c[i][j],
        epsilon_bb[i][j],sigma_bb[i][j],cut_bb_ast[i][j],b_bb[i][j],cut_bb_c[i][j]);
}

/* ---------------------------------------------------------------------- */

void *PairOxdnaExcv::extract(const char *str, int &dim)
{
  dim = 2;

  if (strcmp(str,"epsilon_ss") == 0) return (void *) epsilon_ss;
  if (strcmp(str,"sigma_ss") == 0) return (void *) sigma_ss;
  if (strcmp(str,"cut_ss_ast") == 0) return (void *) cut_ss_ast;
  if (strcmp(str,"b_ss") == 0) return (void *) b_ss;
  if (strcmp(str,"cut_ss_c") == 0) return (void *) cut_ss_c;
  if (strcmp(str,"epsilon_sb") == 0) return (void *) epsilon_sb;
  if (strcmp(str,"sigma_sb") == 0) return (void *) sigma_sb;
  if (strcmp(str,"cut_sb_ast") == 0) return (void *) cut_sb_ast;
  if (strcmp(str,"b_sb") == 0) return (void *) b_sb;
  if (strcmp(str,"cut_sb_c") == 0) return (void *) cut_sb_c;
  if (strcmp(str,"epsilon_bb") == 0) return (void *) epsilon_bb;
  if (strcmp(str,"sigma_bb") == 0) return (void *) sigma_bb;
  if (strcmp(str,"cut_bb_ast") == 0) return (void *) cut_bb_ast;
  if (strcmp(str,"b_bb") == 0) return (void *) b_bb;
  if (strcmp(str,"cut_bb_c") == 0) return (void *) cut_bb_c;

  return NULL;
}
