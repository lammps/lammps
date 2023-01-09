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
   Contributing author: Mark Chaimovich(RSM) mark.chaimovich@russianschool.com
------------------------------------------------------------------------- */

#include "pair_lj_relres.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"

#include <cmath>

using namespace LAMMPS_NS;

static const char cite_relres[] =
  "Pair style lj/relres: doi:10.1021/acs.jctc.0c01003, doi:10.1021/acs.jctc.0c01003\n\n"
  "@Article{Chaimovich1,\n"
  " author = {A. Chaimovich, C. Peter, K. Kremer},\n"
  " title = {Relative Resolution: {A} Hybrid Formalism for Fluid Mixtures},\n"
  " journal = {J.~Chem.\\ Phys.},\n"
  " year =    2015,\n"
  " volume =  143,\n"
  " pages =   {243107}\n"
  "@Article{Chaimovich2,\n"
  " author = {M. Chaimovich and A. Chaimovich},\n"
  " title = {Relative Resolution: A Computationally Efficient Implementation in LAMMPS},\n"
  " journal = {J.~Chem.\\ Theory Comput.},\n"
  " year =    2021,\n"
  " volume =  17,\n"
  " pages =   {1045--1059}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

PairLJRelRes::PairLJRelRes(LAMMPS *lmp) : Pair(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_relres);
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairLJRelRes::~PairLJRelRes()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cutfsq);

    memory->destroy(cut);
    memory->destroy(cut_inner);
    memory->destroy(cut_inner_sq);
    memory->destroy(cutf);
    memory->destroy(cutf_inner);
    memory->destroy(cutf_inner_sq);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(epsilonf);
    memory->destroy(sigmaf);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(ljsw0);
    memory->destroy(ljsw1);
    memory->destroy(ljsw2);
    memory->destroy(ljsw3);
    memory->destroy(ljsw4);
    memory->destroy(ljf1);
    memory->destroy(ljf2);
    memory->destroy(ljf3);
    memory->destroy(ljf4);
    memory->destroy(ljswf0);
    memory->destroy(ljswf1);
    memory->destroy(ljswf2);
    memory->destroy(ljswf3);
    memory->destroy(ljswf4);
    memory->destroy(offset);
    memory->destroy(offsetsm);
    memory->destroy(offsetsp);
  }
}

/* ---------------------------------------------------------------------- */

void PairLJRelRes::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj;
  double r,t,tsq,fskin;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        if (rsq < cutf_inner_sq[itype][jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv*(ljf1[itype][jtype]*r6inv-ljf2[itype][jtype]);
        } else if (rsq < cutfsq[itype][jtype]) {
          r = sqrt(rsq);
          t = r - cutf_inner[itype][jtype];
          tsq = t*t;
          fskin = ljswf1[itype][jtype]+ljswf2[itype][jtype]*t+
            ljswf3[itype][jtype]*tsq+ljswf4[itype][jtype]*tsq*t;
          forcelj = fskin*r;
        } else if (rsq < cut_inner_sq[itype][jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv*(lj1[itype][jtype]*r6inv-lj2[itype][jtype]);
        } else {
          r = sqrt(rsq);
          t = r-cut_inner[itype][jtype];
          tsq = t*t;
          fskin = ljsw1[itype][jtype]+ljsw2[itype][jtype]*t+
            ljsw3[itype][jtype]*tsq+ljsw4[itype][jtype]*tsq*t;
          forcelj = fskin*r;
        }

        fpair = factor_lj*forcelj*r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          if (rsq < cutf_inner_sq[itype][jtype]) {
            evdwl = r6inv*(ljf3[itype][jtype]*r6inv-ljf4[itype][jtype])-offsetsm[itype][jtype];
          } else if (rsq < cutfsq[itype][jtype]) {
            evdwl = ljswf0[itype][jtype]-ljswf1[itype][jtype]*t-
              ljswf2[itype][jtype]*tsq/2.0-ljswf3[itype][jtype]*tsq*t/3.0-
              ljswf4[itype][jtype]*tsq*tsq/4.0-offsetsp[itype][jtype];
          } else if (rsq < cut_inner_sq[itype][jtype]) {
            evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype])-offset[itype][jtype];
          } else {
            evdwl = ljsw0[itype][jtype]-ljsw1[itype][jtype]*t-
              ljsw2[itype][jtype]*tsq/2.0-ljsw3[itype][jtype]*tsq*t/3.0-
              ljsw4[itype][jtype]*tsq*tsq/4.0-offset[itype][jtype];
          }
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJRelRes::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cutfsq,n+1,n+1,"pair:cutfsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(cutf,n+1,n+1,"pair:cutf");
  memory->create(cut_inner,n+1,n+1,"pair:cut_inner");
  memory->create(cutf_inner,n+1,n+1,"pair:cutf_inner");
  memory->create(cut_inner_sq,n+1,n+1,"pair:cut_inner_sq");
  memory->create(cutf_inner_sq,n+1,n+1,"pair:cutf_inner_sq");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(epsilonf,n+1,n+1,"pair:epsilonf");
  memory->create(sigmaf,n+1,n+1,"pair:sigmaf");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(ljf1,n+1,n+1,"pair:ljf1");
  memory->create(ljf2,n+1,n+1,"pair:ljf2");
  memory->create(ljf3,n+1,n+1,"pair:ljf3");
  memory->create(ljf4,n+1,n+1,"pair:ljf4");
  memory->create(ljsw0,n+1,n+1,"pair:ljsw0");
  memory->create(ljsw1,n+1,n+1,"pair:ljsw1");
  memory->create(ljsw2,n+1,n+1,"pair:ljsw2");
  memory->create(ljsw3,n+1,n+1,"pair:ljsw3");
  memory->create(ljsw4,n+1,n+1,"pair:ljsw4");
  memory->create(ljswf0,n+1,n+1,"pair:ljswf0");
  memory->create(ljswf1,n+1,n+1,"pair:ljswf1");
  memory->create(ljswf2,n+1,n+1,"pair:ljswf2");
  memory->create(ljswf3,n+1,n+1,"pair:ljswf3");
  memory->create(ljswf4,n+1,n+1,"pair:ljswf4");
  memory->create(offset,n+1,n+1,"pair:offset");
  memory->create(offsetsp,n+1,n+1,"pair:offsetsp");
  memory->create(offsetsm,n+1,n+1,"pair:offsetsm");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJRelRes::settings(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Illegal pair_style command");

  cutf_inner_global = utils::numeric(FLERR,arg[0],false,lmp);
  cutf_global = utils::numeric(FLERR,arg[1],false,lmp);
  cut_inner_global = utils::numeric(FLERR,arg[2],false,lmp);
  cut_global = utils::numeric(FLERR,arg[3],false,lmp);
  if (cut_inner_global <= 0.0 || cut_inner_global > cut_global)
    error->all(FLERR,"Illegal pair_style command");
  if (cutf_inner_global <= 0.0 || cutf_inner_global > cutf_global)
    error->all(FLERR,"Illegal pair_style command");
  if (cutf_global > cut_inner_global)
    error->all(FLERR,"Illegal pair_style command");

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_inner[i][j] = cut_inner_global;
          cut[i][j] = cut_global;
          cutf_inner[i][j] = cutf_inner_global;
          cutf[i][j] = cutf_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJRelRes::coeff(int narg, char **arg)
{
  if (narg != 6 && narg != 10)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double epsilonf_one = utils::numeric(FLERR,arg[2],false,lmp);
  double sigmaf_one = utils::numeric(FLERR,arg[3],false,lmp);
  double epsilon_one = utils::numeric(FLERR,arg[4],false,lmp);
  double sigma_one = utils::numeric(FLERR,arg[5],false,lmp);

  double cut_inner_one = cut_inner_global;
  double cut_one = cut_global;
  double cutf_inner_one = cutf_inner_global;
  double cutf_one = cutf_global;

  if (narg == 10) {
    cutf_inner_one = utils::numeric(FLERR,arg[6],false,lmp);
    cutf_one = utils::numeric(FLERR,arg[7],false,lmp);
    cut_inner_one = utils::numeric(FLERR,arg[8],false,lmp);
    cut_one = utils::numeric(FLERR,arg[9],false,lmp);
  }

  if (cut_inner_one <= 0.0 || cut_inner_one > cut_one)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (cutf_inner_one <= 0.0 || cutf_inner_one > cutf_one)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (cutf_one > cut_inner_one)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (epsilon_one == 0.0) {  //set cutoff for fg interactions
    cut_inner_one = cutf_one;
    cut_one = cutf_one;
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      epsilonf[i][j] = epsilonf_one;
      sigmaf[i][j] = sigmaf_one;
      cut_inner[i][j] = cut_inner_one;
      cut[i][j] = cut_one;
      cutf_inner[i][j] = cutf_inner_one;
      cutf[i][j] = cutf_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJRelRes::init_one(int i, int j)
{ double ljswc0,ljswc3,ljswc4;
// mixing rules:
//   fg and cg - no mixing;
//   fg and fg or fg anf hybrid - mixing fg parameters only
//   cg and cg of cg and hybrid - mixing cg parameters only
//   hybrid and hybrid - mixing fg and cg parameters

  if (setflag[i][j] == 0) {
    if (((epsilon[i][i] == 0.0) && (epsilonf[j][j] == 0.0))
        || ((epsilonf[i][i] == 0.0) && (epsilon[j][j] == 0.0))) { //no mixing
      epsilon[i][j] = 0.0;
      epsilonf[i][j] = 0.0;
      sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
      sigmaf[i][j] = mix_distance(sigmaf[i][i],sigmaf[j][j]);
      cut_inner[i][j] = cutf[i][j] = cutf_inner[i][j] = cut[i][j] = 0.0;
    } else if ((epsilon[i][i] == 0.0) || (epsilon[j][j] == 0.0)) { // fg only
      epsilon[i][j] = 0.0;
      sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
      epsilonf[i][j] = mix_energy(epsilonf[i][i],epsilonf[j][j],
                                  sigmaf[i][i],sigmaf[j][j]);
      sigmaf[i][j] = mix_distance(sigmaf[i][i],sigmaf[j][j]);
      cutf_inner[i][j] = mix_distance(cutf_inner[i][i],cutf_inner[j][j]);
      cutf[i][j] = mix_distance(cutf[i][i],cutf[j][j]);
      cut_inner[i][j] = cutf[i][j];
      cut[i][j] = cutf[i][j];
    } else if ((epsilonf[i][i] == 0.0) || (epsilonf[j][j] == 0.0)) { // cg only
      epsilonf[i][j] = 0.0;
      epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                                 sigma[i][i],sigma[j][j]);
      sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
      sigmaf[i][j] = mix_distance(sigmaf[i][i],sigmaf[j][j]);
      cut_inner[i][j] = mix_distance(cut_inner[i][i],cut_inner[j][j]);
      cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
      cutf_inner[i][j] = mix_distance(cutf_inner[i][i],cutf_inner[j][j]);
      cutf[i][j] = mix_distance(cutf[i][i],cutf[j][j]);
    } else {                                              // fg and cg
      epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                                 sigma[i][i],sigma[j][j]);
      sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
      epsilonf[i][j] = mix_energy(epsilonf[i][i],epsilonf[j][j],
                                  sigmaf[i][i],sigmaf[j][j]);
      sigmaf[i][j] = mix_distance(sigmaf[i][i],sigmaf[j][j]);
      cut_inner[i][j] = mix_distance(cut_inner[i][i],cut_inner[j][j]);
      cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
      cutf_inner[i][j] = mix_distance(cutf_inner[i][i],cutf_inner[j][j]);
      cutf[i][j] = mix_distance(cutf[i][i],cutf[j][j]);
    }
  }

  cut_inner_sq[i][j] = cut_inner[i][j]*cut_inner[i][j];
  cutf_inner_sq[i][j] = cutf_inner[i][j]*cutf_inner[i][j];
  cutfsq[i][j] = cutf[i][j]*cutf[i][j];

  if (epsilon[i][j] != 0) { // cg or fg+cg (cut coefficients)
    lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
    lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
    lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
    lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
    if (cut_inner[i][j] != cut[i][j]) {
      double r6inv = 1.0/pow(cut_inner[i][j],6.0);
      double t = cut[i][j] - cut_inner[i][j];
      double tsq = t*t;
      double ratio = sigma[i][j] / cut_inner[i][j];
      ljsw0[i][j] = 4.0*epsilon[i][j]*(pow(ratio,12.0) - pow(ratio,6.0));
      ljsw1[i][j] = r6inv*(lj1[i][j]*r6inv-lj2[i][j]) / cut_inner[i][j];
      ljsw2[i][j] = -r6inv * (13.0*lj1[i][j]*r6inv - 7.0*lj2[i][j]) /
        cut_inner_sq[i][j];
      ljsw3[i][j] = -(3.0/tsq) * (ljsw1[i][j] + 2.0/3.0*ljsw2[i][j]*t);
      ljsw4[i][j] = -1.0/(3.0*tsq) * (ljsw2[i][j] + 2.0*ljsw3[i][j]*t);
      if (offset_flag) {
        offset[i][j] = ljsw0[i][j] - ljsw1[i][j]*t - ljsw2[i][j]*tsq/2.0 -
          ljsw3[i][j]*tsq*t/3.0 - ljsw4[i][j]*tsq*tsq/4.0;
      } else offset[i][j] = 0.0;
    } else {
      ljsw0[i][j] = 0.0;
      ljsw1[i][j] = 0.0;
      ljsw2[i][j] = 0.0;
      ljsw3[i][j] = 0.0;
      ljsw4[i][j] = 0.0;
      double ratio = sigma[i][j] / cut_inner[i][j];
      if (offset_flag)
        offset[i][j] = 4.0*epsilon[i][j]*(pow(ratio,12.0) - pow(ratio,6.0));
      else offset[i][j] = 0.0;
    }
  } else {
    ljsw0[i][j] = 0.0;
    ljsw1[i][j] = 0.0;
    ljsw2[i][j] = 0.0;
    ljsw3[i][j] = 0.0;
    ljsw4[i][j] = 0.0;
    lj1[i][j] = 0.0;
    lj2[i][j] = 0.0;
    lj3[i][j] = 0.0;
    lj4[i][j] = 0.0;
    offset[i][j] = 0.0;
  }

  if (epsilonf[i][j] != 0 ) {  // fg (cut=cutf coefficients)
    ljf1[i][j] = 48.0 * epsilonf[i][j] * pow(sigmaf[i][j],12.0);
    ljf2[i][j] = 24.0 * epsilonf[i][j] * pow(sigmaf[i][j],6.0);
    ljf3[i][j] = 4.0 * epsilonf[i][j] * pow(sigmaf[i][j],12.0);
    ljf4[i][j] = 4.0 * epsilonf[i][j] * pow(sigmaf[i][j],6.0);
    if (cutf_inner[i][j] != cutf[i][j]) {
      double r6inv = 1.0/pow(cutf_inner[i][j],6.0);
      double t = cutf[i][j] - cutf_inner[i][j];
      double tsq = t*t;
      double ratio = sigmaf[i][j] / cutf_inner[i][j];
      ljswf0[i][j] = 4.0*epsilonf[i][j]*(pow(ratio,12.0) - pow(ratio,6.0));
      ljswf1[i][j] = r6inv*(ljf1[i][j]*r6inv-ljf2[i][j]) / cutf_inner[i][j];
      ljswf2[i][j] = -r6inv * (13.0*ljf1[i][j]*r6inv - 7.0*ljf2[i][j]) /
        cutf_inner_sq[i][j];
      ljswf3[i][j] = -(3.0/tsq) * (ljswf1[i][j] + 2.0/3.0*ljswf2[i][j]*t);
      ljswf4[i][j] = -1.0/(3.0*tsq) * (ljswf2[i][j] + 2.0*ljswf3[i][j]*t);
      offsetsp[i][j] = ljswf0[i][j] - ljswf1[i][j]*t - ljswf2[i][j]*tsq/2.0-
        ljswf3[i][j]*tsq*t/3.0 - ljswf4[i][j]*tsq*tsq/4.0;
    } else {
      ljswf0[i][j] = 0.0;
      ljswf1[i][j] = 0.0;
      ljswf2[i][j] = 0.0;
      ljswf3[i][j] = 0.0;
      ljswf4[i][j] = 0.0;
      double ratio = sigmaf[i][j] / cutf_inner[i][j];
      offsetsp[i][j] = 4.0*epsilonf[i][j]*(pow(ratio,12.0) - pow(ratio,6.0));
    }
  } else {
    ljswf0[i][j] = 0.0;
    ljswf1[i][j] = 0.0;
    ljswf2[i][j] = 0.0;
    ljswf3[i][j] = 0.0;
    ljswf4[i][j] = 0.0;
    ljf4[i][j] = 0.0;
    ljf1[i][j] = 0.0;
    ljf2[i][j] = 0.0;
    ljf3[i][j] = 0.0;
    offsetsp[i][j] = 0.0;
  }

  if (epsilon[i][j] != 0) {  // cg or fg+cg (cutf coefficients)
    if (cutf_inner[i][j] != cutf[i][j]) {
      double r2inv = 1.0/pow(cutf[i][j],2.0);
      double r6inv = r2inv * r2inv * r2inv;
      double t = cutf[i][j] - cutf_inner[i][j];
      double tsq = t*t;
      double tsqinv = 1.0/tsq;
      double ratio = sigma[i][j] / cutf[i][j];
      double Et = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
      double Ft = r6inv * (lj1[i][j] * r6inv - lj2[i][j]) / cutf[i][j];
      double dFt = -r6inv * (13.0*lj1[i][j]*r6inv - 7.0*lj2[i][j]) * r2inv;
      double A = Ft + dFt * t / 3.0;

      ljswc3 = 3.0 * A * tsqinv;
      ljswc4 = -(2.0 * Ft + dFt * t) * tsqinv / t;
      ljswc0 = Et + ljswc3 * t * tsq /3.0 + ljswc4 * tsq * tsq / 4.0;
      offsetsm[i][j] = ljswc0;
    } else {
      ljswc0 = 0.0;
      ljswc3 = 0.0;
      ljswc4 = 0.0;
      double ratio = sigma[i][j] / cutf_inner[i][j];
      offsetsm[i][j] = 4.0*epsilon[i][j]*(pow(ratio,12.0) - pow(ratio,6.0));
    }
  } else {
    ljswc0 = 0.0;
    ljswc3 = 0.0;
    ljswc4 = 0.0;
    offsetsm[i][j] = 0.0;
  }
  // combine cutf coefficients
  ljswf0[i][j] += ljswc0;
  ljswf3[i][j] += ljswc3;
  ljswf4[i][j] += ljswc4;

  // combine shifting constants
  offsetsp[i][j] += offset[i][j];
  offsetsm[i][j] = offsetsp[i][j] - offsetsm[i][j];

  if (i !=j) {
    cut[j][i] = cut[i][j];
    cutsq[j][i] = cutsq[i][j];
    cutf[j][i] = cutf[i][j];
    cutfsq[j][i] = cutfsq[i][j];
    cut_inner[j][i] = cut_inner[i][j];
    cut_inner_sq[j][i] = cut_inner_sq[i][j];
    cutf_inner[j][i] = cutf_inner[i][j];
    cutf_inner_sq[j][i] = cutf_inner_sq[i][j];
    lj1[j][i] = lj1[i][j];
    lj2[j][i] = lj2[i][j];
    lj3[j][i] = lj3[i][j];
    lj4[j][i] = lj4[i][j];
    ljsw0[j][i] = ljsw0[i][j];
    ljsw1[j][i] = ljsw1[i][j];
    ljsw2[j][i] = ljsw2[i][j];
    ljsw3[j][i] = ljsw3[i][j];
    ljsw4[j][i] = ljsw4[i][j];
    offset[j][i] = offset[i][j];
    ljf1[j][i] = ljf1[i][j];
    ljf2[j][i] = ljf2[i][j];
    ljf3[j][i] = ljf3[i][j];
    ljf4[j][i] = ljf4[i][j];
    ljswf0[j][i] = ljswf0[i][j];
    ljswf1[j][i] = ljswf1[i][j];
    ljswf2[j][i] = ljswf2[i][j];
    ljswf3[j][i] = ljswf3[i][j];
    ljswf4[j][i] = ljswf4[i][j];
    offsetsp[j][i] = offsetsp[i][j];
    offsetsm[j][i] = offsetsm[i][j];
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJRelRes::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilonf[i][j],sizeof(double),1,fp);
        fwrite(&sigmaf[i][j],sizeof(double),1,fp);
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cutf_inner[i][j],sizeof(double),1,fp);
        fwrite(&cutf[i][j],sizeof(double),1,fp);
        fwrite(&cut_inner[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJRelRes::read_restart(FILE *fp)
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
          utils::sfread(FLERR,&epsilonf[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&sigmaf[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&epsilon[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&sigma[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cutf_inner[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cutf[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_inner[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&epsilonf[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigmaf[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cutf_inner[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cutf[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_inner[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJRelRes::write_restart_settings(FILE *fp)
{
  fwrite(&cutf_inner_global,sizeof(double),1,fp);
  fwrite(&cutf_global,sizeof(double),1,fp);
  fwrite(&cut_inner_global,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJRelRes::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR,&cutf_inner_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cutf_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_inner_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cutf_inner_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cutf_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_inner_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJRelRes::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,epsilonf[i][i],sigmaf[i][i],
            epsilon[i][i],sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJRelRes::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g %g %g %g\n",i,j,
              epsilonf[i][j],sigmaf[i][j],epsilon[i][j],sigma[i][j],
              cutf_inner[i][j],cutf[i][j],cut_inner[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */
double PairLJRelRes::single(int /*i*/, int /*j*/, int itype, int jtype,
                            double rsq, double /*factor_coul*/,
                            double factor_lj, double &fforce)
{
  double r2inv,r6inv,forcelj,philj,r,t,tsq,fskin;

  r2inv = 1.0/rsq;
  if (rsq < cutf_inner_sq[itype][jtype]) {
    r6inv = r2inv*r2inv*r2inv;
    forcelj = r6inv*(ljf1[itype][jtype]*r6inv-ljf2[itype][jtype]);
  } else if (rsq < cutfsq[itype][jtype]) {
    r = sqrt(rsq);
    t = r - cutf_inner[itype][jtype];
    tsq = t*t;
    fskin = ljswf1[itype][jtype]+ljswf2[itype][jtype]*t+
      ljswf3[itype][jtype]*tsq+ljswf4[itype][jtype]*tsq*t;
    forcelj = fskin*r;
  } else if (rsq < cut_inner_sq[itype][jtype]) {
    r6inv = r2inv*r2inv*r2inv;
    forcelj = r6inv * (lj1[itype][jtype]*r6inv-lj2[itype][jtype]);
  } else {
    r = sqrt(rsq);
    t = r - cut_inner[itype][jtype];
    tsq = t*t;
    fskin = ljsw1[itype][jtype] + ljsw2[itype][jtype]*t +
      ljsw3[itype][jtype]*tsq + ljsw4[itype][jtype]*tsq*t;
    forcelj = fskin*r;
  }
  fforce = factor_lj*forcelj*r2inv;

  if (rsq < cutf_inner_sq[itype][jtype]) {
    philj = r6inv*(ljf3[itype][jtype]*r6inv-ljf4[itype][jtype])-offsetsm[itype][jtype];
  } else if (rsq < cutfsq[itype][jtype]) {
    philj = ljswf0[itype][jtype]-ljswf1[itype][jtype]*t-
      ljswf2[itype][jtype]*tsq/2.0-ljswf3[itype][jtype]*tsq*t/3.0-
      ljswf4[itype][jtype]*tsq*tsq/4.0-offsetsp[itype][jtype];
  } else if (rsq < cut_inner_sq[itype][jtype]) {
    philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]) - offset[itype][jtype];
  } else {
    philj = ljsw0[itype][jtype] - ljsw1[itype][jtype]*t -
      ljsw2[itype][jtype]*tsq/2.0 - ljsw3[itype][jtype]*tsq*t/3.0 -
      ljsw4[itype][jtype]*tsq*tsq/4.0 - offset[itype][jtype];
  }
  return factor_lj*philj;
}
