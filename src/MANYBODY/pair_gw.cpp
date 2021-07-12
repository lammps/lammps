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
   Contributing author: German Samolyuk (ORNL)
   based on PairTersoff by Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "pair_gw.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "potential_file_reader.h"
#include "tokenizer.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathExtra;

#define DELTA 4

/* ---------------------------------------------------------------------- */

PairGW::PairGW(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  unit_convert_flag = utils::get_supported_conversions(utils::ENERGY);

  params = nullptr;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairGW::~PairGW()
{
  memory->destroy(params);
  memory->destroy(elem3param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairGW::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum;
  int itag,jtag,itype,jtype,ktype,iparam_ij,iparam_ijk;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fi[3],fj[3],fk[3];
  double zeta_ij, prefactor;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // two-body interactions, skip half of them

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtag = tag[j];

      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < x[i][2]) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      jtype = map[type[j]];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      iparam_ij = elem3param[itype][jtype][jtype];
      if (rsq > params[iparam_ij].cutsq) continue;

      repulsive(&params[iparam_ij],rsq,fpair,eflag,evdwl);

      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,fpair,delx,dely,delz);
    }

    // three-body interactions
    // skip immediately if I-J is not within cutoff

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];
      iparam_ij = elem3param[itype][jtype][jtype];

      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 > params[iparam_ij].cutsq) continue;

      // accumulate bondorder zeta for each i-j interaction via loop over k

      zeta_ij = 1.0;

      for (kk = 0; kk < jnum; kk++) {
        if (jj == kk) continue;
        k = jlist[kk];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        iparam_ijk = elem3param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 > params[iparam_ijk].cutsq) continue;

        zeta_ij += zeta(&params[iparam_ijk],rsq1,rsq2,delr1,delr2);
      }

      // pairwise force due to zeta

      force_zeta(&params[iparam_ij],rsq1,zeta_ij,fpair,prefactor,eflag,evdwl);

      f[i][0] += delr1[0]*fpair;
      f[i][1] += delr1[1]*fpair;
      f[i][2] += delr1[2]*fpair;
      f[j][0] -= delr1[0]*fpair;
      f[j][1] -= delr1[1]*fpair;
      f[j][2] -= delr1[2]*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,-fpair,-delr1[0],-delr1[1],-delr1[2]);

      // attractive term via loop over k

      for (kk = 0; kk < jnum; kk++) {
        if (jj == kk) continue;
        k = jlist[kk];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        iparam_ijk = elem3param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 > params[iparam_ijk].cutsq) continue;

        attractive(&params[iparam_ijk],prefactor,
                   rsq1,rsq2,delr1,delr2,fi,fj,fk);

        f[i][0] += fi[0];
        f[i][1] += fi[1];
        f[i][2] += fi[2];
        f[j][0] += fj[0];
        f[j][1] += fj[1];
        f[j][2] += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];

        if (vflag_either) v_tally3(i,j,k,fj,fk,delr1,delr2);
      } // kk
    } // jj
  } // ii

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairGW::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGW::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGW::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  map_element2type(narg-3,arg+3);

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup_params();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGW::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style GW requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style GW requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGW::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairGW::read_file(char *file)
{
  memory->sfree(params);
  params = nullptr;
  nparams = maxparam = 0;

  // open file on proc 0

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, file, "gw", unit_convert_flag);
    char * line;

    // transparently convert units for supported conversions

    int unit_convert = reader.get_unit_convert();
    double conversion_factor = utils::get_conversion_factor(utils::ENERGY,
                                                            unit_convert);
    while ((line = reader.next_line(NPARAMS_PER_LINE))) {
      try {
        ValueTokenizer values(line);

        std::string iname = values.next_string();
        std::string jname = values.next_string();
        std::string kname = values.next_string();

        // ielement,jelement,kelement = 1st args
        // if all 3 args are in element list, then parse this line
        // else skip to next entry in file
        int ielement, jelement, kelement;

        for (ielement = 0; ielement < nelements; ielement++)
          if (iname == elements[ielement]) break;
        if (ielement == nelements) continue;
        for (jelement = 0; jelement < nelements; jelement++)
          if (jname == elements[jelement]) break;
        if (jelement == nelements) continue;
        for (kelement = 0; kelement < nelements; kelement++)
          if (kname == elements[kelement]) break;
        if (kelement == nelements) continue;

        // load up parameter settings and error check their values

        if (nparams == maxparam) {
          maxparam += DELTA;
          params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                              "pair:params");

          // make certain all addional allocated storage is initialized
          // to avoid false positives when checking with valgrind

          memset(params + nparams, 0, DELTA*sizeof(Param));
        }

        params[nparams].ielement = ielement;
        params[nparams].jelement = jelement;
        params[nparams].kelement = kelement;
        params[nparams].powerm = values.next_double();
        params[nparams].gamma  = values.next_double();
        params[nparams].lam3   = values.next_double();
        params[nparams].c      = values.next_double();
        params[nparams].d      = values.next_double();
        params[nparams].h      = values.next_double();
        params[nparams].powern = values.next_double();
        params[nparams].beta   = values.next_double();
        params[nparams].lam2   = values.next_double();
        params[nparams].bigb   = values.next_double();
        params[nparams].bigr   = values.next_double();
        params[nparams].bigd   = values.next_double();
        params[nparams].lam1   = values.next_double();
        params[nparams].biga   = values.next_double();
        params[nparams].powermint = int(params[nparams].powerm);

        if (unit_convert) {
          params[nparams].biga *= conversion_factor;
          params[nparams].bigb *= conversion_factor;
        }
      } catch (TokenizerException &e) {
        error->one(FLERR, e.what());
      }

      // currently only allow m exponent of 1 or 3
      if (params[nparams].c < 0.0 || params[nparams].d < 0.0 ||
          params[nparams].powern < 0.0 || params[nparams].beta < 0.0 ||
          params[nparams].lam2 < 0.0 || params[nparams].bigb < 0.0 ||
          params[nparams].bigr < 0.0 ||params[nparams].bigd < 0.0 ||
          params[nparams].bigd > params[nparams].bigr ||
          params[nparams].lam1 < 0.0 || params[nparams].biga < 0.0 ||
          params[nparams].powerm - params[nparams].powermint != 0.0 ||
          (params[nparams].powermint != 3 && params[nparams].powermint != 1) ||
          params[nparams].gamma < 0.0)
        error->one(FLERR,"Illegal GW parameter");

      nparams++;
    }
  }

  MPI_Bcast(&nparams, 1, MPI_INT, 0, world);
  MPI_Bcast(&maxparam, 1, MPI_INT, 0, world);

  if (comm->me != 0) {
    params = (Param *) memory->srealloc(params,maxparam*sizeof(Param), "pair:params");
  }

  MPI_Bcast(params, maxparam*sizeof(Param), MPI_BYTE, 0, world);
}

/* ---------------------------------------------------------------------- */

void PairGW::setup_params()
{
  int i,j,k,m,n;

  // set elem3param for all element triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem3param);
  memory->create(elem3param,nelements,nelements,nelements,"pair:elem3param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement &&
              k == params[m].kelement) {
            if (n >= 0)
              error->all(FLERR,"Potential file has duplicate entry");
            n = m;
          }
        }
        if (n < 0) error->all(FLERR,"Potential file is missing an entry");
        elem3param[i][j][k] = n;
      }


  // compute parameter values derived from inputs

  for (m = 0; m < nparams; m++) {
    params[m].cut = params[m].bigr + params[m].bigd;
    params[m].cutsq = params[m].cut*params[m].cut;

    params[m].c1 = pow(2.0*params[m].powern*1.0e-16,-1.0/params[m].powern);
    params[m].c2 = pow(2.0*params[m].powern*1.0e-8,-1.0/params[m].powern);
    params[m].c3 = 1.0/params[m].c2;
    params[m].c4 = 1.0/params[m].c1;
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++)
    if (params[m].cut > cutmax) cutmax = params[m].cut;
}

/* ---------------------------------------------------------------------- */

void PairGW::repulsive(Param *param, double rsq, double &fforce,
                            int eflag, double &eng)
{
  double r,tmp_fc,tmp_fc_d,tmp_exp;

  r = sqrt(rsq);
  tmp_fc = gw_fc(r,param);
  tmp_fc_d = gw_fc_d(r,param);
  tmp_exp = exp(-param->lam1 * r);
  fforce = -param->biga * tmp_exp * (tmp_fc_d - tmp_fc*param->lam1) / r;
  if (eflag) eng = tmp_fc * param->biga * tmp_exp;
}

/* ---------------------------------------------------------------------- */

double PairGW::zeta(Param *param, double rsqij, double rsqik,
                         double *delrij, double *delrik)
{
  double rij,rik,costheta,arg,ex_delr;

  rij = sqrt(rsqij);
  rik = sqrt(rsqik);
  costheta = (delrij[0]*delrik[0] + delrij[1]*delrik[1] +
              delrij[2]*delrik[2]) / (rij*rik);

  if (param->powermint == 3) arg = pow(param->lam3 * (rij-rik),3.0);
  else arg = param->lam3 * (rij-rik);

  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(arg);

  return gw_fc(rik,param) * gw_gijk(costheta,param) * ex_delr;
}

/* ---------------------------------------------------------------------- */

void PairGW::force_zeta(Param *param_i, double rsq, double zeta_ij,
                             double &fforce, double &prefactor,
                             int eflag, double &eng)
{
  double r,fa,fa_d,bij;

  r = sqrt(rsq);
  fa = gw_fa(r,param_i);
  fa_d = gw_fa_d(r,param_i);
  bij = gw_bij(zeta_ij,param_i);
  fforce = 0.5*bij*fa_d / r;
  prefactor = -0.5*fa * gw_bij_d(zeta_ij,param_i);
  if (eflag) eng = 0.5*bij*fa;
}

/* ----------------------------------------------------------------------
   attractive term
   use param_ij cutoff for rij test
   use param_ijk cutoff for rik test
------------------------------------------------------------------------- */

void PairGW::attractive(Param *param, double prefactor,
                             double rsqij, double rsqik,
                             double *delrij, double *delrik,
                             double *fi, double *fj, double *fk)
{
  double rij_hat[3],rik_hat[3];
  double rij,rijinv,rik,rikinv;

  rij = sqrt(rsqij);
  rijinv = 1.0/rij;
  scale3(rijinv,delrij,rij_hat);

  rik = sqrt(rsqik);
  rikinv = 1.0/rik;
  scale3(rikinv,delrik,rik_hat);

  gw_zetaterm_d(prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk,param);
}

/* ---------------------------------------------------------------------- */

double PairGW::gw_fc(double r, Param *param)
{
  double gw_R = param->bigr;
  double gw_D = param->bigd;

  if (r < gw_R-gw_D) return 1.0;
  if (r > gw_R+gw_D) return 0.0;
  return 0.5*(1.0 - sin(MY_PI2*(r - gw_R)/gw_D));
}

/* ---------------------------------------------------------------------- */

double PairGW::gw_fc_d(double r, Param *param)
{
  double gw_R = param->bigr;
  double gw_D = param->bigd;

  if (r < gw_R-gw_D) return 0.0;
  if (r > gw_R+gw_D) return 0.0;
  return -(MY_PI4/gw_D) * cos(MY_PI2*(r - gw_R)/gw_D);
}

/* ---------------------------------------------------------------------- */

double PairGW::gw_fa(double r, Param *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return -param->bigb * exp(-param->lam2 * r) * gw_fc(r,param);
}

/* ---------------------------------------------------------------------- */

double PairGW::gw_fa_d(double r, Param *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return param->bigb * exp(-param->lam2 * r) *
    (param->lam2 * gw_fc(r,param) - gw_fc_d(r,param));
}

/* ---------------------------------------------------------------------- */

double PairGW::gw_bij(double zeta_ij, Param *param_i)
{
  double tmp = param_i->beta * zeta_ij;
  return pow(tmp,-param_i->powern);
}

/* ---------------------------------------------------------------------- */

double PairGW::gw_bij_d(double zeta_ij, Param *param_i)
{
  double tmp = param_i->beta * zeta_ij;
  return - param_i->powern * pow(tmp,-param_i->powern-1)*tmp / zeta_ij;
}

/* ---------------------------------------------------------------------- */

void PairGW::gw_zetaterm_d(double prefactor,
                                  double *rij_hat, double rij,
                                  double *rik_hat, double rik,
                                  double *dri, double *drj, double *drk,
                                  Param *param)
{
  double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = gw_fc(rik,param);
  dfc = gw_fc_d(rik,param);
  if (param->powermint == 3) tmp = pow(param->lam3 * (rij-rik),3.0);
  else tmp = param->lam3 * (rij-rik);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (param->powermint == 3)
    ex_delr_d = 3.0*pow(param->lam3,3.0) * pow(rij-rik,2.0)*ex_delr;
  else ex_delr_d = param->lam3 * ex_delr;

  cos_theta = dot3(rij_hat,rik_hat);
  gijk = gw_gijk(cos_theta,param);
  gijk_d = gw_gijk_d(cos_theta,param);
  costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

  // compute the derivative wrt Ri
  // dri = -dfc*gijk*ex_delr*rik_hat;
  // dri += fc*gijk_d*ex_delr*dcosdri;
  // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat);

  scale3(-dfc*gijk*ex_delr,rik_hat,dri);
  scaleadd3(fc*gijk_d*ex_delr,dcosdri,dri,dri);
  scaleadd3(fc*gijk*ex_delr_d,rik_hat,dri,dri);
  scaleadd3(-fc*gijk*ex_delr_d,rij_hat,dri,dri);
  scale3(prefactor,dri);

  // compute the derivative wrt Rj
  // drj = fc*gijk_d*ex_delr*dcosdrj;
  // drj += fc*gijk*ex_delr_d*rij_hat;

  scale3(fc*gijk_d*ex_delr,dcosdrj,drj);
  scaleadd3(fc*gijk*ex_delr_d,rij_hat,drj,drj);
  scale3(prefactor,drj);

  // compute the derivative wrt Rk
  // drk = dfc*gijk*ex_delr*rik_hat;
  // drk += fc*gijk_d*ex_delr*dcosdrk;
  // drk += -fc*gijk*ex_delr_d*rik_hat;

  scale3(dfc*gijk*ex_delr,rik_hat,drk);
  scaleadd3(fc*gijk_d*ex_delr,dcosdrk,drk,drk);
  scaleadd3(-fc*gijk*ex_delr_d,rik_hat,drk,drk);
  scale3(prefactor,drk);
}

/* ---------------------------------------------------------------------- */

void PairGW::costheta_d(double *rij_hat, double rij,
                             double *rik_hat, double rik,
                             double *dri, double *drj, double *drk)
{
  // first element is devative wrt Ri, second wrt Rj, third wrt Rk

  double cos_theta = dot3(rij_hat,rik_hat);

  scaleadd3(-cos_theta,rij_hat,rik_hat,drj);
  scale3(1.0/rij,drj);
  scaleadd3(-cos_theta,rik_hat,rij_hat,drk);
  scale3(1.0/rik,drk);
  add3(drj,drk,dri);
  scale3(-1.0,dri);
}
