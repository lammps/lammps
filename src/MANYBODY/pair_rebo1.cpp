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
   Contributing author: Ray Shan (Materials Design)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_rebo1.h"
#include "atom.h"
#include "update.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include "math_const.h"
using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

PairREBO1::PairREBO1(LAMMPS *lmp) : PairTersoff(lmp)
{
  moliere_aB = 0.5292;
  moliere_c[0] = 0.35;
  moliere_c[1] = 0.55;
  moliere_c[2] = 0.10;
  moliere_d[0] = 0.30;
  moliere_d[1] = 1.20;
  moliere_d[2] = 6.00;

  firsov_const = pow(9.0*MY_PI*MY_PI/128.0, 1.0/3.0);
  NCl = NSi = NULL;
  nmax = 0;

  // set comm size needed by this Pair
  comm_forward = 1;
  comm_reverse = 1;

  // read spline coefficients for coordination term, A12

  read_lib();
}

/* ---------------------------------------------------------------------- */

void PairREBO1::read_file(char *file)
{
  int params_per_line = 27;
  char **words = new char*[params_per_line+1];

  memory->sfree(params);
  params = NULL;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = force->open_potential(file);
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open REBO1 potential file %s",file);
      error->one(FLERR,str);
    }
  }

  // read each line out of file, skipping blank lines or leading '#'
  // store line of params if all 3 element tags are in element list

  int n,nwords,ielement,jelement,kelement;
  char line[MAXLINE],*ptr;
  int eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR,"Incorrect format in REBO1 potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

    // ielement,jelement,kelement = 1st args
    // if all 3 args are in element list, then parse this line
    // else skip to next line

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0],elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1],elements[jelement]) == 0) break;
    if (jelement == nelements) continue;
    for (kelement = 0; kelement < nelements; kelement++)
      if (strcmp(words[2],elements[kelement]) == 0) break;
    if (kelement == nelements) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      maxparam += DELTA;
     params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                          "pair:params");
    }

    params[nparams].ielement = ielement;
    params[nparams].jelement = jelement;
    params[nparams].kelement = kelement;
    params[nparams].powerm = atof(words[3]);
    params[nparams].gamma = atof(words[4]);
    params[nparams].lam3 = atof(words[5]);
    params[nparams].c = atof(words[6]);
    params[nparams].d = atof(words[7]);
    params[nparams].h = atof(words[8]);
    params[nparams].powern = atof(words[9]);
    params[nparams].beta = atof(words[10]);
    params[nparams].lam2 = atof(words[11]);
    params[nparams].bigb = atof(words[12]);
    params[nparams].bigr = atof(words[13]);
    params[nparams].bigd = atof(words[14]);
    params[nparams].lam1 = atof(words[15]);
    params[nparams].biga = atof(words[16]);
    params[nparams].powereta = atof(words[17]);
    params[nparams].Z_i = atof(words[18]);
    params[nparams].Z_j = atof(words[19]);
    params[nparams].spl_ra = atof(words[20]);
    params[nparams].spl_rb = atof(words[21]);
    params[nparams].spl_a = atof(words[22]);
    params[nparams].spl_b = atof(words[23]);
    params[nparams].spl_c = atof(words[24]);
    params[nparams].spl_s = atof(words[25]);
    params[nparams].Re = atof(words[26]);

    // currently only allow m exponent of 1 or 3

    params[nparams].powermint = int(params[nparams].powerm);

    if (
        params[nparams].lam3 < 0.0 || params[nparams].c < 0.0 ||
        params[nparams].d < 0.0 || params[nparams].powern < 0.0 ||
        params[nparams].beta < 0.0 || params[nparams].lam2 < 0.0 ||
        params[nparams].bigb < 0.0 || params[nparams].bigr < 0.0 ||
        params[nparams].bigd < 0.0 ||
        params[nparams].bigd > params[nparams].bigr ||
        params[nparams].lam3 < 0.0 || params[nparams].biga < 0.0 ||
        params[nparams].powerm - params[nparams].powermint != 0.0 ||
        params[nparams].gamma < 0.0 ||
        params[nparams].Z_i < 1.0 || params[nparams].Z_j < 1.0 ||
        params[nparams].spl_ra < 0.0 || params[nparams].spl_rb < 0.0)
      error->all(FLERR,"Illegal REBO1 parameter");

    nparams++;
  }

  delete [] words;
}

/* ---------------------------------------------------------------------- */

void PairREBO1::read_lib()
{
  unsigned int maxlib = 1024;
  int i,j,k,l,nwords,m;
  int ii,jj,kk,ll,mm,iii;
  char s[maxlib];
  char **words = new char*[80];
  int nspl = 401;

  // open libraray file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = force->open_potential("lib.rebo1");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open REBO1 lib.rebo1 file");
      error->one(FLERR,str);
    }

    for (i=0; i<4; i++) {
      for (l=0; l<nspl; l++) {
        fgets(s,maxlib,fp);
        nwords = 0;
        words[nwords++] = strtok(s," \t\n\r\f");
        while ((words[nwords++] = strtok(NULL," \t\n\r\f")))continue;
        coordnumber[i][l] = atof(words[0]);
        coordenergy[i][l] = atof(words[1]);
        coordforce[i][l]  = atof(words[2]);
      }
      for (l=0; l<nspl; l++) {
        coordenergy[4][l] = 0.0;
        coordforce[4][l]  = 0.0;
      }
    }
  }

  MPI_Bcast(&coordnumber[0][0],1604,MPI_DOUBLE,0,world);
  MPI_Bcast(&coordenergy[0][0],1604,MPI_DOUBLE,0,world);
  MPI_Bcast(&coordforce[0][0],1604,MPI_DOUBLE,0,world);

  delete [] words;

}

/* ---------------------------------------------------------------------- */

void PairREBO1::count_neigh()
{
  int i, j, ii, jj, itype, jtype, n;
  int inum, jnum, param, *ilist, *jlist, *numneigh, **firstneigh;
  double r, rsq, delrij[3];
  const double cutshortsq = cutmax*cutmax;

  double **x = atom->x;
  int *type = atom->type;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->grow(NCl,nmax,"pair:NCl");
    memory->grow(NSi,nmax,"pair:NSi");
  }

  inum  = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = map[type[i]];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // skip immediately if center atom is not Si
    if (strcmp(elements[map[itype+1]],"Si") != 0) continue;

    NCl[i] = NSi[i] = 0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj] & NEIGHMASK;
      jtype = map[type[j]];


      delrij[0] = x[i][0] - x[j][0];
      delrij[1] = x[i][1] - x[j][1];
      delrij[2] = x[i][2] - x[j][2];
      rsq = vec3_dot(delrij,delrij);
      param = elem2param[itype][jtype][jtype];

      if (rsq > cutshortsq) continue;

      r = sqrt(rsq);
      if (strcmp(elements[map[jtype+1]],"Cl") == 0) NCl[i] += ters_fc(r,&params[param]);
      if (strcmp(elements[map[jtype+1]],"Si") == 0) NSi[i] += ters_fc(r,&params[param]);
    }
  }

  // communicating coordination number to all nodes
  pack_flag = 1;
  comm->forward_comm_pair(this);
  pack_flag = 2;
  comm->forward_comm_pair(this);

}

/* ---------------------------------------------------------------------- */

void PairREBO1::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum;
  int itype,jtype,ktype,iparam_ij,iparam_ijk,iparam_ik;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fi[3],fj[3],fk[3];
  double zeta_ij,prefactor;
  int nncl, nnsi;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = vflag_atom = 0;

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  const double cutshortsq = cutmax*cutmax;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // count number of nearest neighbors; needs communications
  count_neigh();

  double fxtmp,fytmp,fztmp;

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp = fytmp = fztmp = 0.0;

    // two-body interactions, skip half of them

    jlist = firstneigh[i];
    jnum = numneigh[i];
    int numshort = 0;

    nncl = int((NCl[i]-1.0)*100.0);
    nnsi = int(NSi[i]);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutshortsq) {
        neighshort[numshort++] = j;
        if (numshort >= maxshort) {
          maxshort += maxshort/2;
          memory->grow(neighshort,maxshort,"pair:neighshort");
        }
      }

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
      iparam_ij = elem2param[itype][jtype][jtype];
      if (rsq >= params[iparam_ij].cutsq) continue;

      repulsive(&params[iparam_ij],rsq,fpair,eflag,evdwl);

      fxtmp += delx*fpair;
      fytmp += dely*fpair;
      fztmp += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,fpair,delx,dely,delz);
    }

    // three-body interactions
    // skip immediately if I-J is not within cutoff
    double fjxtmp,fjytmp,fjztmp;

    for (jj = 0; jj < numshort; jj++) {
      j = neighshort[jj];
      jtype = map[type[j]];
      iparam_ij = elem2param[itype][jtype][jtype];

      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 >= params[iparam_ij].cutsq) continue;

      // accumulate bondorder zeta for each i-j interaction via loop over k

      fjxtmp = fjytmp = fjztmp = 0.0;
      zeta_ij = 0.0;

      for (kk = 0; kk < numshort; kk++) {
        if (jj == kk) continue;
        k = neighshort[kk];
        ktype = map[type[k]];
        iparam_ijk = elem2param[itype][jtype][ktype];
	iparam_ik  = elem2param[itype][ktype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 >= params[iparam_ijk].cutsq) continue;

        zeta_ij += zeta(&params[iparam_ijk],&params[iparam_ik],rsq1,rsq2,delr1,delr2);
      }

      // pairwise force due to zeta

      force_zeta(&params[iparam_ij],rsq1,zeta_ij,fpair,prefactor,eflag,evdwl,nncl,nnsi);

      fxtmp += delr1[0]*fpair;
      fytmp += delr1[1]*fpair;
      fztmp += delr1[2]*fpair;
      fjxtmp -= delr1[0]*fpair;
      fjytmp -= delr1[1]*fpair;
      fjztmp -= delr1[2]*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,-fpair,-delr1[0],-delr1[1],-delr1[2]);

      // attractive term via loop over k

      for (kk = 0; kk < numshort; kk++) {
        if (jj == kk) continue;
        k = neighshort[kk];
        ktype = map[type[k]];
        iparam_ijk = elem2param[itype][jtype][ktype];
	iparam_ik  = elem2param[itype][ktype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 >= params[iparam_ijk].cutsq) continue;

        attractive(&params[iparam_ijk],&params[iparam_ik],prefactor,
                   rsq1,rsq2,delr1,delr2,fi,fj,fk);

        fxtmp += fi[0];
        fytmp += fi[1];
        fztmp += fi[2];
        fjxtmp += fj[0];
        fjytmp += fj[1];
        fjztmp += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];

        if (vflag_atom) v_tally3(i,j,k,fj,fk,delr1,delr2);
      }
      f[j][0] += fjxtmp;
      f[j][1] += fjytmp;
      f[j][2] += fjztmp;
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairREBO1::repulsive(Param *param, double rsq, double &fforce,
                               int eflag, double &eng)
{
  double r,rinv,tmp_fc,tmp_fc_d,tmp_exp;
  double energy, forces;
  r = sqrt(rsq);
  rinv = 1.0/r;

  double spl_ra = param->spl_ra;
  double spl_rb = param->spl_rb;

  if (r > spl_rb) {
  // Tersoff repulsive (of Morse form): when r > rb; Eq. A3
 
    tmp_fc = ters_fc(r,param);
    tmp_fc_d = ters_fc_d(r,param);
    tmp_exp = exp(-param->lam1 * r);
    forces = param->biga * tmp_exp * (tmp_fc_d - tmp_fc*param->lam1);
    energy = tmp_fc * param->biga * tmp_exp;

  } else if (r < spl_rb && r > spl_ra) {
  // Spline repulsive: when ra < r < rb; Eq. A9

    double spl_a = param->spl_a;
    double spl_b = param->spl_b;
    double spl_c = param->spl_c;
    energy = spl_c + exp(spl_a * r + spl_b);
    forces = spl_a * exp(spl_a * r + spl_b);

  } else if (r < spl_ra) {
  // Moliere repulsive: when r < ra; Eqa. A7-8

    double zi = param->Z_i;
    double zj = param->Z_j;
    double spl_s = param->spl_s;

    // Eq. A7, first term
    double zizj = zi * zj * force->qqr2e * force->qelectron * force->qelectron * rinv;

    // Eq. A8
    double zizj_half = pow(zi,0.5) + pow(zj,0.5);
    double screen_a = 0.83 * firsov_const * moliere_aB * pow(zizj_half,-2.0/3.0);

    // Eq. A7, second term
    double cexpdra = 0.0;
    for (int i = 0; i < 3; i++)
      cexpdra += moliere_c[i] * exp(-1.0 * moliere_d[i] * r / screen_a);

    energy = zizj * cexpdra + spl_s;

    // forces: derivative of first term, Eq. A7
    double r2inv = 1.0/rsq;
    double dzizj = -1.0 * zizj * r2inv;
  
    // forces: derivative of second term, Eq. A7
    double dcexpdra = 0.0;
    for (int i = 0; i < 3; i++)
      dcexpdra += -moliere_c[i] * moliere_d[i] / screen_a * exp(-moliere_d[i] * r / screen_a);

    forces = dzizj * cexpdra + zizj * dcexpdra;
  }

  fforce = -forces / r;
  if (eflag) eng = energy;
}

/* ---------------------------------------------------------------------- */

double PairREBO1::ters_fa(double r, Param *param)
{
  // Tersoff attraction term; Eq. A4

  if (r > param->bigr + param->bigd) return 0.0;
  return -param->bigb * exp(-param->lam2 * r) * ters_fc(r,param);
}

/* ---------------------------------------------------------------------- */

double PairREBO1::ters_fa_d(double r, Param *param)
{
  // derivative of Tersoff attraction; Eq. A4
  
  if (r > param->bigr + param->bigd) return 0.0;
  return param->bigb * exp(-param->lam2 * r) *
    (param->lam2 * ters_fc(r,param) - ters_fc_d(r,param));
}

/* ---------------------------------------------------------------------- */

double PairREBO1::ters_fc(double r, Param *param)
{
  // Teroff cutoff function; Eq. A6
  
  double ters_R = param->bigr;	// Rmin
  double ters_D = param->bigd;	

  if (r < ters_R) return 1.0;
  if (r > ters_R+ters_D) return 0.0;
  double t = (r - ters_R)/ters_D;
  double t2 = t*t;
  double t3 = t2*t;
  return 1.0 - t3*(6.0*t2-15.0*t+10.0);
}

/* ---------------------------------------------------------------------- */

double PairREBO1::ters_fc_d(double r, Param *param)
{
  // derivative of Tersoff cutoff function; Eq. A6
  
  double ters_R = param->bigr;
  double ters_D = param->bigd;

  if (r < ters_R) return 0.0;
  if (r > ters_R+ters_D) return 0.0;
  double t = (r - ters_R);
  double t2 = t*t;
  return -30.0*t2*(t-ters_D)*(t-ters_D)/pow(ters_D,5);

}

/* ---------------------------------------------------------------------- */

double PairREBO1::zeta(Param *paramij, Param *paramik, double rsqij, double rsqik,
                         double *delrij, double *delrik)
{
  // zeta term inside the bond order term; Eq. A13
  
  double rij,rik,cos_theta,tmp,ex_delr;
  double Re_ij = paramij->Re;
  double Re_ik = paramik->Re;

  rij = sqrt(rsqij);
  rik = sqrt(rsqik);
  cos_theta = (delrij[0]*delrik[0] + delrij[1]*delrik[1] +
              delrij[2]*delrik[2]) / (rij*rik);

  tmp = paramij->lam3 * pow(((rij-Re_ij)-(rik-Re_ik)),paramij->powermint);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  ex_delr = exp(tmp);

  return ters_fc(rik,paramij) * ters_gijk(cos_theta,paramij) * ex_delr;
}

/* ---------------------------------------------------------------------- */

void PairREBO1::force_zeta(Param *param, double rsq, double zeta_ij,
                             double &fforce, double &prefactor,
                             int eflag, double &eng, int nncl, int nnsi)
{
  double r,fa,fa_d,bij;

  r = sqrt(rsq);
  fa = ters_fa(r,param);
  fa_d = ters_fa_d(r,param);
  bij = ters_bij(zeta_ij,param,nncl,nnsi);
  fforce = 0.5*bij*fa_d / r;
  prefactor = -0.5*fa * ters_bij_d(zeta_ij,param,nncl,nnsi);
  if (eflag) eng = 0.5*bij*fa;
}

/* ---------------------------------------------------------------------- */

void PairREBO1::ters_zetaterm_d(Param *paramij, Param *paramik, double prefactor,
                                  double *rij_hat, double rij,
                                  double *rik_hat, double rik,
                                  double *dri, double *drj, double *drk)
{
  // derivate of zeta term inside the bond order term; Eq. A13
  
  double gijk,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];
  double gijkdri[3],gijkdrj[3],gijkdrk[3];
  double Re_ij = paramij->Re;
  double Re_ik = paramik->Re;

  fc = ters_fc(rik,paramij);
  dfc = ters_fc_d(rik,paramij);
  tmp = paramij->lam3 * pow(((rij-Re_ij)-(rik-Re_ik)),paramij->powermint);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  ex_delr_d = paramij->lam3 * paramij->powermint * pow(((rij-Re_ij)-(rik-Re_ik)),paramij->powermint-1) * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk(cos_theta,paramij);
  ters_gijk_d(paramij,rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

  // compute the derivative wrt Ri
  // dri = -dfc*gijk*ex_delr*rik_hat;
  // dri += fc*gijk_d*ex_delr*dcosdri;
  // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat);

  vec3_scale(-dfc*gijk*ex_delr,rik_hat,dri);
  vec3_scaleadd(fc*ex_delr,dcosdri,dri,dri);
  vec3_scaleadd(fc*gijk*ex_delr_d,rik_hat,dri,dri);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rij_hat,dri,dri);
  vec3_scale(prefactor,dri,dri);

  // compute the derivative wrt Rj
  // drj = fc*gijk_d*ex_delr*dcosdrj;
  // drj += fc*gijk*ex_delr_d*rij_hat;

  vec3_scale(fc*ex_delr,dcosdrj,drj);
  vec3_scaleadd(fc*gijk*ex_delr_d,rij_hat,drj,drj);
  vec3_scale(prefactor,drj,drj);

  // compute the derivative wrt Rk
  // drk = dfc*gijk*ex_delr*rik_hat;
  // drk += fc*gijk_d*ex_delr*dcosdrk;
  // drk += -fc*gijk*ex_delr_d*rik_hat;

  vec3_scale(dfc*gijk*ex_delr,rik_hat,drk);
  vec3_scaleadd(fc*ex_delr,dcosdrk,drk,drk);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rik_hat,drk,drk);
  vec3_scale(prefactor,drk,drk);
}

/* ---------------------------------------------------------------------- */

void PairREBO1::ters_gijk_d(Param *param, double *rij_hat, double rij,
                             double *rik_hat, double rik,
                             double *dri, double *drj, double *drk)
{
  // derivative of cos(theta) in gijk term; Eq. A14
  // first element is devative wrt Ri, second wrt Rj, third wrt Rk

  double cos_theta = vec3_dot(rij_hat,rik_hat);

  const double ters_c = param->c;
  const double ters_d = param->d;
  const double hcth = param->h - cos_theta;
  const double tmp = 2.0 * ters_d * hcth;

  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,drj);
  vec3_scale(-tmp/rij,drj,drj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,drk);
  vec3_scale(-tmp/rik,drk,drk);
  vec3_add(drj,drk,dri);
  vec3_scale(-1.0,dri,dri);
}

/* ---------------------------------------------------------------------- */

double PairREBO1::ters_bij(double zeta, Param *param, int nncl, int nnsi)
{
  // REBO1 bond order term; Eq. A12

  double Hcoord;
  if (nnsi > 4) nnsi = 4;
  if (nnsi > 0 && nncl > 0) Hcoord = coordenergy[nnsi][nncl];
  else Hcoord = 0.0;

  double powern = param->powern;
  double powereta = param->powereta;
  double tmp = zeta + Hcoord;
  if (tmp < 0.0) tmp = 1.0e-6;

  return pow(1.0 + pow(tmp,powereta),-powern);
}

/* ---------------------------------------------------------------------- */

double PairREBO1::ters_bij_d(double zeta, Param *param, int nncl, int nnsi)
{
  // derivative of REBO1 bond order term; Eq. A12
  
  double Hcoord;
  if (nnsi > 4) nnsi = 4;
  if (nnsi > 0 && nncl > 0) Hcoord = coordenergy[nnsi][nncl];
  else Hcoord = 0.0;

  double powern = param->powern;
  double powereta = param->powereta;
  double tmp = zeta + Hcoord;
  if (tmp < 0.0) tmp = 1.0e-6;

  return -powereta * powern * pow(tmp,powereta-1) * 
  	pow((1+pow(tmp,powereta)),-powern-1);
}

/* ---------------------------------------------------------------------- */

void PairREBO1::attractive(Param *paramij, Param *paramik, double prefactor,
                             double rsqij, double rsqik,
                             double *delrij, double *delrik,
                             double *fi, double *fj, double *fk)
{
  double rij_hat[3],rik_hat[3];
  double rij,rijinv,rik,rikinv;

  rij = sqrt(rsqij);
  rijinv = 1.0/rij;
  vec3_scale(rijinv,delrij,rij_hat);

  rik = sqrt(rsqik);
  rikinv = 1.0/rik;
  vec3_scale(rikinv,delrik,rik_hat);

  ters_zetaterm_d(paramij, paramik, prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk);
}

/* ---------------------------------------------------------------------- */

int PairREBO1::pack_forward_comm(int n, int *list, double *buf,
                                 int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  if (pack_flag == 1) {
    for (i = 0; i < n; i ++) {
      j = list[i];
      buf[m++] = NCl[j];
    }
  } else if (pack_flag == 2) {
    for (i = 0; i < n; i ++) {
      j = list[i];
      buf[m++] = NSi[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairREBO1::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n ;
  if (pack_flag == 1) {
    for (i = first; i < last; i++)
      NCl[i] = buf[m++];
  } else if (pack_flag == 2) {
    for (i = first; i < last; i++)
      NSi[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int PairREBO1::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (pack_flag == 1) {
    for (i = first; i < last; i++)
      buf[m++] = NCl[i];
  } else if (pack_flag == 2) {
    for (i = first; i < last; i++)
      buf[m++] = NSi[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairREBO1::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  if (pack_flag == 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      NCl[j] += buf[m++];
    }
  } else if (pack_flag == 2) {
    for (i = 0; i < n; i++) {
      j = list[i];
      NSi[j] += buf[m++];
    }
  }
}


