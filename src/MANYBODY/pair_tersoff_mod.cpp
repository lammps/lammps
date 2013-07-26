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
   Contributing author: Aidan Thompson (SNL) - original Tersoff implementation
                        Vitaly Dozhdikov (JIHT of RAS) - MOD addition 
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_tersoff_mod.h"
#include "atom.h"
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

PairTersoffMOD::PairTersoffMOD(LAMMPS *lmp) : PairTersoff(lmp) {}

/* ---------------------------------------------------------------------- */

void PairTersoffMOD::read_file(char *file)
{
  int params_per_line = 20;
  char **words = new char*[params_per_line+1];

  memory->sfree(params);
  params = NULL;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = fopen(file,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open Tersoff potential file %s",file);
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

    if (ptr = strchr(line,'#')) *ptr = '\0';
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
      if (ptr = strchr(line,'#')) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR,"Incorrect format in Tersoff potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while (words[nwords++] = strtok(NULL," \t\n\r\f")) continue;

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
    params[nparams].lam3 = atof(words[4]);
    params[nparams].h = atof(words[5]);
    params[nparams].powern = atof(words[6]);
    params[nparams].beta = atof(words[7]);
    params[nparams].lam2 = atof(words[8]);
    params[nparams].bigb = atof(words[9]);
    params[nparams].bigr = atof(words[10]);
    params[nparams].bigd = atof(words[11]);
    params[nparams].lam1 = atof(words[12]);
    params[nparams].biga = atof(words[13]);
    params[nparams].powern_del = atof(words[14]);
    params[nparams].c1 = atof(words[15]);
    params[nparams].c2 = atof(words[16]);
    params[nparams].c3 = atof(words[17]);
    params[nparams].c4 = atof(words[18]);
    params[nparams].c5 = atof(words[19]);

    // currently only allow m exponent of 1 

    params[nparams].powermint = int(params[nparams].powerm);

    if (
	params[nparams].lam3 < 0.0 || params[nparams].powern < 0.0 || 
	params[nparams].beta < 0.0 || params[nparams].lam2 < 0.0 || 
	params[nparams].bigb < 0.0 || params[nparams].bigr < 0.0 ||
	params[nparams].bigd < 0.0 || 
                               params[nparams].bigd > params[nparams].bigr ||
	params[nparams].lam3 < 0.0 || params[nparams].biga < 0.0 ||
	params[nparams].powerm - params[nparams].powermint != 0.0 ||
    (params[nparams].powermint != 3 && params[nparams].powermint != 1))
      error->all(FLERR,"Illegal Tersoff parameter");

    nparams++;
  }

  delete [] words;
}

/* ---------------------------------------------------------------------- */

void PairTersoffMOD::setup()
{
  int i,j,k,m,n;

  // set elem2param for all element triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem2param);
  memory->create(elem2param,nelements,nelements,nelements,"pair:elem2param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
	n = -1;
	for (m = 0; m < nparams; m++) {
	  if (i == params[m].ielement && j == params[m].jelement && 
	      k == params[m].kelement) {
	    if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
	    n = m;
	  }
	}
	if (n < 0) error->all(FLERR,"Potential file is missing an entry");
	elem2param[i][j][k] = n;
      }


  // compute parameter values derived from inputs

  for (m = 0; m < nparams; m++) {
    params[m].cut = params[m].bigr + params[m].bigd;
    params[m].cutsq = params[m].cut*params[m].cut;

    params[m].ca1 = pow(2.0*params[m].powern_del*1.0e-16,-1.0/params[m].powern);
    params[m].ca4 = 1.0/params[m].ca1;
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++)
    if (params[m].cut > cutmax) cutmax = params[m].cut;
}  

/* ---------------------------------------------------------------------- */

double PairTersoffMOD::zeta(Param *param, double rsqij, double rsqik,
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
  
  return ters_fc(rik,param) * ters_gijk_mod(costheta,param) * ex_delr;
}

/* ---------------------------------------------------------------------- */

double PairTersoffMOD::ters_fc(double r, Param *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;
  
  if (r < ters_R-ters_D) return 1.0;
  if (r > ters_R+ters_D) return 0.0;
  return 0.5*(1.0 - 1.125*sin(MY_PI2*(r - ters_R)/ters_D) - 
              0.125*sin(3*MY_PI2*(r - ters_R)/ters_D));
}

/* ---------------------------------------------------------------------- */

double PairTersoffMOD::ters_fc_d(double r, Param *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;
  
  if (r < ters_R-ters_D) return 0.0;
  if (r > ters_R+ters_D) return 0.0;
  return -(0.375*MY_PI4/ters_D) * (3*cos(MY_PI2*(r - ters_R)/ters_D) + 
                                   cos(3*MY_PI2*(r - ters_R)/ters_D));
}

/* ---------------------------------------------------------------------- */

double PairTersoffMOD::ters_bij(double zeta, Param *param)
{
  double tmp = param->beta * zeta;
  if (tmp > param->ca1) return pow(tmp, -param->powern/(2.0*param->powern_del));
  if (tmp < param->ca4) return 1.0;
  return pow(1.0 + pow(tmp,param->powern), -1.0/(2.0*param->powern_del));
}

/* ---------------------------------------------------------------------- */

double PairTersoffMOD::ters_bij_d(double zeta, Param *param)
{
  double tmp = param->beta * zeta;
  if (tmp > param->ca1) return -0.5*(param->powern/param->powern_del)*
	  pow(tmp,-0.5*(param->powern/param->powern_del)) / zeta;
  if (tmp < param->ca4) return 0.0;
			  
  double tmp_n = pow(tmp,param->powern);
  return -0.5 *(param->powern/param->powern_del)* 
	  pow(1.0+tmp_n, -1.0-(1.0/(2.0*param->powern_del)))*tmp_n / zeta;
}

/* ---------------------------------------------------------------------- */

void PairTersoffMOD::ters_zetaterm_d(double prefactor,
				  double *rij_hat, double rij,
				  double *rik_hat, double rik,
				  double *dri, double *drj, double *drk,
				  Param *param)
{
  double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc(rik,param);
  dfc = ters_fc_d(rik,param);
  if (param->powermint == 3) tmp = pow(param->lam3 * (rij-rik),3.0);
  else tmp = param->lam3 * (rij-rik);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (param->powermint == 3)
    ex_delr_d = 3.0*pow(param->lam3,3.0) * pow(rij-rik,2.0)*ex_delr;
  else ex_delr_d = param->lam3 * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk_mod(cos_theta,param);
  gijk_d = ters_gijk_d_mod(cos_theta,param);
  costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

  // compute the derivative wrt Ri
  // dri = -dfc*gijk*ex_delr*rik_hat;
  // dri += fc*gijk_d*ex_delr*dcosdri;
  // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat);

  vec3_scale(-dfc*gijk*ex_delr,rik_hat,dri);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdri,dri,dri);
  vec3_scaleadd(fc*gijk*ex_delr_d,rik_hat,dri,dri);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rij_hat,dri,dri);
  vec3_scale(prefactor,dri,dri);

  // compute the derivative wrt Rj
  // drj = fc*gijk_d*ex_delr*dcosdrj;
  // drj += fc*gijk*ex_delr_d*rij_hat;

  vec3_scale(fc*gijk_d*ex_delr,dcosdrj,drj);
  vec3_scaleadd(fc*gijk*ex_delr_d,rij_hat,drj,drj);
  vec3_scale(prefactor,drj,drj);

  // compute the derivative wrt Rk
  // drk = dfc*gijk*ex_delr*rik_hat;
  // drk += fc*gijk_d*ex_delr*dcosdrk;
  // drk += -fc*gijk*ex_delr_d*rik_hat;

  vec3_scale(dfc*gijk*ex_delr,rik_hat,drk);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdrk,drk,drk);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rik_hat,drk,drk);
  vec3_scale(prefactor,drk,drk);
}
