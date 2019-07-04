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
   Contributing author: German Samolyuk (ORNL)
   Based on PairTersoffZBL by Aidan Thompson (SNL) and David Farrell (NWU)
------------------------------------------------------------------------- */

#include "pair_gw_zbl.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

PairGWZBL::PairGWZBL(LAMMPS *lmp) : PairGW(lmp)
{
  // hard-wired constants in metal or real units
  // a0 = Bohr radius
  // epsilon0 = permittivity of vacuum = q / energy-distance units
  // e = unit charge
  // 1 Kcal/mole = 0.043365121 eV

  if (strcmp(update->unit_style,"metal") == 0) {
    global_a_0 = 0.529;
    global_epsilon_0 = 0.00552635;
    global_e = 1.0;
  } else if (strcmp(update->unit_style,"real") == 0) {
    global_a_0 = 0.529;
    global_epsilon_0 = 0.00552635 * 0.043365121;
    global_e = 1.0;
  } else error->all(FLERR,"Pair gw/zbl requires metal or real units");
}

/* ---------------------------------------------------------------------- */

void PairGWZBL::read_file(char *file)
{
  int params_per_line = 21;
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
      snprintf(str,128,"Cannot open GW potential file %s",file);
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
      error->all(FLERR,"Incorrect format in GW potential file");

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
    params[nparams].Z_i = atof(words[17]);
    params[nparams].Z_j = atof(words[18]);
    params[nparams].ZBLcut = atof(words[19]);
    params[nparams].ZBLexpscale = atof(words[20]);

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
        (params[nparams].powermint != 3 && params[nparams].powermint != 1) ||
        params[nparams].gamma < 0.0 ||
        params[nparams].Z_i < 1.0 || params[nparams].Z_j < 1.0 ||
        params[nparams].ZBLcut < 0.0 || params[nparams].ZBLexpscale < 0.0)
      error->all(FLERR,"Illegal GW parameter");

    nparams++;
  }

  delete [] words;
}

/* ---------------------------------------------------------------------- */

void PairGWZBL::repulsive(Param *param, double rsq, double &fforce,
                               int eflag, double &eng)
{
  double r,tmp_fc,tmp_fc_d,tmp_exp;

  // GW repulsive portion

  r = sqrt(rsq);
  tmp_fc = gw_fc(r,param);
  tmp_fc_d = gw_fc_d(r,param);
  tmp_exp = exp(-param->lam1 * r);
  double fforce_gw = param->biga * tmp_exp * (tmp_fc_d - tmp_fc*param->lam1);
  double eng_gw = tmp_fc * param->biga * tmp_exp;

  // ZBL repulsive portion

  double esq = pow(global_e,2.0);
  double a_ij = (0.8854*global_a_0) /
    (pow(param->Z_i,0.23) + pow(param->Z_j,0.23));
  double premult = (param->Z_i * param->Z_j * esq)/(4.0*MY_PI*global_epsilon_0);
  double r_ov_a = r/a_ij;
  double phi = 0.1818*exp(-3.2*r_ov_a) + 0.5099*exp(-0.9423*r_ov_a) +
    0.2802*exp(-0.4029*r_ov_a) + 0.02817*exp(-0.2016*r_ov_a);
  double dphi = (1.0/a_ij) * (-3.2*0.1818*exp(-3.2*r_ov_a) -
                              0.9423*0.5099*exp(-0.9423*r_ov_a) -
                              0.4029*0.2802*exp(-0.4029*r_ov_a) -
                              0.2016*0.02817*exp(-0.2016*r_ov_a));
  double fforce_ZBL = premult*-phi/rsq + premult*dphi/r;
  double eng_ZBL = premult*(1.0/r)*phi;

  // combine two parts with smoothing by Fermi-like function

  fforce = -(-F_fermi_d(r,param) * eng_ZBL +
             (1.0 - F_fermi(r,param))*fforce_ZBL +
             F_fermi_d(r,param)*eng_gw + F_fermi(r,param)*fforce_gw) / r;

  if (eflag)
    eng = (1.0 - F_fermi(r,param))*eng_ZBL + F_fermi(r,param)*eng_gw;
}

/* ---------------------------------------------------------------------- */

double PairGWZBL::gw_fa(double r, Param *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return -param->bigb * exp(-param->lam2 * r) * gw_fc(r,param) *
    F_fermi(r,param);
}

/* ---------------------------------------------------------------------- */

double PairGWZBL::gw_fa_d(double r, Param *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return param->bigb * exp(-param->lam2 * r) *
    (param->lam2 * gw_fc(r,param) * F_fermi(r,param) -
     gw_fc_d(r,param) * F_fermi(r,param) - gw_fc(r,param) *
     F_fermi_d(r,param));
}

/* ----------------------------------------------------------------------
   Fermi-like smoothing function
------------------------------------------------------------------------- */

double PairGWZBL::F_fermi(double r, Param *param)
{
  return 1.0 / (1.0 + exp(-param->ZBLexpscale*(r-param->ZBLcut)));
}

/* ----------------------------------------------------------------------
   Fermi-like smoothing function derivative with respect to r
------------------------------------------------------------------------- */

double PairGWZBL::F_fermi_d(double r, Param *param)
{
  return param->ZBLexpscale*exp(-param->ZBLexpscale*(r-param->ZBLcut)) /
    pow(1.0 + exp(-param->ZBLexpscale*(r-param->ZBLcut)),2.0);
}
