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
   Contributing author: Aidan Thompson (SNL) - original Tersoff implementation
                        Vitaly Dozhdikov (JIHT of RAS) - MOD addition
------------------------------------------------------------------------- */

#include "pair_tersoff_mod.h"

#include "comm.h"
#include "error.h"
#include "math_const.h"
#include "math_extra.h"
#include "math_special.h"
#include "memory.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathExtra;
using namespace MathSpecial;

#define DELTA 4

/* ---------------------------------------------------------------------- */

PairTersoffMOD::PairTersoffMOD(LAMMPS *lmp) : PairTersoff(lmp) {}

/* ---------------------------------------------------------------------- */

void PairTersoffMOD::read_file(char *file)
{
  memory->sfree(params);
  params = nullptr;
  nparams = maxparam = 0;

  // open file on proc 0

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, file, "tersoff/mod", unit_convert_flag);
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
        params[nparams].powerm     = values.next_double();
        params[nparams].lam3       = values.next_double();
        params[nparams].h          = values.next_double();
        params[nparams].powern     = values.next_double();
        params[nparams].beta       = values.next_double();
        params[nparams].lam2       = values.next_double();
        params[nparams].bigb       = values.next_double();
        params[nparams].bigr       = values.next_double();
        params[nparams].bigd       = values.next_double();
        params[nparams].lam1       = values.next_double();
        params[nparams].biga       = values.next_double();
        params[nparams].powern_del = values.next_double();
        params[nparams].c1         = values.next_double();
        params[nparams].c2         = values.next_double();
        params[nparams].c3         = values.next_double();
        params[nparams].c4         = values.next_double();
        params[nparams].c5         = values.next_double();
        params[nparams].powermint = int(params[nparams].powerm);

        if (unit_convert) {
          params[nparams].biga *= conversion_factor;
          params[nparams].bigb *= conversion_factor;
        }
      } catch (TokenizerException &e) {
        error->one(FLERR, e.what());
      }

      // currently only allow m exponent of 1 or 3
      if (params[nparams].powern < 0.0 ||
          params[nparams].beta < 0.0 ||
          params[nparams].lam2 < 0.0 ||
          params[nparams].bigb < 0.0 ||
          params[nparams].bigr < 0.0 ||
          params[nparams].bigd < 0.0 ||
          params[nparams].bigd > params[nparams].bigr ||
          params[nparams].lam1 < 0.0 ||
          params[nparams].biga < 0.0 ||
          params[nparams].powerm - params[nparams].powermint != 0.0 ||
          (params[nparams].powermint != 3 &&
          params[nparams].powermint != 1)
          )
        error->one(FLERR,"Illegal Tersoff parameter");

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

void PairTersoffMOD::setup_params()
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
            if (n >= 0) error->all(FLERR,"Potential file has a duplicate entry for: {} {} {}",
                                   elements[i], elements[j], elements[k]);
            n = m;
          }
        }
        if (n < 0) error->all(FLERR,"Potential file is missing an entry for: {} {} {}",
                              elements[i], elements[j], elements[k]);
        elem3param[i][j][k] = n;
      }


  // compute parameter values derived from inputs

  for (m = 0; m < nparams; m++) {
    params[m].cut = params[m].bigr + params[m].bigd;
    params[m].cutsq = params[m].cut*params[m].cut;

    if (params[m].powern > 0.0) {
      params[m].ca1 = pow(2.0*params[m].powern_del*1.0e-16,-1.0/params[m].powern);
      params[m].ca4 = 1.0/params[m].ca1;
    } else params[m].ca1 = params[m].ca4 = 0.0;
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++)
    if (params[m].cut > cutmax) cutmax = params[m].cut;
}

/* ---------------------------------------------------------------------- */

double PairTersoffMOD::zeta(Param *param, double rsqij, double rsqik,
                         double *rij_hat, double *rik_hat)
{
  double rij,rik,costheta,arg,ex_delr;

  rij = sqrt(rsqij);
  rik = sqrt(rsqik);
  costheta = dot3(rij_hat,rik_hat);

  if (param->powermint == 3) arg = cube(param->lam3 * (rij-rik));
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
                                  double *rij_hat, double rij, double rijinv,
                                  double *rik_hat, double rik, double rikinv,
                                  double *dri, double *drj, double *drk,
                                  Param *param)
{
  double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc(rik,param);
  dfc = ters_fc_d(rik,param);
  if (param->powermint == 3) tmp = cube(param->lam3 * (rij-rik));
  else tmp = param->lam3 * (rij-rik);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (param->powermint == 3)
    ex_delr_d = 3.0*cube(param->lam3) * square(rij-rik)*ex_delr;
  else ex_delr_d = param->lam3 * ex_delr;

  cos_theta = dot3(rij_hat,rik_hat);
  gijk = ters_gijk_mod(cos_theta,param);
  gijk_d = ters_gijk_d_mod(cos_theta,param);
  costheta_d(rij_hat,rijinv,rik_hat,rikinv,dcosdri,dcosdrj,dcosdrk);

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
