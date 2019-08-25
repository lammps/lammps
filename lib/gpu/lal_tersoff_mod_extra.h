/// **************************************************************************
//                             tersoff_mod_extra.h
//                             -------------------
//                              Trung Dac Nguyen
//
//  Device code for Tersoff math routines
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : ndactrung@gmail.com
// ***************************************************************************/*

#ifndef LAL_TERSOFF_MOD_EXTRA_H
#define LAL_TERSOFF_MOD_EXTRA_H

#ifdef NV_KERNEL
#include "lal_aux_fun1.h"
#else
#endif

#define MY_PI2 (numtyp)1.57079632679489661923
#define MY_PI4 (numtyp)0.78539816339744830962

/* ---------------------------------------------------------------------- */

ucl_inline numtyp vec3_dot(const numtyp x[3], const numtyp y[3])
{
  return (x[0]*y[0] + x[1]*y[1] + x[2]*y[2]);
}

ucl_inline void vec3_add(const numtyp x[3], const numtyp y[3], numtyp z[3])
{
  z[0] = x[0]+y[0]; z[1] = x[1]+y[1]; z[2] = x[2]+y[2];
}

ucl_inline void vec3_scale(const numtyp k, const numtyp x[3], numtyp y[3])
{
  y[0] = k*x[0]; y[1] = k*x[1]; y[2] = k*x[2];
}

ucl_inline void vec3_scaleadd(const numtyp k, const numtyp x[3],
                              const numtyp y[3], numtyp z[3])
{
  z[0] = k*x[0]+y[0]; z[1] = k*x[1]+y[1]; z[2] = k*x[2]+y[2];
}

/* ---------------------------------------------------------------------- */

ucl_inline numtyp ters_gijk_mod(const numtyp costheta,
                                const numtyp param_c1,
                                const numtyp param_c2,
                                const numtyp param_c3,
                                const numtyp param_c4,
                                const numtyp param_c5,
                                const numtyp param_h)
{
  const numtyp tmp_h = (param_h - costheta)*(param_h - costheta);
  return param_c1 + (param_c2*tmp_h/(param_c3 + tmp_h)) *
    ((numtyp)1.0 + param_c4*ucl_exp(-param_c5*tmp_h));
}

/* ---------------------------------------------------------------------- */

ucl_inline numtyp ters_gijk_d_mod(const numtyp costheta,
                                  const numtyp param_c2,
                                  const numtyp param_c3,
                                  const numtyp param_c4,
                                  const numtyp param_c5,
                                  const numtyp param_h)
{
  const numtyp tmp_h = (param_h - costheta)*(param_h - costheta);
  const numtyp g1 = (param_h - costheta)/(param_c3 + tmp_h);
  const numtyp g2 = ucl_exp(-param_c5*tmp_h);
  return (numtyp)-2.0*param_c2*g1*((1 + param_c4*g2) *
         (1 + g1*(costheta - param_h)) - tmp_h*param_c4*param_c5*g2);
}

/* ---------------------------------------------------------------------- */

ucl_inline void costheta_d(const numtyp rij_hat[3],
                           const numtyp rij,
                           const numtyp rik_hat[3],
                           const numtyp rik,
                           numtyp *dri,
                           numtyp *drj,
                           numtyp *drk)
{
  // first element is derivative wrt Ri, second wrt Rj, third wrt Rk

  numtyp cos_theta = vec3_dot(rij_hat,rik_hat);

  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,drj);
  vec3_scale(ucl_recip(rij),drj,drj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,drk);
  vec3_scale(ucl_recip(rik),drk,drk);
  vec3_add(drj,drk,dri);
  vec3_scale((numtyp)-1.0,dri,dri);
}

/* ---------------------------------------------------------------------- */

ucl_inline numtyp ters_fc(const numtyp r,
                          const numtyp param_bigr,
                          const numtyp param_bigd)
{
  if (r < param_bigr-param_bigd) return (numtyp)1.0;
  if (r > param_bigr+param_bigd) return (numtyp)0.0;
  return (numtyp)0.5*((numtyp)1.0 -
         (numtyp)1.125*sin(MY_PI2*(r - param_bigr)/param_bigd) -
         (numtyp)0.125*sin(3*MY_PI2*(r - param_bigr)/param_bigd));
}

/* ---------------------------------------------------------------------- */

ucl_inline numtyp ters_fc_d(const numtyp r,
                            const numtyp param_bigr,
                            const numtyp param_bigd)
{
  if (r < param_bigr-param_bigd) return (numtyp)0.0;
  if (r > param_bigr+param_bigd) return (numtyp)0.0;
  return -((numtyp)0.375*MY_PI4/param_bigd) *
           ((numtyp)3*cos(MY_PI2*(r - param_bigr)/param_bigd) +
            cos((numtyp)3*MY_PI2*(r - param_bigr)/param_bigd));
}

/* ---------------------------------------------------------------------- */

ucl_inline numtyp ters_fa(const numtyp r,
                          const numtyp param_bigb,
                          const numtyp param_bigr,
                          const numtyp param_bigd,
                          const numtyp param_lam2)
{
  if (r > param_bigr + param_bigd) return (numtyp)0.0;
  return -param_bigb * ucl_exp(-param_lam2 * r) *
    ters_fc(r,param_bigr,param_bigd);
}

/* ---------------------------------------------------------------------- */

ucl_inline numtyp ters_fa_d(const numtyp r,
                            const numtyp param_bigb,
                            const numtyp param_bigr,
                            const numtyp param_bigd,
                            const numtyp param_lam2)
{
  if (r > param_bigr + param_bigd) return (numtyp)0.0;
  return param_bigb * ucl_exp(-param_lam2 * r) * (param_lam2 *
    ters_fc(r,param_bigr,param_bigd) - ters_fc_d(r,param_bigr,param_bigd));
}

/* ---------------------------------------------------------------------- */

ucl_inline numtyp ters_bij(const numtyp zeta,
                           const numtyp param_beta,
                           const numtyp param_powern,
                           const numtyp param_powern_del,
                           const numtyp param_ca1,
                           const numtyp param_ca4)
{
  numtyp tmp = param_beta * zeta;
  if (tmp > param_ca1)
    return ucl_powr(tmp, -param_powern/((numtyp)2.0*param_powern_del));
  if (tmp < param_ca4) return (numtyp)1.0;
  return ucl_powr((numtyp)1.0 + ucl_powr(tmp,param_powern),
    (numtyp)-1.0/((numtyp)2.0*param_powern_del));
}

/* ---------------------------------------------------------------------- */

ucl_inline numtyp ters_bij_d(const numtyp zeta,
                             const numtyp param_beta,
                             const numtyp param_powern,
                             const numtyp param_powern_del,
                             const numtyp param_ca1,
                             const numtyp param_ca4)
{
  numtyp tmp = param_beta * zeta;
  if (tmp > param_ca1) return (numtyp)-0.5*(param_powern/param_powern_del) *
          ucl_powr(tmp,(numtyp)-0.5*(param_powern/param_powern_del)) / zeta;
  if (tmp < param_ca4) return (numtyp)0.0;

  numtyp tmp_n = ucl_powr(tmp,param_powern);
  return (numtyp)-0.5 *(param_powern/param_powern_del) *
          ucl_powr((numtyp)1.0+tmp_n, (numtyp)-1.0-((numtyp)1.0 /
    ((numtyp)2.0*param_powern_del)))*tmp_n / zeta;
}

/* ---------------------------------------------------------------------- */

ucl_inline void ters_zetaterm_d(const numtyp prefactor,
                                const numtyp rij_hat[3],
                                const numtyp rij,
                                const numtyp rik_hat[3],
                                const numtyp rik,
                                const numtyp param_bigr,
                                const numtyp param_bigd,
                                const numtyp param_powermint,
                                const numtyp param_lam3,
                                const numtyp param_h,
                                const numtyp param_c1,
                                const numtyp param_c2,
                                const numtyp param_c3,
                                const numtyp param_c4,
                                const numtyp param_c5,
                                numtyp dri[3],
                                numtyp drj[3],
                                numtyp drk[3])
{
  numtyp gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  numtyp dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc(rik,param_bigr,param_bigd);
  dfc = ters_fc_d(rik,param_bigr,param_bigd);

  numtyp t = param_lam3*(rij-rik);
  if ((int)param_powermint == 3) tmp = t*t*t;
  else tmp = t;

  if (tmp > (numtyp)69.0776) ex_delr = (numtyp)1.e30;
  else if (tmp < (numtyp)-69.0776) ex_delr = (numtyp)0.0;
  else ex_delr = ucl_exp(tmp);

  if ((int)param_powermint == 3)
    ex_delr_d = (numtyp)3.0*param_lam3*t*t*ex_delr;
  else ex_delr_d = param_lam3 * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk_mod(cos_theta,param_c1,param_c2,param_c3,param_c4,param_c5,param_h);
  gijk_d = ters_gijk_d_mod(cos_theta,param_c2,param_c3,param_c4,param_c5,param_h);
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

ucl_inline void ters_zetaterm_d_fi(const numtyp prefactor,
                                   const numtyp rij_hat[3],
                                   const numtyp rij,
                                   const numtyp rik_hat[3],
                                   const numtyp rik,
                                   const numtyp param_bigr,
                                   const numtyp param_bigd,
                                   const numtyp param_powermint,
                                   const numtyp param_lam3,
                                   const numtyp param_h,
                                   const numtyp param_c1,
                                   const numtyp param_c2,
                                   const numtyp param_c3,
                                   const numtyp param_c4,
                                   const numtyp param_c5,
                                   numtyp dri[3])
{
  numtyp gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  numtyp dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc(rik,param_bigr,param_bigd);
  dfc = ters_fc_d(rik,param_bigr,param_bigd);

  numtyp t = param_lam3*(rij-rik);
  if ((int)param_powermint == 3) tmp = t*t*t;
  else tmp = t;

  if (tmp > (numtyp)69.0776) ex_delr = (numtyp)1.e30;
  else if (tmp < (numtyp)-69.0776) ex_delr = (numtyp)0.0;
  else ex_delr = ucl_exp(tmp);

  if ((int)param_powermint == 3)
    ex_delr_d = (numtyp)3.0*param_lam3*t*t*ex_delr;
  else ex_delr_d = param_lam3 * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk_mod(cos_theta,param_c1,param_c2,param_c3,param_c4,param_c5,param_h);
  gijk_d = ters_gijk_d_mod(cos_theta,param_c2,param_c3,param_c4,param_c5,param_h);
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
}

ucl_inline void ters_zetaterm_d_fj(const numtyp prefactor,
                                   const numtyp rij_hat[3],
                                   const numtyp rij,
                                   const numtyp rik_hat[3],
                                   const numtyp rik,
                                   const numtyp param_bigr,
                                   const numtyp param_bigd,
                                   const numtyp param_powermint,
                                   const numtyp param_lam3,
                                   const numtyp param_h,
                                   const numtyp param_c1,
                                   const numtyp param_c2,
                                   const numtyp param_c3,
                                   const numtyp param_c4,
                                   const numtyp param_c5,
                                   numtyp drj[3])
{
  numtyp gijk,gijk_d,ex_delr,ex_delr_d,fc,cos_theta,tmp;
  numtyp dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc(rik,param_bigr,param_bigd);

  numtyp t = param_lam3*(rij-rik);
  if ((int)param_powermint == 3) tmp = t*t*t;
  else tmp = t;

  if (tmp > (numtyp)69.0776) ex_delr = (numtyp)1.e30;
  else if (tmp < (numtyp)-69.0776) ex_delr = (numtyp)0.0;
  else ex_delr = ucl_exp(tmp);

  if ((int)param_powermint == 3)
    ex_delr_d = (numtyp)3.0*param_lam3*t*t*ex_delr;
  else ex_delr_d = param_lam3 * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk_mod(cos_theta,param_c1,param_c2,param_c3,param_c4,param_c5,param_h);
  gijk_d = ters_gijk_d_mod(cos_theta,param_c2,param_c3,param_c4,param_c5,param_h);
  costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

  // compute the derivative wrt Rj
  // drj = fc*gijk_d*ex_delr*dcosdrj;
  // drj += fc*gijk*ex_delr_d*rij_hat;

  vec3_scale(fc*gijk_d*ex_delr,dcosdrj,drj);
  vec3_scaleadd(fc*gijk*ex_delr_d,rij_hat,drj,drj);
  vec3_scale(prefactor,drj,drj);
}

ucl_inline void ters_zetaterm_d_fk(const numtyp prefactor,
                                   const numtyp rij_hat[3],
                                   const numtyp rij,
                                   const numtyp rik_hat[3],
                                   const numtyp rik,
                                   const numtyp param_bigr,
                                   const numtyp param_bigd,
                                   const numtyp param_powermint,
                                   const numtyp param_lam3,
                                   const numtyp param_h,
                                   const numtyp param_c1,
                                   const numtyp param_c2,
                                   const numtyp param_c3,
                                   const numtyp param_c4,
                                   const numtyp param_c5,
                                   numtyp drk[3])
{
  numtyp gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  numtyp dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc(rik,param_bigr,param_bigd);
  dfc = ters_fc_d(rik,param_bigr,param_bigd);

  numtyp t = param_lam3*(rij-rik);
  if ((int)param_powermint == 3) tmp = t*t*t;
  else tmp = t;

  if (tmp > (numtyp)69.0776) ex_delr = (numtyp)1.e30;
  else if (tmp < (numtyp)-69.0776) ex_delr = (numtyp)0.0;
  else ex_delr = ucl_exp(tmp);

  if ((int)param_powermint == 3)
    ex_delr_d = (numtyp)3.0*param_lam3*t*t*ex_delr;
  else ex_delr_d = param_lam3 * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk_mod(cos_theta,param_c1,param_c2,param_c3,param_c4,param_c5,param_h);
  gijk_d = ters_gijk_d_mod(cos_theta,param_c2,param_c3,param_c4,param_c5,param_h);
  costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

  // compute the derivative wrt Rk
  // drk = dfc*gijk*ex_delr*rik_hat;
  // drk += fc*gijk_d*ex_delr*dcosdrk;
  // drk += -fc*gijk*ex_delr_d*rik_hat;

  vec3_scale(dfc*gijk*ex_delr,rik_hat,drk);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdrk,drk,drk);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rik_hat,drk,drk);
  vec3_scale(prefactor,drk,drk);
}

/* ---------------------------------------------------------------------- */

ucl_inline void repulsive(const numtyp param_bigr,
                          const numtyp param_bigd,
                          const numtyp param_lam1,
                          const numtyp param_biga,
                          const numtyp rsq,
                          const int eflag,
                          numtyp *ans)
{
  numtyp r,tmp_fc,tmp_fc_d,tmp_exp;
  r = ucl_sqrt(rsq);
  tmp_fc = ters_fc(r,param_bigr,param_bigd);
  tmp_fc_d = ters_fc_d(r,param_bigr,param_bigd);
  tmp_exp = ucl_exp(-param_lam1 * r);
  // fforce
  ans[0] = -param_biga*tmp_exp*(tmp_fc_d - tmp_fc*param_lam1)*ucl_recip(r);
  // eng
  if (eflag) ans[1] = tmp_fc * param_biga * tmp_exp;
}

/* ---------------------------------------------------------------------- */

ucl_inline numtyp zeta(const numtyp param_powermint,
                       const numtyp param_lam3,
                       const numtyp param_bigr,
                       const numtyp param_bigd,
                       const numtyp param_h,
                       const numtyp param_c1,
                       const numtyp param_c2,
                       const numtyp param_c3,
                       const numtyp param_c4,
                       const numtyp param_c5,
                       const numtyp rsqij,
                       const numtyp rsqik,
                       const numtyp4 delrij,
                       const numtyp4 delrik)
{
  numtyp rij,rik,costheta,arg,ex_delr;

  rij = ucl_sqrt(rsqij);
  rik = ucl_sqrt(rsqik);
  costheta = (delrij.x*delrik.x + delrij.y*delrik.y +
              delrij.z*delrik.z) / (rij*rik);

  numtyp t = param_lam3*(rij-rik);
  if ((int)param_powermint == 3) arg = t*t*t;
  else arg = t;

  if (arg > (numtyp)69.0776) ex_delr = (numtyp)1.e30;
  else if (arg < (numtyp)-69.0776) ex_delr = (numtyp)0.0;
  else ex_delr = ucl_exp(arg);

  return ters_fc(rik,param_bigr,param_bigd) *
    ters_gijk_mod(costheta,param_c1,param_c2,param_c3,param_c4,param_c5,
                  param_h) * ex_delr;
}

/* ---------------------------------------------------------------------- */

ucl_inline void force_zeta(const numtyp param_bigb,
                           const numtyp param_bigr,
                           const numtyp param_bigd,
                           const numtyp param_lam2,
                           const numtyp param_beta,
                           const numtyp param_powern,
                           const numtyp param_powern_del,
                           const numtyp param_ca1,
                           const numtyp param_ca4,
                           const numtyp rsq,
                           const numtyp zeta_ij,
                           const int eflag,
                           numtyp fpfeng[4])
{
  numtyp r,fa,fa_d,bij;

  r = ucl_sqrt(rsq);
  fa = ters_fa(r,param_bigb,param_bigr,param_bigd,param_lam2);
  fa_d = ters_fa_d(r,param_bigb,param_bigr,param_bigd,param_lam2);
  bij = ters_bij(zeta_ij,param_beta,param_powern,
                 param_powern_del,param_ca1,param_ca4);
  fpfeng[0] = (numtyp)0.5*bij*fa_d * ucl_recip(r); // fforce
  fpfeng[1] = (numtyp)-0.5*fa * ters_bij_d(zeta_ij,param_beta, param_powern,
           param_powern_del,param_ca1,param_ca4); // prefactor
  if (eflag) fpfeng[2] = (numtyp)0.5*bij*fa; // eng
}

/* ----------------------------------------------------------------------
   attractive term
   use param_ij cutoff for rij test
   use param_ijk cutoff for rik test
------------------------------------------------------------------------- */

ucl_inline void attractive(const numtyp param_bigr,
                           const numtyp param_bigd,
                           const numtyp param_powermint,
                           const numtyp param_lam3,
                           const numtyp param_h,
                           const numtyp param_c1,
                           const numtyp param_c2,
                           const numtyp param_c3,
                           const numtyp param_c4,
                           const numtyp param_c5,
                           const numtyp prefactor,
                           const numtyp rij,
                           const numtyp rijinv,
                           const numtyp rik,
                           const numtyp rikinv,
                           const numtyp delrij[3],
                           const numtyp delrik[3],
                           numtyp fi[3],
                           numtyp fj[3],
                           numtyp fk[3])
{
  numtyp rij_hat[3],rik_hat[3];
  vec3_scale(rijinv,delrij,rij_hat);
  vec3_scale(rikinv,delrik,rik_hat);
  ters_zetaterm_d(prefactor,rij_hat,rij,rik_hat,rik,
                  param_bigr, param_bigd, param_powermint, param_lam3,
                  param_h, param_c1, param_c2, param_c3, param_c4, param_c5,
                  fi, fj, fk);
}

ucl_inline void attractive_fi(const numtyp param_bigr,
                              const numtyp param_bigd,
                              const numtyp param_powermint,
                              const numtyp param_lam3,
                              const numtyp param_h,
                              const numtyp param_c1,
                              const numtyp param_c2,
                              const numtyp param_c3,
                              const numtyp param_c4,
                              const numtyp param_c5,
                              const numtyp prefactor,
                              const numtyp rij,
                              const numtyp rijinv,
                              const numtyp rik,
                              const numtyp rikinv,
                              const numtyp delrij[3],
                              const numtyp delrik[3],
                              numtyp fi[3])
{
  numtyp rij_hat[3],rik_hat[3];
  vec3_scale(rijinv,delrij,rij_hat);
  vec3_scale(rikinv,delrik,rik_hat);
  ters_zetaterm_d_fi(prefactor,rij_hat,rij,rik_hat,rik,
                  param_bigr, param_bigd, param_powermint, param_lam3,
                  param_h, param_c1, param_c2, param_c3, param_c4, param_c5,
                  fi);
}

ucl_inline void attractive_fj(const numtyp param_bigr,
                              const numtyp param_bigd,
                              const numtyp param_powermint,
                              const numtyp param_lam3,
                              const numtyp param_h,
                              const numtyp param_c1,
                              const numtyp param_c2,
                              const numtyp param_c3,
                              const numtyp param_c4,
                              const numtyp param_c5,
                              const numtyp prefactor,
                              const numtyp rij,
                              const numtyp rijinv,
                              const numtyp rik,
                              const numtyp rikinv,
                              const numtyp delrij[3],
                              const numtyp delrik[3],
                              numtyp fj[3])
{
  numtyp rij_hat[3],rik_hat[3];
  vec3_scale(rijinv,delrij,rij_hat);
  vec3_scale(rikinv,delrik,rik_hat);
  ters_zetaterm_d_fj(prefactor,rij_hat,rij,rik_hat,rik,
                     param_bigr, param_bigd, param_powermint, param_lam3,
                     param_h, param_c1, param_c2, param_c3, param_c4, param_c5,
                     fj);
}

ucl_inline void attractive_fk(const numtyp param_bigr,
                              const numtyp param_bigd,
                              const numtyp param_powermint,
                              const numtyp param_lam3,
                              const numtyp param_h,
                              const numtyp param_c1,
                              const numtyp param_c2,
                              const numtyp param_c3,
                              const numtyp param_c4,
                              const numtyp param_c5,
                              const numtyp prefactor,
                              const numtyp rij,
                              const numtyp rijinv,
                              const numtyp rik,
                              const numtyp rikinv,
                              const numtyp delrij[3],
                              const numtyp delrik[3],
                              numtyp fk[3])
{
  numtyp rij_hat[3],rik_hat[3];
  vec3_scale(rijinv,delrij,rij_hat);
  vec3_scale(rikinv,delrik,rik_hat);
  ters_zetaterm_d_fk(prefactor,rij_hat,rij,rik_hat,rik,
                     param_bigr, param_bigd, param_powermint, param_lam3,
                     param_h, param_c1, param_c2, param_c3, param_c4, param_c5,
                     fk);
}


#endif


