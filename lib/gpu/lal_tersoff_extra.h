/// **************************************************************************
//                              tersoff_extra.h
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

#ifndef LAL_TERSOFF_EXTRA_H
#define LAL_TERSOFF_EXTRA_H

#if defined(NV_KERNEL) || defined(USE_HIP)
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

ucl_inline numtyp ters_gijk(const numtyp costheta,
                            const numtyp param_c,
                            const numtyp param_d,
                            const numtyp param_h,
                            const numtyp param_gamma)
{
  const numtyp hcth = param_h - costheta;
  return param_gamma*((numtyp)1.0 + param_c*ucl_recip(param_d) -
         param_c *ucl_recip(param_d + hcth*hcth));
}

/* ---------------------------------------------------------------------- */

ucl_inline numtyp ters_gijk_d(const numtyp costheta,
                              const numtyp param_c,
                              const numtyp param_d,
                              const numtyp param_h,
                              const numtyp param_gamma,
                              numtyp *ans_d)
{
  const numtyp hcth = param_h - costheta;
  const numtyp idhh=ucl_recip(param_d + hcth*hcth);
  const numtyp numerator = (numtyp)-2.0 * param_c * hcth;
  *ans_d=param_gamma*numerator*idhh*idhh;
  return param_gamma*((numtyp)1.0+param_c*ucl_recip(param_d)-param_c*idhh);
}

/* ---------------------------------------------------------------------- */

ucl_inline void costheta_d(const numtyp cos_theta,
                           const numtyp rij_hat[3],
                           const numtyp rij,
                           const numtyp rik_hat[3],
                           const numtyp rik,
                           numtyp *dri,
                           numtyp *drj,
                           numtyp *drk)
{
  // first element is derivative wrt Ri, second wrt Rj, third wrt Rk
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
  #ifndef ONETYPE
  if (r > param_bigr+param_bigd) return (numtyp)0.0;
  #endif
  return (numtyp)0.5*((numtyp)1.0 - sin(MY_PI2*(r - param_bigr)/param_bigd));
}

/* ---------------------------------------------------------------------- */

ucl_inline numtyp ters_fc_d(const numtyp r,
                            const numtyp param_bigr,
                            const numtyp param_bigd,
                            numtyp *ans_d)
{
  if (r < param_bigr-param_bigd) {
    *ans_d=(numtyp)0.0;
    return (numtyp)1.0;
  }
  #ifndef ONETYPE
  if (r > param_bigr+param_bigd) {
    *ans_d=(numtyp)0.0;
    return (numtyp)0.0;
  }
  #endif
  const numtyp ibigd = ucl_recip(param_bigd);
  const numtyp angle = MY_PI2*(r - param_bigr)*ibigd;
  *ans_d=-(MY_PI4*ibigd) * cos(angle);
  return (numtyp)0.5*((numtyp)1.0 - sin(angle));
}

/* ---------------------------------------------------------------------- */

ucl_inline numtyp ters_fa_d(const numtyp r,
                            const numtyp param_bigb,
                            const numtyp param_bigr,
                            const numtyp param_bigd,
                            const numtyp param_lam2,
                            numtyp *ans_d)
{
  #ifndef ONETYPE
  if (r > param_bigr + param_bigd) {
    *ans_d = (numtyp)0.0;
    return (numtyp)0.0;
  }
  #endif
  numtyp dfc;
  const numtyp fc=ters_fc_d(r,param_bigr,param_bigd,&dfc);
  const numtyp blr = param_bigb * ucl_exp(-param_lam2 * r);
  *ans_d = blr * (param_lam2 * fc - dfc);
  return -blr * fc;
}

/* ---------------------------------------------------------------------- */

ucl_inline numtyp ters_bij_d(const numtyp zeta,
                             const numtyp param_beta,
                             const numtyp param_powern,
                             const numtyp param_c1,
                             const numtyp param_c2,
                             const numtyp param_c3,
                             const numtyp param_c4,
                             numtyp *ans_d)
{
  const numtyp tmp = param_beta * zeta;
  if (tmp > param_c1) {
    *ans_d = param_beta * (numtyp)-0.5*ucl_powr(tmp,(numtyp)-1.5);
    return ucl_rsqrt(tmp);
  }
  if (tmp > param_c2) {
    const numtyp ptmp = ucl_powr(tmp,-param_powern);
    const numtyp i2n = ucl_recip((numtyp)2.0 * param_powern);
    *ans_d = param_beta * ((numtyp)-0.5*ucl_powr(tmp,(numtyp)-1.5) *
                           ((numtyp)1.0 - ((numtyp)1.0 + (numtyp)1.0 * i2n) *
                            ptmp));
    return ((numtyp)1.0 - ptmp * i2n)*ucl_rsqrt(tmp);
  }
  if (tmp < param_c4) {
    *ans_d = (numtyp)0.0;
    return (numtyp)1.0;
  }
  if (tmp < param_c3) {
    *ans_d = (numtyp)-0.5*param_beta * ucl_powr(tmp,param_powern-(numtyp)1.0);
    return (numtyp)1.0 - ucl_powr(tmp,param_powern)/((numtyp)2.0*param_powern);
  }
  const numtyp tmp_n = (numtyp)1.0+ucl_powr(tmp,param_powern);
  const numtyp i2n = -ucl_recip((numtyp)2.0*param_powern);
  *ans_d = (numtyp)-0.5*ucl_powr(tmp_n,(numtyp)-1.0+i2n)*(tmp_n-(numtyp)1.0)/
    zeta;
  return ucl_powr(tmp_n, i2n);
}

/* ---------------------------------------------------------------------- */

ucl_inline void ters_zetaterm_d(const numtyp prefactor,
                                const numtyp rij_hat[3],
                                const numtyp rij,
                                const numtyp rik_hat[3],
                                const numtyp rik,
                                const numtyp param_bigr,
                                const numtyp param_bigd,
                                const int param_powermint,
                                const numtyp param_lam3,
                                const numtyp param_c,
                                const numtyp param_d,
                                const numtyp param_h,
                                const numtyp param_gamma,
                                numtyp dri[3],
                                numtyp drj[3],
                                numtyp drk[3])
{
  numtyp gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  numtyp dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc_d(rik,param_bigr,param_bigd,&dfc);

  numtyp t = param_lam3*(rij-rik);
  if (param_powermint == 3) tmp = t*t*t;
  else tmp = t;

  if (tmp > (numtyp)69.0776) ex_delr = (numtyp)1.e30;
  else if (tmp < (numtyp)-69.0776) ex_delr = (numtyp)0.0;
  else ex_delr = ucl_exp(tmp);

  if (param_powermint == 3)
    ex_delr_d = (numtyp)3.0*param_lam3*t*t*ex_delr;
  else ex_delr_d = param_lam3 * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk_d(cos_theta,param_c,param_d,param_h,param_gamma,&gijk_d);
  costheta_d(cos_theta,rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

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
                                   const int param_powermint,
                                   const numtyp param_lam3,
                                   const numtyp param_c,
                                   const numtyp param_d,
                                   const numtyp param_h,
                                   const numtyp param_gamma,
                                   numtyp dri[3])
{
  numtyp gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  numtyp dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc_d(rik,param_bigr,param_bigd,&dfc);

  numtyp t = param_lam3*(rij-rik);
  if (param_powermint == 3) tmp = t*t*t;
  else tmp = t;

  if (tmp > (numtyp)69.0776) ex_delr = (numtyp)1.e30;
  else if (tmp < (numtyp)-69.0776) ex_delr = (numtyp)0.0;
  else ex_delr = ucl_exp(tmp);

  if (param_powermint == 3)
    ex_delr_d = (numtyp)3.0*param_lam3*t*t*ex_delr;
  else ex_delr_d = param_lam3 * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk_d(cos_theta,param_c,param_d,param_h,param_gamma,&gijk_d);
  costheta_d(cos_theta,rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

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
                                   const int param_powermint,
                                   const numtyp param_lam3,
                                   const numtyp param_c,
                                   const numtyp param_d,
                                   const numtyp param_h,
                                   const numtyp param_gamma,
                                   numtyp drj[3])
{
  numtyp gijk,gijk_d,ex_delr,ex_delr_d,fc,cos_theta,tmp;
  numtyp dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc(rik,param_bigr,param_bigd);

  numtyp t = param_lam3*(rij-rik);
  if (param_powermint == 3) tmp = t*t*t;
  else tmp = t;

  if (tmp > (numtyp)69.0776) ex_delr = (numtyp)1.e30;
  else if (tmp < (numtyp)-69.0776) ex_delr = (numtyp)0.0;
  else ex_delr = ucl_exp(tmp);

  if (param_powermint == 3)
    ex_delr_d = (numtyp)3.0*param_lam3*t*t*ex_delr;
  else ex_delr_d = param_lam3 * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk_d(cos_theta,param_c,param_d,param_h,param_gamma,&gijk_d);
  costheta_d(cos_theta,rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

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
                                   const int param_powermint,
                                   const numtyp param_lam3,
                                   const numtyp param_c,
                                   const numtyp param_d,
                                   const numtyp param_h,
                                   const numtyp param_gamma,
                                   numtyp drk[3])
{
  numtyp gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  numtyp dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc_d(rik,param_bigr,param_bigd,&dfc);

  numtyp t = param_lam3*(rij-rik);
  if (param_powermint == 3) tmp = t*t*t;
  else tmp = t;

  if (tmp > (numtyp)69.0776) ex_delr = (numtyp)1.e30;
  else if (tmp < (numtyp)-69.0776) ex_delr = (numtyp)0.0;
  else ex_delr = ucl_exp(tmp);

  if (param_powermint == 3)
    ex_delr_d = (numtyp)3.0*param_lam3*t*t*ex_delr;
  else ex_delr_d = param_lam3 * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk_d(cos_theta,param_c,param_d,param_h,param_gamma,&gijk_d);
  costheta_d(cos_theta,rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

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
  tmp_fc = ters_fc_d(r,param_bigr,param_bigd,&tmp_fc_d);
  tmp_exp = param_biga * ucl_exp(-param_lam1 * r);
  // fforce
  ans[0] = -tmp_exp*(tmp_fc_d - tmp_fc*param_lam1)*ucl_recip(r);
  // eng
  if (EVFLAG && eflag) ans[1] = tmp_fc * tmp_exp;
}

/* ---------------------------------------------------------------------- */

ucl_inline numtyp zeta(const int param_powermint,
                       const numtyp param_lam3,
                       const numtyp param_bigr,
                       const numtyp param_bigd,
                       const numtyp param_c,
                       const numtyp param_d,
                       const numtyp param_h,
                       const numtyp param_gamma,
                       const numtyp rij,
                       const numtyp rsqik,
                       const numtyp4 delrij,
                       const numtyp4 delrik)
{
  numtyp rik,costheta,arg,ex_delr;

  rik = ucl_sqrt(rsqik);
  costheta = (delrij.x*delrik.x + delrij.y*delrik.y +
              delrij.z*delrik.z) / (rij*rik);

  numtyp t = param_lam3*(rij-rik);
  if (param_powermint == 3) arg = t*t*t;
  else arg = t;

  if (arg > (numtyp)69.0776) ex_delr = (numtyp)1.e30;
  else if (arg < (numtyp)-69.0776) ex_delr = (numtyp)0.0;
  else ex_delr = ucl_exp(arg);

  return ters_fc(rik,param_bigr,param_bigd) *
         ters_gijk(costheta,param_c, param_d, param_h, param_gamma) * ex_delr;
}

/* ---------------------------------------------------------------------- */

ucl_inline void force_zeta(const numtyp param_bigb,
                           const numtyp param_bigr,
                           const numtyp param_bigd,
                           const numtyp param_lam2,
                           const numtyp param_beta,
                           const numtyp param_powern,
                           const numtyp param_c1,
                           const numtyp param_c2,
                           const numtyp param_c3,
                           const numtyp param_c4,
                           const numtyp r,
                           const numtyp zeta_ij,
                           const int eflag,
                           numtyp fpfeng[4])
{
  numtyp fa,fa_d,bij,bij_d;

  fa = ters_fa_d(r,param_bigb,param_bigr,param_bigd,param_lam2,&fa_d);
  bij = ters_bij_d(zeta_ij,param_beta,param_powern,
                   param_c1,param_c2, param_c3, param_c4, &bij_d);
  fpfeng[0] = (numtyp)0.5*bij*fa_d*ucl_recip(r); // fforce
  fpfeng[1] = (numtyp)-0.5*fa*bij_d; // prefactor
  if (EVFLAG && eflag) fpfeng[2] = (numtyp)0.5*bij*fa; // eng
}

/* ----------------------------------------------------------------------
   attractive term
   use param_ij cutoff for rij test
   use param_ijk cutoff for rik test
------------------------------------------------------------------------- */

ucl_inline void attractive(const numtyp param_bigr,
                           const numtyp param_bigd,
                           const int param_powermint,
                           const numtyp param_lam3,
                           const numtyp param_c,
                           const numtyp param_d,
                           const numtyp param_h,
                           const numtyp param_gamma,
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
                  param_c, param_d, param_h, param_gamma, fi, fj, fk);
}

ucl_inline void attractive_fi(const numtyp param_bigr,
                              const numtyp param_bigd,
                              const int param_powermint,
                              const numtyp param_lam3,
                              const numtyp param_c,
                              const numtyp param_d,
                              const numtyp param_h,
                              const numtyp param_gamma,
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
                  param_c, param_d, param_h, param_gamma, fi);
}

ucl_inline void attractive_fj(const numtyp param_bigr,
                              const numtyp param_bigd,
                              const int param_powermint,
                              const numtyp param_lam3,
                              const numtyp param_c,
                              const numtyp param_d,
                              const numtyp param_h,
                              const numtyp param_gamma,
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
                     param_c, param_d, param_h, param_gamma, fj);
}

ucl_inline void attractive_fk(const numtyp param_bigr,
                              const numtyp param_bigd,
                              const int param_powermint,
                              const numtyp param_lam3,
                              const numtyp param_c,
                              const numtyp param_d,
                              const numtyp param_h,
                              const numtyp param_gamma,
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
                     param_c, param_d, param_h, param_gamma, fk);
}


#endif
