/// **************************************************************************
//                              hippo_extra.h
//                             -------------------
//                              Trung Dac Nguyen
//
//  Device code for hippo math routines
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : ndactrung@gmail.com
// ***************************************************************************/*

#ifndef LAL_HIPPO_EXTRA_H
#define LAL_HIPPO_EXTRA_H

#if defined(NV_KERNEL) || defined(USE_HIP)
#include "lal_aux_fun1.h"
#else
#endif

#define MY_PI2 (numtyp)1.57079632679489661923
#define MY_PI4 (numtyp)0.78539816339744830962

/* ----------------------------------------------------------------------
   damprep generates coefficients for the Pauli repulsion
   damping function for powers of the interatomic distance

   literature reference:

   J. A. Rackers and J. W. Ponder, "Classical Pauli Repulsion: An
   Anisotropic, Atomic Multipole Model", Journal of Chemical Physics,
   150, 084104 (2019)
------------------------------------------------------------------------- */

ucl_inline void damprep(const numtyp r, const numtyp r2, const numtyp rr1,
                        const numtyp rr3, const numtyp rr5, const numtyp rr7,
                        const numtyp rr9, const numtyp rr11, const int rorder,
                        const numtyp dmpi, const numtyp dmpk, numtyp dmpik[11])
{
  numtyp r3,r4;
  numtyp r5,r6,r7,r8;
  numtyp s,ds,d2s;
  numtyp d3s,d4s,d5s;
  numtyp dmpi2,dmpk2;
  numtyp dmpi22,dmpi23;
  numtyp dmpi24,dmpi25;
  numtyp dmpi26,dmpi27;
  numtyp dmpk22,dmpk23;
  numtyp dmpk24,dmpk25;
  numtyp dmpk26;
  numtyp eps,diff;
  numtyp expi,expk;
  numtyp dampi,dampk;
  numtyp pre,term,tmp;

  // compute tolerance value for damping exponents

  eps = (numtyp)0.001;
  diff = dmpi-dmpk; // fabs(dmpi-dmpk)
  if (diff < (numtyp)0) diff = -diff;

  // treat the case where alpha damping exponents are equal

  if (diff < eps) {
    r3 = r2 * r;
    r4 = r3 * r;
    r5 = r4 * r;
    r6 = r5 * r;
    r7 = r6 * r;
    dmpi2 = (numtyp)0.5 * dmpi;
    dampi = dmpi2 * r;
    expi = ucl_exp(-dampi);
    dmpi22 = dmpi2 * dmpi2;
    dmpi23 = dmpi22 * dmpi2;
    dmpi24 = dmpi23 * dmpi2;
    dmpi25 = dmpi24 * dmpi2;
    dmpi26 = dmpi25 * dmpi2;
    pre = (numtyp)128.0;
    s = (r + dmpi2*r2 + dmpi22*r3/(numtyp)3.0) * expi;

    ds = (dmpi22*r3 + dmpi23*r4) * expi / (numtyp)3.0;
    d2s = dmpi24 * expi * r5 / (numtyp)9.0;
    d3s = dmpi25 * expi * r6 / (numtyp)45.0;
    d4s = (dmpi25*r6 + dmpi26*r7) * expi / (numtyp)315.0;
    if (rorder >= 11) {
      r8 = r7 * r;
      dmpi27 = dmpi2 * dmpi26;
      d5s = (dmpi25*r6 + dmpi26*r7 + dmpi27*r8/(numtyp)3.0) * expi / (numtyp)945.0;
    }

  // treat the case where alpha damping exponents are unequal

  } else {
    r3 = r2 * r;
    r4 = r3 * r;
    r5 = r4 * r;
    dmpi2 = (numtyp)0.5 * dmpi;
    dmpk2 = (numtyp)0.5 * dmpk;
    dampi = dmpi2 * r;
    dampk = dmpk2 * r;
    expi = ucl_exp(-dampi);
    expk = ucl_exp(-dampk);
    dmpi22 = dmpi2 * dmpi2;
    dmpi23 = dmpi22 * dmpi2;
    dmpi24 = dmpi23 * dmpi2;
    dmpi25 = dmpi24 * dmpi2;
    dmpk22 = dmpk2 * dmpk2;
    dmpk23 = dmpk22 * dmpk2;
    dmpk24 = dmpk23 * dmpk2;
    dmpk25 = dmpk24 * dmpk2;
    term = dmpi22 - dmpk22;
    pre = (numtyp)8192.0 * dmpi23 * dmpk23 / (term*term*term*term); //ucl_powr(term,(numtyp)4.0);
    tmp = (numtyp)4.0 * dmpi2 * dmpk2 / term;
    s = (dampi-tmp)*expk + (dampk+tmp)*expi;

    ds = (dmpi2*dmpk2*r2 - (numtyp)4.0*dmpi2*dmpk22*r/term -
          (numtyp)4.0*dmpi2*dmpk2/term) * expk +
      (dmpi2*dmpk2*r2 + (numtyp)4.0*dmpi22*dmpk2*r/term + (numtyp)4.0*dmpi2*dmpk2/term) * expi;
    d2s = (dmpi2*dmpk2*r2/3.0 + dmpi2*dmpk22*r3/(numtyp)3.0 -
           ((numtyp)4.0/(numtyp)3.0)*dmpi2*dmpk23*r2/term - (numtyp)4.0*dmpi2*dmpk22*r/term -
           (numtyp)4.0*dmpi2*dmpk2/term) * expk +
      (dmpi2*dmpk2*r2/(numtyp)3.0 + dmpi22*dmpk2*r3/(numtyp)3.0 +
       ((numtyp)4.0/(numtyp)3.0)*dmpi23*dmpk2*r2/term + (numtyp)4.0*dmpi22*dmpk2*r/term +
       (numtyp)4.0*dmpi2*dmpk2/term) * expi;
    d3s = (dmpi2*dmpk23*r4/(numtyp)15.0 + dmpi2*dmpk22*r3/(numtyp)5.0 + dmpi2*dmpk2*r2/(numtyp)5.0 -
           ((numtyp)4.0/(numtyp)15.0)*dmpi2*dmpk24*r3/term - ((numtyp)8.0/(numtyp)5.0)*dmpi2*dmpk23*r2/term -
           (numtyp)4.0*dmpi2*dmpk22*r/term - (numtyp)4.0/term*dmpi2*dmpk2) * expk +
      (dmpi23*dmpk2*r4/(numtyp)15.0 + dmpi22*dmpk2*r3/(numtyp)5.0 + dmpi2*dmpk2*r2/(numtyp)5.0 +
       ((numtyp)4.0/(numtyp)15.0)*dmpi24*dmpk2*r3/term + ((numtyp)8.0/(numtyp)5.0)*dmpi23*dmpk2*r2/term +
       (numtyp)4.0*dmpi22*dmpk2*r/term + (numtyp)4.0/term*dmpi2*dmpk2) * expi;
    d4s = (dmpi2*dmpk24*r5/(numtyp)105.0 + ((numtyp)2.0/(numtyp)35.0)*dmpi2*dmpk23*r4 +
           dmpi2*dmpk22*r3/(numtyp)7.0 + dmpi2*dmpk2*r2/(numtyp)7.0 -
           ((numtyp)4.0/(numtyp)105.0)*dmpi2*dmpk25*r4/term - ((numtyp)8.0/21.0)*dmpi2*dmpk24*r3/term -
           ((numtyp)12.0/(numtyp)7.0)*dmpi2*dmpk23*r2/term - (numtyp)4.0*dmpi2*dmpk22*r/term -
           (numtyp)4.0*dmpi2*dmpk2/term) * expk +
      (dmpi24*dmpk2*r5/(numtyp)105.0 + ((numtyp)2.0/(numtyp)35.0)*dmpi23*dmpk2*r4 +
       dmpi22*dmpk2*r3/(numtyp)7.0 + dmpi2*dmpk2*r2/(numtyp)7.0 +
       ((numtyp)4.0/(numtyp)105.0)*dmpi25*dmpk2*r4/term + ((numtyp)8.0/(numtyp)21.0)*dmpi24*dmpk2*r3/term +
       ((numtyp)12.0/(numtyp)7.0)*dmpi23*dmpk2*r2/term + (numtyp)4.0*dmpi22*dmpk2*r/term +
       (numtyp)4.0*dmpi2*dmpk2/term) * expi;

    if (rorder >= 11) {
      r6 = r5 * r;
      dmpi26 = dmpi25 * dmpi2;
      dmpk26 = dmpk25 * dmpk2;
      d5s = (dmpi2*dmpk25*r6/(numtyp)945.0 + ((numtyp)2.0/(numtyp)189.0)*dmpi2*dmpk24*r5 +
             dmpi2*dmpk23*r4/(numtyp)21.0 + dmpi2*dmpk22*r3/(numtyp)9.0 + dmpi2*dmpk2*r2/(numtyp)9.0 -
             ((numtyp)4.0/(numtyp)945.0)*dmpi2*dmpk26*r5/term -
             ((numtyp)4.0/(numtyp)63.0)*dmpi2*dmpk25*r4/term - ((numtyp)4.0/(numtyp)9.0)*dmpi2*dmpk24*r3/term -
             ((numtyp)16.0/(numtyp)9.0)*dmpi2*dmpk23*r2/term - (numtyp)4.0*dmpi2*dmpk22*r/term -
             (numtyp)4.0*dmpi2*dmpk2/term) * expk +
        (dmpi25*dmpk2*r6/(numtyp)945.0 + ((numtyp)2.0/(numtyp)189.0)*dmpi24*dmpk2*r5 +
         dmpi23*dmpk2*r4/(numtyp)21.0 + dmpi22*dmpk2*r3/(numtyp)9.0 + dmpi2*dmpk2*r2/(numtyp)9.0 +
         ((numtyp)4.0/(numtyp)945.0)*dmpi26*dmpk2*r5/term + ((numtyp)4.0/(numtyp)63.0)*dmpi25*dmpk2*r4/term +
         ((numtyp)4.0/(numtyp)9.0)*dmpi24*dmpk2*r3/term + ((numtyp)16.0/(numtyp)9.0)*dmpi23*dmpk2*r2/term +
         (numtyp)4.0*dmpi22*dmpk2*r/term + (numtyp)4.0*dmpi2*dmpk2/term) * expi;
    }
  }

  // convert partial derivatives into full derivatives

  s = s * rr1;
  ds = ds * rr3;
  d2s = d2s * rr5;
  d3s = d3s * rr7;
  d4s = d4s * rr9;
  d5s = d5s * rr11;
  dmpik[0] = (numtyp)0.5 * pre * s * s;
  dmpik[2] = pre * s * ds;
  dmpik[4] = pre * (s*d2s + ds*ds);
  dmpik[6] = pre * (s*d3s + (numtyp)3.0*ds*d2s);
  dmpik[8] = pre * (s*d4s + (numtyp)4.0*ds*d3s + (numtyp)3.0*d2s*d2s);

  if (rorder >= 11) dmpik[10] = pre * (s*d5s + (numtyp)5.0*ds*d4s + (numtyp)10.0*d2s*d3s);
}

/* ----------------------------------------------------------------------
   damppole generates coefficients for the charge penetration
   damping function for powers of the interatomic distance

   literature references:

   L. V. Slipchenko and M. S. Gordon, "Electrostatic Energy in the
   Effective Fragment Potential Method: Theory and Application to
   the Benzene Dimer", Journal of Computational Chemistry, 28,
   276-291 (2007)  [Gordon f1 and f2 models]

   J. A. Rackers, Q. Wang, C. Liu, J.-P. Piquemal, P. Ren and
   J. W. Ponder, "An Optimized Charge Penetration Model for Use with
   the AMOEBA Force Field", Physical Chemistry Chemical Physics, 19,
   276-291 (2017)
------------------------------------------------------------------------- */

ucl_inline void damppole(const numtyp r, const int rorder,
                         const numtyp alphai, const numtyp alphak,
                         numtyp dmpi[9], numtyp dmpk[9], numtyp dmpik[11])
{
  numtyp termi,termk;
  numtyp termi2,termk2;
  numtyp alphai2,alphak2;
  numtyp eps,diff;
  numtyp expi,expk;
  numtyp dampi,dampk;
  numtyp dampi2,dampi3;
  numtyp dampi4,dampi5;
  numtyp dampi6,dampi7;
  numtyp dampi8;
  numtyp dampk2,dampk3;
  numtyp dampk4,dampk5;
  numtyp dampk6;

  // compute tolerance and exponential damping factors

  eps = (numtyp)0.001;
  diff = alphai-alphak;
  if (diff < (numtyp)0) diff = -diff;
  dampi = alphai * r;
  dampk = alphak * r;
  expi = ucl_exp(-dampi);
  expk = ucl_exp(-dampk);

  // core-valence charge penetration damping for Gordon f1

  dampi2 = dampi * dampi;
  dampi3 = dampi * dampi2;
  dampi4 = dampi2 * dampi2;
  dampi5 = dampi2 * dampi3;
  dmpi[0] = (numtyp)1.0 - ((numtyp)1.0 + (numtyp)0.5*dampi)*expi;
  dmpi[2] = (numtyp)1.0 - ((numtyp)1.0 + dampi + (numtyp)0.5*dampi2)*expi;
  dmpi[4] = (numtyp)1.0 - ((numtyp)1.0 + dampi + (numtyp)0.5*dampi2 + dampi3/(numtyp)6.0)*expi;
  dmpi[6] = (numtyp)1.0 - ((numtyp)1.0 + dampi + (numtyp)0.5*dampi2 + dampi3/(numtyp)6.0 + dampi4/(numtyp)30.0)*expi;
  dmpi[8] = (numtyp)1.0 - ((numtyp)1.0 + dampi + (numtyp)0.5*dampi2 + dampi3/(numtyp)6.0 +
                   (numtyp)4.0*dampi4/(numtyp)105.0 + dampi5/(numtyp)210.0)*expi;
  if (diff < eps) {
    dmpk[0] = dmpi[0];
    dmpk[2] = dmpi[2];
    dmpk[4] = dmpi[4];
    dmpk[6] = dmpi[6];
    dmpk[8] = dmpi[8];
  } else {
    dampk2 = dampk * dampk;
    dampk3 = dampk * dampk2;
    dampk4 = dampk2 * dampk2;
    dampk5 = dampk2 * dampk3;
    dmpk[0] = (numtyp)1.0 - ((numtyp)1.0 + (numtyp)0.5*dampk)*expk;
    dmpk[2] = (numtyp)1.0 - ((numtyp)1.0 + dampk + (numtyp)0.5*dampk2)*expk;
    dmpk[4] = (numtyp)1.0 - ((numtyp)1.0 + dampk + (numtyp)0.5*dampk2 + dampk3/(numtyp)6.0)*expk;
    dmpk[6] = (numtyp)1.0 - ((numtyp)1.0 + dampk + (numtyp)0.5*dampk2 + dampk3/(numtyp)6.0 + dampk4/(numtyp)30.0)*expk;
    dmpk[8] = (numtyp)1.0 - ((numtyp)1.0 + dampk + (numtyp)0.5*dampk2 + dampk3/(numtyp)6.0 +
                     (numtyp)4.0*dampk4/(numtyp)105.0 + dampk5/(numtyp)210.0)*expk;
  }

  // valence-valence charge penetration damping for Gordon f1

  if (diff < eps) {
    dampi6 = dampi3 * dampi3;
    dampi7 = dampi3 * dampi4;
    dmpik[0] = (numtyp)1.0 - ((numtyp)1.0 + (numtyp)11.0*dampi/(numtyp)16.0 + (numtyp)3.0*dampi2/(numtyp)16.0 +
                      dampi3/(numtyp)48.0)*expi;
    dmpik[2] = (numtyp)1.0 - ((numtyp)1.0 + dampi + (numtyp)0.5*dampi2 +
                      (numtyp)7.0*dampi3/(numtyp)48.0 + dampi4/(numtyp)48.0)*expi;
    dmpik[4] = (numtyp)1.0 - ((numtyp)1.0 + dampi + (numtyp)0.5*dampi2 + dampi3/(numtyp)6.0 +
                      dampi4/(numtyp)24.0 + dampi5/(numtyp)144.0)*expi;
    dmpik[6] = (numtyp)1.0 - ((numtyp)1.0 + dampi + (numtyp)0.5*dampi2 + dampi3/(numtyp)6.0 +
                      dampi4/(numtyp)24.0 + dampi5/(numtyp)120.0 + dampi6/(numtyp)720.0)*expi;
    dmpik[8] = (numtyp)1.0 - ((numtyp)1.0 + dampi + (numtyp)0.5*dampi2 + dampi3/(numtyp)6.0 +
                      dampi4/(numtyp)24.0 + dampi5/(numtyp)120.0 + dampi6/(numtyp)720.0 +
                      dampi7/(numtyp)5040.0)*expi;
    if (rorder >= 11) {
      dampi8 = dampi4 * dampi4;
      dmpik[10] = (numtyp)1.0 - ((numtyp)1.0 + dampi + (numtyp)0.5*dampi2 + dampi3/(numtyp)6.0 +
                         dampi4/(numtyp)24.0 + dampi5/(numtyp)120.0 + dampi6/(numtyp)720.0 +
                         dampi7/(numtyp)5040.0 + dampi8/(numtyp)45360.0)*expi;
    }

  } else {
    alphai2 = alphai * alphai;
    alphak2 = alphak * alphak;
    termi = alphak2 / (alphak2-alphai2);
    termk = alphai2 / (alphai2-alphak2);
    termi2 = termi * termi;
    termk2 = termk * termk;
    dmpik[0] = (numtyp)1.0 - termi2*(1.0 + (numtyp)2.0*termk + (numtyp)0.5*dampi)*expi -
      termk2*((numtyp)1.0 + (numtyp)2.0*termi + (numtyp)0.5*dampk)*expk;
    dmpik[2] = (numtyp)1.0 - termi2*((numtyp)1.0+dampi+(numtyp)0.5*dampi2)*expi -
      termk2*((numtyp)1.0+dampk+(numtyp)0.5*dampk2)*expk -
      (numtyp)2.0*termi2*termk*((numtyp)1.0+dampi)*expi -
      (numtyp)2.0*termk2*termi*((numtyp)1.0+dampk)*expk;
    dmpik[4] = (numtyp)1.0 - termi2*((numtyp)1.0 + dampi + (numtyp)0.5*dampi2 + dampi3/(numtyp)6.0)*expi -
      termk2*(1.0 + dampk + (numtyp)0.5*dampk2 + dampk3/(numtyp)6.0)*expk -
      (numtyp)2.0*termi2*termk*((numtyp)1.0 + dampi + dampi2/(numtyp)3.0)*expi -
      (numtyp)2.0*termk2*termi*((numtyp)1.0 + dampk + dampk2/(numtyp)3.0)*expk;
    dmpik[6] = (numtyp)1.0 - termi2*((numtyp)1.0 + dampi + (numtyp)0.5*dampi2 +
                             dampi3/(numtyp)6.0 + dampi4/(numtyp)30.0)*expi -
      termk2*((numtyp)1.0 + dampk + 0.5*dampk2 + dampk3/(numtyp)6.0 + dampk4/(numtyp)30.0)*expk -
      (numtyp)2.0*termi2*termk*((numtyp)1.0 + dampi + (numtyp)2.0*dampi2/(numtyp)5.0 + dampi3/(numtyp)15.0)*expi -
      (numtyp)2.0*termk2*termi*((numtyp)1.0 + dampk + (numtyp)2.0*dampk2/(numtyp)5.0 + dampk3/(numtyp)15.0)*expk;
    dmpik[8] = (numtyp)1.0 - termi2*((numtyp)1.0 + dampi + (numtyp)0.5*dampi2 + dampi3/(numtyp)6.0 +
                             (numtyp)4.0*dampi4/(numtyp)105.0 + dampi5/(numtyp)210.0)*expi -
      termk2*((numtyp)1.0 + dampk + (numtyp)0.5*dampk2 + dampk3/(numtyp)6.0 +
              (numtyp)4.0*dampk4/105.0 + dampk5/(numtyp)210.0)*expk -
      (numtyp)2.0*termi2*termk*((numtyp)1.0 + dampi + (numtyp)3.0*dampi2/(numtyp)7.0 +
                        (numtyp)2.0*dampi3/(numtyp)21.0 + dampi4/(numtyp)105.0)*expi -
      (numtyp)2.0*termk2*termi*((numtyp)1.0 + dampk + (numtyp)3.0*dampk2/(numtyp)7.0 +
                        (numtyp)2.0*dampk3/(numtyp)21.0 + dampk4/(numtyp)105.0)*expk;

    if (rorder >= 11) {
      dampi6 = dampi3 * dampi3;
      dampk6 = dampk3 * dampk3;
      dmpik[10] = (numtyp)1.0 - termi2*((numtyp)1.0 + dampi + (numtyp)0.5*dampi2 + dampi3/(numtyp)6.0 +
                                (numtyp)5.0*dampi4/(numtyp)126.0 + (numtyp)2.0*dampi5/(numtyp)315.0 +
                                dampi6/(numtyp)1890.0)*expi -
        termk2*((numtyp)1.0 + dampk + (numtyp)0.5*dampk2 + dampk3/(numtyp)6.0 + (numtyp)5.0*dampk4/(numtyp)126.0 +
                (numtyp)2.0*dampk5/(numtyp)315.0 + dampk6/(numtyp)1890.0)*expk -
        (numtyp)2.0*termi2*termk*((numtyp)1.0 + dampi + (numtyp)4.0*dampi2/(numtyp)9.0 + dampi3/(numtyp)9.0 +
                          dampi4/(numtyp)63.0 + dampi5/(numtyp)945.0)*expi -
        (numtyp)2.0*termk2*termi*((numtyp)1.0 + dampk + 4.0*dampk2/(numtyp)9.0 + dampk3/(numtyp)9.0 +
                          dampk4/(numtyp)63.0 + dampk5/(numtyp)945.0)*expk;
    }
  }
}

/* ----------------------------------------------------------------------
   dampdir = direct field damping coefficents
   dampdir generates coefficients for the direct field damping
   function for powers of the interatomic distance
------------------------------------------------------------------------- */

ucl_inline void dampdir(numtyp r, numtyp alphai, numtyp alphak, numtyp *dmpi, numtyp *dmpk)
{
  numtyp eps,diff;
  numtyp expi,expk;
  numtyp dampi,dampk;
  numtyp dampi2,dampk2;
  numtyp dampi3,dampk3;
  numtyp dampi4,dampk4;

  // compute tolerance and exponential damping factors

  eps = (numtyp)0.001;
  diff = alphai-alphak; // fabs(alphai-alphak);
  if (diff < (numtyp)0) diff = -diff;
  dampi = alphai * r;
  dampk = alphak * r;
  expi = ucl_exp(-dampi);
  expk = ucl_exp(-dampk);

  // core-valence charge penetration damping for Gordon f1 (HIPPO)

  dampi2 = dampi * dampi;
  dampi3 = dampi * dampi2;
  dampi4 = dampi2 * dampi2;
  dmpi[2] = (numtyp)1.0 - ((numtyp)1.0 + dampi + (numtyp)0.5*dampi2)*expi;
  dmpi[4] = (numtyp)1.0 - ((numtyp)1.0 + dampi + (numtyp)0.5*dampi2 + dampi3/(numtyp)6.0)*expi;
  dmpi[6] = (numtyp)1.0 - ((numtyp)1.0 + dampi + (numtyp)0.5*dampi2 + dampi3/(numtyp)6.0 + dampi4/(numtyp)30.0)*expi;
  if (diff < eps) {
    dmpk[2] = dmpi[2];
    dmpk[4] = dmpi[4];
    dmpk[6] = dmpi[6];
  } else {
    dampk2 = dampk * dampk;
    dampk3 = dampk * dampk2;
    dampk4 = dampk2 * dampk2;
    dmpk[2] = (numtyp)1.0 - ((numtyp)1.0 + dampk + (numtyp)0.5*dampk2)*expk;
    dmpk[4] = (numtyp)1.0 - ((numtyp)1.0 + dampk + (numtyp)0.5*dampk2 + dampk3/(numtyp)6.0)*expk;
    dmpk[6] = (numtyp)1.0 - ((numtyp)1.0 + dampk + (numtyp)0.5*dampk2 + dampk3/(numtyp)6.0 + dampk4/30.0)*expk;
  }
}

/* ----------------------------------------------------------------------
   dampmut = mutual field damping coefficents
   dampmut generates coefficients for the mutual field damping
   function for powers of the interatomic distance
------------------------------------------------------------------------- */

ucl_inline void dampmut(numtyp r, numtyp alphai, numtyp alphak, numtyp dmpik[5])
{
  numtyp termi,termk;
  numtyp termi2,termk2;
  numtyp alphai2,alphak2;
  numtyp eps,diff;
  numtyp expi,expk;
  numtyp dampi,dampk;
  numtyp dampi2,dampi3;
  numtyp dampi4,dampi5;
  numtyp dampk2,dampk3;

  // compute tolerance and exponential damping factors

  eps = (numtyp)0.001;
  diff = alphai-alphak; // fabs(alphai-alphak);
  if (diff < (numtyp)0) diff = -diff;
  dampi = alphai * r;
  dampk = alphak * r;
  expi = ucl_exp(-dampi);
  expk = ucl_exp(-dampk);

  // valence-valence charge penetration damping for Gordon f1 (HIPPO)

  dampi2 = dampi * dampi;
  dampi3 = dampi * dampi2;
  if (diff < eps) {
    dampi4 = dampi2 * dampi2;
    dampi5 = dampi2 * dampi3;
    dmpik[2] = (numtyp)1.0 - ((numtyp)1.0 + dampi + (numtyp)0.5*dampi2 +
                      7.0*dampi3/(numtyp)48.0 + dampi4/48.0)*expi;
    dmpik[4] = (numtyp)1.0 - ((numtyp)1.0 + dampi + (numtyp)0.5*dampi2 + dampi3/(numtyp)6.0 +
                      dampi4/(numtyp)24.0 + dampi5/(numtyp)144.0)*expi;
  } else {
    dampk2 = dampk * dampk;
    dampk3 = dampk * dampk2;
    alphai2 = alphai * alphai;
    alphak2 = alphak * alphak;
    termi = alphak2 / (alphak2-alphai2);
    termk = alphai2 / (alphai2-alphak2);
    termi2 = termi * termi;
    termk2 = termk * termk;
    dmpik[2] = (numtyp)1.0 - termi2*((numtyp)1.0+dampi+(numtyp)0.5*dampi2)*expi -
      termk2*((numtyp)1.0+dampk+(numtyp)0.5*dampk2)*expk -
      (numtyp)2.0*termi2*termk*((numtyp)1.0+dampi)*expi - (numtyp)2.0*termk2*termi*((numtyp)1.0+dampk)*expk;
    dmpik[4] = (numtyp)1.0 - termi2*((numtyp)1.0+dampi+(numtyp)0.5*dampi2 + dampi3/(numtyp)6.0)*expi -
      termk2*((numtyp)1.0+dampk+(numtyp)0.5*dampk2 + dampk3/(numtyp)6.00)*expk -
      (numtyp)2.0*termi2*termk *((numtyp)1.0+dampi+dampi2/(numtyp)3.0)*expi -
      (numtyp)2.0*termk2*termi *((numtyp)1.0+dampk+dampk2/(numtyp)3.0)*expk;
  }
}

#endif
