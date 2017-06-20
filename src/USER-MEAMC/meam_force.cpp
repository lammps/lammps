#include "meam.h"
#include <math.h>
#include <algorithm>
#include "math_special.h"

using namespace LAMMPS_NS;
//     Extern "C" declaration has the form:
//
//  void meam_force_(int *, int *, int *, double *, int *, int *, int *, double
//  *,
//		 int *, int *, int *, int *, double *, double *,
//		 double *, double *, double *, double *, double *, double *,
//		 double *, double *, double *, double *, double *, double *,
//		 double *, double *, double *, double *, double *, double *, int
//*);
//
//     Call from pair_meam.cpp has the form:
//
//    meam_force_(&i,&nmax,&eflag_either,&eflag_global,&eflag_atom,&vflag_atom,
//              &eng_vdwl,eatom,&ntype,type,fmap,&x[0][0],
//	       &numneigh[i],firstneigh[i],&numneigh_full[i],firstneigh_full[i],
//	       &scrfcn[offset],&dscrfcn[offset],&fcpair[offset],
//	       dgamma1,dgamma2,dgamma3,rho0,rho1,rho2,rho3,frhop,
//	       &arho1[0][0],&arho2[0][0],arho2b,&arho3[0][0],&arho3b[0][0],
//	       &t_ave[0][0],&tsq_ave[0][0],&f[0][0],&vatom[0][0],&errorflag);
//

void
MEAM::meam_force(int* iptr, int* eflag_either, int* eflag_global,
            int* eflag_atom, int* vflag_atom, double* eng_vdwl, double* eatom,
            int* ntype, int* type, int* fmap, double** x, int* numneigh,
            int* firstneigh, int* numneigh_full, int* firstneigh_full,
            int fnoffset, double** f, double** vatom,
            int* errorflag)
{
  int i, j, jn, k, kn, kk, m, n, p, q;
  int nv2, nv3, elti, eltj, eltk, ind;
  double xitmp, yitmp, zitmp, delij[3 + 1], rij2, rij, rij3;
  double delik[3 + 1], deljk[3 + 1], v[6 + 1], fi[3 + 1], fj[3 + 1];
  double third, sixth;
  double pp, dUdrij, dUdsij, dUdrijm[3 + 1], force, forcem;
  double r, recip, phi, phip;
  double sij;
  double a1, a1i, a1j, a2, a2i, a2j;
  double a3i, a3j;
  double shpi[3 + 1], shpj[3 + 1];
  double ai, aj, ro0i, ro0j, invrei, invrej;
  double rhoa0j, drhoa0j, rhoa0i, drhoa0i;
  double rhoa1j, drhoa1j, rhoa1i, drhoa1i;
  double rhoa2j, drhoa2j, rhoa2i, drhoa2i;
  double a3, a3a, rhoa3j, drhoa3j, rhoa3i, drhoa3i;
  double drho0dr1, drho0dr2, drho0ds1, drho0ds2;
  double drho1dr1, drho1dr2, drho1ds1, drho1ds2;
  double drho1drm1[3 + 1], drho1drm2[3 + 1];
  double drho2dr1, drho2dr2, drho2ds1, drho2ds2;
  double drho2drm1[3 + 1], drho2drm2[3 + 1];
  double drho3dr1, drho3dr2, drho3ds1, drho3ds2;
  double drho3drm1[3 + 1], drho3drm2[3 + 1];
  double dt1dr1, dt1dr2, dt1ds1, dt1ds2;
  double dt2dr1, dt2dr2, dt2ds1, dt2ds2;
  double dt3dr1, dt3dr2, dt3ds1, dt3ds2;
  double drhodr1, drhodr2, drhods1, drhods2, drhodrm1[3 + 1], drhodrm2[3 + 1];
  double arg;
  double arg1i1, arg1j1, arg1i2, arg1j2, arg1i3, arg1j3, arg3i3, arg3j3;
  double dsij1, dsij2, force1, force2;
  double t1i, t2i, t3i, t1j, t2j, t3j;

  *errorflag = 0;
  third = 1.0 / 3.0;
  sixth = 1.0 / 6.0;

  //: aliased
  i = *iptr;

  //     Compute forces atom i

  elti = fmap[arr1v(type, i)];

  if (elti >= 0) {
    xitmp = arr2v(x, 1, i);
    yitmp = arr2v(x, 2, i);
    zitmp = arr2v(x, 3, i);

    //     Treat each pair
    for (jn = 1; jn <= *numneigh; jn++) {
      j = arr1v(firstneigh, jn);
      eltj = fmap[arr1v(type, j)];

      if (!iszero(arr1v(scrfcn, fnoffset + jn)) && eltj >= 0) {

        sij = arr1v(scrfcn, fnoffset + jn) * arr1v(fcpair, fnoffset + jn);
        delij[1] = arr2v(x, 1, j) - xitmp;
        delij[2] = arr2v(x, 2, j) - yitmp;
        delij[3] = arr2v(x, 3, j) - zitmp;
        rij2 = delij[1] * delij[1] + delij[2] * delij[2] + delij[3] * delij[3];
        if (rij2 < this->cutforcesq) {
          rij = sqrt(rij2);
          r = rij;

          //     Compute phi and phip
          ind = this->eltind[elti][eltj];
          pp = rij * this->rdrar;
          kk = (int)pp;
          kk = std::min(kk, this->nrar - 2);
          pp = pp - kk;
          pp = std::min(pp, 1.0);
          phi = ((this->phirar3[ind][kk] * pp +
                  this->phirar2[ind][kk]) *
                   pp +
                 this->phirar1[ind][kk]) *
                  pp +
                this->phirar[ind][kk];
          phip = (this->phirar6[ind][kk] * pp +
                  this->phirar5[ind][kk]) *
                   pp +
                 this->phirar4[ind][kk];
          recip = 1.0 / r;

          if (*eflag_either != 0) {
            if (*eflag_global != 0)
              *eng_vdwl = *eng_vdwl + phi * sij;
            if (*eflag_atom != 0) {
              arr1v(eatom, i) = arr1v(eatom, i) + 0.5 * phi * sij;
              arr1v(eatom, j) = arr1v(eatom, j) + 0.5 * phi * sij;
            }
          }

          //     write(1,*) "force_meamf: phi: ",phi
          //     write(1,*) "force_meamf: phip: ",phip

          //     Compute pair densities and derivatives
          invrei = 1.0 / this->re_meam[elti][elti];
          ai = rij * invrei - 1.0;
          ro0i = this->rho0_meam[elti];
          rhoa0i = ro0i * MathSpecial::fm_exp(-this->beta0_meam[elti] * ai);
          drhoa0i = -this->beta0_meam[elti] * invrei * rhoa0i;
          rhoa1i = ro0i * MathSpecial::fm_exp(-this->beta1_meam[elti] * ai);
          drhoa1i = -this->beta1_meam[elti] * invrei * rhoa1i;
          rhoa2i = ro0i * MathSpecial::fm_exp(-this->beta2_meam[elti] * ai);
          drhoa2i = -this->beta2_meam[elti] * invrei * rhoa2i;
          rhoa3i = ro0i * MathSpecial::fm_exp(-this->beta3_meam[elti] * ai);
          drhoa3i = -this->beta3_meam[elti] * invrei * rhoa3i;

          if (elti != eltj) {
            invrej = 1.0 / this->re_meam[eltj][eltj];
            aj = rij * invrej - 1.0;
            ro0j = this->rho0_meam[eltj];
            rhoa0j = ro0j * MathSpecial::fm_exp(-this->beta0_meam[eltj] * aj);
            drhoa0j = -this->beta0_meam[eltj] * invrej * rhoa0j;
            rhoa1j = ro0j * MathSpecial::fm_exp(-this->beta1_meam[eltj] * aj);
            drhoa1j = -this->beta1_meam[eltj] * invrej * rhoa1j;
            rhoa2j = ro0j * MathSpecial::fm_exp(-this->beta2_meam[eltj] * aj);
            drhoa2j = -this->beta2_meam[eltj] * invrej * rhoa2j;
            rhoa3j = ro0j * MathSpecial::fm_exp(-this->beta3_meam[eltj] * aj);
            drhoa3j = -this->beta3_meam[eltj] * invrej * rhoa3j;
          } else {
            rhoa0j = rhoa0i;
            drhoa0j = drhoa0i;
            rhoa1j = rhoa1i;
            drhoa1j = drhoa1i;
            rhoa2j = rhoa2i;
            drhoa2j = drhoa2i;
            rhoa3j = rhoa3i;
            drhoa3j = drhoa3i;
          }

          if (this->ialloy == 1) {
            rhoa1j = rhoa1j * this->t1_meam[eltj];
            rhoa2j = rhoa2j * this->t2_meam[eltj];
            rhoa3j = rhoa3j * this->t3_meam[eltj];
            rhoa1i = rhoa1i * this->t1_meam[elti];
            rhoa2i = rhoa2i * this->t2_meam[elti];
            rhoa3i = rhoa3i * this->t3_meam[elti];
            drhoa1j = drhoa1j * this->t1_meam[eltj];
            drhoa2j = drhoa2j * this->t2_meam[eltj];
            drhoa3j = drhoa3j * this->t3_meam[eltj];
            drhoa1i = drhoa1i * this->t1_meam[elti];
            drhoa2i = drhoa2i * this->t2_meam[elti];
            drhoa3i = drhoa3i * this->t3_meam[elti];
          }

          nv2 = 1;
          nv3 = 1;
          arg1i1 = 0.0;
          arg1j1 = 0.0;
          arg1i2 = 0.0;
          arg1j2 = 0.0;
          arg1i3 = 0.0;
          arg1j3 = 0.0;
          arg3i3 = 0.0;
          arg3j3 = 0.0;
          for (n = 1; n <= 3; n++) {
            for (p = n; p <= 3; p++) {
              for (q = p; q <= 3; q++) {
                arg = delij[n] * delij[p] * delij[q] * this->v3D[nv3];
                arg1i3 = arg1i3 + arr2v(arho3, nv3, i) * arg;
                arg1j3 = arg1j3 - arr2v(arho3, nv3, j) * arg;
                nv3 = nv3 + 1;
              }
              arg = delij[n] * delij[p] * this->v2D[nv2];
              arg1i2 = arg1i2 + arr2v(arho2, nv2, i) * arg;
              arg1j2 = arg1j2 + arr2v(arho2, nv2, j) * arg;
              nv2 = nv2 + 1;
            }
            arg1i1 = arg1i1 + arr2v(arho1, n, i) * delij[n];
            arg1j1 = arg1j1 - arr2v(arho1, n, j) * delij[n];
            arg3i3 = arg3i3 + arr2v(arho3b, n, i) * delij[n];
            arg3j3 = arg3j3 - arr2v(arho3b, n, j) * delij[n];
          }

          //     rho0 terms
          drho0dr1 = drhoa0j * sij;
          drho0dr2 = drhoa0i * sij;

          //     rho1 terms
          a1 = 2 * sij / rij;
          drho1dr1 = a1 * (drhoa1j - rhoa1j / rij) * arg1i1;
          drho1dr2 = a1 * (drhoa1i - rhoa1i / rij) * arg1j1;
          a1 = 2.0 * sij / rij;
          for (m = 1; m <= 3; m++) {
            drho1drm1[m] = a1 * rhoa1j * arr2v(arho1, m, i);
            drho1drm2[m] = -a1 * rhoa1i * arr2v(arho1, m, j);
          }

          //     rho2 terms
          a2 = 2 * sij / rij2;
          drho2dr1 = a2 * (drhoa2j - 2 * rhoa2j / rij) * arg1i2 -
                     2.0 / 3.0 * arr1v(arho2b, i) * drhoa2j * sij;
          drho2dr2 = a2 * (drhoa2i - 2 * rhoa2i / rij) * arg1j2 -
                     2.0 / 3.0 * arr1v(arho2b, j) * drhoa2i * sij;
          a2 = 4 * sij / rij2;
          for (m = 1; m <= 3; m++) {
            drho2drm1[m] = 0.0;
            drho2drm2[m] = 0.0;
            for (n = 1; n <= 3; n++) {
              drho2drm1[m] = drho2drm1[m] +
                             arr2v(arho2, this->vind2D[m][n], i) * delij[n];
              drho2drm2[m] = drho2drm2[m] -
                             arr2v(arho2, this->vind2D[m][n], j) * delij[n];
            }
            drho2drm1[m] = a2 * rhoa2j * drho2drm1[m];
            drho2drm2[m] = -a2 * rhoa2i * drho2drm2[m];
          }

          //     rho3 terms
          rij3 = rij * rij2;
          a3 = 2 * sij / rij3;
          a3a = 6.0 / 5.0 * sij / rij;
          drho3dr1 = a3 * (drhoa3j - 3 * rhoa3j / rij) * arg1i3 -
                     a3a * (drhoa3j - rhoa3j / rij) * arg3i3;
          drho3dr2 = a3 * (drhoa3i - 3 * rhoa3i / rij) * arg1j3 -
                     a3a * (drhoa3i - rhoa3i / rij) * arg3j3;
          a3 = 6 * sij / rij3;
          a3a = 6 * sij / (5 * rij);
          for (m = 1; m <= 3; m++) {
            drho3drm1[m] = 0.0;
            drho3drm2[m] = 0.0;
            nv2 = 1;
            for (n = 1; n <= 3; n++) {
              for (p = n; p <= 3; p++) {
                arg = delij[n] * delij[p] * this->v2D[nv2];
                drho3drm1[m] = drho3drm1[m] +
                               arr2v(arho3, this->vind3D[m][n][p], i) * arg;
                drho3drm2[m] = drho3drm2[m] +
                               arr2v(arho3, this->vind3D[m][n][p], j) * arg;
                nv2 = nv2 + 1;
              }
            }
            drho3drm1[m] =
              (a3 * drho3drm1[m] - a3a * arr2v(arho3b, m, i)) * rhoa3j;
            drho3drm2[m] =
              (-a3 * drho3drm2[m] + a3a * arr2v(arho3b, m, j)) * rhoa3i;
          }

          //     Compute derivatives of weighting functions t wrt rij
          t1i = arr2v(t_ave, 1, i);
          t2i = arr2v(t_ave, 2, i);
          t3i = arr2v(t_ave, 3, i);
          t1j = arr2v(t_ave, 1, j);
          t2j = arr2v(t_ave, 2, j);
          t3j = arr2v(t_ave, 3, j);

          if (this->ialloy == 1) {

            a1i = 0.0;
            a1j = 0.0;
            a2i = 0.0;
            a2j = 0.0;
            a3i = 0.0;
            a3j = 0.0;
            if (!iszero(arr2v(tsq_ave, 1, i)))
              a1i = drhoa0j * sij / arr2v(tsq_ave, 1, i);
            if (!iszero(arr2v(tsq_ave, 1, j)))
              a1j = drhoa0i * sij / arr2v(tsq_ave, 1, j);
            if (!iszero(arr2v(tsq_ave, 2, i)))
              a2i = drhoa0j * sij / arr2v(tsq_ave, 2, i);
            if (!iszero(arr2v(tsq_ave, 2, j)))
              a2j = drhoa0i * sij / arr2v(tsq_ave, 2, j);
            if (!iszero(arr2v(tsq_ave, 3, i)))
              a3i = drhoa0j * sij / arr2v(tsq_ave, 3, i);
            if (!iszero(arr2v(tsq_ave, 3, j)))
              a3j = drhoa0i * sij / arr2v(tsq_ave, 3, j);

            dt1dr1 = a1i * (this->t1_meam[eltj] -
                            t1i * pow(this->t1_meam[eltj], 2));
            dt1dr2 = a1j * (this->t1_meam[elti] -
                            t1j * pow(this->t1_meam[elti], 2));
            dt2dr1 = a2i * (this->t2_meam[eltj] -
                            t2i * pow(this->t2_meam[eltj], 2));
            dt2dr2 = a2j * (this->t2_meam[elti] -
                            t2j * pow(this->t2_meam[elti], 2));
            dt3dr1 = a3i * (this->t3_meam[eltj] -
                            t3i * pow(this->t3_meam[eltj], 2));
            dt3dr2 = a3j * (this->t3_meam[elti] -
                            t3j * pow(this->t3_meam[elti], 2));

          } else if (this->ialloy == 2) {

            dt1dr1 = 0.0;
            dt1dr2 = 0.0;
            dt2dr1 = 0.0;
            dt2dr2 = 0.0;
            dt3dr1 = 0.0;
            dt3dr2 = 0.0;

          } else {

            ai = 0.0;
            if (!iszero(arr1v(rho0, i)))
              ai = drhoa0j * sij / arr1v(rho0, i);
            aj = 0.0;
            if (!iszero(arr1v(rho0, j)))
              aj = drhoa0i * sij / arr1v(rho0, j);

            dt1dr1 = ai * (this->t1_meam[elti] - t1i);
            dt1dr2 = aj * (this->t1_meam[elti] - t1j);
            dt2dr1 = ai * (this->t2_meam[elti] - t2i);
            dt2dr2 = aj * (this->t2_meam[elti] - t2j);
            dt3dr1 = ai * (this->t3_meam[elti] - t3i);
            dt3dr2 = aj * (this->t3_meam[elti] - t3j);
          }

          //     Compute derivatives of total density wrt rij, sij and rij(3)
          get_shpfcn(shpi, this->lattce_meam[elti][elti]);
          get_shpfcn(shpj, this->lattce_meam[eltj][eltj]);
          drhodr1 =
            arr1v(dgamma1, i) * drho0dr1 +
            arr1v(dgamma2, i) * (dt1dr1 * arr1v(rho1, i) + t1i * drho1dr1 +
                                 dt2dr1 * arr1v(rho2, i) + t2i * drho2dr1 +
                                 dt3dr1 * arr1v(rho3, i) + t3i * drho3dr1) -
            arr1v(dgamma3, i) *
              (shpi[1] * dt1dr1 + shpi[2] * dt2dr1 + shpi[3] * dt3dr1);
          drhodr2 =
            arr1v(dgamma1, j) * drho0dr2 +
            arr1v(dgamma2, j) * (dt1dr2 * arr1v(rho1, j) + t1j * drho1dr2 +
                                 dt2dr2 * arr1v(rho2, j) + t2j * drho2dr2 +
                                 dt3dr2 * arr1v(rho3, j) + t3j * drho3dr2) -
            arr1v(dgamma3, j) *
              (shpj[1] * dt1dr2 + shpj[2] * dt2dr2 + shpj[3] * dt3dr2);
          for (m = 1; m <= 3; m++) {
            drhodrm1[m] = 0.0;
            drhodrm2[m] = 0.0;
            drhodrm1[m] =
              arr1v(dgamma2, i) *
              (t1i * drho1drm1[m] + t2i * drho2drm1[m] + t3i * drho3drm1[m]);
            drhodrm2[m] =
              arr1v(dgamma2, j) *
              (t1j * drho1drm2[m] + t2j * drho2drm2[m] + t3j * drho3drm2[m]);
          }

          //     Compute derivatives wrt sij, but only if necessary
          if (!iszero(arr1v(dscrfcn, fnoffset + jn))) {
            drho0ds1 = rhoa0j;
            drho0ds2 = rhoa0i;
            a1 = 2.0 / rij;
            drho1ds1 = a1 * rhoa1j * arg1i1;
            drho1ds2 = a1 * rhoa1i * arg1j1;
            a2 = 2.0 / rij2;
            drho2ds1 =
              a2 * rhoa2j * arg1i2 - 2.0 / 3.0 * arr1v(arho2b, i) * rhoa2j;
            drho2ds2 =
              a2 * rhoa2i * arg1j2 - 2.0 / 3.0 * arr1v(arho2b, j) * rhoa2i;
            a3 = 2.0 / rij3;
            a3a = 6.0 / (5.0 * rij);
            drho3ds1 = a3 * rhoa3j * arg1i3 - a3a * rhoa3j * arg3i3;
            drho3ds2 = a3 * rhoa3i * arg1j3 - a3a * rhoa3i * arg3j3;

            if (this->ialloy == 1) {

              a1i = 0.0;
              a1j = 0.0;
              a2i = 0.0;
              a2j = 0.0;
              a3i = 0.0;
              a3j = 0.0;
              if (!iszero(arr2v(tsq_ave, 1, i)))
                a1i = rhoa0j / arr2v(tsq_ave, 1, i);
              if (!iszero(arr2v(tsq_ave, 1, j)))
                a1j = rhoa0i / arr2v(tsq_ave, 1, j);
              if (!iszero(arr2v(tsq_ave, 2, i)))
                a2i = rhoa0j / arr2v(tsq_ave, 2, i);
              if (!iszero(arr2v(tsq_ave, 2, j)))
                a2j = rhoa0i / arr2v(tsq_ave, 2, j);
              if (!iszero(arr2v(tsq_ave, 3, i)))
                a3i = rhoa0j / arr2v(tsq_ave, 3, i);
              if (!iszero(arr2v(tsq_ave, 3, j)))
                a3j = rhoa0i / arr2v(tsq_ave, 3, j);

              dt1ds1 = a1i * (this->t1_meam[eltj] -
                              t1i * pow(this->t1_meam[eltj], 2));
              dt1ds2 = a1j * (this->t1_meam[elti] -
                              t1j * pow(this->t1_meam[elti], 2));
              dt2ds1 = a2i * (this->t2_meam[eltj] -
                              t2i * pow(this->t2_meam[eltj], 2));
              dt2ds2 = a2j * (this->t2_meam[elti] -
                              t2j * pow(this->t2_meam[elti], 2));
              dt3ds1 = a3i * (this->t3_meam[eltj] -
                              t3i * pow(this->t3_meam[eltj], 2));
              dt3ds2 = a3j * (this->t3_meam[elti] -
                              t3j * pow(this->t3_meam[elti], 2));

            } else if (this->ialloy == 2) {

              dt1ds1 = 0.0;
              dt1ds2 = 0.0;
              dt2ds1 = 0.0;
              dt2ds2 = 0.0;
              dt3ds1 = 0.0;
              dt3ds2 = 0.0;

            } else {

              ai = 0.0;
              if (!iszero(arr1v(rho0, i)))
                ai = rhoa0j / arr1v(rho0, i);
              aj = 0.0;
              if (!iszero(arr1v(rho0, j)))
                aj = rhoa0i / arr1v(rho0, j);

              dt1ds1 = ai * (this->t1_meam[eltj] - t1i);
              dt1ds2 = aj * (this->t1_meam[elti] - t1j);
              dt2ds1 = ai * (this->t2_meam[eltj] - t2i);
              dt2ds2 = aj * (this->t2_meam[elti] - t2j);
              dt3ds1 = ai * (this->t3_meam[eltj] - t3i);
              dt3ds2 = aj * (this->t3_meam[elti] - t3j);
            }

            drhods1 =
              arr1v(dgamma1, i) * drho0ds1 +
              arr1v(dgamma2, i) * (dt1ds1 * arr1v(rho1, i) + t1i * drho1ds1 +
                                   dt2ds1 * arr1v(rho2, i) + t2i * drho2ds1 +
                                   dt3ds1 * arr1v(rho3, i) + t3i * drho3ds1) -
              arr1v(dgamma3, i) *
                (shpi[1] * dt1ds1 + shpi[2] * dt2ds1 + shpi[3] * dt3ds1);
            drhods2 =
              arr1v(dgamma1, j) * drho0ds2 +
              arr1v(dgamma2, j) * (dt1ds2 * arr1v(rho1, j) + t1j * drho1ds2 +
                                   dt2ds2 * arr1v(rho2, j) + t2j * drho2ds2 +
                                   dt3ds2 * arr1v(rho3, j) + t3j * drho3ds2) -
              arr1v(dgamma3, j) *
                (shpj[1] * dt1ds2 + shpj[2] * dt2ds2 + shpj[3] * dt3ds2);
          }

          //     Compute derivatives of energy wrt rij, sij and rij[3]
          dUdrij = phip * sij + arr1v(frhop, i) * drhodr1 + arr1v(frhop, j) * drhodr2;
          dUdsij = 0.0;
          if (!iszero(arr1v(dscrfcn, fnoffset + jn))) {
            dUdsij = phi + arr1v(frhop, i) * drhods1 + arr1v(frhop, j) * drhods2;
          }
          for (m = 1; m <= 3; m++) {
            dUdrijm[m] =
              arr1v(frhop, i) * drhodrm1[m] + arr1v(frhop, j) * drhodrm2[m];
          }

          //     Add the part of the force due to dUdrij and dUdsij

          force = dUdrij * recip + dUdsij * arr1v(dscrfcn, fnoffset + jn);
          for (m = 1; m <= 3; m++) {
            forcem = delij[m] * force + dUdrijm[m];
            arr2v(f, m, i) = arr2v(f, m, i) + forcem;
            arr2v(f, m, j) = arr2v(f, m, j) - forcem;
          }

          //     Tabulate per-atom virial as symmetrized stress tensor

          if (*vflag_atom != 0) {
            fi[1] = delij[1] * force + dUdrijm[1];
            fi[2] = delij[2] * force + dUdrijm[2];
            fi[3] = delij[3] * force + dUdrijm[3];
            v[1] = -0.5 * (delij[1] * fi[1]);
            v[2] = -0.5 * (delij[2] * fi[2]);
            v[3] = -0.5 * (delij[3] * fi[3]);
            v[4] = -0.25 * (delij[1] * fi[2] + delij[2] * fi[1]);
            v[5] = -0.25 * (delij[1] * fi[3] + delij[3] * fi[1]);
            v[6] = -0.25 * (delij[2] * fi[3] + delij[3] * fi[2]);

            arr2v(vatom, 1, i) = arr2v(vatom, 1, i) + v[1];
            arr2v(vatom, 2, i) = arr2v(vatom, 2, i) + v[2];
            arr2v(vatom, 3, i) = arr2v(vatom, 3, i) + v[3];
            arr2v(vatom, 4, i) = arr2v(vatom, 4, i) + v[4];
            arr2v(vatom, 5, i) = arr2v(vatom, 5, i) + v[5];
            arr2v(vatom, 6, i) = arr2v(vatom, 6, i) + v[6];
            arr2v(vatom, 1, j) = arr2v(vatom, 1, j) + v[1];
            arr2v(vatom, 2, j) = arr2v(vatom, 2, j) + v[2];
            arr2v(vatom, 3, j) = arr2v(vatom, 3, j) + v[3];
            arr2v(vatom, 4, j) = arr2v(vatom, 4, j) + v[4];
            arr2v(vatom, 5, j) = arr2v(vatom, 5, j) + v[5];
            arr2v(vatom, 6, j) = arr2v(vatom, 6, j) + v[6];
          }

          //     Now compute forces on other atoms k due to change in sij

          if (iszero(sij) || iszero(sij - 1.0))
            continue; //: cont jn loop
          for (kn = 1; kn <= *numneigh_full; kn++) {
            k = arr1v(firstneigh_full, kn);
            eltk = fmap[arr1v(type, k)];
            if (k != j && eltk >= 0) {
              dsij(i, j, k, jn, *numneigh, rij2, &dsij1, &dsij2, *ntype,
                   type, fmap, x, &scrfcn[fnoffset], &fcpair[fnoffset]);
              if (!iszero(dsij1) || !iszero(dsij2)) {
                force1 = dUdsij * dsij1;
                force2 = dUdsij * dsij2;
                for (m = 1; m <= 3; m++) {
                  delik[m] = arr2v(x, m, k) - arr2v(x, m, i);
                  deljk[m] = arr2v(x, m, k) - arr2v(x, m, j);
                }
                for (m = 1; m <= 3; m++) {
                  arr2v(f, m, i) = arr2v(f, m, i) + force1 * delik[m];
                  arr2v(f, m, j) = arr2v(f, m, j) + force2 * deljk[m];
                  arr2v(f, m, k) =
                    arr2v(f, m, k) - force1 * delik[m] - force2 * deljk[m];
                }

                //     Tabulate per-atom virial as symmetrized stress tensor

                if (*vflag_atom != 0) {
                  fi[1] = force1 * delik[1];
                  fi[2] = force1 * delik[2];
                  fi[3] = force1 * delik[3];
                  fj[1] = force2 * deljk[1];
                  fj[2] = force2 * deljk[2];
                  fj[3] = force2 * deljk[3];
                  v[1] = -third * (delik[1] * fi[1] + deljk[1] * fj[1]);
                  v[2] = -third * (delik[2] * fi[2] + deljk[2] * fj[2]);
                  v[3] = -third * (delik[3] * fi[3] + deljk[3] * fj[3]);
                  v[4] = -sixth * (delik[1] * fi[2] + deljk[1] * fj[2] +
                                   delik[2] * fi[1] + deljk[2] * fj[1]);
                  v[5] = -sixth * (delik[1] * fi[3] + deljk[1] * fj[3] +
                                   delik[3] * fi[1] + deljk[3] * fj[1]);
                  v[6] = -sixth * (delik[2] * fi[3] + deljk[2] * fj[3] +
                                   delik[3] * fi[2] + deljk[3] * fj[2]);

                  arr2v(vatom, 1, i) = arr2v(vatom, 1, i) + v[1];
                  arr2v(vatom, 2, i) = arr2v(vatom, 2, i) + v[2];
                  arr2v(vatom, 3, i) = arr2v(vatom, 3, i) + v[3];
                  arr2v(vatom, 4, i) = arr2v(vatom, 4, i) + v[4];
                  arr2v(vatom, 5, i) = arr2v(vatom, 5, i) + v[5];
                  arr2v(vatom, 6, i) = arr2v(vatom, 6, i) + v[6];
                  arr2v(vatom, 1, j) = arr2v(vatom, 1, j) + v[1];
                  arr2v(vatom, 2, j) = arr2v(vatom, 2, j) + v[2];
                  arr2v(vatom, 3, j) = arr2v(vatom, 3, j) + v[3];
                  arr2v(vatom, 4, j) = arr2v(vatom, 4, j) + v[4];
                  arr2v(vatom, 5, j) = arr2v(vatom, 5, j) + v[5];
                  arr2v(vatom, 6, j) = arr2v(vatom, 6, j) + v[6];
                  arr2v(vatom, 1, k) = arr2v(vatom, 1, k) + v[1];
                  arr2v(vatom, 2, k) = arr2v(vatom, 2, k) + v[2];
                  arr2v(vatom, 3, k) = arr2v(vatom, 3, k) + v[3];
                  arr2v(vatom, 4, k) = arr2v(vatom, 4, k) + v[4];
                  arr2v(vatom, 5, k) = arr2v(vatom, 5, k) + v[5];
                  arr2v(vatom, 6, k) = arr2v(vatom, 6, k) + v[6];
                }
              }
            }
            //     end of k loop
          }
        }
      }
      //     end of j loop
    }

    //     else if elti<0, this is not a meam atom
  }
}
