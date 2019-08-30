#include "meam.h"
#include <cmath>
#include <algorithm>
#include "math_special.h"

using namespace LAMMPS_NS;


void
MEAM::meam_force(int i, int eflag_either, int eflag_global, int eflag_atom, int vflag_atom, double* eng_vdwl,
                 double* eatom, int /*ntype*/, int* type, int* fmap, double** scale, double** x, int numneigh, int* firstneigh,
                 int numneigh_full, int* firstneigh_full, int fnoffset, double** f, double** vatom)
{
  int j, jn, k, kn, kk, m, n, p, q;
  int nv2, nv3, elti, eltj, eltk, ind;
  double xitmp, yitmp, zitmp, delij[3], rij2, rij, rij3;
  double v[6], fi[3], fj[3];
  double third, sixth;
  double pp, dUdrij, dUdsij, dUdrijm[3], force, forcem;
  double r, recip, phi, phip;
  double sij;
  double a1, a1i, a1j, a2, a2i, a2j;
  double a3i, a3j;
  double shpi[3], shpj[3];
  double ai, aj, ro0i, ro0j, invrei, invrej;
  double rhoa0j, drhoa0j, rhoa0i, drhoa0i;
  double rhoa1j, drhoa1j, rhoa1i, drhoa1i;
  double rhoa2j, drhoa2j, rhoa2i, drhoa2i;
  double a3, a3a, rhoa3j, drhoa3j, rhoa3i, drhoa3i;
  double drho0dr1, drho0dr2, drho0ds1, drho0ds2;
  double drho1dr1, drho1dr2, drho1ds1, drho1ds2;
  double drho1drm1[3], drho1drm2[3];
  double drho2dr1, drho2dr2, drho2ds1, drho2ds2;
  double drho2drm1[3], drho2drm2[3];
  double drho3dr1, drho3dr2, drho3ds1, drho3ds2;
  double drho3drm1[3], drho3drm2[3];
  double dt1dr1, dt1dr2, dt1ds1, dt1ds2;
  double dt2dr1, dt2dr2, dt2ds1, dt2ds2;
  double dt3dr1, dt3dr2, dt3ds1, dt3ds2;
  double drhodr1, drhodr2, drhods1, drhods2, drhodrm1[3], drhodrm2[3];
  double arg;
  double arg1i1, arg1j1, arg1i2, arg1j2, arg1i3, arg1j3, arg3i3, arg3j3;
  double dsij1, dsij2, force1, force2;
  double t1i, t2i, t3i, t1j, t2j, t3j;
  double scaleij;

  third = 1.0 / 3.0;
  sixth = 1.0 / 6.0;

  //     Compute forces atom i

  elti = fmap[type[i]];
  if (elti < 0) return;

  xitmp = x[i][0];
  yitmp = x[i][1];
  zitmp = x[i][2];

  //     Treat each pair
  for (jn = 0; jn < numneigh; jn++) {
    j = firstneigh[jn];
    eltj = fmap[type[j]];
    scaleij = scale[type[i]][type[j]];

    if (!iszero(scrfcn[fnoffset + jn]) && eltj >= 0) {

      sij = scrfcn[fnoffset + jn] * fcpair[fnoffset + jn];
      delij[0] = x[j][0] - xitmp;
      delij[1] = x[j][1] - yitmp;
      delij[2] = x[j][2] - zitmp;
      rij2 = delij[0] * delij[0] + delij[1] * delij[1] + delij[2] * delij[2];
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
        phi = ((this->phirar3[ind][kk] * pp + this->phirar2[ind][kk]) * pp + this->phirar1[ind][kk]) * pp +
          this->phirar[ind][kk];
        phip = (this->phirar6[ind][kk] * pp + this->phirar5[ind][kk]) * pp + this->phirar4[ind][kk];
        recip = 1.0 / r;

        if (eflag_either != 0) {
          double phi_sc = phi * scaleij;
          if (eflag_global != 0)
            *eng_vdwl = *eng_vdwl + phi_sc * sij;
          if (eflag_atom != 0) {
            eatom[i] = eatom[i] + 0.5 * phi_sc * sij;
            eatom[j] = eatom[j] + 0.5 * phi_sc * sij;
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

        const double t1mi = this->t1_meam[elti];
        const double t2mi = this->t2_meam[elti];
        const double t3mi = this->t3_meam[elti];
        const double t1mj = this->t1_meam[eltj];
        const double t2mj = this->t2_meam[eltj];
        const double t3mj = this->t3_meam[eltj];

        if (this->ialloy == 1) {
          rhoa1j  *= t1mj;
          rhoa2j  *= t2mj;
          rhoa3j  *= t3mj;
          rhoa1i  *= t1mi;
          rhoa2i  *= t2mi;
          rhoa3i  *= t3mi;
          drhoa1j *= t1mj;
          drhoa2j *= t2mj;
          drhoa3j *= t3mj;
          drhoa1i *= t1mi;
          drhoa2i *= t2mi;
          drhoa3i *= t3mi;
        }

        nv2 = 0;
        nv3 = 0;
        arg1i1 = 0.0;
        arg1j1 = 0.0;
        arg1i2 = 0.0;
        arg1j2 = 0.0;
        arg1i3 = 0.0;
        arg1j3 = 0.0;
        arg3i3 = 0.0;
        arg3j3 = 0.0;
        for (n = 0; n < 3; n++) {
          for (p = n; p < 3; p++) {
            for (q = p; q < 3; q++) {
              arg = delij[n] * delij[p] * delij[q] * this->v3D[nv3];
              arg1i3 = arg1i3 + arho3[i][nv3] * arg;
              arg1j3 = arg1j3 - arho3[j][nv3] * arg;
              nv3 = nv3 + 1;
            }
            arg = delij[n] * delij[p] * this->v2D[nv2];
            arg1i2 = arg1i2 + arho2[i][nv2] * arg;
            arg1j2 = arg1j2 + arho2[j][nv2] * arg;
            nv2 = nv2 + 1;
          }
          arg1i1 = arg1i1 + arho1[i][n] * delij[n];
          arg1j1 = arg1j1 - arho1[j][n] * delij[n];
          arg3i3 = arg3i3 + arho3b[i][n] * delij[n];
          arg3j3 = arg3j3 - arho3b[j][n] * delij[n];
        }

        //     rho0 terms
        drho0dr1 = drhoa0j * sij;
        drho0dr2 = drhoa0i * sij;

        //     rho1 terms
        a1 = 2 * sij / rij;
        drho1dr1 = a1 * (drhoa1j - rhoa1j / rij) * arg1i1;
        drho1dr2 = a1 * (drhoa1i - rhoa1i / rij) * arg1j1;
        a1 = 2.0 * sij / rij;
        for (m = 0; m < 3; m++) {
          drho1drm1[m] = a1 * rhoa1j * arho1[i][m];
          drho1drm2[m] = -a1 * rhoa1i * arho1[j][m];
        }

        //     rho2 terms
        a2 = 2 * sij / rij2;
        drho2dr1 = a2 * (drhoa2j - 2 * rhoa2j / rij) * arg1i2 - 2.0 / 3.0 * arho2b[i] * drhoa2j * sij;
        drho2dr2 = a2 * (drhoa2i - 2 * rhoa2i / rij) * arg1j2 - 2.0 / 3.0 * arho2b[j] * drhoa2i * sij;
        a2 = 4 * sij / rij2;
        for (m = 0; m < 3; m++) {
          drho2drm1[m] = 0.0;
          drho2drm2[m] = 0.0;
          for (n = 0; n < 3; n++) {
            drho2drm1[m] = drho2drm1[m] + arho2[i][this->vind2D[m][n]] * delij[n];
            drho2drm2[m] = drho2drm2[m] - arho2[j][this->vind2D[m][n]] * delij[n];
          }
          drho2drm1[m] = a2 * rhoa2j * drho2drm1[m];
          drho2drm2[m] = -a2 * rhoa2i * drho2drm2[m];
        }

        //     rho3 terms
        rij3 = rij * rij2;
        a3 = 2 * sij / rij3;
        a3a = 6.0 / 5.0 * sij / rij;
        drho3dr1 = a3 * (drhoa3j - 3 * rhoa3j / rij) * arg1i3 - a3a * (drhoa3j - rhoa3j / rij) * arg3i3;
        drho3dr2 = a3 * (drhoa3i - 3 * rhoa3i / rij) * arg1j3 - a3a * (drhoa3i - rhoa3i / rij) * arg3j3;
        a3 = 6 * sij / rij3;
        a3a = 6 * sij / (5 * rij);
        for (m = 0; m < 3; m++) {
          drho3drm1[m] = 0.0;
          drho3drm2[m] = 0.0;
          nv2 = 0;
          for (n = 0; n < 3; n++) {
            for (p = n; p < 3; p++) {
              arg = delij[n] * delij[p] * this->v2D[nv2];
              drho3drm1[m] = drho3drm1[m] + arho3[i][this->vind3D[m][n][p]] * arg;
              drho3drm2[m] = drho3drm2[m] + arho3[j][this->vind3D[m][n][p]] * arg;
              nv2 = nv2 + 1;
            }
          }
          drho3drm1[m] = (a3 * drho3drm1[m] - a3a * arho3b[i][m]) * rhoa3j;
          drho3drm2[m] = (-a3 * drho3drm2[m] + a3a * arho3b[j][m]) * rhoa3i;
        }

        //     Compute derivatives of weighting functions t wrt rij
        t1i = t_ave[i][0];
        t2i = t_ave[i][1];
        t3i = t_ave[i][2];
        t1j = t_ave[j][0];
        t2j = t_ave[j][1];
        t3j = t_ave[j][2];

        if (this->ialloy == 1) {

          a1i = fdiv_zero(drhoa0j * sij, tsq_ave[i][0]);
          a1j = fdiv_zero(drhoa0i * sij, tsq_ave[j][0]);
          a2i = fdiv_zero(drhoa0j * sij, tsq_ave[i][1]);
          a2j = fdiv_zero(drhoa0i * sij, tsq_ave[j][1]);
          a3i = fdiv_zero(drhoa0j * sij, tsq_ave[i][2]);
          a3j = fdiv_zero(drhoa0i * sij, tsq_ave[j][2]);

          dt1dr1 = a1i * (t1mj - t1i * MathSpecial::square(t1mj));
          dt1dr2 = a1j * (t1mi - t1j * MathSpecial::square(t1mi));
          dt2dr1 = a2i * (t2mj - t2i * MathSpecial::square(t2mj));
          dt2dr2 = a2j * (t2mi - t2j * MathSpecial::square(t2mi));
          dt3dr1 = a3i * (t3mj - t3i * MathSpecial::square(t3mj));
          dt3dr2 = a3j * (t3mi - t3j * MathSpecial::square(t3mi));

        } else if (this->ialloy == 2) {

          dt1dr1 = 0.0;
          dt1dr2 = 0.0;
          dt2dr1 = 0.0;
          dt2dr2 = 0.0;
          dt3dr1 = 0.0;
          dt3dr2 = 0.0;

        } else {

          ai = 0.0;
          if (!iszero(rho0[i]))
            ai = drhoa0j * sij / rho0[i];
          aj = 0.0;
          if (!iszero(rho0[j]))
            aj = drhoa0i * sij / rho0[j];

          dt1dr1 = ai * (t1mj - t1i);
          dt1dr2 = aj * (t1mi - t1j);
          dt2dr1 = ai * (t2mj - t2i);
          dt2dr2 = aj * (t2mi - t2j);
          dt3dr1 = ai * (t3mj - t3i);
          dt3dr2 = aj * (t3mi - t3j);
        }

        //     Compute derivatives of total density wrt rij, sij and rij(3)
        get_shpfcn(this->lattce_meam[elti][elti], this->stheta_meam[elti][elti], this->ctheta_meam[elti][elti], shpi);
        get_shpfcn(this->lattce_meam[eltj][eltj], this->stheta_meam[elti][elti], this->ctheta_meam[elti][elti], shpj);

        drhodr1 = dgamma1[i] * drho0dr1 +
          dgamma2[i] * (dt1dr1 * rho1[i] + t1i * drho1dr1 + dt2dr1 * rho2[i] + t2i * drho2dr1 +
                        dt3dr1 * rho3[i] + t3i * drho3dr1) -
          dgamma3[i] * (shpi[0] * dt1dr1 + shpi[1] * dt2dr1 + shpi[2] * dt3dr1);
        drhodr2 = dgamma1[j] * drho0dr2 +
          dgamma2[j] * (dt1dr2 * rho1[j] + t1j * drho1dr2 + dt2dr2 * rho2[j] + t2j * drho2dr2 +
                        dt3dr2 * rho3[j] + t3j * drho3dr2) -
          dgamma3[j] * (shpj[0] * dt1dr2 + shpj[1] * dt2dr2 + shpj[2] * dt3dr2);
        for (m = 0; m < 3; m++) {
          drhodrm1[m] = 0.0;
          drhodrm2[m] = 0.0;
          drhodrm1[m] = dgamma2[i] * (t1i * drho1drm1[m] + t2i * drho2drm1[m] + t3i * drho3drm1[m]);
          drhodrm2[m] = dgamma2[j] * (t1j * drho1drm2[m] + t2j * drho2drm2[m] + t3j * drho3drm2[m]);
        }

        //     Compute derivatives wrt sij, but only if necessary
        if (!iszero(dscrfcn[fnoffset + jn])) {
          drho0ds1 = rhoa0j;
          drho0ds2 = rhoa0i;
          a1 = 2.0 / rij;
          drho1ds1 = a1 * rhoa1j * arg1i1;
          drho1ds2 = a1 * rhoa1i * arg1j1;
          a2 = 2.0 / rij2;
          drho2ds1 = a2 * rhoa2j * arg1i2 - 2.0 / 3.0 * arho2b[i] * rhoa2j;
          drho2ds2 = a2 * rhoa2i * arg1j2 - 2.0 / 3.0 * arho2b[j] * rhoa2i;
          a3 = 2.0 / rij3;
          a3a = 6.0 / (5.0 * rij);
          drho3ds1 = a3 * rhoa3j * arg1i3 - a3a * rhoa3j * arg3i3;
          drho3ds2 = a3 * rhoa3i * arg1j3 - a3a * rhoa3i * arg3j3;

          if (this->ialloy == 1) {
            a1i = fdiv_zero(rhoa0j, tsq_ave[i][0]);
            a1j = fdiv_zero(rhoa0i, tsq_ave[j][0]);
            a2i = fdiv_zero(rhoa0j, tsq_ave[i][1]);
            a2j = fdiv_zero(rhoa0i, tsq_ave[j][1]);
            a3i = fdiv_zero(rhoa0j, tsq_ave[i][2]);
            a3j = fdiv_zero(rhoa0i, tsq_ave[j][2]);

            dt1ds1 = a1i * (t1mj - t1i * MathSpecial::square(t1mj));
            dt1ds2 = a1j * (t1mi - t1j * MathSpecial::square(t1mi));
            dt2ds1 = a2i * (t2mj - t2i * MathSpecial::square(t2mj));
            dt2ds2 = a2j * (t2mi - t2j * MathSpecial::square(t2mi));
            dt3ds1 = a3i * (t3mj - t3i * MathSpecial::square(t3mj));
            dt3ds2 = a3j * (t3mi - t3j * MathSpecial::square(t3mi));

          } else if (this->ialloy == 2) {

            dt1ds1 = 0.0;
            dt1ds2 = 0.0;
            dt2ds1 = 0.0;
            dt2ds2 = 0.0;
            dt3ds1 = 0.0;
            dt3ds2 = 0.0;

          } else {

            ai = 0.0;
            if (!iszero(rho0[i]))
              ai = rhoa0j / rho0[i];
            aj = 0.0;
            if (!iszero(rho0[j]))
              aj = rhoa0i / rho0[j];

            dt1ds1 = ai * (t1mj - t1i);
            dt1ds2 = aj * (t1mi - t1j);
            dt2ds1 = ai * (t2mj - t2i);
            dt2ds2 = aj * (t2mi - t2j);
            dt3ds1 = ai * (t3mj - t3i);
            dt3ds2 = aj * (t3mi - t3j);
          }

          drhods1 = dgamma1[i] * drho0ds1 +
            dgamma2[i] * (dt1ds1 * rho1[i] + t1i * drho1ds1 + dt2ds1 * rho2[i] + t2i * drho2ds1 +
                          dt3ds1 * rho3[i] + t3i * drho3ds1) -
            dgamma3[i] * (shpi[0] * dt1ds1 + shpi[1] * dt2ds1 + shpi[2] * dt3ds1);
          drhods2 = dgamma1[j] * drho0ds2 +
            dgamma2[j] * (dt1ds2 * rho1[j] + t1j * drho1ds2 + dt2ds2 * rho2[j] + t2j * drho2ds2 +
                          dt3ds2 * rho3[j] + t3j * drho3ds2) -
            dgamma3[j] * (shpj[0] * dt1ds2 + shpj[1] * dt2ds2 + shpj[2] * dt3ds2);
        }

        //     Compute derivatives of energy wrt rij, sij and rij[3]
        dUdrij = phip * sij + frhop[i] * drhodr1 + frhop[j] * drhodr2;
        dUdsij = 0.0;
        if (!iszero(dscrfcn[fnoffset + jn])) {
          dUdsij = phi + frhop[i] * drhods1 + frhop[j] * drhods2;
        }
        for (m = 0; m < 3; m++) {
          dUdrijm[m] = frhop[i] * drhodrm1[m] + frhop[j] * drhodrm2[m];
        }
        if (!isone(scaleij)) {
          dUdrij *= scaleij;
          dUdsij *= scaleij;
          dUdrijm[0] *= scaleij;
          dUdrijm[1] *= scaleij;
          dUdrijm[2] *= scaleij;
        }

        //     Add the part of the force due to dUdrij and dUdsij

        force = dUdrij * recip + dUdsij * dscrfcn[fnoffset + jn];
        for (m = 0; m < 3; m++) {
          forcem = delij[m] * force + dUdrijm[m];
          f[i][m] = f[i][m] + forcem;
          f[j][m] = f[j][m] - forcem;
        }

        //     Tabulate per-atom virial as symmetrized stress tensor

        if (vflag_atom != 0) {
          fi[0] = delij[0] * force + dUdrijm[0];
          fi[1] = delij[1] * force + dUdrijm[1];
          fi[2] = delij[2] * force + dUdrijm[2];
          v[0] = -0.5 * (delij[0] * fi[0]);
          v[1] = -0.5 * (delij[1] * fi[1]);
          v[2] = -0.5 * (delij[2] * fi[2]);
          v[3] = -0.25 * (delij[0] * fi[1] + delij[1] * fi[0]);
          v[4] = -0.25 * (delij[0] * fi[2] + delij[2] * fi[0]);
          v[5] = -0.25 * (delij[1] * fi[2] + delij[2] * fi[1]);

          for (m = 0; m < 6; m++) {
            vatom[i][m] = vatom[i][m] + v[m];
            vatom[j][m] = vatom[j][m] + v[m];
          }
        }

        //     Now compute forces on other atoms k due to change in sij

        if (iszero(sij) || isone(sij)) continue; //: cont jn loop

        double dxik(0), dyik(0), dzik(0);
        double dxjk(0), dyjk(0), dzjk(0);

        for (kn = 0; kn < numneigh_full; kn++) {
          k = firstneigh_full[kn];
          eltk = fmap[type[k]];
          if (k != j && eltk >= 0) {
            double xik, xjk, cikj, sikj, dfc, a;
            double dCikj1, dCikj2;
            double delc, rik2, rjk2;

            sij = scrfcn[jn+fnoffset] * fcpair[jn+fnoffset];
            const double Cmax = this->Cmax_meam[elti][eltj][eltk];
            const double Cmin = this->Cmin_meam[elti][eltj][eltk];

            dsij1 = 0.0;
            dsij2 = 0.0;
            if (!iszero(sij) && !isone(sij)) {
              const double rbound = rij2 * this->ebound_meam[elti][eltj];
              delc = Cmax - Cmin;
              dxjk = x[k][0] - x[j][0];
              dyjk = x[k][1] - x[j][1];
              dzjk = x[k][2] - x[j][2];
              rjk2 = dxjk * dxjk + dyjk * dyjk + dzjk * dzjk;
              if (rjk2 <= rbound) {
                dxik = x[k][0] - x[i][0];
                dyik = x[k][1] - x[i][1];
                dzik = x[k][2] - x[i][2];
                rik2 = dxik * dxik + dyik * dyik + dzik * dzik;
                if (rik2 <= rbound) {
                  xik = rik2 / rij2;
                  xjk = rjk2 / rij2;
                  a = 1 - (xik - xjk) * (xik - xjk);
                  if (!iszero(a)) {
                    cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
                    if (cikj >= Cmin && cikj <= Cmax) {
                      cikj = (cikj - Cmin) / delc;
                      sikj = dfcut(cikj, dfc);
                      dCfunc2(rij2, rik2, rjk2, dCikj1, dCikj2);
                      a = sij / delc * dfc / sikj;
                      dsij1 = a * dCikj1;
                      dsij2 = a * dCikj2;
                    }
                  }
                }
              }
            }

            if (!iszero(dsij1) || !iszero(dsij2)) {
              force1 = dUdsij * dsij1;
              force2 = dUdsij * dsij2;

              f[i][0] += force1 * dxik;
              f[i][1] += force1 * dyik;
              f[i][2] += force1 * dzik;
              f[j][0] += force2 * dxjk;
              f[j][1] += force2 * dyjk;
              f[j][2] += force2 * dzjk;
              f[k][0] -= force1 * dxik + force2 * dxjk;
              f[k][1] -= force1 * dyik + force2 * dyjk;
              f[k][2] -= force1 * dzik + force2 * dzjk;

              //     Tabulate per-atom virial as symmetrized stress tensor

              if (vflag_atom != 0) {
                fi[0] = force1 * dxik;
                fi[1] = force1 * dyik;
                fi[2] = force1 * dzik;
                fj[0] = force2 * dxjk;
                fj[1] = force2 * dyjk;
                fj[2] = force2 * dzjk;
                v[0] = -third * (dxik * fi[0] + dxjk * fj[0]);
                v[1] = -third * (dyik * fi[1] + dyjk * fj[1]);
                v[2] = -third * (dzik * fi[2] + dzjk * fj[2]);
                v[3] = -sixth * (dxik * fi[1] + dxjk * fj[1] + dyik * fi[0] + dyjk * fj[0]);
                v[4] = -sixth * (dxik * fi[2] + dxjk * fj[2] + dzik * fi[0] + dzjk * fj[0]);
                v[5] = -sixth * (dyik * fi[2] + dyjk * fj[2] + dzik * fi[1] + dzjk * fj[1]);

                for (m = 0; m < 6; m++) {
                  vatom[i][m] = vatom[i][m] + v[m];
                  vatom[j][m] = vatom[j][m] + v[m];
                  vatom[k][m] = vatom[k][m] + v[m];
                }
              }
            }
          }
          //     end of k loop
        }
      }
    }
    //     end of j loop
  }
}
