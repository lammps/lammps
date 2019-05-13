#include "meam.h"
#include "math_special.h"

using namespace LAMMPS_NS;

void
MEAM::meam_dens_setup(int atom_nmax, int nall, int n_neigh)
{
  int i, j;

  // grow local arrays if necessary

  if (atom_nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(rho0);
    memory->destroy(rho1);
    memory->destroy(rho2);
    memory->destroy(rho3);
    memory->destroy(frhop);
    memory->destroy(gamma);
    memory->destroy(dgamma1);
    memory->destroy(dgamma2);
    memory->destroy(dgamma3);
    memory->destroy(arho2b);
    memory->destroy(arho1);
    memory->destroy(arho2);
    memory->destroy(arho3);
    memory->destroy(arho3b);
    memory->destroy(t_ave);
    memory->destroy(tsq_ave);

    nmax = atom_nmax;

    memory->create(rho, nmax, "pair:rho");
    memory->create(rho0, nmax, "pair:rho0");
    memory->create(rho1, nmax, "pair:rho1");
    memory->create(rho2, nmax, "pair:rho2");
    memory->create(rho3, nmax, "pair:rho3");
    memory->create(frhop, nmax, "pair:frhop");
    memory->create(gamma, nmax, "pair:gamma");
    memory->create(dgamma1, nmax, "pair:dgamma1");
    memory->create(dgamma2, nmax, "pair:dgamma2");
    memory->create(dgamma3, nmax, "pair:dgamma3");
    memory->create(arho2b, nmax, "pair:arho2b");
    memory->create(arho1, nmax, 3, "pair:arho1");
    memory->create(arho2, nmax, 6, "pair:arho2");
    memory->create(arho3, nmax, 10, "pair:arho3");
    memory->create(arho3b, nmax, 3, "pair:arho3b");
    memory->create(t_ave, nmax, 3, "pair:t_ave");
    memory->create(tsq_ave, nmax, 3, "pair:tsq_ave");
  }

  if (n_neigh > maxneigh) {
    memory->destroy(scrfcn);
    memory->destroy(dscrfcn);
    memory->destroy(fcpair);
    maxneigh = n_neigh;
    memory->create(scrfcn, maxneigh, "pair:scrfcn");
    memory->create(dscrfcn, maxneigh, "pair:dscrfcn");
    memory->create(fcpair, maxneigh, "pair:fcpair");
  }

  // zero out local arrays

  for (i = 0; i < nall; i++) {
    rho0[i] = 0.0;
    arho2b[i] = 0.0;
    arho1[i][0] = arho1[i][1] = arho1[i][2] = 0.0;
    for (j = 0; j < 6; j++)
      arho2[i][j] = 0.0;
    for (j = 0; j < 10; j++)
      arho3[i][j] = 0.0;
    arho3b[i][0] = arho3b[i][1] = arho3b[i][2] = 0.0;
    t_ave[i][0] = t_ave[i][1] = t_ave[i][2] = 0.0;
    tsq_ave[i][0] = tsq_ave[i][1] = tsq_ave[i][2] = 0.0;
  }
}

void
MEAM::meam_dens_init(int i, int ntype, int* type, int* fmap, double** x,
                     int numneigh, int* firstneigh,
                     int numneigh_full, int* firstneigh_full, int fnoffset)
{
  //     Compute screening function and derivatives
  getscreen(i, &scrfcn[fnoffset], &dscrfcn[fnoffset], &fcpair[fnoffset], x, numneigh, firstneigh,
            numneigh_full, firstneigh_full, ntype, type, fmap);

  //     Calculate intermediate density terms to be communicated
  calc_rho1(i, ntype, type, fmap, x, numneigh, firstneigh, &scrfcn[fnoffset], &fcpair[fnoffset]);
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

void
MEAM::getscreen(int i, double* scrfcn, double* dscrfcn, double* fcpair, double** x, int numneigh,
                int* firstneigh, int numneigh_full, int* firstneigh_full, int /*ntype*/, int* type, int* fmap)
{
  int jn, j, kn, k;
  int elti, eltj, eltk;
  double xitmp, yitmp, zitmp, delxij, delyij, delzij, rij2, rij;
  double xjtmp, yjtmp, zjtmp, delxik, delyik, delzik, rik2 /*,rik*/;
  double xktmp, yktmp, zktmp, delxjk, delyjk, delzjk, rjk2 /*,rjk*/;
  double xik, xjk, sij, fcij, sfcij, dfcij, sikj, dfikj, cikj;
  double Cmin, Cmax, delc, /*ebound,*/ a, coef1, coef2;
  double dCikj;
  double rnorm, fc, dfc, drinv;

  drinv = 1.0 / this->delr_meam;
  elti = fmap[type[i]];
  if (elti < 0) return;

  xitmp = x[i][0];
  yitmp = x[i][1];
  zitmp = x[i][2];

  for (jn = 0; jn < numneigh; jn++) {
    j = firstneigh[jn];

    eltj = fmap[type[j]];
    if (eltj < 0) continue;

    //     First compute screening function itself, sij
    xjtmp = x[j][0];
    yjtmp = x[j][1];
    zjtmp = x[j][2];
    delxij = xjtmp - xitmp;
    delyij = yjtmp - yitmp;
    delzij = zjtmp - zitmp;
    rij2 = delxij * delxij + delyij * delyij + delzij * delzij;

    if (rij2 > this->cutforcesq) {
      dscrfcn[jn] = 0.0;
      scrfcn[jn] = 0.0;
      fcpair[jn] = 0.0;
      continue;
    }

    const double rbound = this->ebound_meam[elti][eltj] * rij2;
    rij = sqrt(rij2);
    rnorm = (this->cutforce - rij) * drinv;
    sij = 1.0;

    //     if rjk2 > ebound*rijsq, atom k is definitely outside the ellipse
    for (kn = 0; kn < numneigh_full; kn++) {
      k = firstneigh_full[kn];
      if (k == j) continue;
      eltk = fmap[type[k]];
      if (eltk < 0) continue;

      xktmp = x[k][0];
      yktmp = x[k][1];
      zktmp = x[k][2];

      delxjk = xktmp - xjtmp;
      delyjk = yktmp - yjtmp;
      delzjk = zktmp - zjtmp;
      rjk2 = delxjk * delxjk + delyjk * delyjk + delzjk * delzjk;
      if (rjk2 > rbound) continue;

      delxik = xktmp - xitmp;
      delyik = yktmp - yitmp;
      delzik = zktmp - zitmp;
      rik2 = delxik * delxik + delyik * delyik + delzik * delzik;
      if (rik2 > rbound) continue;

      xik = rik2 / rij2;
      xjk = rjk2 / rij2;
      a = 1 - (xik - xjk) * (xik - xjk);
      //     if a < 0, then ellipse equation doesn't describe this case and
      //     atom k can't possibly screen i-j
      if (a <= 0.0) continue;

      cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
      Cmax = this->Cmax_meam[elti][eltj][eltk];
      Cmin = this->Cmin_meam[elti][eltj][eltk];
      if (cikj >= Cmax) continue;
      //     note that cikj may be slightly negative (within numerical
      //     tolerance) if atoms are colinear, so don't reject that case here
      //     (other negative cikj cases were handled by the test on "a" above)
      else if (cikj <= Cmin) {
        sij = 0.0;
        break;
      } else {
        delc = Cmax - Cmin;
        cikj = (cikj - Cmin) / delc;
        sikj = fcut(cikj);
      }
      sij *= sikj;
    }

    fc = dfcut(rnorm, dfc);
    fcij = fc;
    dfcij = dfc * drinv;

    //     Now compute derivatives
    dscrfcn[jn] = 0.0;
    sfcij = sij * fcij;
    if (!iszero(sfcij) && !iszero(sfcij - 1.0)) {
      for (kn = 0; kn < numneigh_full; kn++) {
        k = firstneigh_full[kn];
        if (k == j) continue;
        eltk = fmap[type[k]];
        if (eltk < 0) continue;

        delxjk = x[k][0] - xjtmp;
        delyjk = x[k][1] - yjtmp;
        delzjk = x[k][2] - zjtmp;
        rjk2 = delxjk * delxjk + delyjk * delyjk + delzjk * delzjk;
        if (rjk2 > rbound) continue;

        delxik = x[k][0] - xitmp;
        delyik = x[k][1] - yitmp;
        delzik = x[k][2] - zitmp;
        rik2 = delxik * delxik + delyik * delyik + delzik * delzik;
        if (rik2 > rbound) continue;

        xik = rik2 / rij2;
        xjk = rjk2 / rij2;
        a = 1 - (xik - xjk) * (xik - xjk);
        //     if a < 0, then ellipse equation doesn't describe this case and
        //     atom k can't possibly screen i-j
        if (a <= 0.0) continue;

        cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
        Cmax = this->Cmax_meam[elti][eltj][eltk];
        Cmin = this->Cmin_meam[elti][eltj][eltk];
        if (cikj >= Cmax) {
          continue;
          //     Note that cikj may be slightly negative (within numerical
          //     tolerance) if atoms are colinear, so don't reject that case
          //     here
          //     (other negative cikj cases were handled by the test on "a"
          //     above)
          //     Note that we never have 0<cikj<Cmin here, else sij=0
          //     (rejected above)
        } else {
          delc = Cmax - Cmin;
          cikj = (cikj - Cmin) / delc;
          sikj = dfcut(cikj, dfikj);
          coef1 = dfikj / (delc * sikj);
          dCikj = dCfunc(rij2, rik2, rjk2);
          dscrfcn[jn] = dscrfcn[jn] + coef1 * dCikj;
        }
      }
      coef1 = sfcij;
      coef2 = sij * dfcij / rij;
      dscrfcn[jn] = dscrfcn[jn] * coef1 - coef2;
    }

    scrfcn[jn] = sij;
    fcpair[jn] = fcij;
  }
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

void
MEAM::calc_rho1(int i, int /*ntype*/, int* type, int* fmap, double** x, int numneigh, int* firstneigh,
                double* scrfcn, double* fcpair)
{
  int jn, j, m, n, p, elti, eltj;
  int nv2, nv3;
  double xtmp, ytmp, ztmp, delij[3], rij2, rij, sij;
  double ai, aj, rhoa0j, rhoa1j, rhoa2j, rhoa3j, A1j, A2j, A3j;
  // double G,Gbar,gam,shp[3+1];
  double ro0i, ro0j;
  double rhoa0i, rhoa1i, rhoa2i, rhoa3i, A1i, A2i, A3i;

  elti = fmap[type[i]];
  xtmp = x[i][0];
  ytmp = x[i][1];
  ztmp = x[i][2];
  for (jn = 0; jn < numneigh; jn++) {
    if (!iszero(scrfcn[jn])) {
      j = firstneigh[jn];
      sij = scrfcn[jn] * fcpair[jn];
      delij[0] = x[j][0] - xtmp;
      delij[1] = x[j][1] - ytmp;
      delij[2] = x[j][2] - ztmp;
      rij2 = delij[0] * delij[0] + delij[1] * delij[1] + delij[2] * delij[2];
      if (rij2 < this->cutforcesq) {
        eltj = fmap[type[j]];
        rij = sqrt(rij2);
        ai = rij / this->re_meam[elti][elti] - 1.0;
        aj = rij / this->re_meam[eltj][eltj] - 1.0;
        ro0i = this->rho0_meam[elti];
        ro0j = this->rho0_meam[eltj];
        rhoa0j = ro0j * MathSpecial::fm_exp(-this->beta0_meam[eltj] * aj) * sij;
        rhoa1j = ro0j * MathSpecial::fm_exp(-this->beta1_meam[eltj] * aj) * sij;
        rhoa2j = ro0j * MathSpecial::fm_exp(-this->beta2_meam[eltj] * aj) * sij;
        rhoa3j = ro0j * MathSpecial::fm_exp(-this->beta3_meam[eltj] * aj) * sij;
        rhoa0i = ro0i * MathSpecial::fm_exp(-this->beta0_meam[elti] * ai) * sij;
        rhoa1i = ro0i * MathSpecial::fm_exp(-this->beta1_meam[elti] * ai) * sij;
        rhoa2i = ro0i * MathSpecial::fm_exp(-this->beta2_meam[elti] * ai) * sij;
        rhoa3i = ro0i * MathSpecial::fm_exp(-this->beta3_meam[elti] * ai) * sij;
        if (this->ialloy == 1) {
          rhoa1j = rhoa1j * this->t1_meam[eltj];
          rhoa2j = rhoa2j * this->t2_meam[eltj];
          rhoa3j = rhoa3j * this->t3_meam[eltj];
          rhoa1i = rhoa1i * this->t1_meam[elti];
          rhoa2i = rhoa2i * this->t2_meam[elti];
          rhoa3i = rhoa3i * this->t3_meam[elti];
        }
        rho0[i] = rho0[i] + rhoa0j;
        rho0[j] = rho0[j] + rhoa0i;
        // For ialloy = 2, use single-element value (not average)
        if (this->ialloy != 2) {
          t_ave[i][0] = t_ave[i][0] + this->t1_meam[eltj] * rhoa0j;
          t_ave[i][1] = t_ave[i][1] + this->t2_meam[eltj] * rhoa0j;
          t_ave[i][2] = t_ave[i][2] + this->t3_meam[eltj] * rhoa0j;
          t_ave[j][0] = t_ave[j][0] + this->t1_meam[elti] * rhoa0i;
          t_ave[j][1] = t_ave[j][1] + this->t2_meam[elti] * rhoa0i;
          t_ave[j][2] = t_ave[j][2] + this->t3_meam[elti] * rhoa0i;
        }
        if (this->ialloy == 1) {
          tsq_ave[i][0] = tsq_ave[i][0] + this->t1_meam[eltj] * this->t1_meam[eltj] * rhoa0j;
          tsq_ave[i][1] = tsq_ave[i][1] + this->t2_meam[eltj] * this->t2_meam[eltj] * rhoa0j;
          tsq_ave[i][2] = tsq_ave[i][2] + this->t3_meam[eltj] * this->t3_meam[eltj] * rhoa0j;
          tsq_ave[j][0] = tsq_ave[j][0] + this->t1_meam[elti] * this->t1_meam[elti] * rhoa0i;
          tsq_ave[j][1] = tsq_ave[j][1] + this->t2_meam[elti] * this->t2_meam[elti] * rhoa0i;
          tsq_ave[j][2] = tsq_ave[j][2] + this->t3_meam[elti] * this->t3_meam[elti] * rhoa0i;
        }
        arho2b[i] = arho2b[i] + rhoa2j;
        arho2b[j] = arho2b[j] + rhoa2i;

        A1j = rhoa1j / rij;
        A2j = rhoa2j / rij2;
        A3j = rhoa3j / (rij2 * rij);
        A1i = rhoa1i / rij;
        A2i = rhoa2i / rij2;
        A3i = rhoa3i / (rij2 * rij);
        nv2 = 0;
        nv3 = 0;
        for (m = 0; m < 3; m++) {
          arho1[i][m] = arho1[i][m] + A1j * delij[m];
          arho1[j][m] = arho1[j][m] - A1i * delij[m];
          arho3b[i][m] = arho3b[i][m] + rhoa3j * delij[m] / rij;
          arho3b[j][m] = arho3b[j][m] - rhoa3i * delij[m] / rij;
          for (n = m; n < 3; n++) {
            arho2[i][nv2] = arho2[i][nv2] + A2j * delij[m] * delij[n];
            arho2[j][nv2] = arho2[j][nv2] + A2i * delij[m] * delij[n];
            nv2 = nv2 + 1;
            for (p = n; p < 3; p++) {
              arho3[i][nv3] = arho3[i][nv3] + A3j * delij[m] * delij[n] * delij[p];
              arho3[j][nv3] = arho3[j][nv3] - A3i * delij[m] * delij[n] * delij[p];
              nv3 = nv3 + 1;
            }
          }
        }
      }
    }
  }
}

