#include "meam.h"
#include <math.h>
#include "math_special.h"

using namespace LAMMPS_NS;
//     Extern "C" declaration has the form:
//
//  void meam_dens_init_(int *, int *, int *, double *, int *, int *, int *,
//  double *,
//		 int *, int *, int *, int *,
//		 double *, double *, double *, double *, double *, double *,
//		 double *, double *, double *, double *, double *, int *);
//
//
//     Call from pair_meam.cpp has the form:
//
//    meam_dens_init_(&i,&nmax,ntype,type,fmap,&x[0][0],
//	       &numneigh[i],firstneigh[i],&numneigh_full[i],firstneigh_full[i],
//	       &scrfcn[offset],&dscrfcn[offset],&fcpair[offset],
//	       rho0,&arho1[0][0],&arho2[0][0],arho2b,
//	       &arho3[0][0],&arho3b[0][0],&t_ave[0][0],&tsq_ave[0][0],&errorflag);
//

void
MEAM::meam_dens_init(int* i, int* nmax, int* ntype, int* type, int* fmap, double* x,
                int* numneigh, int* firstneigh, int* numneigh_full,
                int* firstneigh_full, double* scrfcn, double* dscrfcn,
                double* fcpair, double* rho0, double* arho1, double* arho2,
                double* arho2b, double* arho3, double* arho3b, double* t_ave,
                double* tsq_ave, int* errorflag)
{
  *errorflag = 0;

  //     Compute screening function and derivatives
  getscreen(*i, *nmax, scrfcn, dscrfcn, fcpair, x, *numneigh, firstneigh,
            *numneigh_full, firstneigh_full, *ntype, type, fmap);

  //     Calculate intermediate density terms to be communicated
  calc_rho1(*i, *nmax, *ntype, type, fmap, x, *numneigh, firstneigh, scrfcn,
            fcpair, rho0, arho1, arho2, arho2b, arho3, arho3b, t_ave, tsq_ave);
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

void
MEAM::getscreen(int i, int nmax, double* scrfcn, double* dscrfcn, double* fcpair,
          double* x, int numneigh, int* firstneigh, int numneigh_full,
          int* firstneigh_full, int ntype, int* type, int* fmap)
{

  arrdim2v(x, 3, nmax);

  int jn, j, kn, k;
  int elti, eltj, eltk;
  double xitmp, yitmp, zitmp, delxij, delyij, delzij, rij2, rij;
  double xjtmp, yjtmp, zjtmp, delxik, delyik, delzik, rik2 /*,rik*/;
  double xktmp, yktmp, zktmp, delxjk, delyjk, delzjk, rjk2 /*,rjk*/;
  double xik, xjk, sij, fcij, sfcij, dfcij, sikj, dfikj, cikj;
  double Cmin, Cmax, delc, /*ebound,*/ rbound, a, coef1, coef2;
  // double coef1a,coef1b,coef2a,coef2b;
  double dCikj;
  /*unused:double dC1a,dC1b,dC2a,dC2b;*/
  double rnorm, fc, dfc, drinv;

  drinv = 1.0 / this->delr_meam;
  elti = arr1v(fmap, arr1v(type, i));

  if (elti > 0) {

    xitmp = arr2v(x, 1, i);
    yitmp = arr2v(x, 2, i);
    zitmp = arr2v(x, 3, i);

    for (jn = 1; jn <= numneigh; jn++) {
      j = arr1v(firstneigh, jn);

      eltj = arr1v(fmap, arr1v(type, j));
      if (eltj > 0) {

        //     First compute screening function itself, sij
        xjtmp = arr2v(x, 1, j);
        yjtmp = arr2v(x, 2, j);
        zjtmp = arr2v(x, 3, j);
        delxij = xjtmp - xitmp;
        delyij = yjtmp - yitmp;
        delzij = zjtmp - zitmp;
        rij2 = delxij * delxij + delyij * delyij + delzij * delzij;
        rij = sqrt(rij2);
        if (rij > this->rc_meam) {
          fcij = 0.0;
          dfcij = 0.0;
          sij = 0.0;
        } else {
          rnorm = (this->rc_meam - rij) * drinv;
          screen(i, j, nmax, x, rij2, &sij, numneigh_full, firstneigh_full,
                 ntype, type, fmap);
          dfcut(rnorm, &fc, &dfc);
          fcij = fc;
          dfcij = dfc * drinv;
        }

        //     Now compute derivatives
        arr1v(dscrfcn, jn) = 0.0;
        sfcij = sij * fcij;
        if (iszero(sfcij) || iszero(sfcij - 1.0))
          goto LABEL_100;
        rbound = this->ebound_meam[elti][eltj] * rij2;
        for (kn = 1; kn <= numneigh_full; kn++) {
          k = arr1v(firstneigh_full, kn);
          if (k == j)
            continue;
          eltk = arr1v(fmap, arr1v(type, k));
          if (eltk == 0)
            continue;
          xktmp = arr2v(x, 1, k);
          yktmp = arr2v(x, 2, k);
          zktmp = arr2v(x, 3, k);
          delxjk = xktmp - xjtmp;
          delyjk = yktmp - yjtmp;
          delzjk = zktmp - zjtmp;
          rjk2 = delxjk * delxjk + delyjk * delyjk + delzjk * delzjk;
          if (rjk2 > rbound)
            continue;
          delxik = xktmp - xitmp;
          delyik = yktmp - yitmp;
          delzik = zktmp - zitmp;
          rik2 = delxik * delxik + delyik * delyik + delzik * delzik;
          if (rik2 > rbound)
            continue;
          xik = rik2 / rij2;
          xjk = rjk2 / rij2;
          a = 1 - (xik - xjk) * (xik - xjk);
          //     if a < 0, then ellipse equation doesn't describe this case and
          //     atom k can't possibly screen i-j
          if (a <= 0.0)
            continue;
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
            dfcut(cikj, &sikj, &dfikj);
            coef1 = dfikj / (delc * sikj);
            dCfunc(rij2, rik2, rjk2, &dCikj);
            arr1v(dscrfcn, jn) = arr1v(dscrfcn, jn) + coef1 * dCikj;
          }
        }
        coef1 = sfcij;
        coef2 = sij * dfcij / rij;
        arr1v(dscrfcn, jn) = arr1v(dscrfcn, jn) * coef1 - coef2;

      LABEL_100:
        arr1v(scrfcn, jn) = sij;
        arr1v(fcpair, jn) = fcij;
      }
    }
  }
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

void
MEAM::calc_rho1(int i, int nmax, int ntype, int* type, int* fmap, double* x,
          int numneigh, int* firstneigh, double* scrfcn, double* fcpair,
          double* rho0, double* arho1, double* arho2, double* arho2b,
          double* arho3, double* arho3b, double* t_ave, double* tsq_ave)
{
  arrdim2v(x, 3, nmax);
  arrdim2v(arho1, 3, nmax);
  arrdim2v(arho2, 6, nmax);
  arrdim2v(arho3, 10, nmax);
  arrdim2v(arho3b, 3, nmax);
  arrdim2v(t_ave, 3, nmax);
  arrdim2v(tsq_ave, 3, nmax);

  int jn, j, m, n, p, elti, eltj;
  int nv2, nv3;
  double xtmp, ytmp, ztmp, delij[3 + 1], rij2, rij, sij;
  double ai, aj, rhoa0j, rhoa1j, rhoa2j, rhoa3j, A1j, A2j, A3j;
  // double G,Gbar,gam,shp[3+1];
  double ro0i, ro0j;
  double rhoa0i, rhoa1i, rhoa2i, rhoa3i, A1i, A2i, A3i;

  elti = arr1v(fmap, arr1v(type, i));
  xtmp = arr2v(x, 1, i);
  ytmp = arr2v(x, 2, i);
  ztmp = arr2v(x, 3, i);
  for (jn = 1; jn <= numneigh; jn++) {
    if (!iszero(arr1v(scrfcn, jn))) {
      j = arr1v(firstneigh, jn);
      sij = arr1v(scrfcn, jn) * arr1v(fcpair, jn);
      delij[1] = arr2v(x, 1, j) - xtmp;
      delij[2] = arr2v(x, 2, j) - ytmp;
      delij[3] = arr2v(x, 3, j) - ztmp;
      rij2 = delij[1] * delij[1] + delij[2] * delij[2] + delij[3] * delij[3];
      if (rij2 < this->cutforcesq) {
        eltj = arr1v(fmap, arr1v(type, j));
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
        arr1v(rho0, i) = arr1v(rho0, i) + rhoa0j;
        arr1v(rho0, j) = arr1v(rho0, j) + rhoa0i;
        // For ialloy = 2, use single-element value (not average)
        if (this->ialloy != 2) {
          arr2v(t_ave, 1, i) =
            arr2v(t_ave, 1, i) + this->t1_meam[eltj] * rhoa0j;
          arr2v(t_ave, 2, i) =
            arr2v(t_ave, 2, i) + this->t2_meam[eltj] * rhoa0j;
          arr2v(t_ave, 3, i) =
            arr2v(t_ave, 3, i) + this->t3_meam[eltj] * rhoa0j;
          arr2v(t_ave, 1, j) =
            arr2v(t_ave, 1, j) + this->t1_meam[elti] * rhoa0i;
          arr2v(t_ave, 2, j) =
            arr2v(t_ave, 2, j) + this->t2_meam[elti] * rhoa0i;
          arr2v(t_ave, 3, j) =
            arr2v(t_ave, 3, j) + this->t3_meam[elti] * rhoa0i;
        }
        if (this->ialloy == 1) {
          arr2v(tsq_ave, 1, i) =
            arr2v(tsq_ave, 1, i) +
            this->t1_meam[eltj] * this->t1_meam[eltj] * rhoa0j;
          arr2v(tsq_ave, 2, i) =
            arr2v(tsq_ave, 2, i) +
            this->t2_meam[eltj] * this->t2_meam[eltj] * rhoa0j;
          arr2v(tsq_ave, 3, i) =
            arr2v(tsq_ave, 3, i) +
            this->t3_meam[eltj] * this->t3_meam[eltj] * rhoa0j;
          arr2v(tsq_ave, 1, j) =
            arr2v(tsq_ave, 1, j) +
            this->t1_meam[elti] * this->t1_meam[elti] * rhoa0i;
          arr2v(tsq_ave, 2, j) =
            arr2v(tsq_ave, 2, j) +
            this->t2_meam[elti] * this->t2_meam[elti] * rhoa0i;
          arr2v(tsq_ave, 3, j) =
            arr2v(tsq_ave, 3, j) +
            this->t3_meam[elti] * this->t3_meam[elti] * rhoa0i;
        }
        arr1v(arho2b, i) = arr1v(arho2b, i) + rhoa2j;
        arr1v(arho2b, j) = arr1v(arho2b, j) + rhoa2i;

        A1j = rhoa1j / rij;
        A2j = rhoa2j / rij2;
        A3j = rhoa3j / (rij2 * rij);
        A1i = rhoa1i / rij;
        A2i = rhoa2i / rij2;
        A3i = rhoa3i / (rij2 * rij);
        nv2 = 1;
        nv3 = 1;
        for (m = 1; m <= 3; m++) {
          arr2v(arho1, m, i) = arr2v(arho1, m, i) + A1j * delij[m];
          arr2v(arho1, m, j) = arr2v(arho1, m, j) - A1i * delij[m];
          arr2v(arho3b, m, i) = arr2v(arho3b, m, i) + rhoa3j * delij[m] / rij;
          arr2v(arho3b, m, j) = arr2v(arho3b, m, j) - rhoa3i * delij[m] / rij;
          for (n = m; n <= 3; n++) {
            arr2v(arho2, nv2, i) =
              arr2v(arho2, nv2, i) + A2j * delij[m] * delij[n];
            arr2v(arho2, nv2, j) =
              arr2v(arho2, nv2, j) + A2i * delij[m] * delij[n];
            nv2 = nv2 + 1;
            for (p = n; p <= 3; p++) {
              arr2v(arho3, nv3, i) =
                arr2v(arho3, nv3, i) + A3j * delij[m] * delij[n] * delij[p];
              arr2v(arho3, nv3, j) =
                arr2v(arho3, nv3, j) - A3i * delij[m] * delij[n] * delij[p];
              nv3 = nv3 + 1;
            }
          }
        }
      }
    }
  }
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

void
MEAM::screen(int i, int j, int nmax, double* x, double rijsq, double* sij,
       int numneigh_full, int* firstneigh_full, int ntype, int* type, int* fmap)
//     Screening function
//     Inputs:  i = atom 1 id (integer)
//     j = atom 2 id (integer)
//     rijsq = squared distance between i and j
//     Outputs: sij = screening function
{

  arrdim2v(x, 3, nmax)

    int k,
    nk /*,m*/;
  int elti, eltj, eltk;
  double delxik, delyik, delzik;
  double delxjk, delyjk, delzjk;
  double riksq, rjksq, xik, xjk, cikj, a, delc, sikj /*,fcij,rij*/;
  double Cmax, Cmin, rbound;

  *sij = 1.0;
  eltj = arr1v(fmap, arr1v(type, j));
  elti = arr1v(fmap, arr1v(type, j));

  //     if rjksq > ebound*rijsq, atom k is definitely outside the ellipse
  rbound = this->ebound_meam[elti][eltj] * rijsq;

  for (nk = 1; nk <= numneigh_full; nk++) {
    k = arr1v(firstneigh_full, nk);
    eltk = arr1v(fmap, arr1v(type, k));
    if (k == j)
      continue;
    delxjk = arr2v(x, 1, k) - arr2v(x, 1, j);
    delyjk = arr2v(x, 2, k) - arr2v(x, 2, j);
    delzjk = arr2v(x, 3, k) - arr2v(x, 3, j);
    rjksq = delxjk * delxjk + delyjk * delyjk + delzjk * delzjk;
    if (rjksq > rbound)
      continue;
    delxik = arr2v(x, 1, k) - arr2v(x, 1, i);
    delyik = arr2v(x, 2, k) - arr2v(x, 2, i);
    delzik = arr2v(x, 3, k) - arr2v(x, 3, i);
    riksq = delxik * delxik + delyik * delyik + delzik * delzik;
    if (riksq > rbound)
      continue;
    xik = riksq / rijsq;
    xjk = rjksq / rijsq;
    a = 1 - (xik - xjk) * (xik - xjk);
    //     if a < 0, then ellipse equation doesn't describe this case and
    //     atom k can't possibly screen i-j
    if (a <= 0.0)
      continue;
    cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
    Cmax = this->Cmax_meam[elti][eltj][eltk];
    Cmin = this->Cmin_meam[elti][eltj][eltk];
    if (cikj >= Cmax)
      continue;
    //     note that cikj may be slightly negative (within numerical
    //     tolerance) if atoms are colinear, so don't reject that case here
    //     (other negative cikj cases were handled by the test on "a" above)
    else if (cikj <= Cmin) {
      *sij = 0.0;
      break;
    } else {
      delc = Cmax - Cmin;
      cikj = (cikj - Cmin) / delc;
      fcut(cikj, &sikj);
    }
    *sij = *sij * sikj;
  }
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

void
MEAM::dsij(int i, int j, int k, int jn, int nmax, int numneigh, double rij2,
     double* dsij1, double* dsij2, int ntype, int* type, int* fmap, double* x,
     double* scrfcn, double* fcpair)
{
  //     Inputs: i,j,k = id's of 3 atom triplet
  //     jn = id of i-j pair
  //     rij2 = squared distance between i and j
  //     Outputs: dsij1 = deriv. of sij w.r.t. rik
  //     dsij2 = deriv. of sij w.r.t. rjk

  arrdim2v(x, 3, nmax) int elti, eltj, eltk;
  double rik2, rjk2;

  double dxik, dyik, dzik;
  double dxjk, dyjk, dzjk;
  double rbound, delc, sij, xik, xjk, cikj, sikj, dfc, a;
  double Cmax, Cmin, dCikj1, dCikj2;

  sij = arr1v(scrfcn, jn) * arr1v(fcpair, jn);
  elti = arr1v(fmap, arr1v(type, i));
  eltj = arr1v(fmap, arr1v(type, j));
  eltk = arr1v(fmap, arr1v(type, k));
  Cmax = this->Cmax_meam[elti][eltj][eltk];
  Cmin = this->Cmin_meam[elti][eltj][eltk];

  *dsij1 = 0.0;
  *dsij2 = 0.0;
  if (!iszero(sij) && !iszero(sij - 1.0)) {
    rbound = rij2 * this->ebound_meam[elti][eltj];
    delc = Cmax - Cmin;
    dxjk = arr2v(x, 1, k) - arr2v(x, 1, j);
    dyjk = arr2v(x, 2, k) - arr2v(x, 2, j);
    dzjk = arr2v(x, 3, k) - arr2v(x, 3, j);
    rjk2 = dxjk * dxjk + dyjk * dyjk + dzjk * dzjk;
    if (rjk2 <= rbound) {
      dxik = arr2v(x, 1, k) - arr2v(x, 1, i);
      dyik = arr2v(x, 2, k) - arr2v(x, 2, i);
      dzik = arr2v(x, 3, k) - arr2v(x, 3, i);
      rik2 = dxik * dxik + dyik * dyik + dzik * dzik;
      if (rik2 <= rbound) {
        xik = rik2 / rij2;
        xjk = rjk2 / rij2;
        a = 1 - (xik - xjk) * (xik - xjk);
        if (!iszero(a)) {
          cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
          if (cikj >= Cmin && cikj <= Cmax) {
            cikj = (cikj - Cmin) / delc;
            dfcut(cikj, &sikj, &dfc);
            dCfunc2(rij2, rik2, rjk2, &dCikj1, &dCikj2);
            a = sij / delc * dfc / sikj;
            *dsij1 = a * dCikj1;
            *dsij2 = a * dCikj2;
          }
        }
      }
    }
  }
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

void
MEAM::fcut(double xi, double* fc)
{
  //     cutoff function
  double a;
  if (xi >= 1.0)
    *fc = 1.0;
  else if (xi <= 0.0)
    *fc = 0.0;
  else {
    a = 1.0 - xi;
    a = a * a;
    a = a * a;
    a = 1.0 - a;
    *fc = a * a;
    //     fc = xi
  }
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

void
MEAM::dfcut(double xi, double* fc, double* dfc)
{
  //     cutoff function and its derivative
  double a, a3, a4;
  if (xi >= 1.0) {
    *fc = 1.0;
    *dfc = 0.0;
  } else if (xi <= 0.0) {
    *fc = 0.0;
    *dfc = 0.0;
  } else {
    a = 1.0 - xi;
    a3 = a * a * a;
    a4 = a * a3;
    *fc = pow((1.0 - a4), 2);
    *dfc = 8 * (1.0 - a4) * a3;
    //     fc = xi
    //     dfc = 1.d0
  }
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

void
MEAM::dCfunc(double rij2, double rik2, double rjk2, double* dCikj)
{
  //     Inputs: rij,rij2,rik2,rjk2
  //     Outputs: dCikj = derivative of Cikj w.r.t. rij
  double rij4, a, b, denom;

  rij4 = rij2 * rij2;
  a = rik2 - rjk2;
  b = rik2 + rjk2;
  denom = rij4 - a * a;
  denom = denom * denom;
  *dCikj = -4 * (-2 * rij2 * a * a + rij4 * b + a * a * b) / denom;
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

void
MEAM::dCfunc2(double rij2, double rik2, double rjk2, double* dCikj1, double* dCikj2)
{
  //     Inputs: rij,rij2,rik2,rjk2
  //     Outputs: dCikj1 = derivative of Cikj w.r.t. rik
  //     dCikj2 = derivative of Cikj w.r.t. rjk
  double rij4, rik4, rjk4, a, b, denom;

  rij4 = rij2 * rij2;
  rik4 = rik2 * rik2;
  rjk4 = rjk2 * rjk2;
  a = rik2 - rjk2;
  b = rik2 + rjk2;
  denom = rij4 - a * a;
  denom = denom * denom;
  *dCikj1 = 4 * rij2 *
            (rij4 + rik4 + 2 * rik2 * rjk2 - 3 * rjk4 - 2 * rij2 * a) / denom;
  *dCikj2 = 4 * rij2 *
            (rij4 - 3 * rik4 + 2 * rik2 * rjk2 + rjk4 + 2 * rij2 * a) / denom;

  (void)(b);
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
