// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#include "meam.h"

#include "math_special.h"
#include "memory.h"

#include <cmath>
#include <algorithm>

using namespace LAMMPS_NS;
using namespace MEAM_NS;


/* ----------------------------------------------------------------------
   compute MEAM pair potential for each pair of element types
------------------------------------------------------------------------- */

void MEAM::compute_pair_meam()
{
  double r;
  int j, a, b, nv2;
  double astar, frac, phizbl;
  int Z1, Z2;
  double arat, rarat, scrn, scrn2;
  double phiaa, phibb /*unused:,phitmp*/;
  double C, s111, s112, s221, S11, S22;
  lattice_t lattaa, lattbb, lattab;

  // check for previously allocated arrays and free them
  if (phir != nullptr)
    memory->destroy(phir);
  if (phirar != nullptr)
    memory->destroy(phirar);
  if (phirar1 != nullptr)
    memory->destroy(phirar1);
  if (phirar2 != nullptr)
    memory->destroy(phirar2);
  if (phirar3 != nullptr)
    memory->destroy(phirar3);
  if (phirar4 != nullptr)
    memory->destroy(phirar4);
  if (phirar5 != nullptr)
    memory->destroy(phirar5);
  if (phirar6 != nullptr)
    memory->destroy(phirar6);

  // allocate memory for array that defines the potential
  memory->create(phir, (neltypes * (neltypes + 1)) / 2, nr, "pair:phir");

  // allocate coeff memory

  memory->create(phirar, (neltypes * (neltypes + 1)) / 2, nr, "pair:phirar");
  memory->create(phirar1, (neltypes * (neltypes + 1)) / 2, nr, "pair:phirar1");
  memory->create(phirar2, (neltypes * (neltypes + 1)) / 2, nr, "pair:phirar2");
  memory->create(phirar3, (neltypes * (neltypes + 1)) / 2, nr, "pair:phirar3");
  memory->create(phirar4, (neltypes * (neltypes + 1)) / 2, nr, "pair:phirar4");
  memory->create(phirar5, (neltypes * (neltypes + 1)) / 2, nr, "pair:phirar5");
  memory->create(phirar6, (neltypes * (neltypes + 1)) / 2, nr, "pair:phirar6");

  // loop over pairs of element types
  nv2 = 0;
  for (a = 0; a < neltypes; a++) {
    for (b = a; b < neltypes; b++) {
      // loop over r values and compute
      for (j = 0; j < nr; j++) {
        r = j * dr;
        phir[nv2][j] = phi_meam(r, a, b);

        // if using second-nearest neighbor, solve recursive problem
        // (see Lee and Baskes, PRB 62(13):8564 eqn.(21))
        if (nn2_meam[a][b] == 1) {
          lattaa = lattce_meam[a][a];
          lattbb = lattce_meam[b][b];
          lattab = lattce_meam[a][b];

          Z1 = get_Zij(lattab);
          Z2 = get_Zij2(lattab, Cmin_meam[a][a][b],
                     Cmax_meam[a][a][b], stheta_meam[a][b], arat, scrn);

          //     The B1, B2,  and L12 cases with NN2 have a trick to them; we need to
          //     compute the contributions from second nearest neighbors, like a-a
          //     pairs, but need to include NN2 contributions to those pairs as
          //     well.
          if (lattab == B1 || lattab == B2 || lattab == L12 || lattab == DIA) {
            rarat = r * arat;

            //               phi_aa
            phiaa = phi_meam(rarat, a, a);
            Z1 = get_Zij(lattaa);
            Z2 = get_Zij2(lattaa, Cmin_meam[a][a][a], Cmax_meam[a][a][a], stheta_meam[a][a], arat, scrn);
            phiaa+= phi_2nn_series(scrn, Z1, Z2, a, a, rarat, arat);

            //               phi_bb
            phibb = phi_meam(rarat, b, b);
            Z1 = get_Zij(lattbb);
            Z2 = get_Zij2(lattbb, Cmin_meam[b][b][b], Cmax_meam[b][b][b], stheta_meam[b][b], arat, scrn);
            phibb+= phi_2nn_series(scrn, Z1, Z2, b, b, rarat, arat);

            if (lattab == B1 || lattab == B2 || lattab == DIA) {
              //     Add contributions to the B1 or B2 potential
              Z1 = get_Zij(lattab);
              Z2 = get_Zij2(lattab, Cmin_meam[a][a][b], Cmax_meam[a][a][b], stheta_meam[a][b], arat, scrn);
              phir[nv2][j] -= Z2 * scrn / (2 * Z1) * phiaa;
              Z2 = get_Zij2(lattab, Cmin_meam[b][b][a], Cmax_meam[b][b][a], stheta_meam[a][b], arat, scrn2);
              phir[nv2][j] -= Z2 * scrn2 / (2 * Z1) * phibb;

            } else if (lattab == L12) {
              //     The L12 case has one last trick; we have to be careful to
              //     compute
              //     the correct screening between 2nd-neighbor pairs.  1-1
              //     second-neighbor pairs are screened by 2 type 1 atoms and
              //     two type
              //     2 atoms.  2-2 second-neighbor pairs are screened by 4 type
              //     1
              //     atoms.
              C = 1.0;
              s111 = get_sijk(C, a, a, a);
              s112 = get_sijk(C, a, a, b);
              s221 = get_sijk(C, b, b, a);
              S11 = s111 * s111 * s112 * s112;
              S22 = pow(s221, 4);
              phir[nv2][j] = phir[nv2][j] - 0.75 * S11 * phiaa - 0.25 * S22 * phibb;
            }

          } else {
            phir[nv2][j]+= phi_2nn_series(scrn, Z1, Z2, a, b, r, arat);
          }
        }

        // For Zbl potential:
        // if astar <= -3
        //   potential is zbl potential
        // else if -3 < astar < -1
        //   potential is linear combination with zbl potential
        // endif
        if (zbl_meam[a][b] == 1) {
          astar = alpha_meam[a][b] * (r / re_meam[a][b] - 1.0);
          if (astar <= -3.0)
            phir[nv2][j] = zbl(r, ielt_meam[a], ielt_meam[b]);
          else if (astar > -3.0 && astar < -1.0) {
            frac = fcut(1 - (astar + 1.0) / (-3.0 + 1.0));
            phizbl = zbl(r, ielt_meam[a], ielt_meam[b]);
            phir[nv2][j] = frac * phir[nv2][j] + (1 - frac) * phizbl;
          }
        }
      }

      // call interpolation
      interpolate_meam(nv2);

      nv2 = nv2 + 1;
    }
  }
}

/* ----------------------------------------------------------------------
   Compute \bar{rho} for both a and b at distance r
   (Eqn. 7 in Huang 1995)
------------------------------------------------------------------------- */

bool MEAM::rhobar12(const double r, const int a, const int b, double &rhobar1, double &rhobar2) const
{
  double rho01, rho11, rho21, rho31;
  double rho02, rho12, rho22, rho32;
  double t11av, t21av, t31av, t12av, t22av, t32av;
  // msmeam
  double rho1m1, rho2m1, rho3m1;
  double rho1m2, rho2m2, rho3m2;
  double t1m1av, t2m1av, t3m1av, t1m2av, t2m2av, t3m2av;

  int errorflag, Z12;
  lattice_t latta;
  double G1, G2, s1[3], s2[3];
  double Gam1, Gam2, Z1, Z2;
  double rho_bkgd1, rho_bkgd2;

  // the last 6 arguments are only touched for msmeam
  get_densref(r, a, b, &rho01, &rho11, &rho21, &rho31, &rho02, &rho12, &rho22, &rho32,
              &rho1m1, &rho2m1, &rho3m1,
              &rho1m2, &rho2m2, &rho3m2);

  // if densities are too small, numerical problems may result; let caller return zero
  if (rho01 <= 1e-14 && rho02 <= 1e-14)
    return false;

  // calculate average weighting factors for the reference structure
  // average weighting factors for the reference structure, eqn. I.8
  get_tavref(&t11av, &t21av, &t31av, &t12av, &t22av, &t32av, t1_meam[a], t2_meam[a],
             t3_meam[a], t1_meam[b], t2_meam[b], t3_meam[b], rho01, rho02, r, a, b,
             lattce_meam[a][b]);
  // with msmeam call twice with different sets of variables
  if (msmeamflag) {
    get_tavref(&t1m1av, &t2m1av, &t3m1av, &t1m2av, &t2m2av, &t3m2av, t1m_meam[a], t2m_meam[a],
              t3m_meam[a], t1m_meam[b], t2m_meam[b], t3m_meam[b], rho01, rho02, r, a, b,
              lattce_meam[a][b]);
  }

  switch (lattce_meam[a][b]) {
    case C11:
      // for c11b structure we have analytic solution, swapped depending on which side has the DIA atom
      Z12 = get_Zij(lattce_meam[a][b]);
      latta = lattce_meam[a][a];
      if (latta == DIA) {
        rhobar1 = MathSpecial::square((Z12 / 2) * (rho02 + rho01))
                  + t11av * MathSpecial::square(rho12 - rho11)
                  + t21av / 6.0 * MathSpecial::square(rho22 + rho21)
                  + 121.0 / 40.0 * t31av * MathSpecial::square(rho32 - rho31);
        rhobar1 = sqrt(rhobar1);
        rhobar2 = MathSpecial::square(Z12 * rho01) + 2.0 / 3.0 * t21av * MathSpecial::square(rho21);
        rhobar2 = sqrt(rhobar2);
      } else {
        rhobar2 = MathSpecial::square((Z12 / 2) * (rho01 + rho02))
                  + t12av * MathSpecial::square(rho11 - rho12)
                  + t22av / 6.0 * MathSpecial::square(rho21 + rho22)
                  + 121.0 / 40.0 * t32av * MathSpecial::square(rho31 - rho32);
        rhobar2 = sqrt(rhobar2);
        rhobar1 = MathSpecial::square(Z12 * rho02) + 2.0 / 3.0 * t22av * MathSpecial::square(rho22);
        rhobar1 = sqrt(rhobar1);
      }
      break;
    default:
      // for other structures, use the densities computed before
      Z1 = get_Zij(lattce_meam[a][a]);
      Z2 = get_Zij(lattce_meam[b][b]);

      if (mix_ref_t == 1) {
        // This is the original formalism based on the reference structure developed in Huang's paper, equation I.7
        // WARNING: this is not correct and only provided for compatibility.
        // To get the equations from the paper, bkgd_dyn=1 and ibar=1
        if (ibar_meam[a] <= 0)
          G1 = 1.0;
        else {
          get_shpfcn(lattce_meam[a][a], stheta_meam[a][a], ctheta_meam[a][a], s1);
          Gam1 = (s1[0] * t11av + s1[1] * t21av + s1[2] * t31av) / (Z1 * Z1);
          G1 = G_gam(Gam1, ibar_meam[a], errorflag);
        }
        if (ibar_meam[b] <= 0)
          G2 = 1.0;
        else {
          get_shpfcn(lattce_meam[b][b], stheta_meam[b][b], ctheta_meam[b][b],  s2);
          Gam2 = (s2[0] * t12av + s2[1] * t22av + s2[2] * t32av) / (Z2 * Z2);
          G2 = G_gam(Gam2, ibar_meam[b], errorflag);
        }
        rho_bkgd1 = rho0_meam[a] * Z1 * G1;
        rho_bkgd2 = rho0_meam[b] * Z2 * G2;
      }

      if (msmeamflag) {
        // no additional use of t's here; all included in definitions of rho's for msmeam
        Gam1 = rho11 + rho21 + rho31 - (rho1m1 + rho2m1 + rho3m1);
        Gam2 = rho12 + rho22 + rho32 - (rho1m2 + rho2m2 + rho3m2);
      } else {
        Gam1 = t11av * rho11 + t21av * rho21 + t31av * rho31;
        Gam2 = t12av * rho12 + t22av * rho22 + t32av * rho32;
      }
      Gam1 = (rho01 < 1.0e-14) ? 0.0 : Gam1 / (rho01 * rho01);
      Gam2 = (rho02 < 1.0e-14) ? 0.0 : Gam2 / (rho02 * rho02);

      G1 = G_gam(Gam1, ibar_meam[a], errorflag);
      G2 = G_gam(Gam2, ibar_meam[b], errorflag);
      if (bkgd_dyn == 1) {
        rho_bkgd1 = rho0_meam[a] * Z1;
        rho_bkgd2 = rho0_meam[b] * Z2;
      } else {
        rho_bkgd1 = rho_ref_meam[a];
        rho_bkgd2 = rho_ref_meam[b];
      }
      rhobar1 = rho01 / rho_bkgd1 * G1;
      rhobar2 = rho02 / rho_bkgd2 * G2;
  }

  return true;
}

/* ----------------------------------------------------------------------
   Invert the EAM energy function for a reference structure with known energy and embedding
   (Eqn. 18 in Huang 1995)
------------------------------------------------------------------------- */

double MEAM::invert_eam(const double r, const int a, const int b, const double Eu, const double F1, const double F2) const
{
  int Z1, Z2, Z12;
  int Z1nn, Z2nn;
  lattice_t latta;
  double phiaa, phibb;
  double arat, scrn, scrn2;
  double b11s, b22s;
  double phi_m = 0.0;

  Z12 = get_Zij(lattce_meam[a][b]);

  switch (lattce_meam[a][b]) {
    case C11:
      latta = lattce_meam[a][a];
      if (latta == DIA) {
        phiaa = phi_meam(r, a, a);
        phi_m = (3 * Eu - F2 - 2 * F1 - 5 * phiaa) / Z12;
      } else {
        phibb = phi_meam(r, b, b);
        phi_m = (3 * Eu - F1 - 2 * F2 - 5 * phibb) / Z12;
      }
      break;
    case L12:
      phiaa = phi_meam(r, a, a);
      //       account for second neighbor a-a potential here...
      Z1nn = get_Zij(lattce_meam[a][a]);
      Z2nn = get_Zij2(lattce_meam[a][a], Cmin_meam[a][a][a],
               Cmax_meam[a][a][a], stheta_meam[a][b], arat, scrn);

      phiaa += phi_2nn_series(scrn, Z1nn, Z2nn, a, a, r, arat);
      phi_m = Eu / 3.0 - F1 / 4.0 - F2 / 12.0 - phiaa;
      break;
    case CH4:
      phi_m = (5 * Eu - F1 - 4*F2)/4;
      break;
    case ZIG:
      if (a==b){
        phi_m = (2 * Eu - F1 - F2) / Z12;
      } else{
        Z1 = get_Zij(lattce_meam[a][b]);
        Z2 = get_Zij2_b2nn(lattce_meam[a][b], Cmin_meam[a][a][b], Cmax_meam[a][a][b], scrn);
        b11s = -Z2/Z1*scrn;
        Z2 = get_Zij2_b2nn(lattce_meam[a][b], Cmin_meam[b][b][a], Cmax_meam[b][b][a], scrn2);
        b22s = -Z2/Z1*scrn2;

        phiaa = phi_meam(2.0*stheta_meam[a][b]*r, a, a);
        phibb = phi_meam(2.0*stheta_meam[a][b]*r, b, b);
        phi_m = (2.0*Eu - F1 - F2 + phiaa*b11s + phibb*b22s) / Z12;
      }
      break;
    case TRI:
      if (a==b){
        phi_m = (3.0*Eu - 2.0*F1 - F2) / Z12;
      } else {
        Z1 = get_Zij(lattce_meam[a][b]);
        Z2 = get_Zij2_b2nn(lattce_meam[a][b], Cmin_meam[a][a][b], Cmax_meam[a][a][b], scrn);
        b11s = -Z2/Z1*scrn;
        phiaa = phi_meam(2.0*stheta_meam[a][b]*r, a, a);
        phi_m = (3.0*Eu - 2.0*F1 - F2 + phiaa*b11s) / Z12;
      }
      break;
    default:
      // potential is computed from Rose function and embedding energy
      phi_m = (2 * Eu - F1 - F2) / Z12;
  }

  return phi_m;
}


/* ----------------------------------------------------------------------
   Compute MEAM pair potential for distance r, element types a and b
------------------------------------------------------------------------- */

double MEAM::phi_meam(double r, int a, int b) const
{
  double rhobar1, rhobar2;
  double F1, F2, dF, Eu, phi_m;

  // if r = 0, just return 0
  if (iszero(r))
    return 0.0;

  // Equation numbers below refer to:
  //   I. Huang et.al., Modelling simul. Mater. Sci. Eng. 3:615

  // calculate background electron densities rhobar1 and rhobar2
  if (!rhobar12(r, a, b, rhobar1, rhobar2))
    return 0.0;

  // compute embedding functions, eqn I.5

  F1 = embedding(A_meam[a], Ec_meam[a][a], rhobar1, dF);
  F2 = embedding(A_meam[b], Ec_meam[b][b], rhobar2, dF);

  // compute Rose function, I.16
  Eu = erose(r, re_meam[a][b], alpha_meam[a][b], Ec_meam[a][b], repuls_meam[a][b],
             attrac_meam[a][b], erose_form);

  phi_m = invert_eam(r, a, b, Eu, F1, F2);

  return phi_m;
}

/* ----------------------------------------------------------------------
   Compute 2NN series terms for phi
     To avoid nan values of phir due to rapid decrease of b2nn^n or/and
     argument of phi_meam, i.e. r*arat^n, in some cases (3NN dia with low Cmin value)
------------------------------------------------------------------------- */

double MEAM::phi_2nn_series(const double scrn, const int Z1, const int Z2, const int a, const int b,
                            const double r, const double arat) const
{
  double phi_sum = 0.0;
  double b2nn, phi_val;
  if (scrn > 0.0) {
    b2nn = -Z2*scrn/Z1;
    for (int n = 1; n <= 10; n++) {
      phi_val = MathSpecial::powint(b2nn,n) * phi_meam(r * MathSpecial::powint(arat, n), a, b);
      if (iszero(phi_val)) {
        // once either term becomes zero at some point, all folliwng will also be zero
        // necessary to avoid numerical error (nan or infty) due to exponential decay in phi_meam
        break;
      }
      phi_sum += phi_val;
    }
  }
  return phi_sum;
}

/* ----------------------------------------------------------------------
   Calculate screening ellipse
------------------------------------------------------------------------- */

double MEAM::get_sijk(double C, int i, int j, int k) const
{
  return Csijk(C, Cmin_meam[i][j][k], Cmax_meam[i][j][k]);
}

/* ----------------------------------------------------------------------
   Interpolate coefficients for phi and phi'
------------------------------------------------------------------------- */

void MEAM::interpolate_meam(int ind)
{
  int j;
  double drar;

  // map to coefficient space

  nrar = nr;
  drar = dr;
  rdrar = 1.0 / drar;

  // phir interp

  for (j = 0; j < nrar; j++) {
    phirar[ind][j] = phir[ind][j];
  }
  phirar1[ind][0] = phirar[ind][1] - phirar[ind][0];
  phirar1[ind][1] = 0.5 * (phirar[ind][2] - phirar[ind][0]);
  phirar1[ind][nrar - 2] =
    0.5 * (phirar[ind][nrar - 1] - phirar[ind][nrar - 3]);
  phirar1[ind][nrar - 1] = 0.0;
  for (j = 2; j < nrar - 2; j++) {
    phirar1[ind][j] = ((phirar[ind][j - 2] - phirar[ind][j + 2]) +
                             8.0 * (phirar[ind][j + 1] - phirar[ind][j - 1])) /
                            12.;
  }

  for (j = 0; j < nrar - 1; j++) {
    phirar2[ind][j] = 3.0 * (phirar[ind][j + 1] - phirar[ind][j]) -
                            2.0 * phirar1[ind][j] - phirar1[ind][j + 1];
    phirar3[ind][j] = phirar1[ind][j] + phirar1[ind][j + 1] -
                            2.0 * (phirar[ind][j + 1] - phirar[ind][j]);
  }
  phirar2[ind][nrar - 1] = 0.0;
  phirar3[ind][nrar - 1] = 0.0;

  for (j = 0; j < nrar; j++) {
    phirar4[ind][j] = phirar1[ind][j] / drar;
    phirar5[ind][j] = 2.0 * phirar2[ind][j] / drar;
    phirar6[ind][j] = 3.0 * phirar3[ind][j] / drar;
  }
}
