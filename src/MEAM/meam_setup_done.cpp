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

void MEAM::meam_setup_done(double* cutmax)
{
  int nv2, nv3, m, n, p;

  //     Force cutoff
  cutforce = rc_meam;
  cutforcesq = cutforce * cutforce;

  //     Pass cutoff back to calling program
  *cutmax = cutforce;

  //     Augment t1 term
  for (int i = 0; i < MAXELT; i++)
    t1_meam[i] = t1_meam[i] + augt1 * 3.0 / 5.0 * t3_meam[i];

  //     Compute off-diagonal alloy parameters
  alloyparams();

  // indices and factors for Voight notation
  nv2 = 0;
  nv3 = 0;
  for (m = 0; m < 3; m++) {
    for (n = m; n < 3; n++) {
      vind2D[m][n] = nv2;
      vind2D[n][m] = nv2;
      nv2 = nv2 + 1;
      for (p = n; p < 3; p++) {
        vind3D[m][n][p] = nv3;
        vind3D[m][p][n] = nv3;
        vind3D[n][m][p] = nv3;
        vind3D[n][p][m] = nv3;
        vind3D[p][m][n] = nv3;
        vind3D[p][n][m] = nv3;
        nv3 = nv3 + 1;
      }
    }
  }

  v2D[0] = 1;
  v2D[1] = 2;
  v2D[2] = 2;
  v2D[3] = 1;
  v2D[4] = 2;
  v2D[5] = 1;

  v3D[0] = 1;
  v3D[1] = 3;
  v3D[2] = 3;
  v3D[3] = 3;
  v3D[4] = 6;
  v3D[5] = 3;
  v3D[6] = 1;
  v3D[7] = 3;
  v3D[8] = 3;
  v3D[9] = 1;

  nv2 = 0;
  for (m = 0; m < neltypes; m++) {
    for (n = m; n < neltypes; n++) {
      eltind[m][n] = nv2;
      eltind[n][m] = nv2;
      nv2 = nv2 + 1;
    }
  }

  //     Compute background densities for reference structure
  compute_reference_density();

  //     Compute pair potentials and setup arrays for interpolation
  nr = 1000;
  dr = 1.1 * rc_meam / nr;
  compute_pair_meam();
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// Fill off-diagonal alloy parameters
void MEAM::alloyparams()
{

  int i, j, k;
  double eb;

  // Loop over pairs
  for (i = 0; i < neltypes; i++) {
    for (j = 0; j < neltypes; j++) {
      // Treat off-diagonal pairs
      // If i>j, set all equal to i<j case (which has already been set,
      // here or in the input file)
      if (i > j) {
        re_meam[i][j] = re_meam[j][i];
        Ec_meam[i][j] = Ec_meam[j][i];
        alpha_meam[i][j] = alpha_meam[j][i];
        lattce_meam[i][j] = lattce_meam[j][i];
        nn2_meam[i][j] = nn2_meam[j][i];
        // theta for lin,tri,zig references
        stheta_meam[i][j] = stheta_meam[j][i];
        ctheta_meam[i][j] = ctheta_meam[j][i];
        // If i<j and term is unset, use default values (e.g. mean of i-i and
        // j-j)
      } else if (j > i) {
        if (iszero(Ec_meam[i][j])) {
          if (lattce_meam[i][j] == L12)
            Ec_meam[i][j] =
              (3 * Ec_meam[i][i] + Ec_meam[j][j]) / 4.0 - delta_meam[i][j];
          else if (lattce_meam[i][j] == C11) {
            if (lattce_meam[i][i] == DIA)
              Ec_meam[i][j] =
                (2 * Ec_meam[i][i] + Ec_meam[j][j]) / 3.0 - delta_meam[i][j];
            else
              Ec_meam[i][j] =
                (Ec_meam[i][i] + 2 * Ec_meam[j][j]) / 3.0 - delta_meam[i][j];
          } else
            Ec_meam[i][j] = (Ec_meam[i][i] + Ec_meam[j][j]) / 2.0 - delta_meam[i][j];
        }
        if (iszero(alpha_meam[i][j]))
          alpha_meam[i][j] = (alpha_meam[i][i] + alpha_meam[j][j]) / 2.0;
        if (iszero(re_meam[i][j]))
          re_meam[i][j] = (re_meam[i][i] + re_meam[j][j]) / 2.0;
      }
    }
  }

  // Cmin[i][k][j] is symmetric in i-j, but not k.  For all triplets
  // where i>j, set equal to the i<j element.  Likewise for Cmax.
  for (i = 1; i < neltypes; i++) {
    for (j = 0; j < i; j++) {
      for (k = 0; k < neltypes; k++) {
        Cmin_meam[i][j][k] = Cmin_meam[j][i][k];
        Cmax_meam[i][j][k] = Cmax_meam[j][i][k];
      }
    }
  }

  // ebound gives the squared distance such that, for rik2 or rjk2>ebound,
  // atom k definitely lies outside the screening function ellipse (so
  // there is no need to calculate its effects).  Here, compute it for all
  // triplets [i][j][k] so that ebound[i][j] is the maximized over k
  for (i = 0; i < neltypes; i++) {
    for (j = 0; j < neltypes; j++) {
      for (k = 0; k < neltypes; k++) {
        eb = (Cmax_meam[i][j][k] * Cmax_meam[i][j][k]) / (4.0 * (Cmax_meam[i][j][k] - 1.0));
        ebound_meam[i][j] = std::max(ebound_meam[i][j], eb);
      }
    }
  }
}

//-----------------------------------------------------------------------
// compute MEAM pair potential for each pair of element types
//

void MEAM::compute_pair_meam()
{
  double r;
  int j, a, b, nv2;
  double astar, frac, phizbl;
  int Z1, Z2;
  double arat, rarat, scrn, scrn2;
  double phiaa, phibb /*unused:,phitmp*/;
  double C, s111, s112, s221, S11, S22;

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
          Z1 = get_Zij(lattce_meam[a][b]);
          Z2 = get_Zij2(lattce_meam[a][b], Cmin_meam[a][a][b],
                     Cmax_meam[a][a][b], stheta_meam[a][b], arat, scrn);

          //     The B1, B2,  and L12 cases with NN2 have a trick to them; we need to
          //     compute the contributions from second nearest neighbors, like a-a
          //     pairs, but need to include NN2 contributions to those pairs as
          //     well.
          if (lattce_meam[a][b] == B1 || lattce_meam[a][b] == B2 ||
              lattce_meam[a][b] == L12 || lattce_meam[a][b] == DIA) {
            rarat = r * arat;

            //               phi_aa
            phiaa = phi_meam(rarat, a, a);
            Z1 = get_Zij(lattce_meam[a][a]);
            Z2 = get_Zij2(lattce_meam[a][a], Cmin_meam[a][a][a],
                     Cmax_meam[a][a][a], stheta_meam[a][a], arat, scrn);
            phiaa+= phi_meam_series(scrn, Z1, Z2, a, a, rarat, arat);

            //               phi_bb
            phibb = phi_meam(rarat, b, b);
            Z1 = get_Zij(lattce_meam[b][b]);
            Z2 = get_Zij2(lattce_meam[b][b], Cmin_meam[b][b][b],
                     Cmax_meam[b][b][b], stheta_meam[b][b], arat, scrn);
            phibb+= phi_meam_series(scrn, Z1, Z2, b, b, rarat, arat);

            if (lattce_meam[a][b] == B1 || lattce_meam[a][b] == B2 ||
                lattce_meam[a][b] == DIA) {
              //     Add contributions to the B1 or B2 potential
              Z1 = get_Zij(lattce_meam[a][b]);
              Z2 = get_Zij2(lattce_meam[a][b], Cmin_meam[a][a][b],
                       Cmax_meam[a][a][b], stheta_meam[a][b],  arat, scrn);
              phir[nv2][j] = phir[nv2][j] - Z2 * scrn / (2 * Z1) * phiaa;
              Z2 = get_Zij2(lattce_meam[a][b], Cmin_meam[b][b][a],
                       Cmax_meam[b][b][a], stheta_meam[a][b], arat, scrn2);

              phir[nv2][j] = phir[nv2][j] - Z2 * scrn2 / (2 * Z1) * phibb;

            } else if (lattce_meam[a][b] == L12) {
              //     The L12 case has one last trick; we have to be careful to
              //     compute
              //     the correct screening between 2nd-neighbor pairs.  1-1
              //     second-neighbor pairs are screened by 2 type 1 atoms and
              //     two type
              //     2 atoms.  2-2 second-neighbor pairs are screened by 4 type
              //     1
              //     atoms.
              C = 1.0;
              get_sijk(C, a, a, a, &s111);
              get_sijk(C, a, a, b, &s112);
              get_sijk(C, b, b, a, &s221);
              S11 = s111 * s111 * s112 * s112;
              S22 = pow(s221, 4);
              phir[nv2][j] = phir[nv2][j] - 0.75 * S11 * phiaa - 0.25 * S22 * phibb;
            }

          } else {
            phir[nv2][j]+= phi_meam_series(scrn, Z1, Z2, a, b, r, arat);
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

//----------------------------------------------------------------------c
// Compute MEAM pair potential for distance r, element types a and b
//
double MEAM::phi_meam(double r, int a, int b)
{
  /*unused:double a1,a2,a12;*/
  double t11av, t21av, t31av, t12av, t22av, t32av;
  double G1, G2, s1[3], s2[3], rho0_1, rho0_2;
  double Gam1, Gam2, Z1, Z2;
  double rhobar1, rhobar2, F1, F2, dF;
  double rho01, rho11, rho21, rho31;
  double rho02, rho12, rho22, rho32;
  double scalfac, phiaa, phibb;
  double Eu;
  double arat, scrn, scrn2;
  int Z12, errorflag;
  int Z1nn, Z2nn;
  lattice_t latta /*unused:,lattb*/;
  double rho_bkgd1, rho_bkgd2;
  double b11s, b22s;
  // msmeam
  double t1m1av, t2m1av, t3m1av, t1m2av, t2m2av, t3m2av;
  double rho1m1, rho2m1, rho3m1;
  double rho1m2, rho2m2, rho3m2;

  double phi_m = 0.0;
  // Equation numbers below refer to:
  //   I. Huang et.al., Modelling simul. Mater. Sci. Eng. 3:615

  // get number of neighbors in the reference structure
  //   Nref[i][j] = # of i's neighbors of type j
  Z1 = get_Zij(lattce_meam[a][a]);
  Z2 = get_Zij(lattce_meam[b][b]);
  Z12 = get_Zij(lattce_meam[a][b]);

  // this function has extra args for msmeam
  if (msmeamflag) {
    get_densref(r, a, b, &rho01, &rho11, &rho21, &rho31, &rho02, &rho12, &rho22, &rho32,
                &rho1m1, &rho2m1, &rho3m1,
                &rho1m2, &rho2m2, &rho3m2);
  } else {
    get_densref(r, a, b, &rho01, &rho11, &rho21, &rho31, &rho02, &rho12, &rho22, &rho32,
                nullptr, nullptr, nullptr,
                nullptr, nullptr, nullptr);
  }
  // if densities are too small, numerical problems may result; just return zero
  if (rho01 <= 1e-14 && rho02 <= 1e-14)
    return 0.0;

  // calculate average weighting factors for the reference structure
  if (lattce_meam[a][b] == C11) {
    if (ialloy == 2) {
      t11av = t1_meam[a];
      t12av = t1_meam[b];
      t21av = t2_meam[a];
      t22av = t2_meam[b];
      t31av = t3_meam[a];
      t32av = t3_meam[b];
    } else {
      scalfac = 1.0 / (rho01 + rho02);
      t11av = scalfac * (t1_meam[a] * rho01 + t1_meam[b] * rho02);
      t12av = t11av;
      t21av = scalfac * (t2_meam[a] * rho01 + t2_meam[b] * rho02);
      t22av = t21av;
      t31av = scalfac * (t3_meam[a] * rho01 + t3_meam[b] * rho02);
      t32av = t31av;
    }
  } else {
    // average weighting factors for the reference structure, eqn. I.8
    get_tavref(&t11av, &t21av, &t31av, &t12av, &t22av, &t32av, t1_meam[a], t2_meam[a],
               t3_meam[a], t1_meam[b], t2_meam[b], t3_meam[b], r, a, b,
               lattce_meam[a][b]);
    // with msmeam call twice with different sets of variables
    if (msmeamflag) {
      get_tavref(&t1m1av, &t2m1av, &t3m1av, &t1m2av, &t2m2av, &t3m2av, t1m_meam[a], t2m_meam[a],
                t3m_meam[a], t1m_meam[b], t2m_meam[b], t3m_meam[b], r, a, b,
                lattce_meam[a][b]);
    }
  }

  // for c11b structure, calculate background electron densities
  if (lattce_meam[a][b] == C11) {
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
  } else {
    // for other structures, use formalism developed in Huang's paper
    //
    //     composition-dependent scaling, equation I.7
    //     If using mixing rule for t, apply to reference structure; else
    //     use precomputed values
    if (mix_ref_t == 1) {
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
      rho0_1 = rho0_meam[a] * Z1 * G1;
      rho0_2 = rho0_meam[b] * Z2 * G2;
    }

    if (msmeamflag) {
      // no additional use of t's here; all included in definitions of rho's for msmeam
      Gam1 = rho11 + rho21 + rho31 - (rho1m1 + rho2m1 + rho3m1);
      if (rho01 < 1.0e-14)
        Gam1 = 0.0;
      else
        Gam1 = Gam1 / (rho01 * rho01);
      Gam2 = rho12 + rho22 + rho32 - (rho1m2 + rho2m2 + rho3m2);
      if (rho02 < 1.0e-14)
        Gam2 = 0.0;
      else
        Gam2 = Gam2 / (rho02 * rho02);

    } else {
      Gam1 = (t11av * rho11 + t21av * rho21 + t31av * rho31);
      if (rho01 < 1.0e-14)
        Gam1 = 0.0;
      else
        Gam1 = Gam1 / (rho01 * rho01);

      Gam2 = (t12av * rho12 + t22av * rho22 + t32av * rho32);
      if (rho02 < 1.0e-14)
        Gam2 = 0.0;
      else
        Gam2 = Gam2 / (rho02 * rho02);
    }

    G1 = G_gam(Gam1, ibar_meam[a], errorflag);
    G2 = G_gam(Gam2, ibar_meam[b], errorflag);
    if (mix_ref_t == 1) {
      rho_bkgd1 = rho0_1;
      rho_bkgd2 = rho0_2;
    } else {
      if (bkgd_dyn == 1) {
        rho_bkgd1 = rho0_meam[a] * Z1;
        rho_bkgd2 = rho0_meam[b] * Z2;
      } else {
        rho_bkgd1 = rho_ref_meam[a];
        rho_bkgd2 = rho_ref_meam[b];
      }
    }
    rhobar1 = rho01 / rho_bkgd1 * G1;
    rhobar2 = rho02 / rho_bkgd2 * G2;
  }

  // compute embedding functions, eqn I.5

  F1 = embedding(A_meam[a], Ec_meam[a][a], rhobar1, dF);
  F2 = embedding(A_meam[b], Ec_meam[b][b], rhobar2, dF);


  // compute Rose function, I.16
  Eu = erose(r, re_meam[a][b], alpha_meam[a][b], Ec_meam[a][b], repuls_meam[a][b],
             attrac_meam[a][b], erose_form);

  // calculate the pair energy
  if (lattce_meam[a][b] == C11) {
    latta = lattce_meam[a][a];
    if (latta == DIA) {
      phiaa = phi_meam(r, a, a);
      phi_m = (3 * Eu - F2 - 2 * F1 - 5 * phiaa) / Z12;
    } else {
      phibb = phi_meam(r, b, b);
      phi_m = (3 * Eu - F1 - 2 * F2 - 5 * phibb) / Z12;
    }
  } else if (lattce_meam[a][b] == L12) {
    phiaa = phi_meam(r, a, a);
    //       account for second neighbor a-a potential here...
    Z1nn = get_Zij(lattce_meam[a][a]);
    Z2nn = get_Zij2(lattce_meam[a][a], Cmin_meam[a][a][a],
             Cmax_meam[a][a][a], stheta_meam[a][b], arat, scrn);


    phiaa += phi_meam_series(scrn, Z1nn, Z2nn, a, a, r, arat);
    phi_m = Eu / 3.0 - F1 / 4.0 - F2 / 12.0 - phiaa;

  } else if (lattce_meam[a][b] == CH4) {
    phi_m = (5 * Eu - F1 - 4*F2)/4;


  } else if (lattce_meam[a][b] == ZIG) {
      if (a==b) {
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

  } else if (lattce_meam[a][b] == TRI) {
      if (a==b) {
        phi_m = (3.0*Eu - 2.0*F1 - F2) / Z12;
     } else {
        Z1 = get_Zij(lattce_meam[a][b]);
        Z2 = get_Zij2_b2nn(lattce_meam[a][b], Cmin_meam[a][a][b], Cmax_meam[a][a][b], scrn);
        b11s = -Z2/Z1*scrn;
        phiaa = phi_meam(2.0*stheta_meam[a][b]*r, a, a);
        phi_m = (3.0*Eu - 2.0*F1 - F2 + phiaa*b11s) / Z12;
      }

  } else {
    // potential is computed from Rose function and embedding energy
    phi_m = (2 * Eu - F1 - F2) / Z12;
  }

  // if r = 0, just return 0
  if (iszero(r)) {
    phi_m = 0.0;
  }

  return phi_m;
}

//----------------------------------------------------------------------c
// Compute 2NN series terms for phi
//   To avoid nan values of phir due to rapid decrease of b2nn^n or/and
//   argument of phi_meam, i.e. r*arat^n, in some cases (3NN dia with low Cmin value)
//
double MEAM::phi_meam_series(const double scrn, const int Z1, const int Z2, const int a, const int b,
                             const double r, const double arat)
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

//----------------------------------------------------------------------c
// Compute background density for reference structure of each element
void MEAM::compute_reference_density()
{
  int a, Z, Z2, errorflag;
  double gam, Gbar, shp[3];
  double rho0, rho0_2nn, arat, scrn;

  // loop over element types
  for (a = 0; a < neltypes; a++) {
    Z = get_Zij(lattce_meam[a][a]);
    if (ibar_meam[a] <= 0)
      Gbar = 1.0;
    else {
      get_shpfcn(lattce_meam[a][a], stheta_meam[a][a], ctheta_meam[a][a], shp);
      gam = (t1_meam[a] * shp[0] + t2_meam[a] * shp[1] + t3_meam[a] * shp[2]) / (Z * Z);
      Gbar = G_gam(gam, ibar_meam[a], errorflag);
    }

    //     The zeroth order density in the reference structure, with
    //     equilibrium spacing, is just the number of first neighbors times
    //     the rho0_meam coefficient...
    rho0 = rho0_meam[a] * Z;

    //     ...unless we have unscreened second neighbors, in which case we
    //     add on the contribution from those (accounting for partial
    //     screening)
    if (nn2_meam[a][a] == 1) {
      Z2 = get_Zij2(lattce_meam[a][a], Cmin_meam[a][a][a],
               Cmax_meam[a][a][a], stheta_meam[a][a], arat, scrn);
      rho0_2nn = rho0_meam[a] * MathSpecial::fm_exp(-beta0_meam[a] * (arat - 1));
      rho0 = rho0 + Z2 * rho0_2nn * scrn;
    }

    rho_ref_meam[a] = rho0 * Gbar;
  }
}

//------------------------------------------------------------------------------c
// Average weighting factors for the reference structure
void MEAM::get_tavref(double* t11av, double* t21av, double* t31av, double* t12av, double* t22av, double* t32av,
                      double t11, double t21, double t31, double t12, double t22, double t32, double r, int a,
                      int b, lattice_t latt)
{
  double rhoa01, rhoa02, a1, a2, rho01 /*,rho02*/;

  //     For ialloy = 2, no averaging is done
  if (ialloy == 2) {
    *t11av = t11;
    *t21av = t21;
    *t31av = t31;
    *t12av = t12;
    *t22av = t22;
    *t32av = t32;
  } else switch (latt)  {
    case FCC:
    case BCC:
    case DIA:
    case DIA3:
    case HCP:
    case B1:
    case DIM:
    case B2:
    case CH4:
    case LIN:
    case ZIG:
    case TRI:
    case SC:
      //     all neighbors are of the opposite type
      *t11av = t12;
      *t21av = t22;
      *t31av = t32;
      *t12av = t11;
      *t22av = t21;
      *t32av = t31;
      break;
    default:
      a1 = r / re_meam[a][a] - 1.0;
      a2 = r / re_meam[b][b] - 1.0;
      rhoa01 = rho0_meam[a] * MathSpecial::fm_exp(-beta0_meam[a] * a1);
      rhoa02 = rho0_meam[b] * MathSpecial::fm_exp(-beta0_meam[b] * a2);
      if (latt == L12) {
        rho01 = 8 * rhoa01 + 4 * rhoa02;
        *t11av = (8 * t11 * rhoa01 + 4 * t12 * rhoa02) / rho01;
        *t12av = t11;
        *t21av = (8 * t21 * rhoa01 + 4 * t22 * rhoa02) / rho01;
        *t22av = t21;
        *t31av = (8 * t31 * rhoa01 + 4 * t32 * rhoa02) / rho01;
        *t32av = t31;
      } else {
        //      call error('Lattice not defined in get_tavref.')
      }
  }
}

//------------------------------------------------------------------------------c
void MEAM::get_sijk(double C, int i, int j, int k, double* sijk)
{
  double x;
  x = (C - Cmin_meam[i][j][k]) / (Cmax_meam[i][j][k] - Cmin_meam[i][j][k]);
  *sijk = fcut(x);
}

//------------------------------------------------------------------------------c
// Calculate density functions, assuming reference configuration
void MEAM::get_densref(double r, int a, int b, double* rho01, double* rho11, double* rho21, double* rho31,
                       double* rho02, double* rho12, double* rho22, double* rho32,
                       double* rho1m1, double* rho2m1, double* rho3m1,
                       double* rho1m2, double* rho2m2, double* rho3m2)
{
  double a1, a2;
  double s[3];
  lattice_t lat;
  int Zij,Zij2nn;
  double rhoa01nn, rhoa02nn;
  double rhoa01, rhoa11, rhoa21, rhoa31;
  double rhoa02, rhoa12, rhoa22, rhoa32;
  double arat, scrn, denom;
  double C, s111, s112, s221, S11, S22;
  // msmeam
  double rhoa1m1, rhoa2m1, rhoa3m1, rhoa1m2, rhoa2m2, rhoa3m2;

  a1 = r / re_meam[a][a] - 1.0;
  a2 = r / re_meam[b][b] - 1.0;

  rhoa01 = rho0_meam[a] * MathSpecial::fm_exp(-beta0_meam[a] * a1);

  if (msmeamflag) {
    // the rho variables are multiplied by t here since ialloy not needed in msmeam
    rhoa11 = rho0_meam[a] * t1_meam[a] * MathSpecial::fm_exp(-beta1_meam[a] * a1);
    rhoa21 = rho0_meam[a] * t2_meam[a] * MathSpecial::fm_exp(-beta2_meam[a] * a1);
    rhoa31 = rho0_meam[a] * t3_meam[a] * MathSpecial::fm_exp(-beta3_meam[a] * a1);
    rhoa02 = rho0_meam[b] * MathSpecial::fm_exp(-beta0_meam[b] * a2);
    rhoa12 = rho0_meam[b] * t1_meam[b] * MathSpecial::fm_exp(-beta1_meam[b] * a2);
    rhoa22 = rho0_meam[b] * t2_meam[b] * MathSpecial::fm_exp(-beta2_meam[b] * a2);
    rhoa32 = rho0_meam[b] * t3_meam[b] * MathSpecial::fm_exp(-beta3_meam[b] * a2);
    // msmeam specific rho vars
    rhoa1m1 = rho0_meam[a] * t1m_meam[a] * MathSpecial::fm_exp(-beta1m_meam[a] * a1);
    rhoa2m1 = rho0_meam[a] * t2m_meam[a] * MathSpecial::fm_exp(-beta2m_meam[a] * a1);
    rhoa3m1 = rho0_meam[a] * t3m_meam[a] * MathSpecial::fm_exp(-beta3m_meam[a] * a1);
    rhoa1m2 = rho0_meam[b] * t1m_meam[b] * MathSpecial::fm_exp(-beta1m_meam[b] * a2);
    rhoa2m2 = rho0_meam[b] * t2m_meam[b] * MathSpecial::fm_exp(-beta2m_meam[b] * a2);
    rhoa3m2 = rho0_meam[b] * t3m_meam[b] * MathSpecial::fm_exp(-beta3m_meam[b] * a2);
  } else {
    rhoa11 = rho0_meam[a] * MathSpecial::fm_exp(-beta1_meam[a] * a1);
    rhoa21 = rho0_meam[a] * MathSpecial::fm_exp(-beta2_meam[a] * a1);
    rhoa31 = rho0_meam[a] * MathSpecial::fm_exp(-beta3_meam[a] * a1);
    rhoa02 = rho0_meam[b] * MathSpecial::fm_exp(-beta0_meam[b] * a2);
    rhoa12 = rho0_meam[b] * MathSpecial::fm_exp(-beta1_meam[b] * a2);
    rhoa22 = rho0_meam[b] * MathSpecial::fm_exp(-beta2_meam[b] * a2);
    rhoa32 = rho0_meam[b] * MathSpecial::fm_exp(-beta3_meam[b] * a2);
  }

  lat = lattce_meam[a][b];

  Zij = get_Zij(lat);

  *rho11 = 0.0;
  *rho21 = 0.0;
  *rho31 = 0.0;
  *rho12 = 0.0;
  *rho22 = 0.0;
  *rho32 = 0.0;
  if (msmeamflag) {
    *rho1m1 = 0.0;
    *rho2m1 = 0.0;
    *rho3m1 = 0.0;
    *rho1m2 = 0.0;
    *rho2m2 = 0.0;
    *rho3m2 = 0.0;
  }

  // keep track of density components separately; combine in the calling subroutine
  switch (lat) {
    case FCC:
      *rho01 = 12.0 * rhoa02;
      *rho02 = 12.0 * rhoa01;
      break;
    case BCC:
      *rho01 = 8.0 * rhoa02;
      *rho02 = 8.0 * rhoa01;
      break;
    case B1:
    case SC:
      *rho01 = 6.0 * rhoa02;
      *rho02 = 6.0 * rhoa01;
      break;
    case DIA:
    case DIA3:
      *rho01 = 4.0 * rhoa02;
      *rho02 = 4.0 * rhoa01;
      *rho31 = 32.0 / 9.0 * rhoa32 * rhoa32;
      *rho32 = 32.0 / 9.0 * rhoa31 * rhoa31;
      if (msmeamflag) {
        *rho3m1 = 32.0 / 9.0 * rhoa3m2 * rhoa3m2;
        *rho3m2 = 32.0 / 9.0 * rhoa3m1 * rhoa3m1;
      }
      break;
    case HCP:
      *rho01 = 12 * rhoa02;
      *rho02 = 12 * rhoa01;
      *rho31 = 1.0 / 3.0 * rhoa32 * rhoa32;
      *rho32 = 1.0 / 3.0 * rhoa31 * rhoa31;
      if (msmeamflag) {
        *rho3m1 = 1.0 / 3.0 * rhoa3m2 * rhoa3m2;
        *rho3m2 = 1.0 / 3.0 * rhoa3m1 * rhoa3m1;
      }
      break;
    case DIM:
      get_shpfcn(DIM, 0, 0, s);
      *rho01 = rhoa02;
      *rho02 = rhoa01;
      *rho11 = s[0] * rhoa12 * rhoa12;
      *rho12 = s[0] * rhoa11 * rhoa11;
      *rho21 = s[1] * rhoa22 * rhoa22;
      *rho22 = s[1] * rhoa21 * rhoa21;
      *rho31 = s[2] * rhoa32 * rhoa32;
      *rho32 = s[2] * rhoa31 * rhoa31;
      if (msmeamflag) {
        *rho1m1 = s[0] * rhoa1m2 * rhoa1m2;
        *rho1m2 = s[0] * rhoa1m1 * rhoa1m1;
        *rho2m1 = s[1] * rhoa2m2 * rhoa2m2;
        *rho2m2 = s[1] * rhoa2m1 * rhoa2m1;
        *rho3m1 = s[2] * rhoa3m2 * rhoa3m2;
        *rho3m2 = s[2] * rhoa3m1 * rhoa3m1;
      }
      break;
    case C11:
      *rho01 = rhoa01;
      *rho02 = rhoa02;
      *rho11 = rhoa11;
      *rho12 = rhoa12;
      *rho21 = rhoa21;
      *rho22 = rhoa22;
      *rho31 = rhoa31;
      *rho32 = rhoa32;
      if (msmeamflag) {
        *rho1m1 = rhoa1m1;
        *rho1m2 = rhoa1m2;
        *rho2m1 = rhoa2m1;
        *rho2m2 = rhoa2m2;
        *rho3m1 = rhoa3m1;
        *rho3m2 = rhoa3m2;
      }
      break;
    case L12:
      *rho01 = 8 * rhoa01 + 4 * rhoa02;
      *rho02 = 12 * rhoa01;
      if (ialloy ==1){
        *rho21 = 8. / 3. * MathSpecial::square(rhoa21 * t2_meam[a] - rhoa22 * t2_meam[b]);
        denom = 8 * rhoa01 * MathSpecial::square(t2_meam[a]) + 4 * rhoa02 * MathSpecial::square(t2_meam[b]);
        if (denom > 0.)
          *rho21 = *rho21 / denom * *rho01;
      } else
        *rho21 = 8. / 3. * (rhoa21 - rhoa22) * (rhoa21 - rhoa22);
      if (msmeamflag) {
        *rho2m1 = 8. / 3. * (rhoa2m1 - rhoa2m2) * (rhoa2m1 - rhoa2m2);
      }
      break;
    case B2:
      *rho01 = 8.0 * rhoa02;
      *rho02 = 8.0 * rhoa01;
      break;
    case CH4:
      *rho01 = 4.0 * rhoa02; //in assumption that 'a' represent carbon
      *rho02 = rhoa01;       //in assumption that 'b' represent hydrogen

      get_shpfcn(DIM, 0, 0, s); //H
      *rho12 = s[0] * rhoa11 * rhoa11;
      *rho22 = s[1] * rhoa21 * rhoa21;
      *rho32 = s[2] * rhoa31 * rhoa31;

      get_shpfcn(CH4, 0, 0, s); //C
      *rho11 = s[0] * rhoa12 * rhoa12;
      *rho21 = s[1] * rhoa22 * rhoa22;
      *rho31 = s[2] * rhoa32 * rhoa32;
      break;
    case LIN:
      *rho01 = rhoa02*Zij;
      *rho02 = rhoa01*Zij;

      get_shpfcn(LIN, stheta_meam[a][b], ctheta_meam[a][b], s);
      *rho12 = s[0] * rhoa11 * rhoa11;
      *rho22 = s[1] * rhoa21 * rhoa21;
      *rho32 = s[2] * rhoa31 * rhoa31;
      *rho11 = s[0] * rhoa12 * rhoa12;
      *rho21 = s[1] * rhoa22 * rhoa22;
      *rho31 = s[2] * rhoa32 * rhoa32;
      break;
    case ZIG:
      *rho01 = rhoa02*Zij;
      *rho02 = rhoa01*Zij;

      get_shpfcn(ZIG, stheta_meam[a][b], ctheta_meam[a][b], s);
      *rho12 = s[0] * rhoa11 * rhoa11;
      *rho22 = s[1] * rhoa21 * rhoa21;
      *rho32 = s[2] * rhoa31 * rhoa31;
      *rho11 = s[0] * rhoa12 * rhoa12;
      *rho21 = s[1] * rhoa22 * rhoa22;
      *rho31 = s[2] * rhoa32 * rhoa32;
      break;
    case TRI:
      *rho01 = rhoa02;
      *rho02 = rhoa01*Zij;

      get_shpfcn(TRI, stheta_meam[a][b], ctheta_meam[a][b], s);
      *rho12 = s[0] * rhoa11 * rhoa11;
      *rho22 = s[1] * rhoa21 * rhoa21;
      *rho32 = s[2] * rhoa31 * rhoa31;
      s[0] = 1.0;
      s[1] = 2.0/3.0;
      s[2] = 1.0 - 0.6*s[0];

      *rho11 = s[0] * rhoa12 * rhoa12;
      *rho21 = s[1] * rhoa22 * rhoa22;
      *rho31 = s[2] * rhoa32 * rhoa32;
      break;


    // default:
    //        call error('Lattice not defined in get_densref.')
  }

  if (nn2_meam[a][b] == 1) {


    Zij2nn = get_Zij2(lat, Cmin_meam[a][a][b], Cmax_meam[a][a][b],
                      stheta_meam[a][b], arat, scrn);

    a1 = arat * r / re_meam[a][a] - 1.0;
    a2 = arat * r / re_meam[b][b] - 1.0;

    rhoa01nn = rho0_meam[a] * MathSpecial::fm_exp(-beta0_meam[a] * a1);
    rhoa02nn = rho0_meam[b] * MathSpecial::fm_exp(-beta0_meam[b] * a2);

    if (lat == L12) {
      //     As usual, L12 thinks it's special; we need to be careful computing
      //     the screening functions
      C = 1.0;
      get_sijk(C, a, a, a, &s111);
      get_sijk(C, a, a, b, &s112);
      get_sijk(C, b, b, a, &s221);
      S11 = s111 * s111 * s112 * s112;
      S22 = s221 * s221 * s221 * s221;
      *rho01 = *rho01 + 6 * S11 * rhoa01nn;
      *rho02 = *rho02 + 6 * S22 * rhoa02nn;

    } else {
      //     For other cases, assume that second neighbor is of same type,
      //     first neighbor may be of different type

      *rho01 = *rho01 + Zij2nn * scrn * rhoa01nn;

      //     Assume Zij2nn and arat don't depend on order, but scrn might
      Zij2nn = get_Zij2(lat, Cmin_meam[b][b][a], Cmax_meam[b][b][a],
                        stheta_meam[a][b], arat, scrn);
      *rho02 = *rho02 + Zij2nn * scrn * rhoa02nn;
    }
  }
}


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
