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

#include "meam_kokkos.h"
#include "math_special.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void
MEAMKokkos<DeviceType>::meam_dens_final(int nlocal, int eflag_either, int eflag_global, int eflag_atom,
                      typename ArrayTypes<DeviceType>::t_efloat_1d eatom, int ntype, typename AT::t_int_1d type, typename AT::t_int_1d d_map, typename AT::t_int_2d d_scale, int& errorflag, EV_FLOAT &ev_all)
{
  EV_FLOAT ev;
  this->eflag_either = eflag_either;
  this->eflag_global = eflag_global;
  this->eflag_atom = eflag_atom;
  this->d_eatom = eatom;
  this->ntype = ntype;
  this->type = type;
  this->d_map = d_map;
  this->d_scale = d_scale;

  Kokkos::deep_copy(d_errorflag,0);

  // Complete the calculation of density

  copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagMEAMDensFinal>(0,nlocal),*this,ev);
  ev_all.evdwl += ev.evdwl;
  copymode = 0;

  auto h_errorflag = Kokkos::create_mirror_view_and_copy(LMPHostType(),d_errorflag);
  errorflag = h_errorflag();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void MEAMKokkos<DeviceType>::operator()(TagMEAMDensFinal, const int &i, EV_FLOAT& ev) const {

  F_FLOAT rhob, G, dG, Gbar, dGbar, gam, shp[3], Z;
  F_FLOAT denom, rho_bkgd, Fl;
  double scaleii;

  int elti = d_map[type[i]];
  if (elti >= 0) {
    scaleii = d_scale(type[i],type[i]);
    d_rho1[i] = 0.0;
    if (msmeamflag) {
      d_rho2[i] = -1.0 / 3.0 * (d_arho2b[i] * d_arho2b[i]
                              - d_arho2mb[i] * d_arho2mb[i]);
    } else{
      d_rho2[i] = -1.0 / 3.0 * d_arho2b[i] * d_arho2b[i];
    }
    d_rho3[i] = 0.0;
    for (int m = 0; m < 3; m++) {
      if (msmeamflag) {
        d_rho1[i] = d_rho1[i] + d_arho1(i, m) * d_arho1(i, m)
                             - d_arho1m(i, m) * d_arho1m(i, m);
        d_rho3[i] = d_rho3[i] - 3.0 / 5.0 * (d_arho3b(i, m) * d_arho3b(i, m)
                                          - d_arho3mb(i, m) * d_arho3mb(i, m));
      } else{
        d_rho1[i] += d_arho1(i,m) * d_arho1(i,m);
        d_rho3[i] -= 3.0 / 5.0 * d_arho3b(i,m) * d_arho3b(i,m);
      }
    }
    for (int m = 0; m < 6; m++){
      if (msmeamflag) {
        d_rho2[i] = d_rho2[i] + v2D[m] * (d_arho2(i, m) * d_arho2(i, m)
                                         - d_arho2m(i, m) * d_arho2m(i, m));
      } else{
        d_rho2[i] += v2D[m] * d_arho2(i,m) * d_arho2(i,m);
      }
    }
    for (int m = 0; m < 10; m++)
      if (msmeamflag) {
        d_rho3[i] = d_rho3[i] + v3D[m] * (d_arho3(i, m) * d_arho3(i, m)
                                        - d_arho3m(i, m) * d_arho3m(i, m));
      } else{
        d_rho3[i] += v3D[m] * d_arho3(i,m) * d_arho3(i,m);
      }

    if (msmeamflag) {
      // with msmeam all t weights are already accounted for in rho
      d_gamma[i] = d_rho1[i] + d_rho2[i] + d_rho3[i];
    } else{
      if (d_rho0[i] > 0.0) {
        if (ialloy == 1) {
          d_t_ave(i,0) = fdiv_zero_kk(d_t_ave(i,0), d_tsq_ave(i,0));
          d_t_ave(i,1) = fdiv_zero_kk(d_t_ave(i,1), d_tsq_ave(i,1));
          d_t_ave(i,2) = fdiv_zero_kk(d_t_ave(i,2), d_tsq_ave(i,2));
        } else if (ialloy == 2) {
          d_t_ave(i,0) = t1_meam[elti];
          d_t_ave(i,1) = t2_meam[elti];
          d_t_ave(i,2) = t3_meam[elti];
        } else {
          d_t_ave(i,0) /= d_rho0[i];
          d_t_ave(i,1) /= d_rho0[i];
          d_t_ave(i,2) /= d_rho0[i];
        }
      }
      d_gamma[i] = d_t_ave(i,0) * d_rho1[i] + d_t_ave(i,1) * d_rho2[i] + d_t_ave(i,2) * d_rho3[i];
    }

    if (d_rho0[i] > 0.0)
      d_gamma[i] /= (d_rho0[i] * d_rho0[i]);

    Z = get_Zij(lattce_meam[elti][elti]);

    G = G_gam(d_gamma[i], ibar_meam[elti], d_errorflag());
    if (d_errorflag() != 0)
      return;

    get_shpfcn(lattce_meam[elti][elti], stheta_meam[elti][elti], ctheta_meam[elti][elti], shp);
    if (ibar_meam[elti] <= 0) {
      Gbar = 1.0;
      dGbar = 0.0;
    } else {
      if (mix_ref_t == 1)
        gam = (d_t_ave(i,0) * shp[0] + d_t_ave(i,1) * shp[1] + d_t_ave(i,2) * shp[2]) / (Z * Z);
      else
        gam = (t1_meam[elti] * shp[0] + t2_meam[elti] * shp[1] + t3_meam[elti] * shp[2]) /
              (Z * Z);
      Gbar = G_gam(gam, ibar_meam[elti], d_errorflag());
    }
    d_rho[i] = d_rho0[i] * G;

    if (mix_ref_t == 1) {
      if (ibar_meam[elti] <= 0) {
        Gbar = 1.0;
        dGbar = 0.0;
      } else {
        gam = (d_t_ave(i,0) * shp[0] + d_t_ave(i,1) * shp[1] + d_t_ave(i,2) * shp[2]) / (Z * Z);
        Gbar = dG_gam(gam, ibar_meam[elti], dGbar);
      }
      rho_bkgd = rho0_meam[elti] * Z * Gbar;
    } else {
      if (bkgd_dyn == 1)
        rho_bkgd = rho0_meam[elti] * Z;
      else
        rho_bkgd = rho_ref_meam[elti];
    }
    rhob = d_rho[i] / rho_bkgd;
    denom = 1.0 / rho_bkgd;

    G = dG_gam(d_gamma[i], ibar_meam[elti], dG);

    d_dgamma1[i] = (G - 2 * dG * d_gamma[i]) * denom;

    if (!iszero_kk(d_rho0[i]))
      d_dgamma2[i] = (dG / d_rho0[i]) * denom;
    else
      d_dgamma2[i] = 0.0;

    // dgamma3 is nonzero only if we are using the "mixed" rule for
    // computing t in the reference system (which is not correct, but
    // included for backward compatibility
    if (mix_ref_t == 1)
      d_dgamma3[i] = d_rho0[i] * G * dGbar / (Gbar * Z * Z) * denom;
    else
      d_dgamma3[i] = 0.0;

    Fl = embedding(A_meam[elti], Ec_meam[elti][elti], rhob, d_frhop[i]);

    if (eflag_either) {
      Fl *= scaleii;
      if (eflag_global) {
        ev.evdwl += Fl;
      }
      if (eflag_atom) {
        d_eatom[i] += Fl;
      }
    }
  }
}

