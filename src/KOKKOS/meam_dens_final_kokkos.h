#include "meam_kokkos.h"
#include "math_special.h"

using namespace LAMMPS_NS;

template<class DeviceType>
void
MEAMKokkos<DeviceType>::meam_dens_final(int nlocal, int eflag_either, int eflag_global, int eflag_atom, double* eng_vdwl,
                      typename ArrayTypes<DeviceType>::t_efloat_1d eatom, int ntype, typename AT::t_int_1d_randomread type, typename AT::t_int_1d_randomread fmap, int& errorflag)
{
  EV_FLOAT ev;
  ev.evdwl = *eng_vdwl;
  this->eflag_either = eflag_either;
  this->eflag_global = eflag_global;
  this->eflag_atom = eflag_atom;
  this->d_eatom = eatom;
  this->ntype = ntype;
  this->type = type;
  this->fmap = fmap;

  //     Complete the calculation of density

  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagMEAMDensFinal>(0,nlocal),*this,ev);
  *eng_vdwl = ev.evdwl;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void MEAMKokkos<DeviceType>::operator()(TagMEAMDensFinal, const int &i, EV_FLOAT& ev) const {

  lattice_t elti;
  int m;
  int errorflag;
  F_FLOAT rhob, G, dG, Gbar, dGbar, gam, shp[3], Z;
  F_FLOAT B, denom, rho_bkgd;

  // may or may not be legal
    elti = static_cast<lattice_t>(fmap[type[i]]);
    if (elti >= 0) {
      d_rho1[i] = 0.0;
      d_rho2[i] = -1.0 / 3.0 * d_arho2b[i] * d_arho2b[i];
      d_rho3[i] = 0.0;
      for (m = 0; m < 3; m++) {
        d_rho1[i] = d_rho1[i] + d_arho1(i,m) * d_arho1(i,m);
        d_rho3[i] = d_rho3[i] - 3.0 / 5.0 * d_arho3b(i,m) * d_arho3b(i,m);
      }
      for (m = 0; m < 6; m++) {
        d_rho2[i] = d_rho2[i] + v2D[m] * d_arho2(i,m) * d_arho2(i,m);
      }
      for (m = 0; m < 10; m++) {
        d_rho3[i] = d_rho3[i] + v3D[m] * d_arho3(i,m) * d_arho3(i,m);
      }

      if (d_rho0[i] > 0.0) {
        if (this->ialloy == 1) {
          d_t_ave(i,0) = fdiv_zero_kk(d_t_ave(i,0), d_tsq_ave(i,0));
          d_t_ave(i,1) = fdiv_zero_kk(d_t_ave(i,1), d_tsq_ave(i,1));
          d_t_ave(i,2) = fdiv_zero_kk(d_t_ave(i,2), d_tsq_ave(i,2));
        } else if (this->ialloy == 2) {
          d_t_ave(i,0) = this->t1_meam[elti];
          d_t_ave(i,1) = this->t2_meam[elti];
          d_t_ave(i,2) = this->t3_meam[elti];
        } else {
          d_t_ave(i,0) = d_t_ave(i,0) / d_rho0[i];
          d_t_ave(i,1) = d_t_ave(i,1) / d_rho0[i];
          d_t_ave(i,2) = d_t_ave(i,2) / d_rho0[i];
        }
      }

      d_gamma[i] = d_t_ave(i,0) * d_rho1[i] + d_t_ave(i,1) * d_rho2[i] + d_t_ave(i,2) * d_rho3[i];

      if (d_rho0[i] > 0.0) {
        d_gamma[i] = d_gamma[i] / (d_rho0[i] * d_rho0[i]);
      }

      // need to double check, the analogous function is
      // Z = get_Zij(this->lattce_meam[elti][elti]); in the non-KOKKOS version?
      Z = get_Zij(elti);

      G = G_gam(d_gamma[i], this->ibar_meam[elti], errorflag);
      if (errorflag != 0)
      {
        //char str[128];
        //sprintf(str,"MEAMKokkos library error %d",errorflag);
        //error->one(FLERR,str);
        return;
      }
      get_shpfcn(this->lattce_meam[elti][elti], shp);
      if (this->ibar_meam[elti] <= 0) {
        Gbar = 1.0;
        dGbar = 0.0;
      } else {
        if (this->mix_ref_t == 1) {
          gam = (d_t_ave(i,0) * shp[0] + d_t_ave(i,1) * shp[1] + d_t_ave(i,2) * shp[2]) / (Z * Z);
        } else {
          gam = (this->t1_meam[elti] * shp[0] + this->t2_meam[elti] * shp[1] + this->t3_meam[elti] * shp[2]) /
                (Z * Z);
        }
        Gbar = G_gam(gam, this->ibar_meam[elti], errorflag);
      }
      d_rho[i] = d_rho0[i] * G;

      if (this->mix_ref_t == 1) {
        if (this->ibar_meam[elti] <= 0) {
          Gbar = 1.0;
          dGbar = 0.0;
        } else {
          gam = (d_t_ave(i,0) * shp[0] + d_t_ave(i,1) * shp[1] + d_t_ave(i,2) * shp[2]) / (Z * Z);
          Gbar = dG_gam(gam, this->ibar_meam[elti], dGbar);
        }
        rho_bkgd = this->rho0_meam[elti] * Z * Gbar;
      } else {
        if (this->bkgd_dyn == 1) {
          rho_bkgd = this->rho0_meam[elti] * Z;
        } else {
          rho_bkgd = this->rho_ref_meam[elti];
        }
      }
      rhob = d_rho[i] / rho_bkgd;
      denom = 1.0 / rho_bkgd;

      G = dG_gam(d_gamma[i], this->ibar_meam[elti], dG);

      d_dgamma1[i] = (G - 2 * dG * d_gamma[i]) * denom;

      if (!iszero_kk(d_rho0[i])) {
        d_dgamma2[i] = (dG / d_rho0[i]) * denom;
      } else {
        d_dgamma2[i] = 0.0;
      }

      //     dgamma3 is nonzero only if we are using the "mixed" rule for
      //     computing t in the reference system (which is not correct, but
      //     included for backward compatibility
      if (this->mix_ref_t == 1) {
        d_dgamma3[i] = d_rho0[i] * G * dGbar / (Gbar * Z * Z) * denom;
      } else {
        d_dgamma3[i] = 0.0;
      }

      B = this->A_meam[elti] * this->Ec_meam[elti][elti];

      if (!iszero_kk(rhob)) {
        if (this->emb_lin_neg == 1 && rhob <= 0) {
          d_frhop[i] = -B;
        } else {
          d_frhop[i] = B * (log(rhob) + 1.0);
        }
        if (eflag_either != 0) {
          if (eflag_global != 0) {
            if (this->emb_lin_neg == 1 && rhob <= 0) {
              //*eng_vdwl = *eng_vdwl - B * rhob;
              ev.evdwl = ev.evdwl - B * rhob;
            } else {
              //*eng_vdwl = *eng_vdwl + B * rhob * log(rhob);
              ev.evdwl = ev.evdwl + B * rhob * log(rhob);
            }
          }
          if (eflag_atom != 0) {
            if (this->emb_lin_neg == 1 && rhob <= 0) {
              d_eatom[i] = d_eatom[i] - B * rhob;
            } else {
              d_eatom[i] = d_eatom[i] + B * rhob * log(rhob);
            }
          }
        }
      } else {
        if (this->emb_lin_neg == 1) {
          d_frhop[i] = -B;
        } else {
          d_frhop[i] = B;
        }
      }
    }
    if (errorflag)
    {
        //char str[128];
        //sprintf(str,"MEAMKokkos library error %d",errorflag);
        //error->one(FLERR,str);
    }
}
