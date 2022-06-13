#include "meam_kokkos.h"
#include "math_special_kokkos.h"
#include <algorithm>

using namespace LAMMPS_NS;
using namespace MathSpecialKokkos;


template<class DeviceType>
void
MEAMKokkos<DeviceType>::meam_force(int inum_half, int eflag_either, int eflag_global, int eflag_atom, int vflag_atom, double* eng_vdwl,
                 typename ArrayTypes<DeviceType>::t_efloat_1d eatom, int ntype, typename AT::t_int_1d_randomread type, typename AT::t_int_1d_randomread fmap, typename AT::t_x_array_randomread x, typename AT::t_int_1d_randomread numneigh,
                 typename AT::t_int_1d_randomread numneigh_full, typename AT::t_f_array f, typename ArrayTypes<DeviceType>::t_virial_array vatom, typename AT::t_int_1d_randomread d_ilist_half, typename AT::t_int_1d_randomread d_offset, typename AT::t_neighbors_2d d_neighbors_half, typename AT::t_neighbors_2d d_neighbors_full, int neighflag)
{
    EV_FLOAT ev;
    ev.evdwl = *eng_vdwl;

    this->eflag_either = eflag_either;
    this->eflag_global = eflag_global;
    this->eflag_atom = eflag_atom;
    this->vflag_atom = vflag_atom;
    this->d_eatom = eatom;
    this->ntype = ntype;
    this->type = type;
    this->fmap = fmap;
    this->x = x;
    this->d_numneigh_half = numneigh;
    this->d_numneigh_full = numneigh_full;
    this->d_neighbors_half = d_neighbors_half;
    this->d_neighbors_full = d_neighbors_full;
    this->f = f;
    this->d_vatom = vatom;
    this->d_ilist_half = d_ilist_half;
    this->d_offset = d_offset;
  
    if (neighflag == FULL) 
       Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagMEAMforce<FULL>>(0,inum_half),*this,ev);
    else if (neighflag == HALF)
       Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagMEAMforce<HALF>>(0,inum_half),*this,ev);
    else if (neighflag == HALFTHREAD)
       Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagMEAMforce<HALFTHREAD>>(0,inum_half),*this,ev);
    *eng_vdwl = ev.evdwl;
}

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void
MEAMKokkos<DeviceType>::operator()(TagMEAMforce<NEIGHFLAG>, const int &ii, EV_FLOAT& ev) const {
  int i, j, jn, k, kn, kk, m, n, p, q;
  int nv2, nv3, elti, eltj, eltk, ind;
  X_FLOAT xitmp, yitmp, zitmp, delij[3];
  F_FLOAT rij2, rij, rij3;
  F_FLOAT v[6], fi[3], fj[3];
  F_FLOAT third, sixth;
  F_FLOAT pp, dUdrij, dUdsij, dUdrijm[3], force, forcem;
  F_FLOAT r, recip, phi, phip;
  F_FLOAT sij;
  F_FLOAT a1, a1i, a1j, a2, a2i, a2j;
  F_FLOAT a3i, a3j;
  F_FLOAT shpi[3], shpj[3];
  F_FLOAT ai, aj, ro0i, ro0j, invrei, invrej;
  F_FLOAT rhoa0j, drhoa0j, rhoa0i, drhoa0i;
  F_FLOAT rhoa1j, drhoa1j, rhoa1i, drhoa1i;
  F_FLOAT rhoa2j, drhoa2j, rhoa2i, drhoa2i;
  F_FLOAT a3, a3a, rhoa3j, drhoa3j, rhoa3i, drhoa3i;
  F_FLOAT drho0dr1, drho0dr2, drho0ds1, drho0ds2;
  F_FLOAT drho1dr1, drho1dr2, drho1ds1, drho1ds2;
  F_FLOAT drho1drm1[3], drho1drm2[3];
  F_FLOAT drho2dr1, drho2dr2, drho2ds1, drho2ds2;
  F_FLOAT drho2drm1[3], drho2drm2[3];
  F_FLOAT drho3dr1, drho3dr2, drho3ds1, drho3ds2;
  F_FLOAT drho3drm1[3], drho3drm2[3];
  F_FLOAT dt1dr1, dt1dr2, dt1ds1, dt1ds2;
  F_FLOAT dt2dr1, dt2dr2, dt2ds1, dt2ds2;
  F_FLOAT dt3dr1, dt3dr2, dt3ds1, dt3ds2;
  F_FLOAT drhodr1, drhodr2, drhods1, drhods2, drhodrm1[3], drhodrm2[3];
  F_FLOAT arg;
  F_FLOAT arg1i1, arg1j1, arg1i2, arg1j2, arg1i3, arg1j3, arg3i3, arg3j3;
  F_FLOAT dsij1, dsij2, force1, force2;
  F_FLOAT t1i, t2i, t3i, t1j, t2j, t3j;
  int fnoffset;

  // To do: update with duplication for OpenMP, atomic for GPUs
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_eatom = d_eatom;
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_vatom = d_vatom;
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;


  i = d_ilist_half[ii];
  fnoffset = d_offset[i];
  third = 1.0 / 3.0;
  sixth = 1.0 / 6.0;
  //     Compute forces atom i

  elti = fmap[type[i]];
  if (elti < 0) return;

  xitmp = x(i,0);
  yitmp = x(i,1);
  zitmp = x(i,2);

  //     Treat each pair
  for (jn = 0; jn < d_numneigh_half[i]; jn++) {
    //j = firstneigh[jn];
    j = d_neighbors_half(i,jn);
    eltj = fmap[type[j]];

    if (!iszero_kk(d_scrfcn[fnoffset + jn]) && eltj >= 0) {

      sij = d_scrfcn[fnoffset + jn] * d_fcpair[fnoffset + jn];
      delij[0] = x(j,0) - xitmp;
      delij[1] = x(j,1) - yitmp;
      delij[2] = x(j,2) - zitmp;
      rij2 = delij[0] * delij[0] + delij[1] * delij[1] + delij[2] * delij[2];
      if (rij2 < this->cutforcesq) {
        rij = sqrt(rij2);
        r = rij;

        //     Compute phi and phip
        ind = this->eltind[elti][eltj];
        pp = rij * this->rdrar;
        kk = (int)pp;
        kk = (kk <= (this->nrar - 2))?kk:this->nrar - 2;
        pp = pp - kk;
        pp = (pp <= 1.0)?pp:1.0;
        phi = ((this->d_phirar3(ind,kk) * pp + this->d_phirar2(ind,kk)) * pp + this->d_phirar1(ind,kk)) * pp +
          this->d_phirar(ind,kk);
        phip = (this->d_phirar6(ind,kk) * pp + this->d_phirar5(ind,kk)) * pp + this->d_phirar4(ind,kk);
        recip = 1.0 / r;

        if (eflag_either != 0) {
          if (eflag_global != 0)
            //*eng_vdwl = *eng_vdwl + phi * sij;
            ev.evdwl = ev.evdwl + phi * sij;
          if (eflag_atom != 0) {
            a_eatom[i] += (0.5 * phi * sij);
            a_eatom[j] += (0.5 * phi * sij);
          }
        }

        //     write(1,*) "force_meamf: phi: ",phi
        //     write(1,*) "force_meamf: phip: ",phip

        //     Compute pair densities and derivatives
        invrei = 1.0 / this->re_meam[elti][elti];
        ai = rij * invrei - 1.0;
        ro0i = this->rho0_meam[elti];
        rhoa0i = ro0i * MathSpecialKokkos::fm_exp(-this->beta0_meam[elti] * ai);
        drhoa0i = -this->beta0_meam[elti] * invrei * rhoa0i;
        rhoa1i = ro0i * MathSpecialKokkos::fm_exp(-this->beta1_meam[elti] * ai);
        drhoa1i = -this->beta1_meam[elti] * invrei * rhoa1i;
        rhoa2i = ro0i * MathSpecialKokkos::fm_exp(-this->beta2_meam[elti] * ai);
        drhoa2i = -this->beta2_meam[elti] * invrei * rhoa2i;
        rhoa3i = ro0i * MathSpecialKokkos::fm_exp(-this->beta3_meam[elti] * ai);
        drhoa3i = -this->beta3_meam[elti] * invrei * rhoa3i;

        if (elti != eltj) {
          invrej = 1.0 / this->re_meam[eltj][eltj];
          aj = rij * invrej - 1.0;
          ro0j = this->rho0_meam[eltj];
          rhoa0j = ro0j * MathSpecialKokkos::fm_exp(-this->beta0_meam[eltj] * aj);
          drhoa0j = -this->beta0_meam[eltj] * invrej * rhoa0j;
          rhoa1j = ro0j * MathSpecialKokkos::fm_exp(-this->beta1_meam[eltj] * aj);
          drhoa1j = -this->beta1_meam[eltj] * invrej * rhoa1j;
          rhoa2j = ro0j * MathSpecialKokkos::fm_exp(-this->beta2_meam[eltj] * aj);
          drhoa2j = -this->beta2_meam[eltj] * invrej * rhoa2j;
          rhoa3j = ro0j * MathSpecialKokkos::fm_exp(-this->beta3_meam[eltj] * aj);
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
              arg1i3 = arg1i3 + d_arho3(i,nv3) * arg;
              arg1j3 = arg1j3 - d_arho3(j,nv3) * arg;
              nv3 = nv3 + 1;
            }
            arg = delij[n] * delij[p] * this->v2D[nv2];
            arg1i2 = arg1i2 + d_arho2(i,nv2) * arg;
            arg1j2 = arg1j2 + d_arho2(j,nv2) * arg;
            nv2 = nv2 + 1;
          }
          arg1i1 = arg1i1 + d_arho1(i,n) * delij[n];
          arg1j1 = arg1j1 - d_arho1(j,n) * delij[n];
          arg3i3 = arg3i3 + d_arho3b(i,n) * delij[n];
          arg3j3 = arg3j3 - d_arho3b(j,n) * delij[n];
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
          drho1drm1[m] = a1 * rhoa1j * d_arho1(i,m);
          drho1drm2[m] = -a1 * rhoa1i * d_arho1(j,m);
        }

        //     rho2 terms
        a2 = 2 * sij / rij2;
        drho2dr1 = a2 * (drhoa2j - 2 * rhoa2j / rij) * arg1i2 - 2.0 / 3.0 * d_arho2b[i] * drhoa2j * sij;
        drho2dr2 = a2 * (drhoa2i - 2 * rhoa2i / rij) * arg1j2 - 2.0 / 3.0 * d_arho2b[j] * drhoa2i * sij;
        a2 = 4 * sij / rij2;
        for (m = 0; m < 3; m++) {
          drho2drm1[m] = 0.0;
          drho2drm2[m] = 0.0;
          for (n = 0; n < 3; n++) {
            drho2drm1[m] = drho2drm1[m] + d_arho2(i,this->vind2D[m][n]) * delij[n];
            drho2drm2[m] = drho2drm2[m] - d_arho2(j,this->vind2D[m][n]) * delij[n];
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
              drho3drm1[m] = drho3drm1[m] + d_arho3(i,this->vind3D[m][n][p]) * arg;
              drho3drm2[m] = drho3drm2[m] + d_arho3(j,this->vind3D[m][n][p]) * arg;
              nv2 = nv2 + 1;
            }
          }
          drho3drm1[m] = (a3 * drho3drm1[m] - a3a * d_arho3b(i,m)) * rhoa3j;
          drho3drm2[m] = (-a3 * drho3drm2[m] + a3a * d_arho3b(j,m)) * rhoa3i;
        }

        //     Compute derivatives of weighting functions t wrt rij
        t1i = d_t_ave(i,0);
        t2i = d_t_ave(i,1);
        t3i = d_t_ave(i,2);
        t1j = d_t_ave(j,0);
        t2j = d_t_ave(j,1);
        t3j = d_t_ave(j,2);

        if (this->ialloy == 1) {

          a1i = fdiv_zero_kk(drhoa0j * sij, d_tsq_ave(i,0));
          a1j = fdiv_zero_kk(drhoa0i * sij, d_tsq_ave(j,0));
          a2i = fdiv_zero_kk(drhoa0j * sij, d_tsq_ave(i,1));
          a2j = fdiv_zero_kk(drhoa0i * sij, d_tsq_ave(j,1));
          a3i = fdiv_zero_kk(drhoa0j * sij, d_tsq_ave(i,2));
          a3j = fdiv_zero_kk(drhoa0i * sij, d_tsq_ave(j,2));

          dt1dr1 = a1i * (t1mj - t1i * MathSpecialKokkos::square(t1mj));
          dt1dr2 = a1j * (t1mi - t1j * MathSpecialKokkos::square(t1mi));
          dt2dr1 = a2i * (t2mj - t2i * MathSpecialKokkos::square(t2mj));
          dt2dr2 = a2j * (t2mi - t2j * MathSpecialKokkos::square(t2mi));
          dt3dr1 = a3i * (t3mj - t3i * MathSpecialKokkos::square(t3mj));
          dt3dr2 = a3j * (t3mi - t3j * MathSpecialKokkos::square(t3mi));

        } else if (this->ialloy == 2) {

          dt1dr1 = 0.0;
          dt1dr2 = 0.0;
          dt2dr1 = 0.0;
          dt2dr2 = 0.0;
          dt3dr1 = 0.0;
          dt3dr2 = 0.0;

        } else {

          ai = 0.0;
          if (!iszero_kk(d_rho0[i]))
            ai = drhoa0j * sij / d_rho0[i];
          aj = 0.0;
          if (!iszero_kk(d_rho0[j]))
            aj = drhoa0i * sij / d_rho0[j];

          dt1dr1 = ai * (t1mj - t1i);
          dt1dr2 = aj * (t1mi - t1j);
          dt2dr1 = ai * (t2mj - t2i);
          dt2dr2 = aj * (t2mi - t2j);
          dt3dr1 = ai * (t3mj - t3i);
          dt3dr2 = aj * (t3mi - t3j);
        }

        //     Compute derivatives of total density wrt rij, sij and rij(3)
         get_shpfcn(this->lattce_meam[elti][elti], shpi);
         get_shpfcn(this->lattce_meam[eltj][eltj], shpj);
        drhodr1 = dgamma1[i] * drho0dr1 +
          dgamma2[i] * (dt1dr1 * d_rho1[i] + t1i * drho1dr1 + dt2dr1 * d_rho2[i] + t2i * drho2dr1 +
                        dt3dr1 * d_rho3[i] + t3i * drho3dr1) -
          dgamma3[i] * (shpi[0] * dt1dr1 + shpi[1] * dt2dr1 + shpi[2] * dt3dr1);
        drhodr2 = dgamma1[j] * drho0dr2 +
          dgamma2[j] * (dt1dr2 * d_rho1[j] + t1j * drho1dr2 + dt2dr2 * d_rho2[j] + t2j * drho2dr2 +
                        dt3dr2 * d_rho3[j] + t3j * drho3dr2) -
          dgamma3[j] * (shpj[0] * dt1dr2 + shpj[1] * dt2dr2 + shpj[2] * dt3dr2);
        for (m = 0; m < 3; m++) {
          drhodrm1[m] = 0.0;
          drhodrm2[m] = 0.0;
          drhodrm1[m] = dgamma2[i] * (t1i * drho1drm1[m] + t2i * drho2drm1[m] + t3i * drho3drm1[m]);
          drhodrm2[m] = dgamma2[j] * (t1j * drho1drm2[m] + t2j * drho2drm2[m] + t3j * drho3drm2[m]);
        }

        //     Compute derivatives wrt sij, but only if necessary
        if (!iszero_kk(d_dscrfcn[fnoffset + jn])) {
          drho0ds1 = rhoa0j;
          drho0ds2 = rhoa0i;
          a1 = 2.0 / rij;
          drho1ds1 = a1 * rhoa1j * arg1i1;
          drho1ds2 = a1 * rhoa1i * arg1j1;
          a2 = 2.0 / rij2;
          drho2ds1 = a2 * rhoa2j * arg1i2 - 2.0 / 3.0 * d_arho2b[i] * rhoa2j;
          drho2ds2 = a2 * rhoa2i * arg1j2 - 2.0 / 3.0 * d_arho2b[j] * rhoa2i;
          a3 = 2.0 / rij3;
          a3a = 6.0 / (5.0 * rij);
          drho3ds1 = a3 * rhoa3j * arg1i3 - a3a * rhoa3j * arg3i3;
          drho3ds2 = a3 * rhoa3i * arg1j3 - a3a * rhoa3i * arg3j3;

          if (this->ialloy == 1) {
            a1i = fdiv_zero_kk(rhoa0j, d_tsq_ave(i,0));
            a1j = fdiv_zero_kk(rhoa0i, d_tsq_ave(j,0));
            a2i = fdiv_zero_kk(rhoa0j, d_tsq_ave(i,1));
            a2j = fdiv_zero_kk(rhoa0i, d_tsq_ave(j,1));
            a3i = fdiv_zero_kk(rhoa0j, d_tsq_ave(i,2));
            a3j = fdiv_zero_kk(rhoa0i, d_tsq_ave(j,2));

            dt1ds1 = a1i * (t1mj - t1i * MathSpecialKokkos::square(t1mj));
            dt1ds2 = a1j * (t1mi - t1j * MathSpecialKokkos::square(t1mi));
            dt2ds1 = a2i * (t2mj - t2i * MathSpecialKokkos::square(t2mj));
            dt2ds2 = a2j * (t2mi - t2j * MathSpecialKokkos::square(t2mi));
            dt3ds1 = a3i * (t3mj - t3i * MathSpecialKokkos::square(t3mj));
            dt3ds2 = a3j * (t3mi - t3j * MathSpecialKokkos::square(t3mi));

          } else if (this->ialloy == 2) {

            dt1ds1 = 0.0;
            dt1ds2 = 0.0;
            dt2ds1 = 0.0;
            dt2ds2 = 0.0;
            dt3ds1 = 0.0;
            dt3ds2 = 0.0;

          } else {

            ai = 0.0;
            if (!iszero_kk(d_rho0[i]))
              ai = rhoa0j / d_rho0[i];
            aj = 0.0;
            if (!iszero_kk(d_rho0[j]))
              aj = rhoa0i / d_rho0[j];

            dt1ds1 = ai * (t1mj - t1i);
            dt1ds2 = aj * (t1mi - t1j);
            dt2ds1 = ai * (t2mj - t2i);
            dt2ds2 = aj * (t2mi - t2j);
            dt3ds1 = ai * (t3mj - t3i);
            dt3ds2 = aj * (t3mi - t3j);
          }

          drhods1 = d_dgamma1[i] * drho0ds1 +
            d_dgamma2[i] * (dt1ds1 * d_rho1[i] + t1i * drho1ds1 + dt2ds1 * d_rho2[i] + t2i * drho2ds1 +
                          dt3ds1 * d_rho3[i] + t3i * drho3ds1) -
            d_dgamma3[i] * (shpi[0] * dt1ds1 + shpi[1] * dt2ds1 + shpi[2] * dt3ds1);
          drhods2 = d_dgamma1[j] * drho0ds2 +
            d_dgamma2[j] * (dt1ds2 * d_rho1[j] + t1j * drho1ds2 + dt2ds2 * d_rho2[j] + t2j * drho2ds2 +
                          dt3ds2 * d_rho3[j] + t3j * drho3ds2) -
            d_dgamma3[j] * (shpj[0] * dt1ds2 + shpj[1] * dt2ds2 + shpj[2] * dt3ds2);
        }

        //     Compute derivatives of energy wrt rij, sij and rij[3]
        dUdrij = phip * sij + d_frhop[i] * drhodr1 + d_frhop[j] * drhodr2;
        dUdsij = 0.0;
        if (!iszero_kk(d_dscrfcn[fnoffset + jn])) {
          dUdsij = phi + d_frhop[i] * drhods1 + d_frhop[j] * drhods2;
        }
        for (m = 0; m < 3; m++) {
          dUdrijm[m] = d_frhop[i] * drhodrm1[m] + d_frhop[j] * drhodrm2[m];
        }

        //     Add the part of the force due to dUdrij and dUdsij
        force = dUdrij * recip + dUdsij * d_dscrfcn[fnoffset + jn];
        for (m = 0; m < 3; m++) {
          forcem = delij[m] * force + dUdrijm[m];
          a_f(i,m) += forcem;
          a_f(j,m) -= -forcem;
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
            a_vatom(i,m) += v[m];
            a_vatom(j,m) += v[m];
          }
        }

        //     Now compute forces on other atoms k due to change in sij

        if (iszero_kk(sij) || iszero_kk(sij - 1.0)) continue; //: cont jn loop

        double dxik(0), dyik(0), dzik(0);
        double dxjk(0), dyjk(0), dzjk(0);

        for (kn = 0; kn < d_numneigh_full[i]; kn++) {
          //k = firstneigh_full[kn];
          k = d_neighbors_full(i,kn);
          eltk = fmap[type[k]];
          if (k != j && eltk >= 0) {
            double xik, xjk, cikj, sikj, dfc, a;
            double dCikj1, dCikj2;
            double delc, rik2, rjk2;

            sij = d_scrfcn[jn+fnoffset] * d_fcpair[jn+fnoffset];
            const double Cmax = this->Cmax_meam[elti][eltj][eltk];
            const double Cmin = this->Cmin_meam[elti][eltj][eltk];

            dsij1 = 0.0;
            dsij2 = 0.0;
            if (!iszero_kk(sij) && !iszero_kk(sij - 1.0)) {
              const double rbound = rij2 * this->ebound_meam[elti][eltj];
              delc = Cmax - Cmin;
              dxjk = x(k,0) - x(j,0);
              dyjk = x(k,1) - x(j,1);
              dzjk = x(k,2) - x(j,2);
              rjk2 = dxjk * dxjk + dyjk * dyjk + dzjk * dzjk;
              if (rjk2 <= rbound) {
                dxik = x(k,0) - x(i,0);
                dyik = x(k,1) - x(i,1);
                dzik = x(k,2) - x(i,2);
                rik2 = dxik * dxik + dyik * dyik + dzik * dzik;
                if (rik2 <= rbound) {
                  xik = rik2 / rij2;
                  xjk = rjk2 / rij2;
                  a = 1 - (xik - xjk) * (xik - xjk);
                  if (!iszero_kk(a)) {
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

            if (!iszero_kk(dsij1) || !iszero_kk(dsij2)) {
              force1 = dUdsij * dsij1;
              force2 = dUdsij * dsij2;

              a_f(i,0) += force1 * dxik;
              a_f(i,1) += force1 * dyik;
              a_f(i,2) += force1 * dzik;
              a_f(j,0) += force2 * dxjk;
              a_f(j,1) += force2 * dyjk;
              a_f(j,2) += force2 * dzjk;
              a_f(k,0) -= force1 * dxik + force2 * dxjk;
              a_f(k,1) -= force1 * dyik + force2 * dyjk;
              a_f(k,2) -= force1 * dzik + force2 * dzjk;

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
                  a_vatom(i,m) += v[m];
                  a_vatom(j,m) += v[m];
                  a_vatom(k,m) += v[m];
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
