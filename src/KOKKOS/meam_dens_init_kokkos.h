#include "meam_kokkos.h"
#include "math_special.h"
#include "mpi.h"

using namespace LAMMPS_NS;
using namespace MathSpecialKokkos;

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void MEAMKokkos<DeviceType>::operator()(TagMEAMDensInit<NEIGHFLAG>, const int &i, EV_FLOAT &ev) const {
  int ii, offsetval;
  ii = d_ilist_half[i];
  offsetval = d_offset[i];
  //     Compute screening function and derivatives
  this->template getscreen<NEIGHFLAG>(ii, offsetval, x, d_numneigh_half, 
            d_numneigh_full, ntype, type, fmap);

  //     Calculate intermediate density terms to be communicated
  this->template calc_rho1<NEIGHFLAG>(ii, ntype, type, fmap, x, d_numneigh_half, offsetval);
  ev.evdwl += d_numneigh_half[i];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void MEAMKokkos<DeviceType>::operator()(TagMEAMInitialize, const int &i) const {
    d_rho0[i] = 0.0;
    d_arho2b[i] = 0.0;
    d_arho1(i,0) = d_arho1(i,1) = d_arho1(i,2) = 0.0;
    for (int j = 0; j < 6; j++)
      d_arho2(i,j) = 0.0;
    for (int j = 0; j < 10; j++)
      d_arho3(i,j) = 0.0;
    d_arho3b(i,0) = d_arho3b(i,1) = d_arho3b(i,2) = 0.0;
    d_t_ave(i,0) = d_t_ave(i,1) = d_t_ave(i,2) = 0.0;
    d_tsq_ave(i,0) = d_tsq_ave(i,1) = d_tsq_ave(i,2) = 0.0;
}

template<class DeviceType>
void
MEAMKokkos<DeviceType>::meam_dens_setup(int atom_nmax, int nall, int n_neigh)
{
  MemoryKokkos *memoryKK = (MemoryKokkos *)memory;

  // grow local arrays if necessary

  if (atom_nmax > nmax) {
    memoryKK->destroy_kokkos(k_rho,rho);
    memoryKK->destroy_kokkos(k_rho0,rho0);
    memoryKK->destroy_kokkos(k_rho1,rho1);
    memoryKK->destroy_kokkos(k_rho2,rho2);
    memoryKK->destroy_kokkos(k_rho3,rho3);
    memoryKK->destroy_kokkos(k_frhop,frhop);
    memoryKK->destroy_kokkos(k_gamma,gamma);
    memoryKK->destroy_kokkos(k_dgamma1,dgamma1);
    memoryKK->destroy_kokkos(k_dgamma2,dgamma2);
    memoryKK->destroy_kokkos(k_dgamma3,dgamma3);
    memoryKK->destroy_kokkos(k_arho2b,arho2b);
    memoryKK->destroy_kokkos(k_arho1,arho1);
    memoryKK->destroy_kokkos(k_arho2,arho2);
    memoryKK->destroy_kokkos(k_arho3,arho3);
    memoryKK->destroy_kokkos(k_arho3b,arho3b);
    memoryKK->destroy_kokkos(k_t_ave,t_ave);
    memoryKK->destroy_kokkos(k_tsq_ave,tsq_ave);

    nmax = atom_nmax;
//    memory->create(rho, nmax, "pair:rho");
    k_rho = DAT::tdual_ffloat_1d("pair:rho",nmax);
    d_rho = k_rho.template view<DeviceType>();
    h_rho = k_rho.h_view;
 //   memory->create(rho0, nmax, "pair:rho0");
    k_rho0 = DAT::tdual_ffloat_1d("pair:rho0",nmax);
    d_rho0 = k_rho0.template view<DeviceType>();
    h_rho0 = k_rho0.h_view;
    //memory->create(rho1, nmax, "pair:rho1");
    k_rho1 = DAT::tdual_ffloat_1d("pair:rho1",nmax);
    d_rho1 = k_rho1.template view<DeviceType>();
    h_rho1 = k_rho1.h_view;
    //memory->create(rho2, nmax, "pair:rho2");
    k_rho2 = DAT::tdual_ffloat_1d("pair:rho2",nmax);
    d_rho2 = k_rho2.template view<DeviceType>();
    h_rho2 = k_rho2.h_view;
    //memory->create(rho3, nmax, "pair:rho3");
    k_rho3 = DAT::tdual_ffloat_1d("pair:rho3",nmax);
    d_rho3 = k_rho3.template view<DeviceType>();
    h_rho3 = k_rho3.h_view;
    //memory->create(frhop, nmax, "pair:frhop");
    k_frhop = DAT::tdual_ffloat_1d("pair:frhop",nmax);
    d_frhop = k_frhop.template view<DeviceType>();
    h_frhop = k_frhop.h_view;
    //memory->create(gamma, nmax, "pair:gamma");
    k_gamma = DAT::tdual_ffloat_1d("pair:gamma",nmax);
    d_gamma = k_gamma.template view<DeviceType>();
    h_gamma = k_gamma.h_view;
    //memory->create(dgamma1, nmax, "pair:dgamma1");
    k_dgamma1 = DAT::tdual_ffloat_1d("pair:dgamma1",nmax);
    d_dgamma1 = k_dgamma1.template view<DeviceType>();
    h_dgamma1 = k_dgamma1.h_view;
    //memory->create(dgamma2, nmax, "pair:dgamma2");
    k_dgamma2 = DAT::tdual_ffloat_1d("pair:dgamma2",nmax);
    d_dgamma2 = k_dgamma2.template view<DeviceType>();
    h_dgamma2 = k_dgamma2.h_view;
    //memory->create(dgamma3, nmax, "pair:dgamma3");
    k_dgamma3 = DAT::tdual_ffloat_1d("pair:dgamma3",nmax);
    d_dgamma3 = k_dgamma3.template view<DeviceType>();
    h_dgamma3 = k_dgamma3.h_view;
    //memory->create(arho2b, nmax, "pair:arho2b");
    k_arho2b = DAT::tdual_ffloat_1d("pair:arho2b",nmax);
    d_arho2b = k_arho2b.template view<DeviceType>();
    h_arho2b = k_arho2b.h_view;
    //memory->create(arho1, nmax, 3, "pair:arho1");
    k_arho1 = DAT::tdual_ffloat_2d("pair:arho1",nmax, 3);
    d_arho1 = k_arho1.template view<DeviceType>();
    h_arho1 = k_arho1.h_view;
    //memory->create(arho2, nmax, 6, "pair:arho2");
    k_arho2 = DAT::tdual_ffloat_2d("pair:arho2",nmax, 6);
    d_arho2 = k_arho2.template view<DeviceType>();
    h_arho2 = k_arho2.h_view;
    //memory->create(arho3, nmax, 10, "pair:arho3");
    k_arho3 = DAT::tdual_ffloat_2d("pair:arho3",nmax, 10);
    d_arho3 = k_arho3.template view<DeviceType>();
    h_arho3 = k_arho3.h_view;
    //memory->create(arho3b, nmax, 3, "pair:arho3b");
    k_arho3b = DAT::tdual_ffloat_2d("pair:arho3b",nmax, 3);
    d_arho3b = k_arho3b.template view<DeviceType>();
    h_arho3b = k_arho3b.h_view;
    //memory->create(t_ave, nmax, 3, "pair:t_ave");
    k_t_ave = DAT::tdual_ffloat_2d("pair:t_ave",nmax, 3);
    d_t_ave = k_t_ave.template view<DeviceType>();
    h_t_ave = k_t_ave.h_view;
    //memory->create(tsq_ave, nmax, 3, "pair:tsq_ave");
    k_tsq_ave = DAT::tdual_ffloat_2d("pair:tsq_ave",nmax, 3);
    d_tsq_ave = k_tsq_ave.template view<DeviceType>();
    h_tsq_ave = k_tsq_ave.h_view;
  }

  if (n_neigh > maxneigh) {
    memoryKK->destroy_kokkos(k_scrfcn,scrfcn);
    memoryKK->destroy_kokkos(k_dscrfcn,dscrfcn);
    memoryKK->destroy_kokkos(k_fcpair,fcpair);
    maxneigh = n_neigh;
   // memory->create(scrfcn, maxneigh, "pair:scrfcn");
    k_scrfcn = DAT::tdual_ffloat_1d("pair:scrfcn", maxneigh);
    d_scrfcn = k_scrfcn.template view<DeviceType>();
    h_scrfcn = k_scrfcn.h_view;
    //memory->create(dscrfcn, maxneigh, "pair:dscrfcn");
    k_dscrfcn = DAT::tdual_ffloat_1d("pair:dscrfcn", maxneigh);
    d_dscrfcn = k_dscrfcn.template view<DeviceType>();
    h_dscrfcn = k_dscrfcn.h_view;
    //memory->create(fcpair, maxneigh, "pair:fcpair");
    k_fcpair = DAT::tdual_ffloat_1d("pair:fcpair", maxneigh);
    d_fcpair = k_fcpair.template view<DeviceType>();
    h_fcpair = k_fcpair.h_view;
  }

  // zero out local arrays

   Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagMEAMInitialize>(0, nall),*this);   
}

template<class DeviceType>
void
MEAMKokkos<DeviceType>::meam_dens_init(int inum_half, int ntype, typename AT::t_int_1d_randomread type, typename AT::t_int_1d_randomread fmap, typename AT::t_x_array_randomread x, typename AT::t_int_1d_randomread d_numneigh_half, typename AT::t_int_1d_randomread d_numneigh_full,
                     int *fnoffset, typename AT::t_int_1d_randomread d_ilist_half, typename AT::t_neighbors_2d d_neighbors_half, typename AT::t_neighbors_2d d_neighbors_full, typename AT::t_int_1d_randomread d_offset, int neighflag)
{
  EV_FLOAT ev;
 
  ev.evdwl = 0; 
  this->ntype = ntype;
  this->type = type;
  this->fmap = fmap;
  this->x = x;
  this->d_numneigh_half = d_numneigh_half;
  this->d_numneigh_full = d_numneigh_full;
  this->d_ilist_half = d_ilist_half;
  this->d_neighbors_half = d_neighbors_half;
  this->d_neighbors_full = d_neighbors_full;
  this->d_offset = d_offset;
  if (neighflag == FULL)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagMEAMDensInit<FULL> >(0,inum_half),*this, ev);
  else if (neighflag == HALF)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagMEAMDensInit<HALF> >(0,inum_half),*this, ev);
  else if (neighflag == HALFTHREAD)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagMEAMDensInit<HALFTHREAD> >(0,inum_half),*this, ev);
  *fnoffset = (int)ev.evdwl;
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void
MEAMKokkos<DeviceType>::getscreen(int i, int offset, typename AT::t_x_array_randomread x, typename AT::t_int_1d_randomread numneigh_half,
                typename AT::t_int_1d_randomread numneigh_full, int /*ntype*/, typename AT::t_int_1d_randomread type, typename AT::t_int_1d_randomread fmap)
const {
  int jn, j, kn, k;
  int elti, eltj, eltk;
  X_FLOAT xitmp, yitmp, zitmp, delxij, delyij, delzij;
  F_FLOAT  rij2, rij;
  X_FLOAT xjtmp, yjtmp, zjtmp, delxik, delyik, delzik;
  F_FLOAT  rik2 /*,rik*/;
  X_FLOAT xktmp, yktmp, zktmp, delxjk, delyjk, delzjk;
  F_FLOAT  rjk2 /*,rjk*/;
  X_FLOAT xik, xjk;
  F_FLOAT sij, fcij, sfcij, dfcij, sikj, dfikj, cikj;
  F_FLOAT Cmin, Cmax, delc, /*ebound,*/ a, coef1, coef2;
  F_FLOAT dCikj;
  F_FLOAT rnorm, fc, dfc, drinv;

  drinv = 1.0 / this->delr_meam;
  elti = fmap[type[i]];
  if (elti < 0) return;

  xitmp = x(i,0);
  yitmp = x(i,1);
  zitmp = x(i,2);

  for (jn = 0; jn < numneigh_half[i]; jn++) {
    //j = firstneigh[jn];
    j = d_neighbors_half(i,jn);

    eltj = fmap[type[j]];
    if (eltj < 0) continue;

    //     First compute screening function itself, sij
    xjtmp = x(j,0);
    yjtmp = x(j,1);
    zjtmp = x(j,2);
    delxij = xjtmp - xitmp;
    delyij = yjtmp - yitmp;
    delzij = zjtmp - zitmp;

    rij2 = delxij * delxij + delyij * delyij + delzij * delzij;
    rij = sqrt(rij2);

    double rbound = this->ebound_meam[elti][eltj] * rij2;
    if (rij > this->rc_meam) {
      fcij = 0.0;
      dfcij = 0.0;
      sij = 0.0;
    } else {
      rnorm = (this->rc_meam - rij) * drinv;
      sij = 1.0;

      //     if rjk2 > ebound*rijsq, atom k is definitely outside the ellipse
      for (kn = 0; kn < numneigh_full[i]; kn++) {
        //k = firstneigh_full[kn];
        k = d_neighbors_full(i,kn);
        eltk = fmap[type[k]];
        if (eltk < 0) continue;
        if (k == j) continue;

        delxjk = x(k,0) - xjtmp;
        delyjk = x(k,1) - yjtmp;
        delzjk = x(k,2) - zjtmp;
        rjk2 = delxjk * delxjk + delyjk * delyjk + delzjk * delzjk;
        if (rjk2 > rbound) continue;

        delxik = x(k,0) - xitmp;
        delyik = x(k,1) - yitmp;
        delzik = x(k,2) - zitmp;
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
    }
    //     Now compute derivatives
    d_dscrfcn[offset+jn] = 0.0;
    sfcij = sij * fcij;
    if (iszero_kk(sfcij) || iszero_kk(sfcij - 1.0))
    {
      d_scrfcn[offset+jn] = sij;
      d_fcpair[offset+jn] = fcij;
      continue;
      //goto LABEL_100;
    }

    for (kn = 0; kn < numneigh_full[i]; kn++) {
      //k = firstneigh_full[kn];
      k = d_neighbors_full(i,kn);
      if (k == j) continue;
      eltk = fmap[type[k]];
      if (eltk < 0) continue;

      xktmp = x(k,0);
      yktmp = x(k,1);
      zktmp = x(k,2);
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
        d_dscrfcn[offset+jn] = d_dscrfcn[offset+jn] + coef1 * dCikj;
      }
    }
    coef1 = sfcij;
    coef2 = sij * dfcij / rij;
    d_dscrfcn[offset+jn] = d_dscrfcn[offset+jn] * coef1 - coef2;

    //LABEL_100:
    d_scrfcn[offset+jn] = sij;
    d_fcpair[offset+jn] = fcij;
  }
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void
MEAMKokkos<DeviceType>::calc_rho1(int i, int /*ntype*/, typename AT::t_int_1d_randomread type, typename AT::t_int_1d_randomread fmap, typename AT::t_x_array_randomread x, typename AT::t_int_1d_randomread numneigh,
                int offset) const
{
  int jn, j, m, n, p, elti, eltj;
  int nv2, nv3;
  X_FLOAT xtmp, ytmp, ztmp, delij[3];
  F_FLOAT rij2, rij, sij;
  F_FLOAT ai, aj, rhoa0j, rhoa1j, rhoa2j, rhoa3j, A1j, A2j, A3j;
  // double G,Gbar,gam,shp[3+1];
  F_FLOAT ro0i, ro0j;
  F_FLOAT rhoa0i, rhoa1i, rhoa2i, rhoa3i, A1i, A2i, A3i;

  // Likely to do: replace with duplicated view for OpenMP, atomic view for GPU abstraction
  Kokkos::View<F_FLOAT*, typename DAT::t_ffloat_1d::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_rho0 = d_rho0;
  Kokkos::View<F_FLOAT*, typename DAT::t_ffloat_1d::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_arho2b = d_arho2b;
  Kokkos::View<F_FLOAT**, typename DAT::t_ffloat_2d::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_t_ave = d_t_ave;
  Kokkos::View<F_FLOAT**, typename DAT::t_ffloat_2d::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_tsq_ave = d_tsq_ave;
  Kokkos::View<F_FLOAT**, typename DAT::t_ffloat_2d::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_arho1 = d_arho1;
  Kokkos::View<F_FLOAT**, typename DAT::t_ffloat_2d::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_arho2 = d_arho2;
  Kokkos::View<F_FLOAT**, typename DAT::t_ffloat_2d::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_arho3 = d_arho3;
  Kokkos::View<F_FLOAT**, typename DAT::t_ffloat_2d::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_arho3b = d_arho3b;

  elti = fmap[type[i]];
  xtmp = x(i,0);
  ytmp = x(i,1);
  ztmp = x(i,2);
  for (jn = 0; jn < numneigh[i]; jn++) {
    if (!iszero_kk(d_scrfcn[offset+jn])) {
      //j = firstneigh[jn];
      j = d_neighbors_half(i,jn);
      sij = d_scrfcn[offset+jn] * d_fcpair[offset+jn];
      delij[0] = x(j,0) - xtmp;
      delij[1] = x(j,1) - ytmp;
      delij[2] = x(j,2) - ztmp;
      rij2 = delij[0] * delij[0] + delij[1] * delij[1] + delij[2] * delij[2];
      if (rij2 < this->cutforcesq) {
        eltj = fmap[type[j]];
        rij = sqrt(rij2);
        ai = rij / this->re_meam[elti][elti] - 1.0;
        aj = rij / this->re_meam[eltj][eltj] - 1.0;
        ro0i = this->rho0_meam[elti];
        ro0j = this->rho0_meam[eltj];
        rhoa0j = ro0j * MathSpecialKokkos::fm_exp(-this->beta0_meam[eltj] * aj) * sij;
        rhoa1j = ro0j * MathSpecialKokkos::fm_exp(-this->beta1_meam[eltj] * aj) * sij;
        rhoa2j = ro0j * MathSpecialKokkos::fm_exp(-this->beta2_meam[eltj] * aj) * sij;
        rhoa3j = ro0j * MathSpecialKokkos::fm_exp(-this->beta3_meam[eltj] * aj) * sij;
        rhoa0i = ro0i * MathSpecialKokkos::fm_exp(-this->beta0_meam[elti] * ai) * sij;
        rhoa1i = ro0i * MathSpecialKokkos::fm_exp(-this->beta1_meam[elti] * ai) * sij;
        rhoa2i = ro0i * MathSpecialKokkos::fm_exp(-this->beta2_meam[elti] * ai) * sij;
        rhoa3i = ro0i * MathSpecialKokkos::fm_exp(-this->beta3_meam[elti] * ai) * sij;
        if (this->ialloy == 1) {
          rhoa1j = rhoa1j * this->t1_meam[eltj];
          rhoa2j = rhoa2j * this->t2_meam[eltj];
          rhoa3j = rhoa3j * this->t3_meam[eltj];
          rhoa1i = rhoa1i * this->t1_meam[elti];
          rhoa2i = rhoa2i * this->t2_meam[elti];
          rhoa3i = rhoa3i * this->t3_meam[elti];
        }
        a_rho0[i] += rhoa0j;
        a_rho0[j] += rhoa0i;
        // For ialloy = 2, use single-element value (not average)
        if (this->ialloy != 2) {
          a_t_ave(i,0) = a_t_ave(i,0) + this->t1_meam[eltj] * rhoa0j;
          a_t_ave(i,1) = a_t_ave(i,1) + this->t2_meam[eltj] * rhoa0j;
          a_t_ave(i,2) = a_t_ave(i,2) + this->t3_meam[eltj] * rhoa0j;
          a_t_ave(j,0) = a_t_ave(j,0) + this->t1_meam[elti] * rhoa0i;
          a_t_ave(j,1) = a_t_ave(j,1) + this->t2_meam[elti] * rhoa0i;
          a_t_ave(j,2) = a_t_ave(j,2) + this->t3_meam[elti] * rhoa0i;
        }
        if (this->ialloy == 1) {
          a_tsq_ave(i,0) = a_tsq_ave(i,0) + this->t1_meam[eltj] * this->t1_meam[eltj] * rhoa0j;
          a_tsq_ave(i,1) = a_tsq_ave(i,1) + this->t2_meam[eltj] * this->t2_meam[eltj] * rhoa0j;
          a_tsq_ave(i,2) = a_tsq_ave(i,2) + this->t3_meam[eltj] * this->t3_meam[eltj] * rhoa0j;
          a_tsq_ave(j,0) = a_tsq_ave(j,0) + this->t1_meam[elti] * this->t1_meam[elti] * rhoa0i;
          a_tsq_ave(j,1) = a_tsq_ave(j,1) + this->t2_meam[elti] * this->t2_meam[elti] * rhoa0i;
          a_tsq_ave(j,2) = a_tsq_ave(j,2) + this->t3_meam[elti] * this->t3_meam[elti] * rhoa0i;
        }
        a_arho2b[i] += rhoa2j;
        a_arho2b[j] += rhoa2i;

        A1j = rhoa1j / rij;
        A2j = rhoa2j / rij2;
        A3j = rhoa3j / (rij2 * rij);
        A1i = rhoa1i / rij;
        A2i = rhoa2i / rij2;
        A3i = rhoa3i / (rij2 * rij);
        nv2 = 0;
        nv3 = 0;
        for (m = 0; m < 3; m++) {
          a_arho1(i,m) += A1j * delij[m];
          a_arho1(j,m) += (-A1i * delij[m]);
          a_arho3b(i,m) += rhoa3j * delij[m] / rij;
          a_arho3b(j,m) += (- rhoa3i * delij[m] / rij);
          for (n = m; n < 3; n++) {
            a_arho2(i,nv2) += A2j * delij[m] * delij[n];
            a_arho2(j,nv2) += A2i * delij[m] * delij[n];
            nv2 = nv2 + 1;
            for (p = n; p < 3; p++) {
              a_arho3(i,nv3) += A3j * delij[m] * delij[n] * delij[p];
              a_arho3(j,nv3) += (- A3i * delij[m] * delij[n] * delij[p]);
              nv3 = nv3 + 1;
            }
          }
        }
      }
    }
  }
}

//Cutoff function and derivative

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double MEAMKokkos<DeviceType>::dfcut(const double xi, double& dfc) const {
    double a, a3, a4, a1m4;
    if (xi >= 1.0) {
      dfc = 0.0;
      return 1.0;
    } else if (xi <= 0.0) {
      dfc = 0.0;
      return 0.0;
    } else {
      a = 1.0 - xi;
      a3 = a * a * a;
      a4 = a * a3;
      a1m4 = 1.0-a4;

      dfc = 8 * a1m4 * a3;
      return a1m4*a1m4;
    }
  }

  //-----------------------------------------------------------------------------
  //  // Derivative of Cikj w.r.t. rij
  //    //     Inputs: rij,rij2,rik2,rjk2
  //      //
  //

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double MEAMKokkos<DeviceType>::dCfunc(const double rij2, const double rik2, const double rjk2) const
{
    double rij4, a, asq, b,denom;

    rij4 = rij2 * rij2;
    a = rik2 - rjk2;
    b = rik2 + rjk2;
    asq = a*a;
    denom = rij4 - asq;
    denom = denom * denom;
    return -4 * (-2 * rij2 * asq + rij4 * b + asq * b) / denom;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void MEAMKokkos<DeviceType>::dCfunc2(const double rij2, const double rik2, const double rjk2, double& dCikj1, double& dCikj2) const
{
    double rij4, rik4, rjk4, a, denom;

    rij4 = rij2 * rij2;
    rik4 = rik2 * rik2;
    rjk4 = rjk2 * rjk2;
    a = rik2 - rjk2;
    denom = rij4 - a * a;
    denom = denom * denom;
    dCikj1 = 4 * rij2 * (rij4 + rik4 + 2 * rik2 * rjk2 - 3 * rjk4 - 2 * rij2 * a) / denom;
    dCikj2 = 4 * rij2 * (rij4 - 3 * rik4 + 2 * rik2 * rjk2 + rjk4 + 2 * rij2 * a) / denom;

}
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double MEAMKokkos<DeviceType>::fcut(const double xi) const{
    double a;
    if (xi >= 1.0)
      return 1.0;
    else if (xi <= 0.0)
      return 0.0;
    else {
      a = 1.0 - xi;
      a *= a; a *= a;
      a = 1.0 - a;
      return a * a;
    }
  }

