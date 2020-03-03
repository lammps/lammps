/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Christian Trott (SNL), Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "sna_kokkos.h"
#include <cmath>
#include <cstring>
#include <cstdlib>

namespace LAMMPS_NS {

static const double MY_PI  = 3.14159265358979323846; // pi

template<class DeviceType>
inline
SNAKokkos<DeviceType>::SNAKokkos(double rfac0_in,
         int twojmax_in,
         double rmin0_in, int switch_flag_in, int bzero_flag_in)
{
  wself = 1.0;

  rfac0 = rfac0_in;
  rmin0 = rmin0_in;
  switch_flag = switch_flag_in;
  bzero_flag = bzero_flag_in;

  twojmax = twojmax_in;

  ncoeff = compute_ncoeff();

  nmax = 0;

  build_indexlist();

  int jdimpq = twojmax + 2;
  rootpqarray = t_sna_2d("SNAKokkos::rootpqarray",jdimpq,jdimpq);

  cglist = t_sna_1d("SNAKokkos::cglist",idxcg_max);

  if (bzero_flag) {
    bzero = Kokkos::View<double*, Kokkos::LayoutRight, DeviceType>("sna:bzero",twojmax+1);
    auto h_bzero = Kokkos::create_mirror_view(bzero);

    double www = wself*wself*wself;
    for(int j = 0; j <= twojmax; j++)
      h_bzero[j] = www*(j+1);
    Kokkos::deep_copy(bzero,h_bzero);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
SNAKokkos<DeviceType>::~SNAKokkos()
{
}

template<class DeviceType>
inline
void SNAKokkos<DeviceType>::build_indexlist()
{
  // index list for cglist

  int jdim = twojmax + 1;
  idxcg_block = Kokkos::View<int***, DeviceType>("SNAKokkos::idxcg_block",jdim,jdim,jdim);
  auto h_idxcg_block = Kokkos::create_mirror_view(idxcg_block);

  int idxcg_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        h_idxcg_block(j1,j2,j) = idxcg_count; 
        for (int m1 = 0; m1 <= j1; m1++)
          for (int m2 = 0; m2 <= j2; m2++)
            idxcg_count++;
      }
  idxcg_max = idxcg_count;
  Kokkos::deep_copy(idxcg_block,h_idxcg_block);

  // index list for uarray
  // need to include both halves

  idxu_block = Kokkos::View<int*, DeviceType>("SNAKokkos::idxu_block",jdim);
  auto h_idxu_block = Kokkos::create_mirror_view(idxu_block);

  int idxu_count = 0;
  
  for(int j = 0; j <= twojmax; j++) {
    h_idxu_block[j] = idxu_count; 
    for(int mb = 0; mb <= j; mb++)
      for(int ma = 0; ma <= j; ma++)
        idxu_count++;
  }
  idxu_max = idxu_count;
  Kokkos::deep_copy(idxu_block,h_idxu_block);

  // index list for beta and B

  int idxb_count = 0;  
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) idxb_count++;
  
  idxb_max = idxb_count;
  idxb = Kokkos::View<int*[3], DeviceType>("SNAKokkos::idxb",idxb_max);
  auto h_idxb = Kokkos::create_mirror_view(idxb);
  
  idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) {
          h_idxb(idxb_count,0) = j1;
          h_idxb(idxb_count,1) = j2;
          h_idxb(idxb_count,2) = j;
          idxb_count++;
        }
  Kokkos::deep_copy(idxb,h_idxb);

  // reverse index list for beta and b

  idxb_block = Kokkos::View<int***, DeviceType>("SNAKokkos::idxb_block",jdim,jdim,jdim);
  auto h_idxb_block = Kokkos::create_mirror_view(idxb_block);

  idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        if (j >= j1) {
          h_idxb_block(j1,j2,j) = idxb_count; 
          idxb_count++;
        }
      }
  Kokkos::deep_copy(idxb_block,h_idxb_block);

  // index list for zlist

  int idxz_count = 0;

  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++)
            idxz_count++;
  
  idxz_max = idxz_count;
  idxz = Kokkos::View<int*[10], DeviceType>("SNAKokkos::idxz",idxz_max);
  auto h_idxz = Kokkos::create_mirror_view(idxz);

  idxz_block = Kokkos::View<int***, DeviceType>("SNAKokkos::idxz_block", jdim,jdim,jdim);
  auto h_idxz_block = Kokkos::create_mirror_view(idxz_block);
  
  idxz_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        h_idxz_block(j1,j2,j) = idxz_count; 

        // find right beta(ii,jjb) entry
        // multiply and divide by j+1 factors
        // account for multiplicity of 1, 2, or 3

        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {
            h_idxz(idxz_count,0) = j1;
            h_idxz(idxz_count,1) = j2;
            h_idxz(idxz_count,2) = j;
            h_idxz(idxz_count,3) = MAX(0, (2 * ma - j - j2 + j1) / 2);
            h_idxz(idxz_count,4) = (2 * ma - j - (2 * h_idxz(idxz_count,3) - j1) + j2) / 2;
            h_idxz(idxz_count,5) = MAX(0, (2 * mb - j - j2 + j1) / 2);
            h_idxz(idxz_count,6) = (2 * mb - j - (2 * h_idxz(idxz_count,5) - j1) + j2) / 2;
            h_idxz(idxz_count,7) = MIN(j1, (2 * ma - j + j2 + j1) / 2) - h_idxz(idxz_count,3) + 1;
            h_idxz(idxz_count,8) = MIN(j1, (2 * mb - j + j2 + j1) / 2) - h_idxz(idxz_count,5) + 1;

            // apply to z(j1,j2,j,ma,mb) to unique element of y(j)

            const int jju = h_idxu_block[j] + (j+1)*mb + ma;
            h_idxz(idxz_count,9) = jju;

            idxz_count++;
          }
      }
  Kokkos::deep_copy(idxz,h_idxz);
  Kokkos::deep_copy(idxz_block,h_idxz_block);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
inline
void SNAKokkos<DeviceType>::init()
{
  init_clebsch_gordan();
  init_rootpqarray();
}

template<class DeviceType>
inline
void SNAKokkos<DeviceType>::grow_rij(int newnatom, int newnmax)
{
  if(newnatom <= natom && newnmax <= nmax) return;
  natom = newnatom;
  nmax = newnmax;

  rij = t_sna_3d("sna:rij",natom,nmax,3);
  inside = t_sna_2i("sna:inside",natom,nmax);
  wj = t_sna_2d("sna:wj",natom,nmax);
  rcutij = t_sna_2d("sna:rcutij",natom,nmax);
  dedr = t_sna_3d("sna:dedr",natom,nmax,3);

  blist = t_sna_2d_ll("sna:blist",idxb_max,natom);
  //ulisttot = t_sna_2c("sna:ulisttot",natom,idxu_max);
  ulisttot = t_sna_2c_ll("sna:ulisttot",idxu_max,natom);
  
  zlist = t_sna_2c_ll("sna:zlist",idxz_max,natom);

  //ulist = t_sna_3c("sna:ulist",natom,nmax,idxu_max);
  ulist = t_sna_3c_ll("sna:ulist",idxu_max,natom,nmax);
  //ylist = t_sna_2c_lr("sna:ylist",natom,idxu_max);
  ylist = t_sna_2c_ll("sna:ylist",idxu_max,natom);

  //dulist = t_sna_4c("sna:dulist",natom,nmax,idxu_max);
  dulist = t_sna_4c_ll("sna:dulist",idxu_max,natom,nmax);
}

/* ----------------------------------------------------------------------
 *    compute Ui by summing over neighbors j
 *    ------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::pre_ui(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, const int& iatom)
{
  for (int j = 0; j <= twojmax; j++) {
    const int jju = idxu_block(j);

    // Only diagonal elements get initialized
    // for (int m = 0; m < (j+1)*(j+1); m++)
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, (j+1)*(j+1)),
      [&] (const int m) {

      const int jjup = jju + m;

      // if m is on the "diagonal", initialize it with the self energy.
      // Otherwise zero it out
      SNAcomplex init = {0., 0.};
      if (m % (j+2) == 0) { init = {wself, 0.0}; }

      ulisttot(jjup, iatom) = init;
    });
  }

}

/* ----------------------------------------------------------------------
   compute Ui by summing over bispectrum components
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_ui(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom, int jnbor)
{
  double rsq, r, x, y, z, z0, theta0;

  // utot(j,ma,mb) = 0 for all j,ma,ma
  // utot(j,ma,ma) = 1 for all j,ma
  // for j in neighbors of i:
  //   compute r0 = (x,y,z,z0)
  //   utot(j,ma,mb) += u(r0;j,ma,mb) for all j,ma,mb

  x = rij(iatom,jnbor,0);
  y = rij(iatom,jnbor,1);
  z = rij(iatom,jnbor,2);
  rsq = x * x + y * y + z * z;
  r = sqrt(rsq);

  theta0 = (r - rmin0) * rfac0 * MY_PI / (rcutij(iatom,jnbor) - rmin0);
  //    theta0 = (r - rmin0) * rscale0;
  z0 = r / tan(theta0);

  compute_uarray(team, iatom, jnbor, x, y, z, z0, r);

  // if we're on the GPU, accumulating into uarraytot is done in a separate kernel.
  // if we're not, it's more efficient to include it in compute_uarray.
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_ui_cpu(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom, int jnbor)
{
  double rsq, r, x, y, z, z0, theta0;

  // utot(j,ma,mb) = 0 for all j,ma,ma
  // utot(j,ma,ma) = 1 for all j,ma
  // for j in neighbors of i:
  //   compute r0 = (x,y,z,z0)
  //   utot(j,ma,mb) += u(r0;j,ma,mb) for all j,ma,mb

  x = rij(iatom,jnbor,0);
  y = rij(iatom,jnbor,1);
  z = rij(iatom,jnbor,2);
  rsq = x * x + y * y + z * z;
  r = sqrt(rsq);

  theta0 = (r - rmin0) * rfac0 * MY_PI / (rcutij(iatom,jnbor) - rmin0);
  //    theta0 = (r - rmin0) * rscale0;
  z0 = r / tan(theta0);

  compute_uarray_cpu(team, iatom, jnbor, x, y, z, z0, r);
  add_uarraytot(team, iatom, jnbor, r, wj(iatom,jnbor), rcutij(iatom,jnbor));
  
}

/* ----------------------------------------------------------------------
   compute UiTot by summing over neighbors
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_uitot(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int idx, int iatom, int ninside)
{
  // fuse initialize in, avoid this load?
  SNAcomplex utot = ulisttot(idx, iatom);
  for (int jnbor = 0; jnbor < ninside; jnbor++) {

    const auto x = rij(iatom,jnbor,0);
    const auto y = rij(iatom,jnbor,1);
    const auto z = rij(iatom,jnbor,2);
    const auto rsq = x * x + y * y + z * z;
    const auto r = sqrt(rsq);

    const double wj_local = wj(iatom, jnbor);
    const double rcut = rcutij(iatom, jnbor);
    const double sfac = compute_sfac(r, rcut) * wj_local;

    auto ulist_local = ulist(idx, iatom, jnbor);
    utot.re += sfac * ulist_local.re;
    utot.im += sfac * ulist_local.im;
  }

  ulisttot(idx, iatom) = utot;

}

/* ----------------------------------------------------------------------
   compute Zi by summing over products of Ui
   not updated yet
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_zi(const int& iter)
{
  const int iatom = iter / idxz_max;
  const int jjz = iter % idxz_max;

  const int j1 = idxz(jjz,0);
  const int j2 = idxz(jjz,1);
  const int j = idxz(jjz,2);
  const int ma1min = idxz(jjz,3);
  const int ma2max = idxz(jjz,4);
  const int mb1min = idxz(jjz,5);
  const int mb2max = idxz(jjz,6);
  const int na = idxz(jjz,7);
  const int nb = idxz(jjz,8);

  const double* cgblock = cglist.data() + idxcg_block(j1,j2,j);

  zlist(jjz,iatom).re = 0.0; 
  zlist(jjz,iatom).im = 0.0;

  int jju1 = idxu_block[j1] + (j1+1)*mb1min;
  int jju2 = idxu_block[j2] + (j2+1)*mb2max;
  int icgb = mb1min*(j2+1) + mb2max;
  for(int ib = 0; ib < nb; ib++) {

    double suma1_r = 0.0;
    double suma1_i = 0.0;

    int ma1 = ma1min;
    int ma2 = ma2max;
    int icga = ma1min*(j2+1) + ma2max;
    for(int ia = 0; ia < na; ia++) {
      suma1_r += cgblock[icga] * (ulisttot(jju1+ma1,iatom).re * ulisttot(jju2+ma2,iatom).re - ulisttot(jju1+ma1,iatom).im * ulisttot(jju2+ma2,iatom).im);
      suma1_i += cgblock[icga] * (ulisttot(jju1+ma1,iatom).re * ulisttot(jju2+ma2,iatom).im + ulisttot(jju1+ma1,iatom).im * ulisttot(jju2+ma2,iatom).re);
      ma1++;
      ma2--;
      icga += j2;
    } // end loop over ia

    zlist(jjz,iatom).re += cgblock[icgb] * suma1_r;
    zlist(jjz,iatom).im += cgblock[icgb] * suma1_i;

    jju1 += j1+1;
    jju2 -= j2+1;
    icgb += j2;
  } // end loop over ib
}

/* ----------------------------------------------------------------------
   compute Yi from Ui without storing Zi, looping over zlist indices
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::zero_yi(const int& idx, const int& iatom)
{
  ylist(idx,iatom) = {0.0, 0.0};
}

/* ----------------------------------------------------------------------
   compute Yi from Ui without storing Zi, looping over zlist indices
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_yi(int iter,
 const Kokkos::View<F_FLOAT**, DeviceType> &beta)
{
  double betaj;
  const int iatom = iter / idxz_max;
  const int jjz = iter % idxz_max;

  const int j1 = idxz(jjz,0);
  const int j2 = idxz(jjz,1);
  const int j = idxz(jjz,2);
  const int ma1min = idxz(jjz,3);
  const int ma2max = idxz(jjz,4);
  const int mb1min = idxz(jjz,5);
  const int mb2max = idxz(jjz,6);
  const int na = idxz(jjz,7);
  const int nb = idxz(jjz,8);
  const int jju = idxz(jjz,9);

  const double* cgblock = cglist.data() + idxcg_block(j1,j2,j);
  //int mb = (2 * (mb1min+mb2max) - j1 - j2 + j) / 2;
  //int ma = (2 * (ma1min+ma2max) - j1 - j2 + j) / 2;

  double ztmp_r = 0.0;
  double ztmp_i = 0.0;

  int jju1 = idxu_block[j1] + (j1+1)*mb1min;
  int jju2 = idxu_block[j2] + (j2+1)*mb2max;
  int icgb = mb1min*(j2+1) + mb2max;
  for(int ib = 0; ib < nb; ib++) {

    double suma1_r = 0.0;
    double suma1_i = 0.0;

    int ma1 = ma1min;
    int ma2 = ma2max;
    int icga = ma1min*(j2+1) + ma2max;

    for(int ia = 0; ia < na; ia++) {
      suma1_r += cgblock[icga] * (ulisttot(jju1+ma1,iatom).re * ulisttot(jju2+ma2,iatom).re - ulisttot(jju1+ma1,iatom).im * ulisttot(jju2+ma2,iatom).im);
      suma1_i += cgblock[icga] * (ulisttot(jju1+ma1,iatom).re * ulisttot(jju2+ma2,iatom).im + ulisttot(jju1+ma1,iatom).im * ulisttot(jju2+ma2,iatom).re);
      ma1++;
      ma2--;
      icga += j2;
    } // end loop over ia

    ztmp_r += cgblock[icgb] * suma1_r;
    ztmp_i += cgblock[icgb] * suma1_i;
    jju1 += j1+1;
    jju2 -= j2+1;
    icgb += j2;
  } // end loop over ib

  // apply to z(j1,j2,j,ma,mb) to unique element of y(j)
  // find right y_list[jju] and beta(iatom,jjb) entries
  // multiply and divide by j+1 factors
  // account for multiplicity of 1, 2, or 3

  // pick out right beta value

  if (j >= j1) {
    const int jjb = idxb_block(j1,j2,j);
    if (j1 == j) {
      if (j2 == j) betaj = 3*beta(jjb,iatom);
      else betaj = 2*beta(jjb,iatom);
    } else betaj = beta(jjb,iatom); 
  } else if (j >= j2) {
    const int jjb = idxb_block(j,j2,j1);
    if (j2 == j) betaj = 2*beta(jjb,iatom)*(j1+1)/(j+1.0);
    else betaj = beta(jjb,iatom)*(j1+1)/(j+1.0);
  } else {
    const int jjb = idxb_block(j2,j,j1);
    betaj = beta(jjb,iatom)*(j1+1)/(j+1.0);
  }

  Kokkos::atomic_add(&(ylist(jju,iatom).re), betaj*ztmp_r);
  Kokkos::atomic_add(&(ylist(jju,iatom).im), betaj*ztmp_i);
}

/* ----------------------------------------------------------------------
   compute dEidRj
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_deidrj(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom, int jnbor)
{
  t_scalar3<double> final_sum;

  // Like in ComputeUi/ComputeDuidrj, regular loop over j.
  for (int j = 0; j <= twojmax; j++) {
    int jju = idxu_block(j);

    // Flatten loop over ma, mb, reduce w/in

    const int n_ma = j+1;
    // for (int mb = 0; 2*mb <= j; mb++)
    const int n_mb = j/2+1;

    const int total_iters = n_ma * n_mb;

    t_scalar3<double> sum;

    //for (int m = 0; m < total_iters; m++) {
    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, total_iters),
      [&] (const int m, t_scalar3<double>& sum_tmp) {

      // ma fast, mb slow
      int ma = m % n_ma;
      int mb = m / n_ma;

      // get index
      const int jju_index = jju+mb+mb*j+ma;

      // get ylist, rescale last element by 0.5
      SNAcomplex y_local = ylist(jju_index,iatom);

      const SNAcomplex du_x = dulist(jju_index,iatom,jnbor,0);
      const SNAcomplex du_y = dulist(jju_index,iatom,jnbor,1);
      const SNAcomplex du_z = dulist(jju_index,iatom,jnbor,2);

      if (j % 2 == 0 && 2*mb == j) {
        if (ma == mb) { y_local = 0.5*y_local; }
        else if (ma > mb) { y_local = { 0., 0. }; } 
        // else the ma < mb gets "double counted", cancelling the 0.5.
      }

      sum_tmp.x += du_x.re * y_local.re + du_x.im * y_local.im;
      sum_tmp.y += du_y.re * y_local.re + du_y.im * y_local.im;
      sum_tmp.z += du_z.re * y_local.re + du_z.im * y_local.im;

    }, sum); // end loop over flattened ma,mb

    final_sum.x += sum.x;
    final_sum.y += sum.y;
    final_sum.z += sum.z;
  }

  Kokkos::single(Kokkos::PerThread(team), [&] () {
    dedr(iatom,jnbor,0) = final_sum.x*2.0;
    dedr(iatom,jnbor,1) = final_sum.y*2.0;
    dedr(iatom,jnbor,2) = final_sum.z*2.0;
  });

}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_deidrj_cpu(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom, int jnbor)
{
  t_scalar3<double> final_sum;

  //for(int j = 0; j <= twojmax; j++) {
  Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,twojmax+1),
      [&] (const int& j, t_scalar3<double>& sum_tmp) {
    int jju = idxu_block[j];

    for(int mb = 0; 2*mb < j; mb++)
      for(int ma = 0; ma <= j; ma++) {
        sum_tmp.x += dulist(jju,iatom,jnbor,0).re * ylist(jju,iatom).re + dulist(jju,iatom,jnbor,0).im * ylist(jju,iatom).im;
        sum_tmp.y += dulist(jju,iatom,jnbor,1).re * ylist(jju,iatom).re + dulist(jju,iatom,jnbor,1).im * ylist(jju,iatom).im;
        sum_tmp.z += dulist(jju,iatom,jnbor,2).re * ylist(jju,iatom).re + dulist(jju,iatom,jnbor,2).im * ylist(jju,iatom).im;
        jju++;
      } //end loop over ma mb

    // For j even, handle middle column

    if (j%2 == 0) {

      int mb = j/2;
      for(int ma = 0; ma < mb; ma++) {
        sum_tmp.x += dulist(jju,iatom,jnbor,0).re * ylist(jju,iatom).re + dulist(jju,iatom,jnbor,0).im * ylist(jju,iatom).im;
        sum_tmp.y += dulist(jju,iatom,jnbor,1).re * ylist(jju,iatom).re + dulist(jju,iatom,jnbor,1).im * ylist(jju,iatom).im;
        sum_tmp.z += dulist(jju,iatom,jnbor,2).re * ylist(jju,iatom).re + dulist(jju,iatom,jnbor,2).im * ylist(jju,iatom).im;
        jju++;
      }

      //int ma = mb;
      sum_tmp.x += (dulist(jju,iatom,jnbor,0).re * ylist(jju,iatom).re + dulist(jju,iatom,jnbor,0).im * ylist(jju,iatom).im)*0.5;
      sum_tmp.y += (dulist(jju,iatom,jnbor,1).re * ylist(jju,iatom).re + dulist(jju,iatom,jnbor,1).im * ylist(jju,iatom).im)*0.5;
      sum_tmp.z += (dulist(jju,iatom,jnbor,2).re * ylist(jju,iatom).re + dulist(jju,iatom,jnbor,2).im * ylist(jju,iatom).im)*0.5;
    } // end if jeven

  },final_sum); // end loop over j

  Kokkos::single(Kokkos::PerThread(team), [&] () {
    dedr(iatom,jnbor,0) = final_sum.x*2.0;
    dedr(iatom,jnbor,1) = final_sum.y*2.0;
    dedr(iatom,jnbor,2) = final_sum.z*2.0;
  });

}

/* ----------------------------------------------------------------------
   compute Bi by summing conj(Ui)*Zi
   not updated yet
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_bi(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom)
{
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        b(j1,j2,j) = 0
  //        for mb = 0,...,jmid
  //          for ma = 0,...,j
  //            b(j1,j2,j) +=
  //              2*Conj(u(j,ma,mb))*z(j1,j2,j,ma,mb)

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,idxb_max),
      [&] (const int& jjb) {
  //for(int jjb = 0; jjb < idxb_max; jjb++) {
    const int j1 = idxb(jjb,0);
    const int j2 = idxb(jjb,1);
    const int j = idxb(jjb,2);

    int jjz = idxz_block(j1,j2,j);
    int jju = idxu_block[j];
    double sumzu = 0.0;
    double sumzu_temp = 0.0;
    const int bound = (j+2)/2;
    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,(j+1)*bound),
        [&] (const int mbma, double& sum) {
        //for(int mb = 0; 2*mb < j; mb++)
          //for(int ma = 0; ma <= j; ma++) {
        const int ma = mbma%(j+1);
        const int mb = mbma/(j+1);
        const int jju_index = jju+mb*(j+1)+ma;
        const int jjz_index = jjz+mb*(j+1)+ma;
        if (2*mb == j) return;
        sum +=
          ulisttot(jju_index,iatom).re * zlist(jjz_index,iatom).re +
          ulisttot(jju_index,iatom).im * zlist(jjz_index,iatom).im;
      },sumzu_temp); // end loop over ma, mb
      sumzu += sumzu_temp;

    // For j even, special treatment for middle column

    if (j%2 == 0) {
      const int mb = j/2;
      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, mb),
          [&] (const int ma, double& sum) {
      //for(int ma = 0; ma < mb; ma++) {
        const int jju_index = jju+(mb-1)*(j+1)+(j+1)+ma;
        const int jjz_index = jjz+(mb-1)*(j+1)+(j+1)+ma;
        sum +=
          ulisttot(jju_index,iatom).re * zlist(jjz_index,iatom).re +
          ulisttot(jju_index,iatom).im * zlist(jjz_index,iatom).im;
      },sumzu_temp); // end loop over ma
      sumzu += sumzu_temp;

      const int ma = mb;
      const int jju_index = jju+(mb-1)*(j+1)+(j+1)+ma;
      const int jjz_index = jjz+(mb-1)*(j+1)+(j+1)+ma;
      sumzu += 0.5*
        (ulisttot(jju_index,iatom).re * zlist(jjz_index,iatom).re +
         ulisttot(jju_index,iatom).im * zlist(jjz_index,iatom).im);
    } // end if jeven

    Kokkos::single(Kokkos::PerThread(team), [&] () {
      sumzu *= 2.0;

      // apply bzero shift

      if (bzero_flag)
        sumzu -= bzero[j];

      blist(jjb,iatom) = sumzu;
    });
  });
      //} // end loop over j
    //} // end loop over j1, j2
}

/* ----------------------------------------------------------------------
   calculate derivative of Ui w.r.t. atom j
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_duidrj(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom, int jnbor)
{
  double rsq, r, x, y, z, z0, theta0, cs, sn;
  double dz0dr;

  x = rij(iatom,jnbor,0);
  y = rij(iatom,jnbor,1);
  z = rij(iatom,jnbor,2);
  rsq = x * x + y * y + z * z;
  r = sqrt(rsq);
  double rscale0 = rfac0 * MY_PI / (rcutij(iatom,jnbor) - rmin0);
  theta0 = (r - rmin0) * rscale0;
  cs = cos(theta0);
  sn = sin(theta0);
  z0 = r * cs / sn;
  dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;

  compute_duarray(team, iatom, jnbor, x, y, z, z0, r, dz0dr, wj(iatom,jnbor), rcutij(iatom,jnbor));
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_duidrj_cpu(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom, int jnbor)
{
  double rsq, r, x, y, z, z0, theta0, cs, sn;
  double dz0dr;

  x = rij(iatom,jnbor,0);
  y = rij(iatom,jnbor,1);
  z = rij(iatom,jnbor,2);
  rsq = x * x + y * y + z * z;
  r = sqrt(rsq);
  double rscale0 = rfac0 * MY_PI / (rcutij(iatom,jnbor) - rmin0);
  theta0 = (r - rmin0) * rscale0;
  cs = cos(theta0);
  sn = sin(theta0);
  z0 = r * cs / sn;
  dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;

  compute_duarray_cpu(team, iatom, jnbor, x, y, z, z0, r, dz0dr, wj(iatom,jnbor), rcutij(iatom,jnbor));
}

/* ----------------------------------------------------------------------
   add Wigner U-functions for one neighbor to the total
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::add_uarraytot(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom, int jnbor,
                                          double r, double wj, double rcut)
{
  const double sfac = compute_sfac(r, rcut) * wj;

  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,ulisttot.extent(0)),
      [&] (const int& i) {
    Kokkos::atomic_add(&(ulisttot(i,iatom).re), sfac * ulist(i,iatom,jnbor).re);
    Kokkos::atomic_add(&(ulisttot(i,iatom).im), sfac * ulist(i,iatom,jnbor).im);
  });
}

/* ----------------------------------------------------------------------
   compute Wigner U-functions for one neighbor
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_uarray(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom, int jnbor,
                         double x, double y, double z,
                         double z0, double r)
{
  // define size of scratch memory buffer
  const int max_m_tile = (twojmax+1)*(twojmax+1);
  const int team_rank = team.team_rank();

  // get scratch memory double buffer
  SNAcomplex* buf1 = (SNAcomplex*)team.team_shmem( ).get_shmem(team.team_size()*max_m_tile*sizeof(SNAcomplex), 0);
  SNAcomplex* buf2 = (SNAcomplex*)team.team_shmem( ).get_shmem(team.team_size()*max_m_tile*sizeof(SNAcomplex), 0);

  // compute Cayley-Klein parameters for unit quaternion,
  // pack into complex number
  double r0inv = 1.0 / sqrt(r * r + z0 * z0);
  SNAcomplex a = { r0inv * z0, -r0inv * z };
  SNAcomplex b = { r0inv * y, -r0inv * x };

  // VMK Section 4.8.2

  // All writes go to global memory and shared memory
  // so we can avoid all global memory reads!
  Kokkos::single(Kokkos::PerThread(team), [=]() {
    ulist(0,iatom,jnbor) = { 1.0, 0.0 };
    buf1[max_m_tile*team_rank] = {1.,0.};
  });

  for (int j = 1; j <= twojmax; j++) {
    const int jju = idxu_block[j];
    int jjup = idxu_block[j-1];

    // fill in left side of matrix layer from previous layer

    // Flatten loop over ma, mb, need to figure out total
    // number of iterations

    // for (int ma = 0; ma <= j; ma++)
    const int n_ma = j+1;
    // for (int mb = 0; 2*mb <= j; mb++)
    const int n_mb = j/2+1;

    const int total_iters = n_ma * n_mb;

    //for (int m = 0; m < total_iters; m++) {
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, total_iters),
      [&] (const int m) {

      // ma fast, mb slow
      int ma = m % n_ma;
      int mb = m / n_ma;

      // index into global memory array
      const int jju_index = jju+mb+mb*j+ma;

      // index into shared memory buffer for previous level
      const int jju_shared_idx = max_m_tile*team_rank+mb+mb*j+ma;

      // index into shared memory buffer for next level
      const int jjup_shared_idx = max_m_tile*team_rank+mb*j+ma;

      SNAcomplex u_accum = {0., 0.};
      
      // VMK recursion relation: grab contribution which is multiplied by a*
      const double rootpq1 = rootpqarray(j - ma, j - mb);
      const SNAcomplex u_up1 = (ma < j)?rootpq1*buf1[jjup_shared_idx]:SNAcomplex(0.,0.);
      caconjxpy(a, u_up1, u_accum);

      // VMK recursion relation: grab contribution which is multiplied by b*
      const double rootpq2 = -rootpqarray(ma, j - mb);
      const SNAcomplex u_up2 = (ma > 0)?rootpq2*buf1[jjup_shared_idx-1]:SNAcomplex(0.,0.);
      caconjxpy(b, u_up2, u_accum);

      ulist(jju_index,iatom,jnbor) = u_accum;

      // We no longer accumulate into ulisttot in this kernel.
      // Instead, we have a separate kernel which avoids atomics.
      // Running two separate kernels is net faster.

      // back up into shared memory for next iter
      if (j != twojmax) buf2[jju_shared_idx] = u_accum;

      // copy left side to right side with inversion symmetry VMK 4.4(2)
      // u[ma-j,mb-j] = (-1)^(ma-mb)*Conj([u[ma,mb))
      // We can avoid this if we're on the last row for an integer j
      if (!(n_ma % 2 == 1 && (mb+1) == n_mb)) {

        int sign_factor = ((ma%2==0)?1:-1)*(mb%2==0?1:-1);
        const int jjup_flip = jju+(j+1-mb)*(j+1)-(ma+1);
        const int jju_shared_flip = max_m_tile*team_rank+(j+1-mb)*(j+1)-(ma+1);

        if (sign_factor == 1) {
          u_accum.im = -u_accum.im;
        } else {
          u_accum.re = -u_accum.re;
        }
        ulist(jjup_flip,iatom,jnbor) = u_accum;
        if (j != twojmax) buf2[jju_shared_flip] = u_accum;
      }
    });
    // In CUDA backend,
    // ThreadVectorRange has a __syncwarp (appropriately masked for 
    // vector lengths < 32) implict at the end

    // swap double buffers
    auto tmp = buf1; buf1 = buf2; buf2 = tmp;
    //std::swap(buf1, buf2); // throws warnings
    

  }
}

// CPU version
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_uarray_cpu(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom, int jnbor,
                         double x, double y, double z,
                         double z0, double r)
{
  double r0inv;
  double a_r, b_r, a_i, b_i;
  double rootpq;

  // compute Cayley-Klein parameters for unit quaternion

  r0inv = 1.0 / sqrt(r * r + z0 * z0);
  a_r = r0inv * z0;
  a_i = -r0inv * z;
  b_r = r0inv * y;
  b_i = -r0inv * x;

  // VMK Section 4.8.2

  ulist(0,iatom,jnbor).re = 1.0;
  ulist(0,iatom,jnbor).im = 0.0;

  for (int j = 1; j <= twojmax; j++) {
    int jju = idxu_block[j];
    int jjup = idxu_block[j-1];

    // fill in left side of matrix layer from previous layer

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,(j+2)/2),
        [&] (const int& mb) {
    //for (int mb = 0; 2*mb <= j; mb++) {
      const int jju_index = jju+mb+mb*j;
      ulist(jju_index,iatom,jnbor).re = 0.0;
      ulist(jju_index,iatom,jnbor).im = 0.0;

      for (int ma = 0; ma < j; ma++) {
        const int jju_index = jju+mb+mb*j+ma;
        const int jjup_index = jjup+mb*j+ma;
        rootpq = rootpqarray(j - ma,j - mb);
        ulist(jju_index,iatom,jnbor,jju).re +=
          rootpq *
          (a_r * ulist(jjup_index,iatom,jnbor).re +
           a_i * ulist(jjup_index,iatom,jnbor).im);
        ulist(jju_index,iatom,jnbor).im +=
          rootpq *
          (a_r * ulist(jjup_index,iatom,jnbor).im -
           a_i * ulist(jjup_index,iatom,jnbor).re);

        rootpq = rootpqarray(ma + 1,j - mb);
        ulist(jju_index+1,iatom,jnbor).re =
          -rootpq *
          (b_r * ulist(jjup_index,iatom,jnbor).re +
           b_i * ulist(jjup_index,iatom,jnbor).im);
        ulist(jju_index+1,iatom,jnbor).im =
          -rootpq *
          (b_r * ulist(jjup_index,iatom,jnbor).im -
           b_i * ulist(jjup_index,iatom,jnbor).re);
      }
    });

    // copy left side to right side with inversion symmetry VMK 4.4(2)
    // u[ma-j,mb-j] = (-1)^(ma-mb)*Conj([u[ma,mb))

    jju = idxu_block[j];
    jjup = jju+(j+1)*(j+1)-1;
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,(j+2)/2),
        [&] (const int& mb) {
//    for (int mb = 0; 2*mb <= j; mb++) {
      int mbpar = (mb)%2==0?1:-1;
      int mapar = mbpar;
      for (int ma = 0; ma <= j; ma++) {
        const int jju_index = jju+mb*(j+1)+ma;
        const int jjup_index = jjup-mb*(j+1)-ma;
        if (mapar == 1) {
          ulist(jjup_index,iatom,jnbor).re = ulist(jju_index,iatom,jnbor).re;
          ulist(jjup_index,iatom,jnbor).im = -ulist(jju_index,iatom,jnbor).im;
        } else {
          ulist(jjup_index,iatom,jnbor).re = -ulist(jju_index,iatom,jnbor).re;
          ulist(jjup_index,iatom,jnbor).im = ulist(jju_index,iatom,jnbor).im;
        }
        mapar = -mapar;
      }
    });
  }
}

/* ----------------------------------------------------------------------
   compute derivatives of Wigner U-functions for one neighbor
   see comments in compute_uarray()
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_duarray(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom, int jnbor,
                          double x, double y, double z,
                          double z0, double r, double dz0dr,
                          double wj, double rcut)
{

  // get shared memory offset
  const int max_m_tile = (twojmax+1)*(twojmax+1);
  const int team_rank = team.team_rank();

  // double buffer for ulist
  SNAcomplex* ulist_buf1 = (SNAcomplex*)team.team_shmem( ).get_shmem(team.team_size()*max_m_tile*sizeof(SNAcomplex), 0);
  SNAcomplex* ulist_buf2 = (SNAcomplex*)team.team_shmem( ).get_shmem(team.team_size()*max_m_tile*sizeof(SNAcomplex), 0);

  // double buffer for dulist
  SNAcomplex* dulist_buf1 = (SNAcomplex*)team.team_shmem( ).get_shmem(team.team_size()*max_m_tile*sizeof(SNAcomplex), 0);
  SNAcomplex* dulist_buf2 = (SNAcomplex*)team.team_shmem( ).get_shmem(team.team_size()*max_m_tile*sizeof(SNAcomplex), 0);

  const double sfac = wj * compute_sfac(r, rcut);
  const double dsfac = wj * compute_dsfac(r, rcut);

  const double rinv = 1.0 / r;

  // extract a single unit vector
  const double u = (dir == 0 ? x * rinv : dir == 1 ? y * rinv : z * rinv);

  // Compute Cayley-Klein parameters for unit quaternion

  const double r0inv = 1.0 / sqrt(r * r + z0 * z0);

  const SNAcomplex a = { r0inv * z0, -r0inv * z };
  const SNAcomplex b = { r0inv * y, -r0inv * x };

  const double dr0invdr = -r0inv * r0inv * r0inv * (r + z0 * dz0dr);
  const double dr0inv = dr0invdr * u;
  const double dz0 = dz0dr * u;

  const SNAcomplex da = { dz0 * r0inv + z0 * dr0inv,
                              - z * dr0inv + (dir == 2 ? - r0inv : 0.) };

  const SNAcomplex db = { y * dr0inv + (dir==1?r0inv:0.),
                              -x * dr0inv + (dir==0?-r0inv:0.) };

  // single has a warp barrier at the end
  Kokkos::single(Kokkos::PerThread(team), [=]() {
    dulist(0,iatom,jnbor,dir) = { dsfac * u, 0. }; // fold in chain rule here
    ulist_buf1[max_m_tile*team_rank] = {1., 0.};
    dulist_buf1[max_m_tile*team_rank] = {0., 0.};
  });


  for (int j = 1; j <= twojmax; j++) {
    int jju = idxu_block[j];
    int jjup = idxu_block[j-1];

    // flatten the loop over ma,mb

    // for (int ma = 0; ma <= j; ma++)
    const int n_ma = j+1;
    // for (int mb = 0; 2*mb <= j; mb++)
    const int n_mb = j/2+1;

    const int total_iters = n_ma * n_mb;

    //for (int m = 0; m < total_iters; m++) {
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, total_iters),
      [&] (const int m) {

      // ma fast, mb slow
      int ma = m % n_ma;
      int mb = m / n_ma;

      const int jju_index = jju+mb+mb*j+ma;

      // index into shared memory
      const int jju_shared_idx = max_m_tile*team_rank+mb+mb*j+ma;
      const int jjup_shared_idx = max_m_tile*team_rank+mb*j+ma;

      // Need to compute and accumulate both u and du (mayhaps, we could probably
      // balance some read and compute by reading u each time).
      SNAcomplex u_accum = { 0., 0. };
      SNAcomplex du_accum = { 0., 0. };

      const double rootpq1 = rootpqarray(j - ma, j - mb);
      const SNAcomplex u_up1 = (ma < j)?rootpq1*ulist_buf1[jjup_shared_idx]:SNAcomplex(0.,0.);
      caconjxpy(a, u_up1, u_accum);

      const double rootpq2 = -rootpqarray(ma, j - mb);
      const SNAcomplex u_up2 = (ma > 0)?rootpq2*ulist_buf1[jjup_shared_idx-1]:SNAcomplex(0.,0.);
      caconjxpy(b, u_up2, u_accum);

      // No need to save u_accum to global memory
      if (j != twojmax) ulist_buf2[jju_shared_idx] = u_accum;

      // Next, spin up du_accum
      const SNAcomplex du_up1 = (ma < j) ? rootpq1*dulist_buf1[jjup_shared_idx] : SNAcomplex(0.,0.);
      caconjxpy(da, u_up1, du_accum);
      caconjxpy(a, du_up1, du_accum);

      const SNAcomplex du_up2 = (ma > 0) ? rootpq2*dulist_buf1[jjup_shared_idx-1] : SNAcomplex(0.,0.);
      caconjxpy(db, u_up2, du_accum);
      caconjxpy(b, du_up2, du_accum);

      dulist(jju_index,iatom,jnbor,dir) = ((dsfac * u)*u_accum) + (sfac*du_accum);

      if (j != twojmax) dulist_buf2[jju_shared_idx] = du_accum;

      // copy left side to right side with inversion symmetry VMK 4.4(2)
      // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])

      int sign_factor = ((ma%2==0)?1:-1)*(mb%2==0?1:-1);
      const int jjup_flip = jju+(j+1-mb)*(j+1)-(ma+1);
      const int jju_shared_flip = max_m_tile*team_rank+(j+1-mb)*(j+1)-(ma+1);

      if (sign_factor == 1) {
        //ulist_alt(iatom,jnbor,jjup_flip).re = u_accum.re;
        //ulist_alt(iatom,jnbor,jjup_flip).im = -u_accum.im;
        u_accum.im = -u_accum.im;
        du_accum.im = -du_accum.im;
      } else {
        //ulist_alt(iatom,jnbor,jjup_flip).re = -u_accum.re;
        //ulist_alt(iatom,jnbor,jjup_flip).im = u_accum.im;
        u_accum.re = -u_accum.re;
        du_accum.re = -du_accum.re;
      }

      dulist(jjup_flip,iatom,jnbor,dir) = ((dsfac * u)*u_accum) + (sfac*du_accum);

      if (j != twojmax) {
        ulist_buf2[jju_shared_flip] = u_accum;
        dulist_buf2[jju_shared_flip] = du_accum;
      }

    });

    // swap buffers
    auto tmp = ulist_buf1; ulist_buf1 = ulist_buf2; ulist_buf2 = tmp;
    tmp = dulist_buf1; dulist_buf1 = dulist_buf2; dulist_buf2 = tmp;
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_duarray_cpu(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom, int jnbor,
                          double x, double y, double z,
                          double z0, double r, double dz0dr,
                          double wj, double rcut)
{
double r0inv;
  double a_r, a_i, b_r, b_i;
  double da_r[3], da_i[3], db_r[3], db_i[3];
  double dz0[3], dr0inv[3], dr0invdr;
  double rootpq;

  double rinv = 1.0 / r;
  double ux = x * rinv;
  double uy = y * rinv;
  double uz = z * rinv;

  r0inv = 1.0 / sqrt(r * r + z0 * z0);
  a_r = z0 * r0inv;
  a_i = -z * r0inv;
  b_r = y * r0inv;
  b_i = -x * r0inv;

  dr0invdr = -r0inv * r0inv * r0inv * (r + z0 * dz0dr);

  dr0inv[0] = dr0invdr * ux;
  dr0inv[1] = dr0invdr * uy;
  dr0inv[2] = dr0invdr * uz;

  dz0[0] = dz0dr * ux;
  dz0[1] = dz0dr * uy;
  dz0[2] = dz0dr * uz;

  for (int k = 0; k < 3; k++) {
    da_r[k] = dz0[k] * r0inv + z0 * dr0inv[k];
    da_i[k] = -z * dr0inv[k];
  }

  da_i[2] += -r0inv;

  for (int k = 0; k < 3; k++) {
    db_r[k] = y * dr0inv[k];
    db_i[k] = -x * dr0inv[k];
  }

  db_i[0] += -r0inv;
  db_r[1] += r0inv;

  dulist(0,iatom,jnbor,0).re = 0.0;
  dulist(0,iatom,jnbor,1).re = 0.0;
  dulist(0,iatom,jnbor,2).re = 0.0;
  dulist(0,iatom,jnbor,0).im = 0.0;
  dulist(0,iatom,jnbor,1).im = 0.0;
  dulist(0,iatom,jnbor,2).im = 0.0;

  for (int j = 1; j <= twojmax; j++) {
    int jju = idxu_block[j];
    int jjup = idxu_block[j-1];
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,(j+2)/2),
        [&] (const int& mb) {
    //for (int mb = 0; 2*mb <= j; mb++) {
      const int jju_index = jju+mb+mb*j;
      dulist(jju_index,iatom,jnbor,0).re = 0.0;
      dulist(jju_index,iatom,jnbor,1).re = 0.0;
      dulist(jju_index,iatom,jnbor,2).re = 0.0;
      dulist(jju_index,iatom,jnbor,0).im = 0.0;
      dulist(jju_index,iatom,jnbor,1).im = 0.0;
      dulist(jju_index,iatom,jnbor,2).im = 0.0;

      for (int ma = 0; ma < j; ma++) {
        const int jju_index = jju+mb+mb*j+ma;
        const int jjup_index = jjup+mb*j+ma;
        rootpq = rootpqarray(j - ma,j - mb);
        for (int k = 0; k < 3; k++) {
          dulist(jju_index,iatom,jnbor,k).re +=
            rootpq * (da_r[k] * ulist(jjup_index,iatom,jnbor).re +
                      da_i[k] * ulist(jjup_index,iatom,jnbor).im +
                      a_r * dulist(jjup_index,iatom,jnbor,k).re +
                      a_i * dulist(jjup_index,iatom,jnbor,k).im);
          dulist(jju_index,iatom,jnbor,k).im +=
            rootpq * (da_r[k] * ulist(jjup_index,iatom,jnbor).im -
                      da_i[k] * ulist(jjup_index,iatom,jnbor).re +
                      a_r * dulist(jjup_index,iatom,jnbor,k).im -
                      a_i * dulist(jjup_index,iatom,jnbor,k).re);
        }

        rootpq = rootpqarray(ma + 1,j - mb);
        for (int k = 0; k < 3; k++) {
          dulist(jju_index+1,iatom,jnbor,k).re =
            -rootpq * (db_r[k] * ulist(jjup_index,iatom,jnbor).re +
                       db_i[k] * ulist(jjup_index,iatom,jnbor).im +
                       b_r * dulist(jjup_index,iatom,jnbor,k).re +
                       b_i * dulist(jjup_index,iatom,jnbor,k).im);
          dulist(jju_index+1,iatom,jnbor,k).im =
            -rootpq * (db_r[k] * ulist(jjup_index,iatom,jnbor).im -
                       db_i[k] * ulist(jjup_index,iatom,jnbor).re +
                       b_r * dulist(jjup_index,iatom,jnbor,k).im -
                       b_i * dulist(jjup_index,iatom,jnbor,k).re);
        }
      }
    });

    // copy left side to right side with inversion symmetry VMK 4.4(2)
    // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])

    jju = idxu_block[j];
    jjup = jju+(j+1)*(j+1)-1;
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,(j+2)/2),
        [&] (const int& mb) {
//    for (int mb = 0; 2*mb <= j; mb++) {
      int mbpar = (mb)%2==0?1:-1;
      int mapar = mbpar;
      for (int ma = 0; ma <= j; ma++) {
        const int jju_index = jju+mb*(j+1)+ma;
        const int jjup_index = jjup-mb*(j+1)-ma;
        if (mapar == 1) {
          for (int k = 0; k < 3; k++) {
            dulist(jjup_index,iatom,jnbor,k).re = dulist(jju_index,iatom,jnbor,k).re;
            dulist(jjup_index,iatom,jnbor,k).im = -dulist(jju_index,iatom,jnbor,k).im;
          }
        } else {
          for (int k = 0; k < 3; k++) {
            dulist(jjup_index,iatom,jnbor,k).re = -dulist(jju_index,iatom,jnbor,k).re;
            dulist(jjup_index,iatom,jnbor,k).im = dulist(jju_index,iatom,jnbor,k).im;
          }
        }
        mapar = -mapar;
      }
    });
  }

  double sfac = compute_sfac(r, rcut);
  double dsfac = compute_dsfac(r, rcut);

  sfac *= wj;
  dsfac *= wj;

  for (int j = 0; j <= twojmax; j++) {
    int jju = idxu_block[j];
    for (int mb = 0; 2*mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++) {
        dulist(jju,iatom,jnbor,0).re = dsfac * ulist(jju,iatom,jnbor).re * ux +
                                  sfac * dulist(jju,iatom,jnbor,0).re;
        dulist(jju,iatom,jnbor,0).im = dsfac * ulist(jju,iatom,jnbor).im * ux +
                                  sfac * dulist(jju,iatom,jnbor,0).im;
        dulist(jju,iatom,jnbor,1).re = dsfac * ulist(jju,iatom,jnbor).re * uy +
                                  sfac * dulist(jju,iatom,jnbor,1).re;
        dulist(jju,iatom,jnbor,1).im = dsfac * ulist(jju,iatom,jnbor).im * uy +
                                  sfac * dulist(jju,iatom,jnbor,1).im;
        dulist(jju,iatom,jnbor,2).re = dsfac * ulist(jju,iatom,jnbor).re * uz +
                                  sfac * dulist(jju,iatom,jnbor,2).re;
        dulist(jju,iatom,jnbor,2).im = dsfac * ulist(jju,iatom,jnbor).im * uz +
                                  sfac * dulist(jju,iatom,jnbor,2).im;

        jju++;
      }
  }
}

/* ----------------------------------------------------------------------
   factorial n, wrapper for precomputed table
------------------------------------------------------------------------- */

template<class DeviceType>
inline
double SNAKokkos<DeviceType>::factorial(int n)
{
  //if (n < 0 || n > nmaxfactorial) {
  //  char str[128];
  //  sprintf(str, "Invalid argument to factorial %d", n);
  //  error->all(FLERR, str);
  //}

  return nfac_table[n];
}

/* ----------------------------------------------------------------------
   factorial n table, size SNA::nmaxfactorial+1
------------------------------------------------------------------------- */

template<class DeviceType>
const double SNAKokkos<DeviceType>::nfac_table[] = {
  1,
  1,
  2,
  6,
  24,
  120,
  720,
  5040,
  40320,
  362880,
  3628800,
  39916800,
  479001600,
  6227020800,
  87178291200,
  1307674368000,
  20922789888000,
  355687428096000,
  6.402373705728e+15,
  1.21645100408832e+17,
  2.43290200817664e+18,
  5.10909421717094e+19,
  1.12400072777761e+21,
  2.5852016738885e+22,
  6.20448401733239e+23,
  1.5511210043331e+25,
  4.03291461126606e+26,
  1.08888694504184e+28,
  3.04888344611714e+29,
  8.8417619937397e+30,
  2.65252859812191e+32,
  8.22283865417792e+33,
  2.63130836933694e+35,
  8.68331761881189e+36,
  2.95232799039604e+38,
  1.03331479663861e+40,
  3.71993326789901e+41,
  1.37637530912263e+43,
  5.23022617466601e+44,
  2.03978820811974e+46,
  8.15915283247898e+47,
  3.34525266131638e+49,
  1.40500611775288e+51,
  6.04152630633738e+52,
  2.65827157478845e+54,
  1.1962222086548e+56,
  5.50262215981209e+57,
  2.58623241511168e+59,
  1.24139155925361e+61,
  6.08281864034268e+62,
  3.04140932017134e+64,
  1.55111875328738e+66,
  8.06581751709439e+67,
  4.27488328406003e+69,
  2.30843697339241e+71,
  1.26964033536583e+73,
  7.10998587804863e+74,
  4.05269195048772e+76,
  2.35056133128288e+78,
  1.3868311854569e+80,
  8.32098711274139e+81,
  5.07580213877225e+83,
  3.14699732603879e+85,
  1.98260831540444e+87,
  1.26886932185884e+89,
  8.24765059208247e+90,
  5.44344939077443e+92,
  3.64711109181887e+94,
  2.48003554243683e+96,
  1.71122452428141e+98,
  1.19785716699699e+100,
  8.50478588567862e+101,
  6.12344583768861e+103,
  4.47011546151268e+105,
  3.30788544151939e+107,
  2.48091408113954e+109,
  1.88549470166605e+111,
  1.45183092028286e+113,
  1.13242811782063e+115,
  8.94618213078297e+116,
  7.15694570462638e+118,
  5.79712602074737e+120,
  4.75364333701284e+122,
  3.94552396972066e+124,
  3.31424013456535e+126,
  2.81710411438055e+128,
  2.42270953836727e+130,
  2.10775729837953e+132,
  1.85482642257398e+134,
  1.65079551609085e+136,
  1.48571596448176e+138,
  1.3520015276784e+140,
  1.24384140546413e+142,
  1.15677250708164e+144,
  1.08736615665674e+146,
  1.03299784882391e+148,
  9.91677934870949e+149,
  9.61927596824821e+151,
  9.42689044888324e+153,
  9.33262154439441e+155,
  9.33262154439441e+157,
  9.42594775983835e+159,
  9.61446671503512e+161,
  9.90290071648618e+163,
  1.02990167451456e+166,
  1.08139675824029e+168,
  1.14628056373471e+170,
  1.22652020319614e+172,
  1.32464181945183e+174,
  1.44385958320249e+176,
  1.58824554152274e+178,
  1.76295255109024e+180,
  1.97450685722107e+182,
  2.23119274865981e+184,
  2.54355973347219e+186,
  2.92509369349301e+188,
  3.3931086844519e+190,
  3.96993716080872e+192,
  4.68452584975429e+194,
  5.5745857612076e+196,
  6.68950291344912e+198,
  8.09429852527344e+200,
  9.8750442008336e+202,
  1.21463043670253e+205,
  1.50614174151114e+207,
  1.88267717688893e+209,
  2.37217324288005e+211,
  3.01266001845766e+213,
  3.8562048236258e+215,
  4.97450422247729e+217,
  6.46685548922047e+219,
  8.47158069087882e+221,
  1.118248651196e+224,
  1.48727070609069e+226,
  1.99294274616152e+228,
  2.69047270731805e+230,
  3.65904288195255e+232,
  5.01288874827499e+234,
  6.91778647261949e+236,
  9.61572319694109e+238,
  1.34620124757175e+241,
  1.89814375907617e+243,
  2.69536413788816e+245,
  3.85437071718007e+247,
  5.5502938327393e+249,
  8.04792605747199e+251,
  1.17499720439091e+254,
  1.72724589045464e+256,
  2.55632391787286e+258,
  3.80892263763057e+260,
  5.71338395644585e+262,
  8.62720977423323e+264,
  1.31133588568345e+267,
  2.00634390509568e+269,
  3.08976961384735e+271,
  4.78914290146339e+273,
  7.47106292628289e+275,
  1.17295687942641e+278,
  1.85327186949373e+280,
  2.94670227249504e+282,
  4.71472363599206e+284,
  7.59070505394721e+286,
  1.22969421873945e+289,
  2.0044015765453e+291,
  3.28721858553429e+293,
  5.42391066613159e+295,
  9.00369170577843e+297,
  1.503616514865e+300, // nmaxfactorial = 167
};

/* ----------------------------------------------------------------------
   the function delta given by VMK Eq. 8.2(1)
------------------------------------------------------------------------- */

template<class DeviceType>
inline
double SNAKokkos<DeviceType>::deltacg(int j1, int j2, int j)
{
  double sfaccg = factorial((j1 + j2 + j) / 2 + 1);
  return sqrt(factorial((j1 + j2 - j) / 2) *
              factorial((j1 - j2 + j) / 2) *
              factorial((-j1 + j2 + j) / 2) / sfaccg);
}

/* ----------------------------------------------------------------------
   assign Clebsch-Gordan coefficients using
   the quasi-binomial formula VMK 8.2.1(3)
------------------------------------------------------------------------- */

template<class DeviceType>
inline
void SNAKokkos<DeviceType>::init_clebsch_gordan()
{
  auto h_cglist = Kokkos::create_mirror_view(cglist);

  double sum,dcg,sfaccg;
  int m, aa2, bb2, cc2;
  int ifac;

  int idxcg_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        for (int m1 = 0; m1 <= j1; m1++) {
          aa2 = 2 * m1 - j1;

          for (int m2 = 0; m2 <= j2; m2++) {

            // -c <= cc <= c

            bb2 = 2 * m2 - j2;
            m = (aa2 + bb2 + j) / 2;

            if(m < 0 || m > j) {
              h_cglist[idxcg_count] = 0.0;
              idxcg_count++;
              continue;
            }

            sum = 0.0;

            for (int z = MAX(0, MAX(-(j - j2 + aa2)
                                    / 2, -(j - j1 - bb2) / 2));
                 z <= MIN((j1 + j2 - j) / 2,
                          MIN((j1 - aa2) / 2, (j2 + bb2) / 2));
                 z++) {
              ifac = z % 2 ? -1 : 1;
              sum += ifac /
                (factorial(z) *
                 factorial((j1 + j2 - j) / 2 - z) *
                 factorial((j1 - aa2) / 2 - z) *
                 factorial((j2 + bb2) / 2 - z) *
                 factorial((j - j2 + aa2) / 2 + z) *
                 factorial((j - j1 - bb2) / 2 + z));
            }

            cc2 = 2 * m - j;
            dcg = deltacg(j1, j2, j);
            sfaccg = sqrt(factorial((j1 + aa2) / 2) *
                          factorial((j1 - aa2) / 2) *
                          factorial((j2 + bb2) / 2) *
                          factorial((j2 - bb2) / 2) *
                          factorial((j  + cc2) / 2) *
                          factorial((j  - cc2) / 2) *
                          (j + 1));
            
            h_cglist[idxcg_count] = sum * dcg * sfaccg;
            idxcg_count++;
          }
        }
      }
  Kokkos::deep_copy(cglist,h_cglist);
}

/* ----------------------------------------------------------------------
   pre-compute table of sqrt[p/m2], p, q = 1,twojmax
   the p = 0, q = 0 entries are allocated and skipped for convenience.
------------------------------------------------------------------------- */

template<class DeviceType>
inline
void SNAKokkos<DeviceType>::init_rootpqarray()
{
  auto h_rootpqarray = Kokkos::create_mirror_view(rootpqarray);
  for (int p = 1; p <= twojmax; p++)
    for (int q = 1; q <= twojmax; q++)
      h_rootpqarray(p,q) = sqrt(static_cast<double>(p)/q);
  Kokkos::deep_copy(rootpqarray,h_rootpqarray);
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
inline
int SNAKokkos<DeviceType>::compute_ncoeff()
{
  int ncount;

  ncount = 0;

  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2;
           j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) ncount++;

  return ncount;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double SNAKokkos<DeviceType>::compute_sfac(double r, double rcut)
{
  if (switch_flag == 0) return 1.0;
  if (switch_flag == 1) {
    if(r <= rmin0) return 1.0;
    else if(r > rcut) return 0.0;
    else {
      double rcutfac = MY_PI / (rcut - rmin0);
      return 0.5 * (cos((r - rmin0) * rcutfac) + 1.0);
    }
  }
  return 0.0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double SNAKokkos<DeviceType>::compute_dsfac(double r, double rcut)
{
  if (switch_flag == 0) return 0.0;
  if (switch_flag == 1) {
    if(r <= rmin0) return 0.0;
    else if(r > rcut) return 0.0;
    else {
      double rcutfac = MY_PI / (rcut - rmin0);
      return -0.5 * sin((r - rmin0) * rcutfac) * rcutfac;
    }
  }
  return 0.0;
}

/* ---------------------------------------------------------------------- */

// efficient complex FMA (i.e., y += a x)
template<class DeviceType>
KOKKOS_FORCEINLINE_FUNCTION
void SNAKokkos<DeviceType>::caxpy(const SNAcomplex& a, const SNAcomplex& x, SNAcomplex& y) {
  y.re += a.re * x.re;
  y.re -= a.im * x.im;
  y.im += a.im * x.re;
  y.im += a.re * x.im;
}

/* ---------------------------------------------------------------------- */

// efficient complex FMA, conjugate of scalar (i.e.) y += (a.re - i a.im) x)
template<class DeviceType>
KOKKOS_FORCEINLINE_FUNCTION
void SNAKokkos<DeviceType>::caconjxpy(const SNAcomplex& a, const SNAcomplex& x, SNAcomplex& y) {
  y.re += a.re * x.re;
  y.re += a.im * x.im;
  y.im -= a.im * x.re;
  y.im += a.re * x.im;
}

/* ---------------------------------------------------------------------- */

// set direction of batched Duidrj
template<class DeviceType>
KOKKOS_FORCEINLINE_FUNCTION
void SNAKokkos<DeviceType>::set_dir(int dir_) {
  dir = dir_;
}

/* ----------------------------------------------------------------------
   memory usage of arrays
------------------------------------------------------------------------- */

template<class DeviceType>
double SNAKokkos<DeviceType>::memory_usage()
{
  int jdimpq = twojmax + 2;
  int jdim = twojmax + 1;
  double bytes;

  bytes = 0;

  bytes += jdimpq*jdimpq * sizeof(double);               // pqarray
  bytes += idxcg_max * sizeof(double);                   // cglist

  bytes += natom * idxu_max * sizeof(double) * 2;        // ulist
  bytes += natom * idxu_max * sizeof(double) * 2;        // ulisttot
  if (!Kokkos::Impl::is_same<typename DeviceType::array_layout,Kokkos::LayoutRight>::value)
    bytes += natom * idxu_max * sizeof(double) * 2;        // ulisttot_lr
  bytes += natom * idxu_max * 3 * sizeof(double) * 2;    // dulist
                                                       
  bytes += natom * idxz_max * sizeof(double) * 2;        // zlist
  bytes += natom * idxb_max * sizeof(double);            // blist
  bytes += natom * idxu_max * sizeof(double) * 2;        // ylist

  bytes += jdim * jdim * jdim * sizeof(int);             // idxcg_block
  bytes += jdim * sizeof(int);                           // idxu_block
  bytes += jdim * jdim * jdim * sizeof(int);             // idxz_block
  bytes += jdim * jdim * jdim * sizeof(int);             // idxb_block

  bytes += idxz_max * 10 * sizeof(int);                  // idxz
  bytes += idxb_max * 3 * sizeof(int);                   // idxb

  bytes += jdim * sizeof(double);                        // bzero

  bytes += natom * nmax * 3 * sizeof(double);            // rij
  bytes += natom * nmax * sizeof(int);                   // inside
  bytes += natom * nmax * sizeof(double);                // wj
  bytes += natom * nmax * sizeof(double);                // rcutij
  bytes += natom * nmax * idxu_max * sizeof(double) * 2; // ulist_ij

  return bytes;
}

} // namespace LAMMPS_NS
