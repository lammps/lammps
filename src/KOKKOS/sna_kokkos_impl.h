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

  blist = t_sna_2d("sna:blist",natom,idxb_max);
  ulisttot = t_sna_2c("sna:ulisttot",natom,idxu_max);
  if (!Kokkos::Impl::is_same<typename DeviceType::array_layout,Kokkos::LayoutRight>::value) 
    ulisttot_lr = t_sna_2c_lr("sna:ulisttot_lr",natom,idxu_max);
  zlist = t_sna_2c("sna:zlist",natom,idxz_max);

  ulist = t_sna_3c("sna:ulist",natom,nmax,idxu_max);
  ylist = t_sna_2c_lr("sna:ylist",natom,idxu_max);

  dulist = t_sna_4c("sna:dulist",natom,nmax,idxu_max);
}

/* ----------------------------------------------------------------------
 *    compute Ui by summing over neighbors j
 *    ------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::pre_ui(const int& iatom)
{
  //if(team.team_rank() == 0) {
    zero_uarraytot(iatom);
    //Kokkos::single(Kokkos::PerThread(team), [&] (){
    addself_uarraytot(iatom,wself);
    //});
  //}
  //team.team_barrier();
}

/* ----------------------------------------------------------------------
   compute Ui by summing over neighbors j
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
  add_uarraytot(team, iatom, jnbor, r, wj(iatom,jnbor), rcutij(iatom,jnbor));
}

/* ----------------------------------------------------------------------
   compute Zi by summing over products of Ui
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_zi(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,idxz_max),
      [&] (const int& jjz) {
  //for(int jjz = 0; jjz < idxz_max; jjz++) {
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

    zlist(iatom,jjz).re = 0.0; 
    zlist(iatom,jjz).im = 0.0;

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
        suma1_r += cgblock[icga] * (ulisttot(iatom,jju1+ma1).re * ulisttot(iatom,jju2+ma2).re - ulisttot(iatom,jju1+ma1).im * ulisttot(iatom,jju2+ma2).im);
        suma1_i += cgblock[icga] * (ulisttot(iatom,jju1+ma1).re * ulisttot(iatom,jju2+ma2).im + ulisttot(iatom,jju1+ma1).im * ulisttot(iatom,jju2+ma2).re);
        ma1++;
        ma2--;
        icga += j2;
      } // end loop over ia

      zlist(iatom,jjz).re += cgblock[icgb] * suma1_r;
      zlist(iatom,jjz).im += cgblock[icgb] * suma1_i;

      jju1 += j1+1;
      jju2 -= j2+1;
      icgb += j2;
    } // end loop over ib

  }); // end loop over jjz
}

/* ----------------------------------------------------------------------
   compute Yi from Ui without storing Zi, looping over zlist indices
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::zero_yi(const int& iatom)
{
    for (int j = 0; j < idxu_max; j++)
      ylist(iatom,j) = {0.0,0.0};
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
      suma1_r += cgblock[icga] * (ulisttot_lr(iatom,jju1+ma1).re * ulisttot_lr(iatom,jju2+ma2).re - ulisttot_lr(iatom,jju1+ma1).im * ulisttot_lr(iatom,jju2+ma2).im);
      suma1_i += cgblock[icga] * (ulisttot_lr(iatom,jju1+ma1).re * ulisttot_lr(iatom,jju2+ma2).im + ulisttot_lr(iatom,jju1+ma1).im * ulisttot_lr(iatom,jju2+ma2).re);
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
      if (j2 == j) betaj = 3*beta(iatom,jjb);
      else betaj = 2*beta(iatom,jjb);
    } else betaj = beta(iatom,jjb); 
  } else if (j >= j2) {
    const int jjb = idxb_block(j,j2,j1);
    if (j2 == j) betaj = 2*beta(iatom,jjb)*(j1+1)/(j+1.0);
    else betaj = beta(iatom,jjb)*(j1+1)/(j+1.0);
  } else {
    const int jjb = idxb_block(j2,j,j1);
    betaj = beta(iatom,jjb)*(j1+1)/(j+1.0);
  }

  Kokkos::atomic_add(&(ylist(iatom,jju).re), betaj*ztmp_r);
  Kokkos::atomic_add(&(ylist(iatom,jju).im), betaj*ztmp_i);
}

/* ----------------------------------------------------------------------
   compute dEidRj
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_deidrj(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom, int jnbor)
{
  t_scalar3<double> sum;

  Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,twojmax+1),
      [&] (const int& j, t_scalar3<double>& sum_tmp) {
  //for(int j = 0; j <= twojmax; j++) {
    int jju = idxu_block[j];

    for(int mb = 0; 2*mb < j; mb++)
      for(int ma = 0; ma <= j; ma++) {
        sum_tmp.x += dulist(iatom,jnbor,jju,0).re * ylist(iatom,jju).re + dulist(iatom,jnbor,jju,0).im * ylist(iatom,jju).im;
        sum_tmp.y += dulist(iatom,jnbor,jju,1).re * ylist(iatom,jju).re + dulist(iatom,jnbor,jju,1).im * ylist(iatom,jju).im;
        sum_tmp.z += dulist(iatom,jnbor,jju,2).re * ylist(iatom,jju).re + dulist(iatom,jnbor,jju,2).im * ylist(iatom,jju).im;
        jju++;
      } //end loop over ma mb

    // For j even, handle middle column

    if (j%2 == 0) {

      int mb = j/2;
      for(int ma = 0; ma < mb; ma++) {
        sum_tmp.x += dulist(iatom,jnbor,jju,0).re * ylist(iatom,jju).re + dulist(iatom,jnbor,jju,0).im * ylist(iatom,jju).im;
        sum_tmp.y += dulist(iatom,jnbor,jju,1).re * ylist(iatom,jju).re + dulist(iatom,jnbor,jju,1).im * ylist(iatom,jju).im;
        sum_tmp.z += dulist(iatom,jnbor,jju,2).re * ylist(iatom,jju).re + dulist(iatom,jnbor,jju,2).im * ylist(iatom,jju).im;
        jju++;
      }

      //int ma = mb;
      sum_tmp.x += (dulist(iatom,jnbor,jju,0).re * ylist(iatom,jju).re + dulist(iatom,jnbor,jju,0).im * ylist(iatom,jju).im)*0.5;
      sum_tmp.y += (dulist(iatom,jnbor,jju,1).re * ylist(iatom,jju).re + dulist(iatom,jnbor,jju,1).im * ylist(iatom,jju).im)*0.5;
      sum_tmp.z += (dulist(iatom,jnbor,jju,2).re * ylist(iatom,jju).re + dulist(iatom,jnbor,jju,2).im * ylist(iatom,jju).im)*0.5;
    } // end if jeven

  },sum); // end loop over j
  //}

  Kokkos::single(Kokkos::PerThread(team), [&] () {
    dedr(iatom,jnbor,0) = sum.x*2.0;
    dedr(iatom,jnbor,1) = sum.y*2.0;
    dedr(iatom,jnbor,2) = sum.z*2.0;
  });

}

/* ----------------------------------------------------------------------
   compute Bi by summing conj(Ui)*Zi
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
          ulisttot(iatom,jju_index).re * zlist(iatom,jjz_index).re +
          ulisttot(iatom,jju_index).im * zlist(iatom,jjz_index).im;
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
          ulisttot(iatom,jju_index).re * zlist(iatom,jjz_index).re +
          ulisttot(iatom,jju_index).im * zlist(iatom,jjz_index).im;
      },sumzu_temp); // end loop over ma
      sumzu += sumzu_temp;

      const int ma = mb;
      const int jju_index = jju+(mb-1)*(j+1)+(j+1)+ma;
      const int jjz_index = jjz+(mb-1)*(j+1)+(j+1)+ma;
      sumzu += 0.5*
        (ulisttot(iatom,jju_index).re * zlist(iatom,jjz_index).re +
         ulisttot(iatom,jju_index).im * zlist(iatom,jjz_index).im);
    } // end if jeven

    Kokkos::single(Kokkos::PerThread(team), [&] () {
      sumzu *= 2.0;

      // apply bzero shift

      if (bzero_flag)
        sumzu -= bzero[j];

      blist(iatom,jjb) = sumzu;
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

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::zero_uarraytot(const int& iatom)
{
  {
    //Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,ulisttot.extent(1)),
    //    [&] (const int& i) {
    for (int i = 0; i < ulisttot.extent(1); i++) {
      ulisttot(iatom,i).re = 0.0;
      ulisttot(iatom,i).im = 0.0;
    }
    //});
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::addself_uarraytot(const int& iatom, const double& wself_in)
{
  //Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,twojmax+1),
  //  [&] (const int& j) {
  for (int j = 0; j <= twojmax; j++) {
    int jju = idxu_block[j];
    for (int ma = 0; ma <= j; ma++) {
      ulisttot(iatom,jju).re = wself_in;
      ulisttot(iatom,jju).im = 0.0;
      jju += j+2;
    }
  }//});
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

  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,ulisttot.extent(1)),
      [&] (const int& i) {
    Kokkos::atomic_add(&(ulisttot(iatom,i).re), sfac * ulist(iatom,jnbor,i).re);
    Kokkos::atomic_add(&(ulisttot(iatom,i).im), sfac * ulist(iatom,jnbor,i).im);
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

  ulist(iatom,jnbor,0).re = 1.0;
  ulist(iatom,jnbor,0).im = 0.0;

  for (int j = 1; j <= twojmax; j++) {
    int jju = idxu_block[j];
    int jjup = idxu_block[j-1];

    // fill in left side of matrix layer from previous layer

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,(j+2)/2),
        [&] (const int& mb) {
    //for (int mb = 0; 2*mb <= j; mb++) {
      const int jju_index = jju+mb+mb*j;
      ulist(iatom,jnbor,jju_index).re = 0.0;
      ulist(iatom,jnbor,jju_index).im = 0.0;

      for (int ma = 0; ma < j; ma++) {
        const int jju_index = jju+mb+mb*j+ma;
        const int jjup_index = jjup+mb*j+ma;
        rootpq = rootpqarray(j - ma,j - mb);
        ulist(iatom,jnbor,jju_index).re +=
          rootpq *
          (a_r * ulist(iatom,jnbor,jjup_index).re +
           a_i * ulist(iatom,jnbor,jjup_index).im);
        ulist(iatom,jnbor,jju_index).im +=
          rootpq *
          (a_r * ulist(iatom,jnbor,jjup_index).im -
           a_i * ulist(iatom,jnbor,jjup_index).re);

        rootpq = rootpqarray(ma + 1,j - mb);
        ulist(iatom,jnbor,jju_index+1).re =
          -rootpq *
          (b_r * ulist(iatom,jnbor,jjup_index).re +
           b_i * ulist(iatom,jnbor,jjup_index).im);
        ulist(iatom,jnbor,jju_index+1).im =
          -rootpq *
          (b_r * ulist(iatom,jnbor,jjup_index).im -
           b_i * ulist(iatom,jnbor,jjup_index).re);
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
          ulist(iatom,jnbor,jjup_index).re = ulist(iatom,jnbor,jju_index).re;
          ulist(iatom,jnbor,jjup_index).im = -ulist(iatom,jnbor,jju_index).im;
        } else {
          ulist(iatom,jnbor,jjup_index).re = -ulist(iatom,jnbor,jju_index).re;
          ulist(iatom,jnbor,jjup_index).im = ulist(iatom,jnbor,jju_index).im;
        }
        mapar = -mapar;
      }
    });
  }
}

template<class DeviceType>
void SNAKokkos<DeviceType>::transpose_ulisttot()
{
  UlisttotHelper<typename DeviceType::array_layout,decltype(ulisttot_lr),decltype(ulisttot)>::transpose(ulisttot_lr,ulisttot);
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

  dr0invdr = -pow(r0inv, 3.0) * (r + z0 * dz0dr);

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

  dulist(iatom,jnbor,0,0).re = 0.0;
  dulist(iatom,jnbor,0,1).re = 0.0;
  dulist(iatom,jnbor,0,2).re = 0.0;
  dulist(iatom,jnbor,0,0).im = 0.0;
  dulist(iatom,jnbor,0,1).im = 0.0;
  dulist(iatom,jnbor,0,2).im = 0.0;

  for (int j = 1; j <= twojmax; j++) {
    int jju = idxu_block[j];
    int jjup = idxu_block[j-1];
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,(j+2)/2),
        [&] (const int& mb) {
    //for (int mb = 0; 2*mb <= j; mb++) {
      const int jju_index = jju+mb+mb*j;
      dulist(iatom,jnbor,jju_index,0).re = 0.0;
      dulist(iatom,jnbor,jju_index,1).re = 0.0;
      dulist(iatom,jnbor,jju_index,2).re = 0.0;
      dulist(iatom,jnbor,jju_index,0).im = 0.0;
      dulist(iatom,jnbor,jju_index,1).im = 0.0;
      dulist(iatom,jnbor,jju_index,2).im = 0.0;

      for (int ma = 0; ma < j; ma++) {
        const int jju_index = jju+mb+mb*j+ma;
        const int jjup_index = jjup+mb*j+ma;
        rootpq = rootpqarray(j - ma,j - mb);
        for (int k = 0; k < 3; k++) {
          dulist(iatom,jnbor,jju_index,k).re +=
            rootpq * (da_r[k] * ulist(iatom,jnbor,jjup_index).re +
                      da_i[k] * ulist(iatom,jnbor,jjup_index).im +
                      a_r * dulist(iatom,jnbor,jjup_index,k).re +
                      a_i * dulist(iatom,jnbor,jjup_index,k).im);
          dulist(iatom,jnbor,jju_index,k).im +=
            rootpq * (da_r[k] * ulist(iatom,jnbor,jjup_index).im -
                      da_i[k] * ulist(iatom,jnbor,jjup_index).re +
                      a_r * dulist(iatom,jnbor,jjup_index,k).im -
                      a_i * dulist(iatom,jnbor,jjup_index,k).re);
        }

        rootpq = rootpqarray(ma + 1,j - mb);
        for (int k = 0; k < 3; k++) {
          dulist(iatom,jnbor,jju_index+1,k).re =
            -rootpq * (db_r[k] * ulist(iatom,jnbor,jjup_index).re +
                       db_i[k] * ulist(iatom,jnbor,jjup_index).im +
                       b_r * dulist(iatom,jnbor,jjup_index,k).re +
                       b_i * dulist(iatom,jnbor,jjup_index,k).im);
          dulist(iatom,jnbor,jju_index+1,k).im =
            -rootpq * (db_r[k] * ulist(iatom,jnbor,jjup_index).im -
                       db_i[k] * ulist(iatom,jnbor,jjup_index).re +
                       b_r * dulist(iatom,jnbor,jjup_index,k).im -
                       b_i * dulist(iatom,jnbor,jjup_index,k).re);
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
            dulist(iatom,jnbor,jjup_index,k).re = dulist(iatom,jnbor,jju_index,k).re;
            dulist(iatom,jnbor,jjup_index,k).im = -dulist(iatom,jnbor,jju_index,k).im;
          }
        } else {
          for (int k = 0; k < 3; k++) {
            dulist(iatom,jnbor,jjup_index,k).re = -dulist(iatom,jnbor,jju_index,k).re;
            dulist(iatom,jnbor,jjup_index,k).im = dulist(iatom,jnbor,jju_index,k).im;
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
        dulist(iatom,jnbor,jju,0).re = dsfac * ulist(iatom,jnbor,jju).re * ux +
                                  sfac * dulist(iatom,jnbor,jju,0).re;
        dulist(iatom,jnbor,jju,0).im = dsfac * ulist(iatom,jnbor,jju).im * ux +
                                  sfac * dulist(iatom,jnbor,jju,0).im;
        dulist(iatom,jnbor,jju,1).re = dsfac * ulist(iatom,jnbor,jju).re * uy +
                                  sfac * dulist(iatom,jnbor,jju,1).re;
        dulist(iatom,jnbor,jju,1).im = dsfac * ulist(iatom,jnbor,jju).im * uy +
                                  sfac * dulist(iatom,jnbor,jju,1).im;
        dulist(iatom,jnbor,jju,2).re = dsfac * ulist(iatom,jnbor,jju).re * uz +
                                  sfac * dulist(iatom,jnbor,jju,2).re;
        dulist(iatom,jnbor,jju,2).im = dsfac * ulist(iatom,jnbor,jju).im * uz +
                                  sfac * dulist(iatom,jnbor,jju,2).im;

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
