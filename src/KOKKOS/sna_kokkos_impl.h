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

  //create_twojmax_arrays();

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

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
SNAKokkos<DeviceType>::SNAKokkos(const SNAKokkos<DeviceType>& sna, const typename Kokkos::TeamPolicy<DeviceType>::member_type& team) {
  wself = sna.wself;

  rfac0 = sna.rfac0;
  rmin0 = sna.rmin0;
  switch_flag = sna.switch_flag;
  bzero_flag = sna.bzero_flag;

  twojmax = sna.twojmax;

  ncoeff = sna.ncoeff;
  nmax = sna.nmax;
  idxz = sna.idxz;
  idxb = sna.idxb;
  idxcg_max = sna.idxcg_max;
  idxu_max = sna.idxu_max;
  idxz_max = sna.idxz_max;
  idxb_max = sna.idxb_max;
  idxcg_block = sna.idxcg_block;
  idxu_block = sna.idxu_block;
  idxz_block = sna.idxz_block;
  idxb_block = sna.idxb_block;
  cglist = sna.cglist;
  rootpqarray = sna.rootpqarray;
  bzero = sna.bzero;
  create_team_scratch_arrays(team);
  create_thread_scratch_arrays(team);
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
  idxb = Kokkos::View<SNAKK_BINDICES*, DeviceType>("SNAKokkos::idxb",idxb_max);
  auto h_idxb = Kokkos::create_mirror_view(idxb);
  
  idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) {
          h_idxb[idxb_count].j1 = j1;
          h_idxb[idxb_count].j2 = j2;
          h_idxb[idxb_count].j = j;
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
  idxz = Kokkos::View<SNAKK_ZINDICES*, DeviceType>("SNAKokkos::idxz",idxz_max);
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
            h_idxz[idxz_count].j1 = j1;
            h_idxz[idxz_count].j2 = j2;
            h_idxz[idxz_count].j = j;
            h_idxz[idxz_count].ma1min = MAX(0, (2 * ma - j - j2 + j1) / 2);
            h_idxz[idxz_count].ma2max = (2 * ma - j - (2 * h_idxz[idxz_count].ma1min - j1) + j2) / 2;
            h_idxz[idxz_count].na = MIN(j1, (2 * ma - j + j2 + j1) / 2) - h_idxz[idxz_count].ma1min + 1;
            h_idxz[idxz_count].mb1min = MAX(0, (2 * mb - j - j2 + j1) / 2);
            h_idxz[idxz_count].mb2max = (2 * mb - j - (2 * h_idxz[idxz_count].mb1min - j1) + j2) / 2;
            h_idxz[idxz_count].nb = MIN(j1, (2 * mb - j + j2 + j1) / 2) - h_idxz[idxz_count].mb1min + 1;

            // apply to z(j1,j2,j,ma,mb) to unique element of y(j)

            const int jju = h_idxu_block[j] + (j+1)*mb + ma;
            h_idxz[idxz_count].jju = jju;

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
void SNAKokkos<DeviceType>::grow_rij(int newnmax)
{
  if(newnmax <= nmax) return;
  nmax = newnmax;
}
/* ----------------------------------------------------------------------
   compute Ui by summing over neighbors j
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_ui(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int jnum)
{
  double rsq, r, x, y, z, z0, theta0;

  // utot(j,ma,mb) = 0 for all j,ma,ma
  // utot(j,ma,ma) = 1 for all j,ma
  // for j in neighbors of i:
  //   compute r0 = (x,y,z,z0)
  //   utot(j,ma,mb) += u(r0;j,ma,mb) for all j,ma,mb

  if(team.team_rank() == 0) {
    zero_uarraytot(team);
    //Kokkos::single(Kokkos::PerThread(team), [&] (){
    addself_uarraytot(team,wself);
    //});
  }
  team.team_barrier();

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,jnum),
      [&] (const int& j) {
  //for(int j = 0; j < jnum; j++) {
    x = rij(j,0);
    y = rij(j,1);
    z = rij(j,2);
    rsq = x * x + y * y + z * z;
    r = sqrt(rsq);

    theta0 = (r - rmin0) * rfac0 * MY_PI / (rcutij[j] - rmin0);
    //    theta0 = (r - rmin0) * rscale0;
    z0 = r / tan(theta0);

    compute_uarray(team,x, y, z, z0, r);
    //Kokkos::single(Kokkos::PerThread(team), [&] (){
    add_uarraytot(team,r, wj[j], rcutij[j], j);
    //});
  });

}

/* ----------------------------------------------------------------------
   compute Zi by summing over products of Ui
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_zi(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,idxz_max),
      [&] (const int& jjz) {
  //for(int jjz = 0; jjz < idxz_max; jjz++) {
    const int j1 = idxz[jjz].j1;
    const int j2 = idxz[jjz].j2;
    const int j = idxz[jjz].j;
    const int ma1min = idxz[jjz].ma1min;
    const int ma2max = idxz[jjz].ma2max;
    const int na = idxz[jjz].na;
    const int mb1min = idxz[jjz].mb1min;
    const int mb2max = idxz[jjz].mb2max;
    const int nb = idxz[jjz].nb;

    const double* cgblock = cglist.data() + idxcg_block(j1,j2,j);

    zlist_r[jjz] = 0.0; 
    zlist_i[jjz] = 0.0;

    int jju1 = idxu_block[j1] + (j1+1)*mb1min;
    int jju2 = idxu_block[j2] + (j2+1)*mb2max;
    int icgb = mb1min*(j2+1) + mb2max;
    for(int ib = 0; ib < nb; ib++) {

      double suma1_r = 0.0;
      double suma1_i = 0.0;

      const double* u1_r = ulisttot_r.data() + jju1;
      const double* u1_i = ulisttot_i.data() + jju1;
      const double* u2_r = ulisttot_r.data() + jju2;
      const double* u2_i = ulisttot_i.data() + jju2;

      int ma1 = ma1min;
      int ma2 = ma2max;
      int icga = ma1min*(j2+1) + ma2max;
      for(int ia = 0; ia < na; ia++) {
        suma1_r += cgblock[icga] * (u1_r[ma1] * u2_r[ma2] - u1_i[ma1] * u2_i[ma2]);
        suma1_i += cgblock[icga] * (u1_r[ma1] * u2_i[ma2] + u1_i[ma1] * u2_r[ma2]);
        ma1++;
        ma2--;
        icga += j2;
      } // end loop over ia

      zlist_r[jjz] += cgblock[icgb] * suma1_r;
      zlist_i[jjz] += cgblock[icgb] * suma1_i;

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
void SNAKokkos<DeviceType>::compute_yi(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team,
 const Kokkos::View<F_FLOAT**, DeviceType> &beta, const int ii)
{
  double betaj;

  {
    double* const ptr = ylist_r.data();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,ylist_r.span()),
        [&] (const int& i) {
      ptr[i] = 0.0;
    });
  }
  {
    double* const ptr = ylist_i.data();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,ylist_i.span()),
        [&] (const int& i) {
      ptr[i] = 0.0;
    });
  }

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,idxz_max),
      [&] (const int& jjz) {
  //for(int jjz = 0; jjz < idxz_max; jjz++) {
    const int j1 = idxz[jjz].j1;
    const int j2 = idxz[jjz].j2;
    const int j = idxz[jjz].j;
    const int ma1min = idxz[jjz].ma1min;
    const int ma2max = idxz[jjz].ma2max;
    const int na = idxz[jjz].na;
    const int mb1min = idxz[jjz].mb1min;
    const int mb2max = idxz[jjz].mb2max;
    const int nb = idxz[jjz].nb;

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

      const double* u1_r = ulisttot_r.data() + jju1;
      const double* u1_i = ulisttot_i.data() + jju1;
      const double* u2_r = ulisttot_r.data() + jju2;
      const double* u2_i = ulisttot_i.data() + jju2;

      int ma1 = ma1min;
      int ma2 = ma2max;
      int icga = ma1min*(j2+1) + ma2max;

      for(int ia = 0; ia < na; ia++) {
        suma1_r += cgblock[icga] * (u1_r[ma1] * u2_r[ma2] - u1_i[ma1] * u2_i[ma2]);
        suma1_i += cgblock[icga] * (u1_r[ma1] * u2_i[ma2] + u1_i[ma1] * u2_r[ma2]);        ma1++;
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
    // find right y_list[jju] and beta(ii,jjb) entries
    // multiply and divide by j+1 factors
    // account for multiplicity of 1, 2, or 3

    const int jju = idxz[jjz].jju;

  // pick out right beta value

    if (j >= j1) {
      const int jjb = idxb_block(j1,j2,j);
      if (j1 == j) {
        if (j2 == j) betaj = 3*beta(ii,jjb);
        else betaj = 2*beta(ii,jjb);
      } else betaj = beta(ii,jjb); 
    } else if (j >= j2) {
      const int jjb = idxb_block(j,j2,j1);
      if (j2 == j) betaj = 2*beta(ii,jjb)*(j1+1)/(j+1.0);
      else betaj = beta(ii,jjb)*(j1+1)/(j+1.0);
    } else {
      const int jjb = idxb_block(j2,j,j1);
      betaj = beta(ii,jjb)*(j1+1)/(j+1.0);
    }

  Kokkos::single(Kokkos::PerThread(team), [&] () {
    ylist_r[jju] += betaj*ztmp_r;
    ylist_i[jju] += betaj*ztmp_i;
  });

  }); // end loop over jjz
}

/* ----------------------------------------------------------------------
   compute dEidRj
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_deidrj(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, double* dedr)
{
  t_scalar3<double> sum;

  // TODO: which loop is faster to parallelize?
  Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,twojmax+1),
      [&] (const int& j, t_scalar3<double>& sum_tmp) {
  //for(int j = 0; j <= twojmax; j++) {
    int jju = idxu_block[j];

    for(int mb = 0; 2*mb < j; mb++)
      for(int ma = 0; ma <= j; ma++) {
        sum_tmp.x += dulist_r(jju,0) * ylist_r[jju] + dulist_i(jju,0) * ylist_i[jju];
        sum_tmp.y += dulist_r(jju,1) * ylist_r[jju] + dulist_i(jju,1) * ylist_i[jju];
        sum_tmp.z += dulist_r(jju,2) * ylist_r[jju] + dulist_i(jju,2) * ylist_i[jju];
        jju++;
      } //end loop over ma mb

    // For j even, handle middle column

    if (j%2 == 0) {

      int mb = j/2;
      for(int ma = 0; ma < mb; ma++) {
        sum_tmp.x += dulist_r(jju,0) * ylist_r[jju] + dulist_i(jju,0) * ylist_i[jju];
        sum_tmp.y += dulist_r(jju,1) * ylist_r[jju] + dulist_i(jju,1) * ylist_i[jju];
        sum_tmp.z += dulist_r(jju,2) * ylist_r[jju] + dulist_i(jju,2) * ylist_i[jju];
        jju++;
      }

      //int ma = mb;
      sum_tmp.x += (dulist_r(jju,0) * ylist_r[jju] + dulist_i(jju,0) * ylist_i[jju])*0.5;
      sum_tmp.y += (dulist_r(jju,1) * ylist_r[jju] + dulist_i(jju,1) * ylist_i[jju])*0.5;
      sum_tmp.z += (dulist_r(jju,2) * ylist_r[jju] + dulist_i(jju,2) * ylist_i[jju])*0.5;
    } // end if jeven

  },sum); // end loop over j

  Kokkos::single(Kokkos::PerThread(team), [&] () {
    dedr[0] = sum.x*2.0;
    dedr[1] = sum.y*2.0;
    dedr[2] = sum.z*2.0;
  });

}

/* ----------------------------------------------------------------------
   compute Bi by summing conj(Ui)*Zi
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_bi(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team)
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
    const int j1 = idxb[jjb].j1;
    const int j2 = idxb[jjb].j2;
    const int j = idxb[jjb].j;

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
          ulisttot_r(jju_index) * zlist_r(jjz_index) +
          ulisttot_i(jju_index) * zlist_i(jjz_index);
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
          ulisttot_r(jju_index) * zlist_r(jjz_index) +
          ulisttot_i(jju_index) * zlist_i(jjz_index);
      },sumzu_temp); // end loop over ma
      sumzu += sumzu_temp;

      const int ma = mb;
      const int jju_index = jju+(mb-1)*(j+1)+(j+1)+ma;
      const int jjz_index = jjz+(mb-1)*(j+1)+(j+1)+ma;
      sumzu += 0.5*
        (ulisttot_r(jju_index) * zlist_r(jjz_index) +
         ulisttot_i(jju_index) * zlist_i(jjz_index));
    } // end if jeven

    Kokkos::single(Kokkos::PerThread(team), [&] () {
      sumzu *= 2.0;

      // apply bzero shift

      if (bzero_flag)
        sumzu -= bzero[j];

      blist(jjb) = sumzu;
    });
  });
      //} // end loop over j
    //} // end loop over j1, j2
}

/* ----------------------------------------------------------------------
   calculate derivative of Bi w.r.t. atom j
   variant using indexlist for j1,j2,j
   variant using symmetry relation
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_dbidrj(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team)
{
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        zdb = 0
  //        for mb = 0,...,jmid
  //          for ma = 0,...,j
  //            zdb +=
  //              Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)
  //        dbdr(j1,j2,j) += 2*zdb
  //        zdb = 0
  //        for mb1 = 0,...,j1mid
  //          for ma1 = 0,...,j1
  //            zdb +=
  //              Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)
  //        dbdr(j1,j2,j) += 2*zdb*(j+1)/(j1+1)
  //        zdb = 0
  //        for mb2 = 0,...,j2mid
  //          for ma2 = 0,...,j2
  //            zdb +=
  //              Conj(dudr(j2,ma2,mb2))*z(j1,j,j2,ma2,mb2)
  //        dbdr(j1,j2,j) += 2*zdb*(j+1)/(j2+1)


  // TODO:   double* dudr_r, *dudr_i;


  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,idxb_max),
          [&] (const int& jjb) {
  //for(int jjb = 0; jjb < idxb_max; jjb++) {
    const int j1 = idxb[jjb].j1;
    const int j2 = idxb[jjb].j2;
    const int j = idxb[jjb].j;

//    dbdr = dblist(jjb);
//    dbdr[0] = 0.0;
//    dbdr[1] = 0.0;
//    dbdr[2] = 0.0;

    t_scalar3<double> dbdr,sumzdu_r;
    // Sum terms Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)

    int jjz = idxz_block(j1,j2,j);
    int jju = idxu_block[j];

    for(int mb = 0; 2*mb < j; mb++)
      for(int ma = 0; ma <= j; ma++) {
        const int jju_index = jju+mb*(j+1)+ma;
        const int jjz_index = jjz+mb*(j+1)+ma;
        sumzdu_r.x += (dulist_r(jju_index,0) * zlist_r[jjz_index] + dulist_i(jju_index,0) * zlist_i[jjz_index]);
        sumzdu_r.y += (dulist_r(jju_index,1) * zlist_r[jjz_index] + dulist_i(jju_index,1) * zlist_i[jjz_index]);
        sumzdu_r.z += (dulist_r(jju_index,2) * zlist_r[jjz_index] + dulist_i(jju_index,2) * zlist_i[jjz_index]);
      } //end loop over ma mb

    // For j even, handle middle column

    if (j%2 == 0) {
      int mb = j/2;
      for(int ma = 0; ma <= mb; ma++) {
        const int jju_index = jju+(mb-1)*(j+1)+(j+1)+ma;
        const int jjz_index = jjz+(mb-1)*(j+1)+(j+1)+ma;
        sumzdu_r.x += (dulist_r(jju_index,0) * zlist_r[jjz_index] + dulist_i(jju_index,0) * zlist_i[jjz_index]);
        sumzdu_r.y += (dulist_r(jju_index,1) * zlist_r[jjz_index] + dulist_i(jju_index,1) * zlist_i[jjz_index]);
        sumzdu_r.z += (dulist_r(jju_index,2) * zlist_r[jjz_index] + dulist_i(jju_index,2) * zlist_i[jjz_index]);
      }
      int ma = mb;
      const int jju_index = jju+(mb-1)*(j+1)+(j+1)+ma;
      const int jjz_index = jjz+(mb-1)*(j+1)+(j+1)+ma;
      for(int k = 0; k < 3; k++) {
        sumzdu_r.x += (dulist_r(jju_index,0) * zlist_r[jjz] + dulist_i(jju_index,0) * zlist_i[jjz_index])*0.5;
        sumzdu_r.y += (dulist_r(jju_index,1) * zlist_r[jjz] + dulist_i(jju_index,1) * zlist_i[jjz_index])*0.5;
        sumzdu_r.z += (dulist_r(jju_index,2) * zlist_r[jjz] + dulist_i(jju_index,2) * zlist_i[jjz_index])*0.5;
      }
    } // end if jeven

      dbdr += 2.0*sumzdu_r;

    // Sum over Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)

    double j1fac = (j+1)/(j1+1.0);

    jjz = idxz_block(j,j2,j1);
    jju = idxu_block[j1];

    sumzdu_r.x = 0.0; sumzdu_r.y = 0.0; sumzdu_r.z = 0.0;

    for(int mb = 0; 2*mb < j1; mb++)
      for(int ma = 0; ma <= j1; ma++) {
        const int jju_index = jju+mb*(j1+1)+ma;
        const int jjz_index = jjz+mb*(j1+1)+ma;
        sumzdu_r.x += (dulist_r(jju_index,0) * zlist_r[jjz_index] + dulist_i(jju_index,0) * zlist_i[jjz_index]);
        sumzdu_r.y += (dulist_r(jju_index,1) * zlist_r[jjz_index] + dulist_i(jju_index,1) * zlist_i[jjz_index]);
        sumzdu_r.z += (dulist_r(jju_index,2) * zlist_r[jjz_index] + dulist_i(jju_index,2) * zlist_i[jjz_index]);
      } //end loop over ma1 mb1

    // For j1 even, handle middle column

    if (j1%2 == 0) {
      const int mb = j1/2;
      for(int ma = 0; ma <= mb; ma++) {
        const int jju_index = jju+(mb-1)*(j1+1)+(j1+1)+ma;
        const int jjz_index = jjz+(mb-1)*(j1+1)+(j1+1)+ma;
        sumzdu_r.x += (dulist_r(jju_index,0) * zlist_r[jjz_index] + dulist_i(jju_index,0) * zlist_i[jjz_index]);
        sumzdu_r.y += (dulist_r(jju_index,1) * zlist_r[jjz_index] + dulist_i(jju_index,1) * zlist_i[jjz_index]);
        sumzdu_r.z += (dulist_r(jju_index,2) * zlist_r[jjz_index] + dulist_i(jju_index,2) * zlist_i[jjz_index]);
      }
      int ma = mb;
      const int jju_index = jju+(mb-1)*(j1+1)+(j1+1)+ma;
      const int jjz_index = jjz+(mb-1)*(j1+1)+(j1+1)+ma;
      for(int k = 0; k < 3; k++) {
        sumzdu_r.x += (dulist_r(jju_index,0) * zlist_r[jjz] + dulist_i(jju_index,0) * zlist_i[jjz_index])*0.5;
        sumzdu_r.y += (dulist_r(jju_index,1) * zlist_r[jjz] + dulist_i(jju_index,1) * zlist_i[jjz_index])*0.5;
        sumzdu_r.z += (dulist_r(jju_index,2) * zlist_r[jjz] + dulist_i(jju_index,2) * zlist_i[jjz_index])*0.5;
      }
    } // end if j1even

      dbdr += 2.0*sumzdu_r*j1fac;

    // Sum over Conj(dudr(j2,ma2,mb2))*z(j1,j,j2,ma2,mb2)

    double j2fac = (j+1)/(j2+1.0);

    jjz = idxz_block(j,j1,j2);
    jju = idxu_block[j2];

    sumzdu_r.x = 0.0; sumzdu_r.y = 0.0; sumzdu_r.z = 0.0;

    for(int mb = 0; 2*mb < j2; mb++)
      for(int ma = 0; ma <= j2; ma++) {
        const int jju_index = jju+mb*(j2+1)+ma;
        const int jjz_index = jjz+mb*(j2+1)+ma;
        sumzdu_r.x += (dulist_r(jju_index,0) * zlist_r[jjz_index] + dulist_i(jju_index,0) * zlist_i[jjz_index]);
        sumzdu_r.y += (dulist_r(jju_index,1) * zlist_r[jjz_index] + dulist_i(jju_index,1) * zlist_i[jjz_index]);
        sumzdu_r.z += (dulist_r(jju_index,2) * zlist_r[jjz_index] + dulist_i(jju_index,2) * zlist_i[jjz_index]);
      } //end loop over ma2 mb2

    // For j2 even, handle middle column

    if (j2%2 == 0) {
      const int mb = j2/2;
      for(int ma = 0; ma <= mb; ma++) {
        const int jju_index = jju+(mb-1)*(j2+1)+(j2+1)+ma;
        const int jjz_index = jjz+(mb-1)*(j2+1)+(j2+1)+ma;
        sumzdu_r.x += (dulist_r(jju_index,0) * zlist_r[jjz_index] + dulist_i(jju_index,0) * zlist_i[jjz_index]);
        sumzdu_r.y += (dulist_r(jju_index,1) * zlist_r[jjz_index] + dulist_i(jju_index,1) * zlist_i[jjz_index]);
        sumzdu_r.z += (dulist_r(jju_index,2) * zlist_r[jjz_index] + dulist_i(jju_index,2) * zlist_i[jjz_index]);
      }
      int ma = mb;
      const int jju_index = jju+(mb-1)*(j2+1)+(j2+1)+ma;
      const int jjz_index = jjz+(mb-1)*(j2+1)+(j2+1)+ma;
      for(int k = 0; k < 3; k++) {
        sumzdu_r.x += (dulist_r(jju_index,0) * zlist_r[jjz] + dulist_i(jju_index,0) * zlist_i[jjz_index])*0.5;
        sumzdu_r.y += (dulist_r(jju_index,1) * zlist_r[jjz] + dulist_i(jju_index,1) * zlist_i[jjz_index])*0.5;
        sumzdu_r.z += (dulist_r(jju_index,2) * zlist_r[jjz] + dulist_i(jju_index,2) * zlist_i[jjz_index])*0.5;
      }
    } // end if j2even

    dbdr += 2.0*sumzdu_r*j2fac;
    dblist(jjb,0) = dbdr.x;
    dblist(jjb,1) = dbdr.y;
    dblist(jjb,2) = dbdr.z;

  }); //end loop over j1 j2 j
}

/* ----------------------------------------------------------------------
   calculate derivative of Ui w.r.t. atom j
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_duidrj(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team,
                         double* rij, double wj, double rcut, int jj)
{
  double rsq, r, x, y, z, z0, theta0, cs, sn;
  double dz0dr;

  x = rij[0];
  y = rij[1];
  z = rij[2];
  rsq = x * x + y * y + z * z;
  r = sqrt(rsq);
  double rscale0 = rfac0 * MY_PI / (rcut - rmin0);
  theta0 = (r - rmin0) * rscale0;
  cs = cos(theta0);
  sn = sin(theta0);
  z0 = r * cs / sn;
  dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;

  compute_duarray(team, x, y, z, z0, r, dz0dr, wj, rcut, jj);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::zero_uarraytot(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team)
{
  {
    double* const ptr = ulisttot_r.data();
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,ulisttot_r.span()),
        [&] (const int& i) {
      ptr[i] = 0.0;
    });
  }
  {
    double* const ptr = ulisttot_i.data();
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,ulisttot_i.span()),
        [&] (const int& i) {
      ptr[i] = 0.0;
    });
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::addself_uarraytot(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, double wself_in)
{
  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,twojmax+1),
    [&] (const int& j) {
  //for (int j = 0; j <= twojmax; j++)
    int jju = idxu_block[j];
    for (int ma = 0; ma <= j; ma++) {
      ulisttot_r[jju] = wself_in;
      ulisttot_i[jju] = 0.0;
      jju += j+2;
    }
  });
}

/* ----------------------------------------------------------------------
   add Wigner U-functions for one neighbor to the total
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::add_uarraytot(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team,
                                          double r, double wj, double rcut, int j)
{
  const double sfac = compute_sfac(r, rcut) * wj;

  const double* const ptr_r = ulist_r.data();
  const double* const ptr_i = ulist_i.data();
  double* const ptrtot_r = ulisttot_r.data();
  double* const ptrtot_i = ulisttot_i.data();

  Kokkos::View<double*,Kokkos::LayoutRight,DeviceType,Kokkos::MemoryTraits<Kokkos::Unmanaged>>
    ulist_r_j(ulist_r_ij,j,Kokkos::ALL);
  Kokkos::View<double*,Kokkos::LayoutRight,DeviceType,Kokkos::MemoryTraits<Kokkos::Unmanaged>>
    ulist_i_j(ulist_i_ij,j,Kokkos::ALL);

  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,ulisttot_r.span()),
      [&] (const int& i) {
    Kokkos::atomic_add(ptrtot_r+i, sfac * ptr_r[i]);
    Kokkos::atomic_add(ptrtot_i+i, sfac * ptr_i[i]);

    ulist_r_j(i) = ulist_r(i);
    ulist_i_j(i) = ulist_i(i);
  });
}

/* ----------------------------------------------------------------------
   compute Wigner U-functions for one neighbor
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_uarray(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team,
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

  ulist_r[0] = 1.0;
  ulist_i[0] = 0.0;

  for (int j = 1; j <= twojmax; j++) {
    int jju = idxu_block[j];
    int jjup = idxu_block[j-1];

    // fill in left side of matrix layer from previous layer

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,(j+2)/2),
        [&] (const int& mb) {
    //for (int mb = 0; 2*mb <= j; mb++) {
      const int jju_index = jju+mb+mb*j;
      ulist_r[jju_index] = 0.0;
      ulist_i[jju_index] = 0.0;

      for (int ma = 0; ma < j; ma++) {
        const int jju_index = jju+mb+mb*j+ma;
        const int jjup_index = jjup+mb*j+ma;
        rootpq = rootpqarray(j - ma,j - mb);
        ulist_r[jju_index] +=
          rootpq *
          (a_r * ulist_r[jjup_index] +
           a_i * ulist_i[jjup_index]);
        ulist_i[jju_index] +=
          rootpq *
          (a_r * ulist_i[jjup_index] -
           a_i * ulist_r[jjup_index]);

        rootpq = rootpqarray(ma + 1,j - mb);
        ulist_r[jju_index+1] =
          -rootpq *
          (b_r * ulist_r[jjup_index] +
           b_i * ulist_i[jjup_index]);
        ulist_i[jju_index+1] =
          -rootpq *
          (b_r * ulist_i[jjup_index] -
           b_i * ulist_r[jjup_index]);
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
          ulist_r[jjup_index] = ulist_r[jju_index];
          ulist_i[jjup_index] = -ulist_i[jju_index];
        } else {
          ulist_r[jjup_index] = -ulist_r[jju_index];
          ulist_i[jjup_index] = ulist_i[jju_index];
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
void SNAKokkos<DeviceType>::compute_duarray(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team,
                          double x, double y, double z,
                          double z0, double r, double dz0dr,
                          double wj, double rcut, int jj)
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

  Kokkos::View<double*,Kokkos::LayoutRight,DeviceType,Kokkos::MemoryTraits<Kokkos::Unmanaged>>
    ulist_r(ulist_r_ij,jj,Kokkos::ALL);
  Kokkos::View<double*,Kokkos::LayoutRight,DeviceType,Kokkos::MemoryTraits<Kokkos::Unmanaged>>
    ulist_i(ulist_i_ij,jj,Kokkos::ALL);

  dulist_r(0,0) = 0.0;
  dulist_r(0,1) = 0.0;
  dulist_r(0,2) = 0.0;
  dulist_i(0,0) = 0.0;
  dulist_i(0,1) = 0.0;
  dulist_i(0,2) = 0.0;

  for (int j = 1; j <= twojmax; j++) {
    int jju = idxu_block[j];
    int jjup = idxu_block[j-1];
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,(j+2)/2),
        [&] (const int& mb) {
    //for (int mb = 0; 2*mb <= j; mb++) {
      const int jju_index = jju+mb+mb*j;
      dulist_r(jju_index,0) = 0.0;
      dulist_r(jju_index,1) = 0.0;
      dulist_r(jju_index,2) = 0.0;
      dulist_i(jju_index,0) = 0.0;
      dulist_i(jju_index,1) = 0.0;
      dulist_i(jju_index,2) = 0.0;

      for (int ma = 0; ma < j; ma++) {
        const int jju_index = jju+mb+mb*j+ma;
        const int jjup_index = jjup+mb*j+ma;
        rootpq = rootpqarray(j - ma,j - mb);
        for (int k = 0; k < 3; k++) {
          dulist_r(jju_index,k) +=
            rootpq * (da_r[k] * ulist_r[jjup_index] +
                      da_i[k] * ulist_i[jjup_index] +
                      a_r * dulist_r(jjup_index,k) +
                      a_i * dulist_i(jjup_index,k));
          dulist_i(jju_index,k) +=
            rootpq * (da_r[k] * ulist_i[jjup_index] -
                      da_i[k] * ulist_r[jjup_index] +
                      a_r * dulist_i(jjup_index,k) -
                      a_i * dulist_r(jjup_index,k));
        }

        rootpq = rootpqarray(ma + 1,j - mb);
        for (int k = 0; k < 3; k++) {
          dulist_r(jju_index+1,k) =
            -rootpq * (db_r[k] * ulist_r[jjup_index] +
                       db_i[k] * ulist_i[jjup_index] +
                       b_r * dulist_r(jjup_index,k) +
                       b_i * dulist_i(jjup_index,k));
          dulist_i(jju_index+1,k) =
            -rootpq * (db_r[k] * ulist_i[jjup_index] -
                       db_i[k] * ulist_r[jjup_index] +
                       b_r * dulist_i(jjup_index,k) -
                       b_i * dulist_r(jjup_index,k));
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
            dulist_r(jjup_index,k) = dulist_r(jju_index,k);
            dulist_i(jjup_index,k) = -dulist_i(jju_index,k);
          }
        } else {
          for (int k = 0; k < 3; k++) {
            dulist_r(jjup_index,k) = -dulist_r(jju_index,k);
            dulist_i(jjup_index,k) = dulist_i(jju_index,k);
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

  for (int j = 0; j <= twojmax; j++) { //TODO: parallelize one of these loops
    int jju = idxu_block[j];
    for (int mb = 0; 2*mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++) {
        dulist_r(jju,0) = dsfac * ulist_r[jju] * ux +
                                  sfac * dulist_r(jju,0);
        dulist_i(jju,0) = dsfac * ulist_i[jju] * ux +
                                  sfac * dulist_i(jju,0);
        dulist_r(jju,1) = dsfac * ulist_r[jju] * uy +
                                  sfac * dulist_r(jju,1);
        dulist_i(jju,1) = dsfac * ulist_i[jju] * uy +
                                  sfac * dulist_i(jju,1);
        dulist_r(jju,2) = dsfac * ulist_r[jju] * uz +
                                  sfac * dulist_r(jju,2);
        dulist_i(jju,2) = dsfac * ulist_i[jju] * uz +
                                  sfac * dulist_i(jju,2);

        jju++;
      }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::create_team_scratch_arrays(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team)
{
  ulisttot_r_a = ulisttot_r = t_sna_1d(team.team_scratch(1),idxu_max);
  ulisttot_i_a = ulisttot_i = t_sna_1d(team.team_scratch(1),idxu_max);
  ylist_r = t_sna_1d(team.team_scratch(1),idxu_max);
  ylist_i = t_sna_1d(team.team_scratch(1),idxu_max);
  zlist_r = t_sna_1d(team.team_scratch(1),idxz_max);
  zlist_i = t_sna_1d(team.team_scratch(1),idxz_max);
  blist = t_sna_1d(team.team_scratch(1),idxb_max);

  rij = t_sna_2d(team.team_scratch(1),nmax,3);
  rcutij = t_sna_1d(team.team_scratch(1),nmax);
  wj = t_sna_1d(team.team_scratch(1),nmax);
  inside = t_sna_1i(team.team_scratch(1),nmax);
  ulist_r_ij = t_sna_2d(team.team_scratch(1),nmax,idxu_max);
  ulist_i_ij = t_sna_2d(team.team_scratch(1),nmax,idxu_max);
}

template<class DeviceType>
inline
T_INT SNAKokkos<DeviceType>::size_team_scratch_arrays() {
  T_INT size = 0;

  size += t_sna_1d::shmem_size(idxu_max)*2; // ulisttot
  size += t_sna_1d::shmem_size(idxu_max)*2; // ylist
  size += t_sna_1d::shmem_size(idxz_max)*2; // zlist
  size += t_sna_1d::shmem_size(idxb_max); // blist

  size += t_sna_2d::shmem_size(nmax,3); // rij
  size += t_sna_1d::shmem_size(nmax); // rcutij
  size += t_sna_1d::shmem_size(nmax); // wj
  size += t_sna_1i::shmem_size(nmax); // inside
  size += t_sna_2d::shmem_size(nmax,idxu_max)*2; // ulist_ij

  return size;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::create_thread_scratch_arrays(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team)
{
  dblist = t_sna_2d(team.thread_scratch(1),idxb_max,3);
  ulist_r = t_sna_1d(team.thread_scratch(1),idxu_max);
  ulist_i = t_sna_1d(team.thread_scratch(1),idxu_max);
  dulist_r = t_sna_2d(team.thread_scratch(1),idxu_max,3);
  dulist_i = t_sna_2d(team.thread_scratch(1),idxu_max,3);
}

template<class DeviceType>
inline
T_INT SNAKokkos<DeviceType>::size_thread_scratch_arrays() {
  T_INT size = 0;

  size += t_sna_2d::shmem_size(idxb_max,3); // dblist
  size += t_sna_1d::shmem_size(idxu_max)*2; // ulist
  size += t_sna_2d::shmem_size(idxu_max,3)*2; // dulist
  return size;
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

  bytes += idxu_max * sizeof(double) * 2;                // ulist
  bytes += idxu_max * sizeof(double) * 2;                // ulisttot
  bytes += idxu_max * 3 * sizeof(double) * 2;            // dulist

  bytes += idxz_max * sizeof(double) * 2;                // zlist
  bytes += idxb_max * sizeof(double);                    // blist
  bytes += idxb_max * 3 * sizeof(double);                // dblist
  bytes += idxu_max * sizeof(double) * 2;                // ylist

  bytes += jdim * jdim * jdim * sizeof(int);             // idxcg_block
  bytes += jdim * sizeof(int);                           // idxu_block
  bytes += jdim * jdim * jdim * sizeof(int);             // idxz_block
  bytes += jdim * jdim * jdim * sizeof(int);             // idxb_block

  bytes += idxz_max * sizeof(SNAKK_ZINDICES);            // idxz
  bytes += idxb_max * sizeof(SNAKK_BINDICES);            // idxb

  bytes += jdim * sizeof(double);                        // bzero

  bytes += nmax * 3 * sizeof(double);                    // rij
  bytes += nmax * sizeof(int);                           // inside
  bytes += nmax * sizeof(double);                        // wj
  bytes += nmax * sizeof(double);                        // rcutij
  bytes += nmax * idxu_max * sizeof(double) * 2;         // ulist_ij

  return bytes;
}

} // namespace LAMMPS_NS
