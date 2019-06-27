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

  int jdim = twojmax + 1;

  cgarray = t_sna_5d("SNAKokkos::cgarray",jdim,jdim,jdim,jdim,jdim);
  rootpqarray = t_sna_2d("SNAKokkos::rootpqarray",jdim+1,jdim+1);

  if (bzero_flag) {
    bzero = Kokkos::View<double*, Kokkos::LayoutRight, DeviceType>("sna:bzero",jdim);
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
  idxj = sna.idxj;
  idxj_max = sna.idxj_max;
  idxj_full = sna.idxj_full;
  idxj_full_max = sna.idxj_full_max;
  cgarray = sna.cgarray;
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
  int idxj_count = 0;
  int idxj_full_count = 0;

  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2) {
	if (j >= j1) idxj_count++;
	idxj_full_count++;
      }

  // indexList can be changed here

  idxj = Kokkos::View<SNAKK_LOOPINDICES*, DeviceType>("SNAKokkos::idxj",idxj_count);
  idxj_full = Kokkos::View<SNAKK_LOOPINDICES*, DeviceType>("SNAKokkos::idxj_full",idxj_full_count);
  auto h_idxj = Kokkos::create_mirror_view(idxj);
  auto h_idxj_full = Kokkos::create_mirror_view(idxj_full);

  idxj_max = idxj_count;
  idxj_full_max = idxj_full_count;

  idxj_count = 0;
  idxj_full_count = 0;

  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2) {
	if (j >= j1) {
	  h_idxj[idxj_count].j1 = j1;
	  h_idxj[idxj_count].j2 = j2;
	  h_idxj[idxj_count].j = j;
	  idxj_count++;
	}
	h_idxj_full[idxj_full_count].j1 = j1;
	h_idxj_full[idxj_full_count].j2 = j2;
	h_idxj_full[idxj_full_count].j = j;
	idxj_full_count++;
      }
  Kokkos::deep_copy(idxj,h_idxj);
  Kokkos::deep_copy(idxj_full,h_idxj_full);

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
    add_uarraytot(team,r, wj[j], rcutij[j]);
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
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        for ma = 0,...,j
  //          for mb = 0,...,jmid
  //            z(j1,j2,j,ma,mb) = 0
  //            for ma1 = Max(0,ma+(j1-j2-j)/2),Min(j1,ma+(j1+j2-j)/2)
  //              sumb1 = 0
  //              ma2 = ma-ma1+(j1+j2-j)/2;
  //              for mb1 = Max(0,mb+(j1-j2-j)/2),Min(j1,mb+(j1+j2-j)/2)
  //                mb2 = mb-mb1+(j1+j2-j)/2;
  //                sumb1 += cg(j1,mb1,j2,mb2,j) *
  //                  u(j1,ma1,mb1) * u(j2,ma2,mb2)
  //              z(j1,j2,j,ma,mb) += sumb1*cg(j1,ma1,j2,ma2,j)

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &starttime);
#endif

  // compute_dbidrj() requires full j1/j2/j chunk of z elements
  // use zarray j1/j2 symmetry

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,idxj_full_max),
      [&] (const int& idx) {
    const int j1 = idxj_full(idx).j1;
    const int j2 = idxj_full(idx).j2;
    const int j =  idxj_full(idx).j;

    const int bound = (j+2)/2;
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,(j+1)*bound),
        [&] (const int mbma ) {
        //for(int mb = 0; 2*mb <= j; mb++)
          //for(int ma = 0; ma <= j; ma++) {
      const int ma = mbma%(j+1);
      const int mb = mbma/(j+1);

            //zarray_r(j1,j2,j,ma,mb) = 0.0;
            //zarray_i(j1,j2,j,ma,mb) = 0.0;
      double z_r = 0.0;
      double z_i = 0.0;

            for(int ma1 = MAX(0, (2 * ma - j - j2 + j1) / 2);
                ma1 <= MIN(j1, (2 * ma - j + j2 + j1) / 2); ma1++) {
              double sumb1_r = 0.0;
              double sumb1_i = 0.0;

              const int ma2 = (2 * ma - j - (2 * ma1 - j1) + j2) / 2;

              for(int mb1  = MAX( 0, (2 * mb - j - j2 + j1) / 2);
                mb1 <= MIN(j1, (2 * mb - j + j2 + j1) / 2); mb1++) {

                const int mb2 = (2 * mb - j - (2 * mb1 - j1) + j2) / 2;
    const double cga = cgarray(j1,j2,j,mb1,mb2);
    const double uat1_r = uarraytot_r(j1,ma1,mb1);
    const double uat1_i = uarraytot_i(j1,ma1,mb1);
    const double uat2_r = uarraytot_r(j2,ma2,mb2);
    const double uat2_i = uarraytot_i(j2,ma2,mb2);
    sumb1_r += cga * (uat1_r * uat2_r - uat1_i * uat2_i);
    sumb1_i += cga * (uat1_r * uat2_i + uat1_i * uat2_r);
                /*sumb1_r += cgarray(j1,j2,j,mb1,mb2) *
                  (uarraytot_r(j1,ma1,mb1) * uarraytot_r(j2,ma2,mb2) -
                   uarraytot_i(j1,ma1,mb1) * uarraytot_i(j2,ma2,mb2));
                sumb1_i += cgarray(j1,j2,j,mb1,mb2) *
                  (uarraytot_r(j1,ma1,mb1) * uarraytot_i(j2,ma2,mb2) +
                   uarraytot_i(j1,ma1,mb1) * uarraytot_r(j2,ma2,mb2));*/
              } // end loop over mb1

        const double cga = cgarray(j1,j2,j,ma1,ma2);
              z_r += sumb1_r * cga;//rray(j1,j2,j,ma1,ma2);
              z_i += sumb1_i * cga;//rray(j1,j2,j,ma1,ma2);
            } // end loop over ma1
      zarray_r(j1,j2,j,mb,ma) = z_r;
            zarray_i(j1,j2,j,mb,ma) = z_i;
          }); // end loop over ma, mb
    //  }
    //}
  });
      //} // end loop over j
    //} // end loop over j1, j2

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &endtime);
  timers[1] += (endtime.tv_sec - starttime.tv_sec + 1.0 *
                (endtime.tv_nsec - starttime.tv_nsec) / 1000000000);
#endif
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

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &starttime);
#endif

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,idxj_full_max),
      [&] (const int& idx) {
    const int j1 = idxj_full(idx).j1;
    const int j2 = idxj_full(idx).j2;
    const int j =  idxj_full(idx).j;

    const int bound = (j+2)/2;
    double b_j1_j2_j = 0.0;
    double b_j1_j2_j_temp = 0.0;
    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,(j+1)*bound),
        [&] (const int mbma, double& sum) {
        //for(int mb = 0; 2*mb <= j; mb++)
          //for(int ma = 0; ma <= j; ma++) {
        const int ma = mbma%(j+1);
        const int mb = mbma/(j+1);
        if (2*mb == j) return;
        sum +=
          uarraytot_r(j,ma,mb) * zarray_r(j1,j2,j,mb,ma) +
          uarraytot_i(j,ma,mb) * zarray_i(j1,j2,j,mb,ma);
      },b_j1_j2_j_temp); // end loop over ma, mb
      b_j1_j2_j += b_j1_j2_j_temp;

    // For j even, special treatment for middle column

    if (j%2 == 0) {
      const int mb = j/2;
      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, mb),
          [&] (const int ma, double& sum) {
      //for(int ma = 0; ma < mb; ma++) {
        sum +=
          uarraytot_r(j,ma,mb) * zarray_r(j1,j2,j,mb,ma) +
          uarraytot_i(j,ma,mb) * zarray_i(j1,j2,j,mb,ma);
      },b_j1_j2_j_temp); // end loop over ma
      b_j1_j2_j += b_j1_j2_j_temp;

      const int ma = mb;
      b_j1_j2_j +=
        (uarraytot_r(j,ma,mb) * zarray_r(j1,j2,j,mb,ma) +
         uarraytot_i(j,ma,mb) * zarray_i(j1,j2,j,mb,ma))*0.5;
    }

    Kokkos::single(Kokkos::PerThread(team), [&] () {
      b_j1_j2_j *= 2.0;
      if (bzero_flag)
        b_j1_j2_j -= bzero[j];

      barray(j1,j2,j) = b_j1_j2_j;
    });
  });
      //} // end loop over j
    //} // end loop over j1, j2

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &endtime);
  timers[2] += (endtime.tv_sec - starttime.tv_sec + 1.0 *
                (endtime.tv_nsec - starttime.tv_nsec) / 1000000000);
#endif

}

/* ----------------------------------------------------------------------
   copy Bi derivatives into a vector
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::copy_bi2bvec(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team)
{
 /* int ncount, j1, j2, j;

  ncount = 0;

  for(j1 = 0; j1 <= twojmax; j1++) {
      for(j2 = 0; j2 <= j1; j2++)
        for(j = abs(j1 - j2);
            j <= MIN(twojmax, j1 + j2); j += 2)
            if (j >= j1) {*/
  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,idxj_max),
          [&] (const int& JJ) {
  //for(int JJ = 0; JJ < idxj_max; JJ++) {
    const int j1 = idxj[JJ].j1;
    const int j2 = idxj[JJ].j2;
    const int j =  idxj[JJ].j;
    bvec(JJ) = barray(j1,j2,j);
    //ncount++;
  });
}

/* ----------------------------------------------------------------------
   calculate derivative of Ui w.r.t. atom j
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_duidrj(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team,
                         double* rij, double wj, double rcut)
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

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &starttime);
#endif

  compute_duarray(team, x, y, z, z0, r, dz0dr, wj, rcut);

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &endtime);
  timers[3] += (endtime.tv_sec - starttime.tv_sec + 1.0 *
                (endtime.tv_nsec - starttime.tv_nsec) / 1000000000);
#endif

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

  double* dudr_r, *dudr_i;
  double jjjmambzarray_r;
  double jjjmambzarray_i;

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &starttime);
#endif
  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,idxj_max),
          [&] (const int& JJ) {
  //for(int JJ = 0; JJ < idxj_max; JJ++) {
    const int j1 = idxj[JJ].j1;
    const int j2 = idxj[JJ].j2;
    const int j = idxj[JJ].j;

//    dbdr = &dbarray(j1,j2,j,0);
//    dbdr[0] = 0.0;
//    dbdr[1] = 0.0;
//    dbdr[2] = 0.0;

    t_scalar3<double> dbdr,sumzdu_r;
    // Sum terms Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)

    // use zarray j1/j2 symmetry (optional)

    int j_,j1_,j2_;
    if (j1 >= j2) {
      //jjjzarray_r = &zarray_r(j1,j2,j);
      //jjjzarray_i = &zarray_i(j1,j2,j);
      j1_ = j1;
      j2_ = j2;
      j_ = j;
    } else {
      j1_ = j2;
      j2_ = j1;
      j_ = j;
      //jjjzarray_r = &zarray_r(j2,j1,j);
      //jjjzarray_i = &zarray_i(j2,j1,j);
    }

    for(int mb = 0; 2*mb < j; mb++)
      for(int ma = 0; ma <= j; ma++) {

        dudr_r = &duarray_r(j,mb,ma,0);
        dudr_i = &duarray_i(j,mb,ma,0);
        jjjmambzarray_r = zarray_r(j1_,j2_,j_,mb,ma);
        jjjmambzarray_i = zarray_i(j1_,j2_,j_,mb,ma);
        sumzdu_r.x += (dudr_r[0] * jjjmambzarray_r + dudr_i[0] * jjjmambzarray_i);
        sumzdu_r.y += (dudr_r[1] * jjjmambzarray_r + dudr_i[1] * jjjmambzarray_i);
        sumzdu_r.z += (dudr_r[2] * jjjmambzarray_r + dudr_i[2] * jjjmambzarray_i);

      } //end loop over ma mb

    // For j even, handle middle column

    if (j%2 == 0) {
      int mb = j/2;
      for(int ma = 0; ma <= mb; ma++) {
        dudr_r = &duarray_r(j,mb,ma,0);
        dudr_i = &duarray_i(j,mb,ma,0);
        const double factor = ma==mb?0.5:1.0;
        jjjmambzarray_r = zarray_r(j1_,j2_,j_,mb,ma) * factor;
        jjjmambzarray_i = zarray_i(j1_,j2_,j_,mb,ma) * factor;
        sumzdu_r.x += (dudr_r[0] * jjjmambzarray_r + dudr_i[0] * jjjmambzarray_i);
        sumzdu_r.y += (dudr_r[1] * jjjmambzarray_r + dudr_i[1] * jjjmambzarray_i);
        sumzdu_r.z += (dudr_r[2] * jjjmambzarray_r + dudr_i[2] * jjjmambzarray_i);
      }
    } // end if jeven

      dbdr += 2.0*sumzdu_r;

    // Sum over Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)

    double j1fac = (j+1)/(j1+1.0);

    sumzdu_r.x = 0.0; sumzdu_r.y = 0.0; sumzdu_r.z = 0.0;

    // use zarray j1/j2 symmetry (optional)

    if (j >= j2) {
      j1_ = j;
      j2_ = j2;
      j_ = j1;

      //jjjzarray_r = zarray_r(j,j2,j1);
      //jjjzarray_i = zarray_i(j,j2,j1);
    } else {
      j1_ = j2;
      j2_ = j;
      j_ = j1;
      //jjjzarray_r = zarray_r(j2,j,j1);
      //jjjzarray_i = zarray_i(j2,j,j1);
    }

    for(int mb1 = 0; 2*mb1 < j1; mb1++)
      for(int ma1 = 0; ma1 <= j1; ma1++) {

        dudr_r = &duarray_r(j1,mb1,ma1,0);
        dudr_i = &duarray_i(j1,mb1,ma1,0);
        jjjmambzarray_r = zarray_r(j1_,j2_,j_,mb1,ma1);
        jjjmambzarray_i = zarray_i(j1_,j2_,j_,mb1,ma1);
        sumzdu_r.x += (dudr_r[0] * jjjmambzarray_r + dudr_i[0] * jjjmambzarray_i);
        sumzdu_r.y += (dudr_r[1] * jjjmambzarray_r + dudr_i[1] * jjjmambzarray_i);
        sumzdu_r.z += (dudr_r[2] * jjjmambzarray_r + dudr_i[2] * jjjmambzarray_i);
      } //end loop over ma1 mb1

    // For j1 even, handle middle column

    if (j1%2 == 0) {
      const int mb1 = j1/2;
      for(int ma1 = 0; ma1 <= mb1; ma1++) {
        dudr_r = &duarray_r(j1,mb1,ma1,0);
        dudr_i = &duarray_i(j1,mb1,ma1,0);
        const double factor = ma1==mb1?0.5:1.0;
        jjjmambzarray_r = zarray_r(j1_,j2_,j_,mb1,ma1) * factor;
        jjjmambzarray_i = zarray_i(j1_,j2_,j_,mb1,ma1) * factor;
        sumzdu_r.x += (dudr_r[0] * jjjmambzarray_r + dudr_i[0] * jjjmambzarray_i);
        sumzdu_r.y += (dudr_r[1] * jjjmambzarray_r + dudr_i[1] * jjjmambzarray_i);
        sumzdu_r.z += (dudr_r[2] * jjjmambzarray_r + dudr_i[2] * jjjmambzarray_i);
      }
    } // end if j1even

      dbdr += 2.0*sumzdu_r*j1fac;

    // Sum over Conj(dudr(j2,ma2,mb2))*z(j1,j,j2,ma2,mb2)

    double j2fac = (j+1)/(j2+1.0);

    sumzdu_r.x = 0.0; sumzdu_r.y = 0.0; sumzdu_r.z = 0.0;

    // use zarray j1/j2 symmetry (optional)

    if (j1 >= j) {
      j1_ = j1;
      j2_ = j;
      j_ = j2;
      //jjjzarray_r = zarray_r(j1,j,j2);
      //jjjzarray_i = zarray_i(j1,j,j2);
    } else {
      j1_ = j;
      j2_ = j1;
      j_ = j2;
      //jjjzarray_r = zarray_r(j,j1,j2);
      //jjjzarray_i = zarray_i(j,j1,j2);
    }

    for(int mb2 = 0; 2*mb2 < j2; mb2++)
      for(int ma2 = 0; ma2 <= j2; ma2++) {

        dudr_r = &duarray_r(j2,mb2,ma2,0);
        dudr_i = &duarray_i(j2,mb2,ma2,0);
        jjjmambzarray_r = zarray_r(j1_,j2_,j_,mb2,ma2);
        jjjmambzarray_i = zarray_i(j1_,j2_,j_,mb2,ma2);
        sumzdu_r.x += (dudr_r[0] * jjjmambzarray_r + dudr_i[0] * jjjmambzarray_i);
        sumzdu_r.y += (dudr_r[1] * jjjmambzarray_r + dudr_i[1] * jjjmambzarray_i);
        sumzdu_r.z += (dudr_r[2] * jjjmambzarray_r + dudr_i[2] * jjjmambzarray_i);
      } //end loop over ma2 mb2

    // For j2 even, handle middle column

    if (j2%2 == 0) {
      const int mb2 = j2/2;
      for(int ma2 = 0; ma2 <= mb2; ma2++) {
        dudr_r = &duarray_r(j2,mb2,ma2,0);
        dudr_i = &duarray_i(j2,mb2,ma2,0);
        const double factor = ma2==mb2?0.5:1.0;
        jjjmambzarray_r = zarray_r(j1_,j2_,j_,mb2,ma2) * factor;
        jjjmambzarray_i = zarray_i(j1_,j2_,j_,mb2,ma2) * factor;
        sumzdu_r.x += (dudr_r[0] * jjjmambzarray_r + dudr_i[0] * jjjmambzarray_i);
        sumzdu_r.y += (dudr_r[1] * jjjmambzarray_r + dudr_i[1] * jjjmambzarray_i);
        sumzdu_r.z += (dudr_r[2] * jjjmambzarray_r + dudr_i[2] * jjjmambzarray_i);
      }
    } // end if j2even

    dbdr += 2.0*sumzdu_r*j2fac;
    dbarray(j1,j2,j,0) = dbdr.x;
    dbarray(j1,j2,j,1) = dbdr.y;
    dbarray(j1,j2,j,2) = dbdr.z;
  }); //end loop over j1 j2 j

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &endtime);
  timers[4] += (endtime.tv_sec - starttime.tv_sec + 1.0 *
                (endtime.tv_nsec - starttime.tv_nsec) / 1000000000);
#endif

}

/* ----------------------------------------------------------------------
   copy Bi derivatives into a vector
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::copy_dbi2dbvec(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team)
{
 /* int ncount, j1, j2, j;

  ncount = 0;

  for(j1 = 0; j1 <= twojmax; j1++) {
      for(j2 = 0; j2 <= j1; j2++)
        for(j = abs(j1 - j2);
            j <= MIN(twojmax, j1 + j2); j += 2)
          if (j >= j1) {*/
  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,idxj_max),
          [&] (const int& JJ) {
  //for(int JJ = 0; JJ < idxj_max; JJ++) {
    const int j1 = idxj[JJ].j1;
    const int j2 = idxj[JJ].j2;
    const int j = idxj[JJ].j;
    dbvec(JJ,0) = dbarray(j1,j2,j,0);
    dbvec(JJ,1) = dbarray(j1,j2,j,1);
    dbvec(JJ,2) = dbarray(j1,j2,j,2);
    //ncount++;
  });
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::zero_uarraytot(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team)
{
  {
    double* const ptr = uarraytot_r.data();
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,uarraytot_r.span()),
        [&] (const int& i) {
      ptr[i] = 0.0;
    });
  }
  {
    double* const ptr = uarraytot_i.data();
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,uarraytot_r.span()),
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
  //for (int j = 0; j <= twojmax; j++)
  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,twojmax+1),
    [&] (const int& j) {
    for (int ma = 0; ma <= j; ma++) {
      uarraytot_r(j,ma,ma) = wself_in;
      uarraytot_i(j,ma,ma) = 0.0;
    }
  });
}

/* ----------------------------------------------------------------------
   add Wigner U-functions for one neighbor to the total
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::add_uarraytot(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, double r, double wj, double rcut)
{
  const double sfac = compute_sfac(r, rcut) * wj;

/*
  for (int j = 0; j <= twojmax; j++)
    for (int ma = 0; ma <= j; ma++)
      for (int mb = 0; mb <= j; mb++) {
        uarraytot_r_a(j,ma,mb) +=
          sfac * uarray_r(j,ma,mb);
        uarraytot_i_a(j,ma,mb) +=
          sfac * uarray_i(j,ma,mb);
      }*/
  const double* const ptr_r = uarray_r.data();
  const double* const ptr_i = uarray_i.data();
  double* const ptrtot_r = uarraytot_r.data();
  double* const ptrtot_i = uarraytot_i.data();
  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,uarraytot_r.span()),
      [&] (const int& i) {
    Kokkos::atomic_add(ptrtot_r+i, sfac * ptr_r[i]);
    Kokkos::atomic_add(ptrtot_i+i, sfac * ptr_i[i]);
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

  uarray_r(0,0,0) = 1.0;
  uarray_i(0,0,0) = 0.0;

  for (int j = 1; j <= twojmax; j++) {

    // fill in left side of matrix layer from previous layer

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,(j+2)/2),
        [&] (const int& mb) {
      //const int mb = 2*mb_2;
    //for (int mb = 0; 2*mb <= j; mb++) {
      uarray_r(j,0,mb) = 0.0;
      uarray_i(j,0,mb) = 0.0;

      for (int ma = 0; ma < j; ma++) {
        rootpq = rootpqarray(j - ma,j - mb);
        uarray_r(j,ma,mb) +=
          rootpq *
          (a_r * uarray_r(j - 1,ma,mb) +
           a_i * uarray_i(j - 1,ma,mb));
        uarray_i(j,ma,mb) +=
          rootpq *
          (a_r * uarray_i(j - 1,ma,mb) -
           a_i * uarray_r(j - 1,ma,mb));

        rootpq = rootpqarray(ma + 1,j - mb);
        uarray_r(j,ma + 1,mb) =
          -rootpq *
          (b_r * uarray_r(j - 1,ma,mb) +
           b_i * uarray_i(j - 1,ma,mb));
        uarray_i(j,ma + 1,mb) =
          -rootpq *
          (b_r * uarray_i(j - 1,ma,mb) -
           b_i * uarray_r(j - 1,ma,mb));
      }
    });

    // copy left side to right side with inversion symmetry VMK 4.4(2)
    // u[ma-j,mb-j] = (-1)^(ma-mb)*Conj([u[ma,mb))

    //int mbpar = -1;
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,(j+2)/2),
        [&] (const int& mb) {
//    for (int mb = 0; 2*mb <= j; mb++) {
      int mbpar = (mb)%2==0?1:-1;
      int mapar = -mbpar;
      for (int ma = 0; ma <= j; ma++) {
        mapar = -mapar;
        if (mapar == 1) {
          uarray_r(j,j-ma,j-mb) = uarray_r(j,ma,mb);
          uarray_i(j,j-ma,j-mb) = -uarray_i(j,ma,mb);
        } else {
          uarray_r(j,j-ma,j-mb) = -uarray_r(j,ma,mb);
          uarray_i(j,j-ma,j-mb) = uarray_i(j,ma,mb);
        }
        //OK
        //printf("%lf %lf %lf %lf %lf %lf %lf SNAP-COMPARE: UARRAY\n",x,y,z,z0,r,uarray_r(j,ma,mb),uarray_i(j,ma,mb));
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

  uarray_r(0,0,0) = 1.0;
  duarray_r(0,0,0,0) = 0.0;
  duarray_r(0,0,0,1) = 0.0;
  duarray_r(0,0,0,2) = 0.0;
  uarray_i(0,0,0) = 0.0;
  duarray_i(0,0,0,0) = 0.0;
  duarray_i(0,0,0,1) = 0.0;
  duarray_i(0,0,0,2) = 0.0;

  for (int j = 1; j <= twojmax; j++) {
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,(j+2)/2),
        [&] (const int& mb) {

    //for (int mb = 0; 2*mb <= j; mb++) {
      uarray_r(j,0,mb) = 0.0;
      duarray_r(j,mb,0,0) = 0.0;
      duarray_r(j,mb,0,1) = 0.0;
      duarray_r(j,mb,0,2) = 0.0;
      uarray_i(j,0,mb) = 0.0;
      duarray_i(j,mb,0,0) = 0.0;
      duarray_i(j,mb,0,1) = 0.0;
      duarray_i(j,mb,0,2) = 0.0;

      for (int ma = 0; ma < j; ma++) {
        rootpq = rootpqarray(j - ma,j - mb);
        uarray_r(j,ma,mb) += rootpq *
                               (a_r *  uarray_r(j - 1,ma,mb) +
                                a_i *  uarray_i(j - 1,ma,mb));
        uarray_i(j,ma,mb) += rootpq *
                               (a_r *  uarray_i(j - 1,ma,mb) -
                                a_i *  uarray_r(j - 1,ma,mb));

        for (int k = 0; k < 3; k++) {
          duarray_r(j,mb,ma,k) +=
            rootpq * (da_r[k] * uarray_r(j - 1,ma,mb) +
                      da_i[k] * uarray_i(j - 1,ma,mb) +
                      a_r * duarray_r(j - 1,mb,ma,k) +
                      a_i * duarray_i(j - 1,mb,ma,k));
          duarray_i(j,mb,ma,k) +=
            rootpq * (da_r[k] * uarray_i(j - 1,ma,mb) -
                      da_i[k] * uarray_r(j - 1,ma,mb) +
                      a_r * duarray_i(j - 1,mb,ma,k) -
                      a_i * duarray_r(j - 1,mb,ma,k));
        }

        rootpq = rootpqarray(ma + 1,j - mb);
        uarray_r(j,ma + 1,mb) =
          -rootpq * (b_r *  uarray_r(j - 1,ma,mb) +
                     b_i *  uarray_i(j - 1,ma,mb));
        uarray_i(j,ma + 1,mb) =
          -rootpq * (b_r *  uarray_i(j - 1,ma,mb) -
                     b_i *  uarray_r(j - 1,ma,mb));

        for (int k = 0; k < 3; k++) {
          duarray_r(j,mb,ma + 1,k) =
            -rootpq * (db_r[k] * uarray_r(j - 1,ma,mb) +
                       db_i[k] * uarray_i(j - 1,ma,mb) +
                       b_r * duarray_r(j - 1,mb,ma,k) +
                       b_i * duarray_i(j - 1,mb,ma,k));
          duarray_i(j,mb,ma + 1,k) =
            -rootpq * (db_r[k] * uarray_i(j - 1,ma,mb) -
                       db_i[k] * uarray_r(j - 1,ma,mb) +
                       b_r * duarray_i(j - 1,mb,ma,k) -
                       b_i * duarray_r(j - 1,mb,ma,k));
        }
      }
    });

    //int mbpar = -1;
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,(j+2)/2),
        [&] (const int& mb) {
//    for (int mb = 0; 2*mb <= j; mb++) {
      int mbpar = (mb)%2==0?1:-1;
      int mapar = -mbpar;
      for (int ma = 0; ma <= j; ma++) {
        mapar = -mapar;
        if (mapar == 1) {
          uarray_r(j,j-ma,j-mb) = uarray_r(j,ma,mb);
          uarray_i(j,j-ma,j-mb) = -uarray_i(j,ma,mb);
          for (int k = 0; k < 3; k++) {
            duarray_r(j,j-mb,j-ma,k) = duarray_r(j,mb,ma,k);
            duarray_i(j,j-mb,j-ma,k) = -duarray_i(j,mb,ma,k);
          }
        } else {
          uarray_r(j,j-ma,j-mb) = -uarray_r(j,ma,mb);
          uarray_i(j,j-ma,j-mb) = uarray_i(j,ma,mb);
          for (int k = 0; k < 3; k++) {
            duarray_r(j,j-mb,j-ma,k) = -duarray_r(j,mb,ma,k);
            duarray_i(j,j-mb,j-ma,k) = duarray_i(j,mb,ma,k);
          }
        }
      }
    });
  }

  double sfac = compute_sfac(r, rcut);
  double dsfac = compute_dsfac(r, rcut);

  sfac *= wj;
  dsfac *= wj;

  for (int j = 0; j <= twojmax; j++)
    for (int mb = 0; mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++) {
        duarray_r(j,mb,ma,0) = dsfac * uarray_r(j,ma,mb) * ux +
                                  sfac * duarray_r(j,mb,ma,0);
        duarray_i(j,mb,ma,0) = dsfac * uarray_i(j,ma,mb) * ux +
                                  sfac * duarray_i(j,mb,ma,0);
        duarray_r(j,mb,ma,1) = dsfac * uarray_r(j,ma,mb) * uy +
                                  sfac * duarray_r(j,mb,ma,1);
        duarray_i(j,mb,ma,1) = dsfac * uarray_i(j,ma,mb) * uy +
                                  sfac * duarray_i(j,mb,ma,1);
        duarray_r(j,mb,ma,2) = dsfac * uarray_r(j,ma,mb) * uz +
                                  sfac * duarray_r(j,mb,ma,2);
        duarray_i(j,mb,ma,2) = dsfac * uarray_i(j,ma,mb) * uz +
                                  sfac * duarray_i(j,mb,ma,2);
      }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::create_team_scratch_arrays(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team)
{
  int jdim = twojmax + 1;
  uarraytot_r_a = uarraytot_r = t_sna_3d(team.team_scratch(1),jdim,jdim,jdim);
  uarraytot_i_a = uarraytot_i = t_sna_3d(team.team_scratch(1),jdim,jdim,jdim);
  zarray_r = t_sna_5d(team.team_scratch(1),jdim,jdim,jdim,jdim,jdim);
  zarray_i = t_sna_5d(team.team_scratch(1),jdim,jdim,jdim,jdim,jdim);
  bvec = Kokkos::View<double*, Kokkos::LayoutRight, DeviceType>(team.team_scratch(1),ncoeff);
  barray = t_sna_3d(team.team_scratch(1),jdim,jdim,jdim);

  rij = t_sna_2d(team.team_scratch(1),nmax,3);
  rcutij = t_sna_1d(team.team_scratch(1),nmax);
  wj = t_sna_1d(team.team_scratch(1),nmax);
  inside = t_sna_1i(team.team_scratch(1),nmax);
}


template<class DeviceType>
inline
T_INT SNAKokkos<DeviceType>::size_team_scratch_arrays() {
  T_INT size = 0;
  int jdim = twojmax + 1;

  size += t_sna_3d::shmem_size(jdim,jdim,jdim); // uarraytot_r_a
  size += t_sna_3d::shmem_size(jdim,jdim,jdim); // uarraytot_i_a
  size += t_sna_5d::shmem_size(jdim,jdim,jdim,jdim,jdim); // zarray_r
  size += t_sna_5d::shmem_size(jdim,jdim,jdim,jdim,jdim); // zarray_i
  size += Kokkos::View<double*, Kokkos::LayoutRight, DeviceType>::shmem_size(ncoeff); // bvec
  size += t_sna_3d::shmem_size(jdim,jdim,jdim); // barray

  size += t_sna_2d::shmem_size(nmax,3); // rij
  size += t_sna_1d::shmem_size(nmax); // rcutij
  size += t_sna_1d::shmem_size(nmax); // wj
  size += t_sna_1i::shmem_size(nmax); // inside

  return size;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::create_thread_scratch_arrays(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team)
{
  int jdim = twojmax + 1;

  dbvec = Kokkos::View<double*[3], Kokkos::LayoutRight, DeviceType>(team.thread_scratch(1),ncoeff);
  dbarray = t_sna_4d(team.thread_scratch(1),jdim,jdim,jdim);

  uarray_r = t_sna_3d(team.thread_scratch(1),jdim,jdim,jdim);
  uarray_i = t_sna_3d(team.thread_scratch(1),jdim,jdim,jdim);
  duarray_r = t_sna_4d(team.thread_scratch(1),jdim,jdim,jdim);
  duarray_i = t_sna_4d(team.thread_scratch(1),jdim,jdim,jdim);
}

template<class DeviceType>
inline
T_INT SNAKokkos<DeviceType>::size_thread_scratch_arrays() {
  T_INT size = 0;
  int jdim = twojmax + 1;

  size += Kokkos::View<double*[3], Kokkos::LayoutRight, DeviceType>::shmem_size(ncoeff); // dbvec
  size += t_sna_4d::shmem_size(jdim,jdim,jdim); // dbarray

  size += t_sna_3d::shmem_size(jdim,jdim,jdim); // uarray_r
  size += t_sna_3d::shmem_size(jdim,jdim,jdim); // uarray_i
  size += t_sna_4d::shmem_size(jdim,jdim,jdim); // duarray_r
  size += t_sna_4d::shmem_size(jdim,jdim,jdim); // duarray_i
  return size;
}

/* ----------------------------------------------------------------------
   factorial n
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double SNAKokkos<DeviceType>::factorial(int n)
{
  double result = 1.0;
  for(int i=1; i<=n; i++)
    result *= 1.0*i;
  return result;
}

/* ----------------------------------------------------------------------
   the function delta given by VMK Eq. 8.2(1)
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
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
  double sum,dcg,sfaccg;
  int m, aa2, bb2, cc2;
  int ifac;
  auto h_cgarray = Kokkos::create_mirror_view(cgarray);

  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= twojmax; j2++)
      for (int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2)
        for (int m1 = 0; m1 <= j1; m1 += 1) {
          aa2 = 2 * m1 - j1;

          for (int m2 = 0; m2 <= j2; m2 += 1) {

            // -c <= cc <= c

            bb2 = 2 * m2 - j2;
            m = (aa2 + bb2 + j) / 2;

            if(m < 0 || m > j) continue;

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

            h_cgarray(j1,j2,j,m1,m2) = sum * dcg * sfaccg;
            //printf("SNAP-COMPARE: CG: %i %i %i %i %i %e\n",j1,j2,j,m1,m2,cgarray(j1,j2,j,m1,m2));
          }
        }
  Kokkos::deep_copy(cgarray,h_cgarray);
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
      for (int j = abs(j1 - j2);
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
  int jdim = twojmax + 1;
  double bytes;
  bytes = jdim * jdim * jdim * jdim * jdim * sizeof(double);
  bytes += 2 * jdim * jdim * jdim * sizeof(std::complex<double>);
  bytes += 2 * jdim * jdim * jdim * sizeof(double);
  bytes += jdim * jdim * jdim * 3 * sizeof(std::complex<double>);
  bytes += jdim * jdim * jdim * 3 * sizeof(double);
  bytes += ncoeff * sizeof(double);
  bytes += jdim * jdim * jdim * jdim * jdim * sizeof(std::complex<double>);
  return bytes;
}

} // namespace LAMMPS_NS
