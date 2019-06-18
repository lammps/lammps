/* -*- c++ -*- -------------------------------------------------------------
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

#ifndef LMP_SNA_KOKKOS_H
#define LMP_SNA_KOKKOS_H

#include <complex>
#include <ctime>
#include <Kokkos_Core.hpp>
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct SNAKK_LOOPINDICES {
  int j1, j2, j;
};

template<class DeviceType>
class SNAKokkos {

public:
  typedef Kokkos::View<int*, DeviceType> t_sna_1i;
  typedef Kokkos::View<double*, DeviceType> t_sna_1d;
  typedef Kokkos::View<double**, Kokkos::LayoutRight, DeviceType> t_sna_2d;
  typedef Kokkos::View<double***, Kokkos::LayoutRight, DeviceType> t_sna_3d;
  typedef Kokkos::View<double***, Kokkos::LayoutRight, DeviceType, Kokkos::MemoryTraits<Kokkos::Atomic> > t_sna_3d_atomic;
  typedef Kokkos::View<double***[3], Kokkos::LayoutRight, DeviceType> t_sna_4d;
  typedef Kokkos::View<double**[3], Kokkos::LayoutRight, DeviceType> t_sna_3d3;
  typedef Kokkos::View<double*****, Kokkos::LayoutRight, DeviceType> t_sna_5d;

inline
  SNAKokkos() {};
  KOKKOS_INLINE_FUNCTION
  SNAKokkos(const SNAKokkos<DeviceType>& sna, const typename Kokkos::TeamPolicy<DeviceType>::member_type& team);

inline
  SNAKokkos(double, int, double, int, int);

  KOKKOS_INLINE_FUNCTION
  ~SNAKokkos();

inline
  void build_indexlist(); // SNAKokkos()

inline
  void init();            //

inline
  T_INT size_team_scratch_arrays();

inline
  T_INT size_thread_scratch_arrays();

  double memory_usage();

  int ncoeff;

  // functions for bispectrum coefficients

  KOKKOS_INLINE_FUNCTION
  void compute_ui(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int); // ForceSNAP
  KOKKOS_INLINE_FUNCTION
  void compute_zi(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team);    // ForceSNAP
  KOKKOS_INLINE_FUNCTION
  void compute_bi(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team);    // ForceSNAP
  KOKKOS_INLINE_FUNCTION
  void copy_bi2bvec(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team); //ForceSNAP

  // functions for derivatives

  KOKKOS_INLINE_FUNCTION
  void compute_duidrj(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, double*, double, double); //ForceSNAP
  KOKKOS_INLINE_FUNCTION
  void compute_dbidrj(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team); //ForceSNAP
  KOKKOS_INLINE_FUNCTION
  void copy_dbi2dbvec(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team); //ForceSNAP
  KOKKOS_INLINE_FUNCTION
  double compute_sfac(double, double); // add_uarraytot, compute_duarray
  KOKKOS_INLINE_FUNCTION
  double compute_dsfac(double, double); // compute_duarray

#ifdef TIMING_INFO
  double* timers;
  timespec starttime, endtime;
  int print;
  int counter;
#endif

  //per sna class instance for OMP use


  // Per InFlight Particle
  t_sna_2d rij;
  t_sna_1i inside;
  t_sna_1d wj;
  t_sna_1d rcutij;
  int nmax;

  void grow_rij(int);

  int twojmax, diagonalstyle;
  // Per InFlight Particle
  t_sna_3d barray;
  t_sna_3d uarraytot_r, uarraytot_i;
  t_sna_3d_atomic uarraytot_r_a, uarraytot_i_a;
  t_sna_5d zarray_r, zarray_i;

  // Per InFlight Interaction
  t_sna_3d uarray_r, uarray_i;

  Kokkos::View<double*, Kokkos::LayoutRight, DeviceType> bvec;

  // derivatives of data
  Kokkos::View<double*[3], Kokkos::LayoutRight, DeviceType> dbvec;
  t_sna_4d duarray_r, duarray_i;
  t_sna_4d dbarray;

private:
  double rmin0, rfac0;

  //use indexlist instead of loops, constructor generates these
  // Same accross all SNAKokkos
  Kokkos::View<SNAKK_LOOPINDICES*, DeviceType> idxj,idxj_full;
  int idxj_max,idxj_full_max;
  // data for bispectrum coefficients

  // Same accross all SNAKokkos
  t_sna_5d cgarray;
  t_sna_2d rootpqarray;


  static const int nmaxfactorial = 167;
  KOKKOS_INLINE_FUNCTION
  double factorial(int);

  KOKKOS_INLINE_FUNCTION
  void create_team_scratch_arrays(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team); // SNAKokkos()
  KOKKOS_INLINE_FUNCTION
  void create_thread_scratch_arrays(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team); // SNAKokkos()

inline
  void init_clebsch_gordan(); // init()

inline
  void init_rootpqarray();    // init()
  KOKKOS_INLINE_FUNCTION
  void zero_uarraytot(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team);      // compute_ui
  KOKKOS_INLINE_FUNCTION
  void addself_uarraytot(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, double); // compute_ui
  KOKKOS_INLINE_FUNCTION
  void add_uarraytot(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, double, double, double); // compute_ui

  KOKKOS_INLINE_FUNCTION
  void compute_uarray(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team,
                      double, double, double,
                      double, double); // compute_ui
  KOKKOS_INLINE_FUNCTION
  double deltacg(int, int, int);  // init_clebsch_gordan

inline
  int compute_ncoeff();           // SNAKokkos()
  KOKKOS_INLINE_FUNCTION
  void compute_duarray(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team,
                       double, double, double, // compute_duidrj
                       double, double, double, double, double);

  // Sets the style for the switching function
  // 0 = none
  // 1 = cosine
  int switch_flag;

  // Self-weight
  double wself;

  int bzero_flag; // 1 if bzero subtracted from barray
  Kokkos::View<double*, Kokkos::LayoutRight, DeviceType> bzero; // array of B values for isolated atoms
};

}

#include "sna_kokkos_impl.h"
#endif

/* ERROR/WARNING messages:

E: Invalid argument to factorial %d

N must be >= 0 and <= 167, otherwise the factorial result is too
large.

*/
