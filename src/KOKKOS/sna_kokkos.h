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

typedef double SNAreal;

//typedef struct { SNAreal re, im; } SNAcomplex;
template <typename real>
struct alignas(2*sizeof(real)) SNAComplex
{
  real re,im;

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex() = default;

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex(real re)
   : re(re), im(static_cast<real>(0.)) { ; }

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex(real re, real im)
   : re(re), im(im) { ; }

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex(const SNAComplex& other)
   : re(other.re), im(other.im) { ; }

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex& operator=(const SNAComplex& other) {
    re = other.re; im = other.im;
    return *this;
  }

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex(SNAComplex&& other)
   : re(other.re), im(other.im) { ; }

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex& operator=(SNAComplex&& other) {
    re = other.re; im = other.im;
    return *this;
  }

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex operator+(SNAComplex const& other) {
    return SNAComplex(re + other.re, im + other.im);
  }

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex& operator+=(SNAComplex const& other) {
    re += other.re; im += other.im;
    return *this;
  }

};

template <typename real>
KOKKOS_FORCEINLINE_FUNCTION SNAComplex<real> operator*(const real& r, const SNAComplex<real>& self) {
  return SNAComplex<real>(r*self.re, r*self.im);
}

typedef SNAComplex<SNAreal> SNAcomplex;

//struct SNAKK_ZINDICES {
//  int j1, j2, j, ma1min, ma2max, mb1min, mb2max, na, nb, jju;
//};
//
//struct SNAKK_BINDICES {
//  int j1, j2, j;
//};

template<class DeviceType>
class SNAKokkos {

public:
  typedef Kokkos::View<int*, DeviceType> t_sna_1i;
  typedef Kokkos::View<double*, DeviceType> t_sna_1d;
  typedef Kokkos::View<double*, DeviceType, Kokkos::MemoryTraits<Kokkos::Atomic> > t_sna_1d_atomic;
  typedef Kokkos::View<int**, DeviceType> t_sna_2i;
  typedef Kokkos::View<double**, DeviceType> t_sna_2d;
  typedef Kokkos::View<double**, Kokkos::LayoutLeft, DeviceType> t_sna_2d_ll;
  typedef Kokkos::View<double***, DeviceType> t_sna_3d;
  typedef Kokkos::View<double***[3], DeviceType> t_sna_4d;
  typedef Kokkos::View<double**[3], DeviceType> t_sna_3d3;
  typedef Kokkos::View<double*****, DeviceType> t_sna_5d;

  typedef Kokkos::View<SNAcomplex*, DeviceType> t_sna_1c;
  typedef Kokkos::View<SNAcomplex*, DeviceType, Kokkos::MemoryTraits<Kokkos::Atomic> > t_sna_1c_atomic;
  typedef Kokkos::View<SNAcomplex**, DeviceType> t_sna_2c;
  typedef Kokkos::View<SNAcomplex**, Kokkos::LayoutLeft, DeviceType> t_sna_2c_ll;
  typedef Kokkos::View<SNAcomplex**, Kokkos::LayoutRight, DeviceType> t_sna_2c_lr;
  typedef Kokkos::View<SNAcomplex***, DeviceType> t_sna_3c;
  typedef Kokkos::View<SNAcomplex***, Kokkos::LayoutLeft, DeviceType> t_sna_3c_ll;
  typedef Kokkos::View<SNAcomplex***[3], DeviceType> t_sna_4c;
  typedef Kokkos::View<SNAcomplex***[3], Kokkos::LayoutLeft, DeviceType> t_sna_4c_ll;
  typedef Kokkos::View<SNAcomplex**[3], DeviceType> t_sna_3c3;
  typedef Kokkos::View<SNAcomplex*****, DeviceType> t_sna_5c;

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

  double memory_usage();

  int ncoeff;

  // functions for bispectrum coefficients
  KOKKOS_INLINE_FUNCTION
  void pre_ui(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team,const int&); // ForceSNAP
  KOKKOS_INLINE_FUNCTION
  void compute_ui(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, const int, const int); // ForceSNAP
  KOKKOS_INLINE_FUNCTION
  void compute_ui_cpu(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int, int); // ForceSNAP
  KOKKOS_INLINE_FUNCTION
  void compute_zi(const int&);    // ForceSNAP
  KOKKOS_INLINE_FUNCTION
  void zero_yi(const int&,const int&); // ForceSNAP
  KOKKOS_INLINE_FUNCTION
  void compute_yi(int,
   const Kokkos::View<F_FLOAT**, DeviceType> &beta); // ForceSNAP
  KOKKOS_INLINE_FUNCTION
  void compute_bi(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int);    // ForceSNAP

  // functions for derivatives

  KOKKOS_INLINE_FUNCTION
  void compute_fused_deidrj(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, const int, const int); //ForceSNAP
  KOKKOS_INLINE_FUNCTION
  void compute_duidrj_cpu(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int, int); //ForceSNAP
  KOKKOS_INLINE_FUNCTION
  void compute_deidrj_cpu(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int, int); // ForceSNAP
  KOKKOS_INLINE_FUNCTION
  double compute_sfac(double, double); // add_uarraytot, compute_duarray
  KOKKOS_INLINE_FUNCTION
  double compute_dsfac(double, double); // compute_duarray

  // efficient complex FMA
  // efficient caxpy (i.e., y += a x)
  static KOKKOS_FORCEINLINE_FUNCTION
  void caxpy(const SNAcomplex& a, const SNAcomplex& x, SNAcomplex& y);

  // efficient complex FMA, conjugate of scalar
  static KOKKOS_FORCEINLINE_FUNCTION
  void caconjxpy(const SNAcomplex& a, const SNAcomplex& x, SNAcomplex& y);

  // Set the direction for split ComputeDuidrj
  KOKKOS_INLINE_FUNCTION
  void set_dir(int);

#ifdef TIMING_INFO
  double* timers;
  timespec starttime, endtime;
  int print;
  int counter;
#endif

  //per sna class instance for OMP use


  // Per InFlight Particle
  t_sna_3d rij;
  t_sna_2i inside;
  t_sna_2d wj;
  t_sna_2d rcutij;
  t_sna_3d dedr;
  int natom, nmax;

  void grow_rij(int, int);

  int twojmax, diagonalstyle;

  t_sna_2d_ll blist;
  t_sna_2c_ll ulisttot;
  t_sna_2c_ll zlist;

  t_sna_3c_ll ulist;
  t_sna_2c_ll ylist;

  // derivatives of data
  t_sna_4c_ll dulist;

  int idxcg_max, idxu_max, idxz_max, idxb_max;

private:
  double rmin0, rfac0;

  //use indexlist instead of loops, constructor generates these
  // Same across all SNAKokkos
  Kokkos::View<int*[10], DeviceType> idxz;
  Kokkos::View<int*[3], DeviceType> idxb;
  Kokkos::View<int***, DeviceType> idxcg_block;
  Kokkos::View<int*, DeviceType> idxu_block;
  Kokkos::View<int***, DeviceType> idxz_block;
  Kokkos::View<int***, DeviceType> idxb_block;

  // data for bispectrum coefficients

  // Same across all SNAKokkos
  t_sna_1d cglist;
  t_sna_2d rootpqarray;

  static const int nmaxfactorial = 167;
  static const double nfac_table[];
  inline
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
  void add_uarraytot(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int, int, double, double, double); // compute_ui

  KOKKOS_INLINE_FUNCTION
  void compute_uarray_cpu(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int, int,
                      double, double, double,
                      double, double); // compute_ui_cpu


  inline
  double deltacg(int, int, int);  // init_clebsch_gordan

inline
  int compute_ncoeff();           // SNAKokkos()
  KOKKOS_INLINE_FUNCTION
  void compute_duarray_cpu(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int, int,
                       double, double, double, // compute_duidrj_cpu
                       double, double, double, double, double);

  // Sets the style for the switching function
  // 0 = none
  // 1 = cosine
  int switch_flag;

  // Self-weight
  double wself;

  int bzero_flag; // 1 if bzero subtracted from barray
  Kokkos::View<double*, DeviceType> bzero; // array of B values for isolated atoms

  // for per-direction dulist calculation, specify the direction.
  int dir;
};

}

#include "sna_kokkos_impl.h"
#endif

/* ERROR/WARNING messages:

E: Invalid argument to factorial %d

N must be >= 0 and <= 167, otherwise the factorial result is too
large.

*/
