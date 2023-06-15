/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(pppm/kk,PPPMKokkos<LMPDeviceType>);
KSpaceStyle(pppm/kk/device,PPPMKokkos<LMPDeviceType>);
KSpaceStyle(pppm/kk/host,PPPMKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PPPM_KOKKOS_H
#define LMP_PPPM_KOKKOS_H

#include "grid3d_kokkos.h"
#include "remap_kokkos.h"
#include "fft3d_kokkos.h"
#include "kokkos_base_fft.h"
#include "fftdata_kokkos.h"
#include "kokkos_type.h"
#include "kokkos_few.h"

// clang-format off

// fix up FFT defines for KOKKOS with CUDA

#ifdef KOKKOS_ENABLE_CUDA
# if defined(FFT_FFTW)
#  undef FFT_FFTW
# endif
# if defined(FFT_FFTW3)
#  undef FFT_FFTW3
# endif
# if defined(FFT_MKL)
#  undef FFT_MKL
# endif
# if !defined(FFT_CUFFT) && !defined(FFT_KISSFFT)
#  define FFT_KISSFFT
# endif
#endif

#include "pppm.h"

namespace LAMMPS_NS {

struct TagPPPM_setup1{};
struct TagPPPM_setup2{};
struct TagPPPM_setup3{};
struct TagPPPM_setup4{};
struct TagPPPM_setup_triclinic1{};
struct TagPPPM_setup_triclinic2{};
struct TagPPPM_compute_gf_ik{};
struct TagPPPM_compute_gf_ik_triclinic{};
struct TagPPPM_self1{};
struct TagPPPM_self2{};
struct TagPPPM_brick2fft{};
struct TagPPPM_particle_map{};
struct TagPPPM_make_rho_zero{};
struct TagPPPM_make_rho_atomic{};
struct TagPPPM_make_rho{};
struct TagPPPM_poisson_ik1{};
struct TagPPPM_poisson_ik2{};
struct TagPPPM_poisson_ik3{};
struct TagPPPM_poisson_ik4{};
struct TagPPPM_poisson_ik5{};
struct TagPPPM_poisson_ik6{};
struct TagPPPM_poisson_ik7{};
struct TagPPPM_poisson_ik8{};
struct TagPPPM_poisson_ik9{};
struct TagPPPM_poisson_ik10{};
struct TagPPPM_poisson_ik_triclinic1{};
struct TagPPPM_poisson_ik_triclinic2{};
struct TagPPPM_poisson_ik_triclinic3{};
struct TagPPPM_poisson_ik_triclinic4{};
struct TagPPPM_poisson_ik_triclinic5{};
struct TagPPPM_poisson_ik_triclinic6{};
struct TagPPPM_poisson_peratom1{};
struct TagPPPM_poisson_peratom2{};
struct TagPPPM_poisson_peratom3{};
struct TagPPPM_poisson_peratom4{};
struct TagPPPM_poisson_peratom5{};
struct TagPPPM_poisson_peratom6{};
struct TagPPPM_poisson_peratom7{};
struct TagPPPM_poisson_peratom8{};
struct TagPPPM_poisson_peratom9{};
struct TagPPPM_poisson_peratom10{};
struct TagPPPM_poisson_peratom11{};
struct TagPPPM_poisson_peratom12{};
struct TagPPPM_poisson_peratom13{};
struct TagPPPM_poisson_peratom14{};
struct TagPPPM_fieldforce_ik{};
struct TagPPPM_fieldforce_peratom{};
struct TagPPPM_pack_forward1{};
struct TagPPPM_pack_forward2{};
struct TagPPPM_unpack_forward1{};
struct TagPPPM_unpack_forward2{};
struct TagPPPM_pack_reverse{};
struct TagPPPM_unpack_reverse{};
struct TagPPPM_slabcorr1{};
struct TagPPPM_slabcorr2{};
struct TagPPPM_slabcorr3{};
struct TagPPPM_slabcorr4{};
struct TagPPPM_timing_zero{};

template<class DeviceType>
class PPPMKokkos : public PPPM, public KokkosBaseFFT {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef FFTArrayTypes<DeviceType> FFT_AT;

  PPPMKokkos(class LAMMPS *);
  ~PPPMKokkos() override;
  void init() override;
  void setup() override;
  void compute(int, int) override;
  int timing_1d(int, double &) override;
  int timing_3d(int, double &) override;
  double memory_usage() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_setup1, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_setup2, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_setup3, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_setup4, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_setup_triclinic1, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_setup_triclinic2, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_compute_gf_ik, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_compute_gf_ik_triclinic, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_self1, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_self2, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_brick2fft, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_particle_map, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_make_rho_zero, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_make_rho_atomic, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_make_rho, typename Kokkos::TeamPolicy<DeviceType, TagPPPM_make_rho>::member_type) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_ik1, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_ik2, const int&, EV_FLOAT &) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_ik3, const int&, EV_FLOAT &) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_ik4, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_ik5, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_ik6, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_ik7, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_ik8, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_ik9, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_ik10, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_ik_triclinic1, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_ik_triclinic2, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_ik_triclinic3, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_ik_triclinic4, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_ik_triclinic5, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_ik_triclinic6, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_peratom1, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_peratom2, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_peratom3, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_peratom4, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_peratom5, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_peratom6, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_peratom7, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_peratom8, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_peratom9, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_peratom10, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_peratom11, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_peratom12, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_peratom13, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_poisson_peratom14, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_fieldforce_ik, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_fieldforce_peratom, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_pack_forward1, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_pack_forward2, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_unpack_forward1, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_unpack_forward2, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_pack_reverse, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_unpack_reverse, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_slabcorr1, const int&, double&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_slabcorr2, const int&, double&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_slabcorr3, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_slabcorr4, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_timing_zero, const int&) const;

  //template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  //KOKKOS_INLINE_FUNCTION
  //void operator()(TagPPPMKernelA<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

 protected:
  double unitkx,unitky,unitkz;
  double scaleinv,s2;
  double qscale,efact,ffact,dipole_all,dipole_r2,zprd;
  double xprd,yprd,zprd_slab;
  int nbx,nby,nbz,twoorder;
  int numx_fft,numy_fft,numz_fft;
  int numx_inout,numy_inout,numz_inout;
  int numx_out,numy_out,numz_out;
  int ix,iy,nlocal;

  // Local copies of the domain box tilt etc.
  Few<double,6> h, h_inv;

  KOKKOS_INLINE_FUNCTION
  void x2lamdaT_kokkos(double* v, double* lamda) const
  {
    double lamda_tmp[3];

    lamda_tmp[0] = h_inv[0]*v[0];
    lamda_tmp[1] = h_inv[5]*v[0] + h_inv[1]*v[1];
    lamda_tmp[2] = h_inv[4]*v[0] + h_inv[3]*v[1] + h_inv[2]*v[2];

    lamda[0] = lamda_tmp[0];
    lamda[1] = lamda_tmp[1];
    lamda[2] = lamda_tmp[2];
  }

  int nx,ny,nz;
  typename AT::t_int_1d_um d_list_index;
  typename FFT_AT::t_FFT_SCALAR_1d_um d_buf;
  int unpack_offset;

  DAT::tdual_int_scalar k_flag;

  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_float_1d_randomread q;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename ArrayTypes<DeviceType>::t_efloat_1d d_eatom;
  typename ArrayTypes<DeviceType>::t_virial_array d_vatom;

  typename FFT_AT::t_FFT_SCALAR_3d d_density_brick;
  typename FFT_AT::t_FFT_SCALAR_3d d_vdx_brick,d_vdy_brick,d_vdz_brick;
  typename FFT_AT::t_FFT_SCALAR_3d d_u_brick;
  typename FFT_AT::t_FFT_SCALAR_3d d_v0_brick,d_v1_brick,d_v2_brick;
  typename FFT_AT::t_FFT_SCALAR_3d d_v3_brick,d_v4_brick,d_v5_brick;
  typename AT::t_float_1d d_greensfn;
  typename AT::t_virial_array d_vg;
  typename AT::t_float_1d d_fkx;
  typename AT::t_float_1d d_fky;
  typename AT::t_float_1d d_fkz;
  FFT_DAT::tdual_FFT_SCALAR_1d k_density_fft;
  FFT_DAT::tdual_FFT_SCALAR_1d k_work1;
  FFT_DAT::tdual_FFT_SCALAR_1d k_work2;
  typename FFT_AT::t_FFT_SCALAR_1d d_density_fft;
  typename FFT_AT::t_FFT_SCALAR_1d d_work1;
  typename FFT_AT::t_FFT_SCALAR_1d d_work2;

  DAT::tdual_float_1d k_gf_b;
  typename AT::t_float_1d d_gf_b;

  //FFT_SCALAR **rho1d,**rho_coeff,**drho1d,**drho_coeff;
  typename FFT_AT::t_FFT_SCALAR_2d_3 d_rho1d;
  FFT_DAT::tdual_FFT_SCALAR_2d k_rho_coeff;
  typename FFT_AT::t_FFT_SCALAR_2d d_rho_coeff;
  FFT_HAT::t_FFT_SCALAR_2d h_rho_coeff;
  //double **acons;
  typename Kokkos::DualView<F_FLOAT[8][7],Kokkos::LayoutRight,DeviceType>::t_host acons;

  // FFTs and grid communication

  FFT3dKokkos<DeviceType> *fft1,*fft2;
  RemapKokkos<DeviceType> *remap;
  Grid3dKokkos<DeviceType> *gc;

  FFT_DAT::tdual_FFT_SCALAR_1d k_gc_buf1,k_gc_buf2;
  int ngc_buf1,ngc_buf2,npergrid;

  //int **part2grid;             // storage for particle -> grid mapping
  typename AT::t_int_1d_3 d_part2grid;

  double boxlo[3];

  void set_grid_local() override;

  void allocate() override;
  void allocate_peratom() override;
  void deallocate() override;
  void deallocate_peratom() override;
  double estimate_ik_error(double, double, bigint) override;
  void compute_gf_denom() override;
  void compute_gf_ik() override;

  void particle_map() override;
  void make_rho() override;
  void brick2fft() override;

  void poisson_ik() override;

  void fieldforce() override;
  void fieldforce_ik() override;

  void poisson_peratom() override;
  void fieldforce_peratom() override;

  KOKKOS_INLINE_FUNCTION
  void compute_rho1d(const int i, const FFT_SCALAR &, const FFT_SCALAR &,
                     const FFT_SCALAR &) const;
  void compute_rho_coeff();
  void slabcorr() override;

  // grid communication

  void pack_forward_grid_kokkos(int, FFT_DAT::tdual_FFT_SCALAR_1d &, int, DAT::tdual_int_2d &, int) override;
  void unpack_forward_grid_kokkos(int, FFT_DAT::tdual_FFT_SCALAR_1d &, int, int, DAT::tdual_int_2d &, int) override;
  void pack_reverse_grid_kokkos(int, FFT_DAT::tdual_FFT_SCALAR_1d &, int, DAT::tdual_int_2d &, int) override;
  void unpack_reverse_grid_kokkos(int, FFT_DAT::tdual_FFT_SCALAR_1d &, int, int, DAT::tdual_int_2d &, int) override;

  // triclinic

  int triclinic;               // domain settings, orthog or triclinic
  void setup_triclinic();
  void compute_gf_ik_triclinic();
  void poisson_ik_triclinic();

/* ----------------------------------------------------------------------
   denominator for Hockney-Eastwood Green's function
     of x,y,z = sin(kx*deltax/2), etc

            inf                 n-1
   S(n,k) = Sum  W(k+pi*j)**2 = Sum b(l)*(z*z)**l
           j=-inf               l=0

          = -(z*z)**n /(2n-1)! * (d/dx)**(2n-1) cot(x)  at z = sin(x)
   gf_b = denominator expansion coeffs
------------------------------------------------------------------------- */

  KOKKOS_INLINE_FUNCTION
  double gf_denom(const double &x, const double &y,
                         const double &z) const {
    double sx,sy,sz;
    sz = sy = sx = 0.0;
    for (int l = order-1; l >= 0; l--) {
      sx = d_gf_b[l] + sx*x;
      sy = d_gf_b[l] + sy*y;
      sz = d_gf_b[l] + sz*z;
    }
    double s = sx*sy*sz;
    return s*s;
  };
};

}

#endif
#endif


