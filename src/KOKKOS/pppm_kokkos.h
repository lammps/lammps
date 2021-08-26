/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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

#include "gridcomm_kokkos.h"
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
  virtual ~PPPMKokkos();
  virtual void init();
  virtual void setup();
  void setup_grid();
  virtual void settings(int, char **);
  virtual void compute(int, int);
  virtual int timing_1d(int, double &);
  virtual int timing_3d(int, double &);
  virtual double memory_usage();

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
  void x2lamdaT(double* v, double* lamda) const
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

  int factors[3];

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
  GridCommKokkos<DeviceType> *gc;

  FFT_DAT::tdual_FFT_SCALAR_1d k_gc_buf1,k_gc_buf2;
  int ngc_buf1,ngc_buf2,npergrid;

  //int **part2grid;             // storage for particle -> grid mapping
  typename AT::t_int_1d_3 d_part2grid;

  //double *boxlo;
  double boxlo[3];

  void set_grid_global();
  void set_grid_local();
  void adjust_gewald();
  double newton_raphson_f();
  double derivf();
  double final_accuracy();

  virtual void allocate();
  virtual void allocate_peratom();
  virtual void deallocate();
  virtual void deallocate_peratom();
  int factorable(int);
  double compute_df_kspace();
  double estimate_ik_error(double, double, bigint);
  virtual void compute_gf_denom();
  virtual void compute_gf_ik();

  virtual void particle_map();
  virtual void make_rho();
  virtual void brick2fft();

  virtual void poisson();
  virtual void poisson_ik();

  virtual void fieldforce();
  virtual void fieldforce_ik();

  virtual void poisson_peratom();
  virtual void fieldforce_peratom();
  void procs2grid2d(int,int,int,int *, int*);

  KOKKOS_INLINE_FUNCTION
  void compute_rho1d(const int i, const FFT_SCALAR &, const FFT_SCALAR &,
                     const FFT_SCALAR &) const;
  void compute_rho_coeff();
  void slabcorr();

  // grid communication

  void pack_forward_grid_kokkos(int, FFT_DAT::tdual_FFT_SCALAR_1d &, int, DAT::tdual_int_2d &, int);
  void unpack_forward_grid_kokkos(int, FFT_DAT::tdual_FFT_SCALAR_1d &, int, int, DAT::tdual_int_2d &, int);
  void pack_reverse_grid_kokkos(int, FFT_DAT::tdual_FFT_SCALAR_1d &, int, DAT::tdual_int_2d &, int);
  void unpack_reverse_grid_kokkos(int, FFT_DAT::tdual_FFT_SCALAR_1d &, int, int, DAT::tdual_int_2d &, int);

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

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot (yet) use PPPM Kokkos with 'kspace_modify diff ad'

UNDOCUMENTED

E: Cannot (yet) use PPPM with triclinic box and slab correction

This feature is not yet supported.

E: Cannot use PPPM with 2d simulation

The kspace style pppm cannot be used in 2d simulations.  You can use
2d PPPM in a 3d simulation; see the kspace_modify command.

E: PPPM can only currently be used with comm_style brick

This is a current restriction in LAMMPS.

E: Kspace style requires atomKK attribute q

UNDOCUMENTED

E: Cannot use non-periodic boundaries with PPPM

For kspace style pppm, all 3 dimensions must have periodic boundaries
unless you use the kspace_modify command to define a 2d slab with a
non-periodic z dimension.

E: Incorrect boundaries with slab PPPM

Must have periodic x,y dimensions and non-periodic z dimension to use
2d slab option with PPPM.

E: PPPM order cannot be < 2 or > than %d

This is a limitation of the PPPM implementation in LAMMPS.

E: KSpace style is incompatible with Pair style

Setting a kspace style requires that a pair style with matching
long-range Coulombic or dispersion components be used.

E: Cannot (yet) use PPPM Kokkos TIP4P

UNDOCUMENTED

W: Reducing PPPM order b/c stencil extends beyond nearest neighbor processor

This may lead to a larger grid than desired.  See the kspace_modify overlap
command to prevent changing of the PPPM order.

E: PPPM order < minimum allowed order

The default minimum order is 2.  This can be reset by the
kspace_modify minorder command.

E: PPPM grid stencil extends beyond nearest neighbor processor

This is not allowed if the kspace_modify overlap setting is no.

E: KSpace accuracy must be > 0

The kspace accuracy designated in the input must be greater than zero.

E: Must use 'kspace_modify gewald' for uncharged system

UNDOCUMENTED

E: PPPM grid is too large

The global PPPM grid is larger than OFFSET in one or more dimensions.
OFFSET is currently set to 4096.  You likely need to decrease the
requested accuracy.

E: Could not compute g_ewald

The Newton-Raphson solver failed to converge to a good value for
g_ewald.  This error should not occur for typical problems.  Please
send an email to the developers.

E: Non-numeric box dimensions - simulation unstable

The box size has apparently blown up.

E: Out of range atoms - cannot compute PPPM

One or more atoms are attempting to map their charge to a PPPM grid
point that is not owned by a processor.  This is likely for one of two
reasons, both of them bad.  First, it may mean that an atom near the
boundary of a processor's sub-domain has moved more than 1/2 the
"neighbor skin distance"_neighbor.html without neighbor lists being
rebuilt and atoms being migrated to new processors.  This also means
you may be missing pairwise interactions that need to be computed.
The solution is to change the re-neighboring criteria via the
"neigh_modify"_neigh_modify command.  The safest settings are "delay 0
every 1 check yes".  Second, it may mean that an atom has moved far
outside a processor's sub-domain or even the entire simulation box.
This indicates bad physics, e.g. due to highly overlapping atoms, too
large a timestep, etc.

U: Cannot (yet) use PPPM with triclinic box and kspace_modify diff ad

This feature is not yet supported.

U: Kspace style requires atom attribute q

The atom style defined does not have these attributes.

U: Pair style is incompatible with TIP4P KSpace style

The pair style does not have the requires TIP4P settings.

U: Bond and angle potentials must be defined for TIP4P

Cannot use TIP4P pair potential unless bond and angle potentials
are defined.

U: Bad TIP4P angle type for PPPM/TIP4P

Specified angle type is not valid.

U: Bad TIP4P bond type for PPPM/TIP4P

Specified bond type is not valid.

U: Cannot (yet) use PPPM with triclinic box and TIP4P

This feature is not yet supported.

U: Could not compute grid size

The code is unable to compute a grid size consistent with the desired
accuracy.  This error should not occur for typical problems.  Please
send an email to the developers.

*/

