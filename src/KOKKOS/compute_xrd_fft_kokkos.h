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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(xrd/fft/kk,ComputeXRDFFTKokkos<LMPDeviceType>);
ComputeStyle(xrd/fft/kk/device,ComputeXRDFFTKokkos<LMPDeviceType>);
ComputeStyle(xrd/fft/kk/host,ComputeXRDFFTKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_XRD_FFT_KOKKOS_H
#define LMP_COMPUTE_XRD_FFT_KOKKOS_H

#include "compute_xrd_fft.h"
#include "kokkos_type.h"
#include "fftdata_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType> class FFT3dKokkos;

struct TagXRDFFTZeroFoot{};
struct TagXRDFFTZeroSlab{};
struct TagXRDFFTSpreadAtomic{};
struct TagXRDFFTSpreadTiled{};
struct TagXRDFFTPack{};
struct TagXRDFFTUnpack{};
struct TagXRDFFTWork{};
struct TagXRDFFTModes{};

// the extent of a rank's atoms in mesh coordinates, reduced over its atoms.
// the default constructor is the identity of the merge below, which is what
// Kokkos initializes a reduction with

struct XRDFFTBounds {
  double umin[3], umax[3];
  int nany;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  XRDFFTBounds() {
    umin[0] = umin[1] = umin[2] = 1.0e300;
    umax[0] = umax[1] = umax[2] = -1.0e300;
    nany = 0;
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  XRDFFTBounds &operator+=(const XRDFFTBounds &rhs) {
    for (int d = 0; d < 3; d++) {
      if (rhs.umin[d] < umin[d]) umin[d] = rhs.umin[d];
      if (rhs.umax[d] > umax[d]) umax[d] = rhs.umax[d];
    }
    nany += rhs.nany;
    return *this;
  }
};

// the two reductions below are their own functors rather than tags on the
// compute, because a functor carries a single value_type and these two differ

template<class DeviceType>
struct XRDFFTBoundsFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef XRDFFTBounds value_type;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_int_1d_randomread mask;
  double mv[3][3];
  int groupbit;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i, XRDFFTBounds &b) const {
    if (!(mask(i) & groupbit)) return;
    const double xi = x(i,0), yi = x(i,1), zi = x(i,2);
    for (int d = 0; d < 3; d++) {
      const double u = xi*mv[d][0] + yi*mv[d][1] + zi*mv[d][2];
      if (u < b.umin[d]) b.umin[d] = u;
      if (u > b.umax[d]) b.umax[d] = u;
    }
    b.nany = 1;
  }
};

template<class DeviceType>
struct XRDFFTBucketFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef int value_type;

  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_int_1d slot_of_type;
  typename AT::t_int_1d slot_atoms;
  int groupbit, slot, offset;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i, int &update, const bool &final) const {
    const int in = ((mask(i) & groupbit) && (slot_of_type(type(i)-1) == slot)) ? 1 : 0;
    if (final && in) slot_atoms(offset + update) = i;
    update += in;
  }
};

template<class DeviceType>
class ComputeXRDFFTKokkos : public ComputeXRDFFT {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef FFTArrayTypes<DeviceType> FFT_AT;

  ComputeXRDFFTKokkos(class LAMMPS *, int, char **);
  ~ComputeXRDFFTKokkos() override;
  void init() override;
  void compute_array() override;
  double memory_usage() override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagXRDFFTZeroFoot, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagXRDFFTZeroSlab, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagXRDFFTSpreadAtomic, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagXRDFFTSpreadTiled,
                  typename Kokkos::TeamPolicy<DeviceType, TagXRDFFTSpreadTiled>::member_type) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagXRDFFTPack, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagXRDFFTUnpack, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagXRDFFTWork, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagXRDFFTModes, const int &) const;

  // the stencil arrays of the spreading kernel live on the stack, so that the
  // kernel needs no per atom scratch in memory.  that fixes the widest stencil
  // this style can spread; wider ones are still available from compute xrd/fft
  // without the accelerator suffix, and are far beyond the width at which the
  // window already reproduces the direct sum to round-off

  static constexpr int MAXORDER = 25;

 protected:
  void allocate_mesh() override;
  void grow_density_own() override;
  int minmax_u(double *, double *) override;
  void refresh_scaling() override;
  void deallocate() override;
  void bucket_atoms() override;
  void fold_reduce(int) override;
  void spread(int) override;

  void setup_device();
  void copy_scaling();
  void set_kernel_state();
  void unpack_kernel(int);

  // the mesh index of one atom along each dimension, and the window there

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int stencil_base(const int, double *, int *) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void stencil_rho(const int, const double, double *) const;

  FFT3dKokkos<DeviceType> *fftkk;

  typename AT::t_kkfloat_1d_3_lr_randomread d_x;
  typename AT::t_int_1d_randomread d_type;
  typename AT::t_int_1d_randomread d_mask;

  DAT::tdual_int_1d k_slot_of_type, k_slot_atoms;
  typename AT::t_int_1d d_slot_of_type, d_slot_atoms;

  DAT::tdual_double_1d k_kb_cheb;
  typename AT::t_double_1d d_kb_cheb;

  DAT::tdual_int_1d k_own_idx;
  typename AT::t_int_1d d_own_idx;

  DAT::tdual_double_1d k_own_deconv, k_own_asf, k_Fre, k_Fim;
  typename AT::t_double_1d d_own_deconv, d_own_asf, d_Fre, d_Fim;

  FFT_DAT::tdual_FFT_SCALAR_1d k_density_own, k_density_slab, k_work1;
  FFT_DAT::tdual_FFT_SCALAR_1d k_sendbuf, k_recvbuf;
  typename FFT_AT::t_FFT_SCALAR_1d d_density_own, d_density_slab, d_work1;
  typename FFT_AT::t_FFT_SCALAR_1d d_sendbuf, d_recvbuf;

  // the buffer an unpack reads from: a message just received, or this rank's
  // own piece, which never became one

  typename FFT_AT::t_FFT_SCALAR_1d d_unpacksrc;

  // an atom whose stencil leaves the footprint of its rank is an internal
  // error, and a kernel cannot report one itself

  DAT::tdual_int_scalar k_overrun;
  typename AT::t_int_scalar d_overrun;

  int device_ready;
  int gpu_aware;      // MPI may read device memory directly
  int mpi_direct;     // ... or the two spaces are the same to begin with

  // state the kernels read, copied with the compute when it is handed to one

  int spread_lo, spread_hi;      // range of slot_atoms this spreading covers
  int bucket_maxatoms;           // allocated length of slot_atoms
  int nmesh_kk[3];               // mesh dimensions
  int foot_lo_kk[3], foot_n_kk[3];
  int fftn_kk[3];
  int nfoot_kk;
  int order_kk, nlower_kk;
  double mesh_vec_kk[3][3];

  // the overlap of a footprint and a brick, as one run of mesh points per
  // dimension, or two once the footprint wraps.  set on the host before each
  // pack or unpack, since deciding what to move is integer work on a handful
  // of numbers rather than something to give a kernel

  int seg_kk[3][2][3];
  int seglen_kk[3];
  int segsrc_n[3];      // dimensions of the buffer the segments index into
  int segdst_lo[3];     // first mesh point of the brick, when unpacking
  int bufoff;           // where this piece starts in the buffer it is packed to
  int slot_off;         // where this element's scattering factors start
};

}    // namespace LAMMPS_NS

#endif
#endif
