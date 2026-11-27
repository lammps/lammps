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

/* ----------------------------------------------------------------------
   lj/cut/coul/long2/kk -- experimental GPU variant of lj/cut/coul/long/kk
   with the same physics.  Force-only steps bypass the generic pair_kokkos.h
   kernel in favor of hand-written CUDA kernels: a warp-compacted kernel
   over the half neighbor list, or (opt-in, single MPI rank) a cluster-pair
   list.  Energy/virial steps and unsupported configurations fall back to
   PairLJCutCoulLongKokkos::compute().  See the comments in the .cpp file.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/coul/long2/kk,PairLJCutCoulLong2Kokkos<LMPDeviceType>);
PairStyle(lj/cut/coul/long2/kk/device,PairLJCutCoulLong2Kokkos<LMPDeviceType>);
PairStyle(lj/cut/coul/long2/kk/host,PairLJCutCoulLong2Kokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_CUT_COUL_LONG2_KOKKOS_H
#define LMP_PAIR_LJ_CUT_COUL_LONG2_KOKKOS_H

#include "pair_lj_cut_coul_long_kokkos.h"
#include "kokkos_base.h"

namespace LAMMPS_NS {

// packed j-side gather record: position and charge of one atom, one aligned load
struct alignas(16) LJCL2XQ { KK_FLOAT x, y, z, q; };

// cluster-pair list entry (8 i x 4 j tile), see the cluster block comment in the .cpp
//   cj    j-cluster index (4 consecutive cluster slots)
//   mask  bit kj*8+ki set iff slot pair (ki,kj) is in the list
//   sidx  index into the special-code array, or -1 if the tile has no special pair
//   shift periodic image of the j-cluster, (sx+1) + 3*(sy+1) + 9*(sz+1)
struct alignas(16) LJCL2CTile { int cj; unsigned int mask; int sidx; int shift; };
// axis-aligned bounding box of a cluster (at list-build positions)
struct LJCL2BB { float lo[3], hi[3]; };

template<class DeviceType>
class PairLJCutCoulLong2Kokkos : public PairLJCutCoulLongKokkos<DeviceType>, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  PairLJCutCoulLong2Kokkos(class LAMMPS *);
  ~PairLJCutCoulLong2Kokkos() override;

  void compute(int, int) override;
  void init_style() override;
  void finish_async_compute() override;
  bool async_compute_pending() const override { return cl_pending; }

  // public only because nvcc forbids extended device lambdas inside
  // private/protected member functions
  void cl_build();
  void cl_step(int eflag, int vflag);
  template<bool EV> void cl_force();
  template<bool EV, bool UNICUT> void cl_force_launch();

 protected:
  // ---- compact kernel: packed (x,y,z,q) record and byte-wide type per atom
  Kokkos::View<LJCL2XQ*, DeviceType> d_xq;
  Kokkos::View<unsigned char*, DeviceType> d_type8;

  template<class FunctorT> void fill_common(FunctorT &ff);
  void pack_refresh();
  void compact_compute();

  // ---- cluster-pair list path (LMP_LJCL2_CLUSTER=1, CUDA, one MPI rank,
  // orthogonal periodic box).  Atoms are binned into x-y columns, sorted by
  // z and cut into 8-atom i-clusters / 4-atom j-clusters; the tile list is
  // searched directly with bounding boxes (no flat list is built: the pair's
  // neighbor request is made occasional), periodic images use shift vectors
  // instead of ghost atoms, and an inner list pruned to rc + delta is
  // refreshed on the device whenever an atom moved more than delta/2.
  int cl_on = -1;                     // -1 = not probed
  bigint cl_built_step = -1;
  int cl_ncx = 0, cl_ncy = 0, cl_ncol = 0, cl_nslot = 0, cl_nicl = 0;
  int cl_maxout = 0, cl_maxsc = 0;
  double cl_colw[2] = {0.0, 0.0};
  double cl_delta = 2.0;              // inner-list buffer (LMP_CL_DELTA)
  Kokkos::View<unsigned int*, DeviceType> d_cl_key;
  Kokkos::View<int*, DeviceType> d_cl_perm, d_cl_colcnt, d_cl_colfirst, d_cl_colslot,
    d_cl_atom, d_cl_icol, d_cl_nout, d_cl_nin, d_cl_ctl;
  Kokkos::View<LJCL2XQ*, DeviceType> d_cl_xq;
  Kokkos::View<float*[4], Kokkos::LayoutRight, DeviceType> d_cl_xp;   // x at last prune
  Kokkos::View<unsigned char*, DeviceType> d_cl_t8;
  Kokkos::View<KK_ACC_FLOAT*[3], Kokkos::LayoutRight, DeviceType> d_cl_f;
  Kokkos::View<LJCL2BB*, DeviceType> d_cl_jbb;
  Kokkos::View<LJCL2CTile**, Kokkos::LayoutRight, DeviceType> d_cl_out, d_cl_in;
  Kokkos::View<unsigned int*[2], Kokkos::LayoutRight, DeviceType> d_cl_sc;
  Kokkos::View<double*, DeviceType> d_cl_ev;

  // LMP_CL_ASYNC: on force-only steps the cluster kernels run on a second
  // stream, overlapping bonded + kspace; finish_async_compute() joins them.
  // While they are pending, VerletKokkos runs kspace before the bonded terms.
  int cl_async = 0;
  int cl_persist = 0;         // LMP_CL_PERSIST: force-kernel blocks per SM (0 = one warp per i-cluster)
  int cl_nsm = 0;
  bool cl_pending = false;
  typename DeviceType::execution_space cl_space;      // the cl_stream instance
  typename DeviceType::execution_space cl_launch_space; // instance for cl_force_launch
#if defined(KOKKOS_ENABLE_CUDA)
  cudaStream_t cl_stream = nullptr;
  cudaEvent_t cl_ev_in = nullptr, cl_ev_out = nullptr;
#endif
  void cl_scatter();
};

}

#endif
#endif
