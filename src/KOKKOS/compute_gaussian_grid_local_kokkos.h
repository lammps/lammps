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
ComputeStyle(gaussian/grid/local/kk,ComputeGaussianGridLocalKokkos<LMPDeviceType>);
ComputeStyle(gaussian/grid/local/kk/device,ComputeGaussianGridLocalKokkos<LMPDeviceType>);
ComputeStyle(gaussian/grid/local/kk/host,ComputeGaussianGridLocalKokkos<LMPHostType>);
// clang-format on

#else

#ifndef LMP_COMPUTE_GAUSSIAN_GRID_LOCAL_KOKKOS_H
#define LMP_COMPUTE_GAUSSIAN_GRID_LOCAL_KOKKOS_H

#include "compute_gaussian_grid_local.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

// clang-format off
struct TagComputeGaussianGridLocalNeigh{};
// clang-format on

template <class DeviceType> class ComputeGaussianGridLocalKokkos : public ComputeGaussianGridLocal {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  // Static team/tile sizes for device offload

#ifdef KOKKOS_ENABLE_HIP
  static constexpr int team_size_compute_neigh = 2;
#else
  static constexpr int team_size_compute_neigh = 4;
#endif

  ComputeGaussianGridLocalKokkos(class LAMMPS *, int, char **);
  ~ComputeGaussianGridLocalKokkos() override;
  void setup() override;
  void init() override;
  void compute_local() override;

  template<class TagStyle>
  void check_team_size_for(int, int&, int);

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeGaussianGridLocalNeigh, const typename Kokkos::TeamPolicy<DeviceType, TagComputeGaussianGridLocalNeigh>::member_type& team) const;

 private:
  Kokkos::View<double*, DeviceType> d_radelem;              // element radii
  Kokkos::View<double*, DeviceType> d_sigmaelem;
  Kokkos::View<double*, DeviceType> d_prefacelem;
  Kokkos::View<double*, DeviceType> d_argfacelem;
  Kokkos::View<int*, DeviceType> d_ninside;                // ninside for all atoms in list
  Kokkos::View<int*, DeviceType> d_map;                    // mapping from atom types to elements

  typedef Kokkos::DualView<F_FLOAT**, DeviceType> tdual_fparams;
  tdual_fparams k_cutsq;
  typedef Kokkos::View<const F_FLOAT**, DeviceType,
      Kokkos::MemoryTraits<Kokkos::RandomAccess> > t_fparams_rnd;
  t_fparams_rnd rnd_cutsq;


  int max_neighs, inum, chunk_size, chunk_offset;
  int host_flag;
  int total_range; // total number of loop iterations in grid
  int xlen, ylen, zlen;
  int chunksize;
  int ntotal;

  typename AT::t_x_array_randomread x;
  typename AT::t_int_1d_randomread type;

  DAT::tdual_float_2d k_alocal;
  typename AT::t_float_2d d_alocal;

  // triclinic vars
  double h0, h1, h2, h3, h4, h5;
  double lo0, lo1, lo2;
};

}    // namespace LAMMPS_NS

#endif
#endif
