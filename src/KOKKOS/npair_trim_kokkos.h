/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Trimright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef NPAIR_CLASS
// clang-format off
NPairStyle(trim/kk/device,
           NPairTrimKokkos<LMPDeviceType>,
           NP_COPY | NP_TRIM | NP_KOKKOS_DEVICE);

NPairStyle(trim/kk/host,
           NPairTrimKokkos<LMPHostType>,
           NP_COPY | NP_TRIM | NP_KOKKOS_HOST);
// clang-format on
#else

// clang-format off
#ifndef LMP_NPAIR_TRIM_KOKKOS_H
#define LMP_NPAIR_TRIM_KOKKOS_H

#include "npair.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagNPairTrim{};

template<class DeviceType>
class NPairTrimKokkos : public NPair {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  NPairTrimKokkos(class LAMMPS *);
  void build(class NeighList *) override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagNPairTrim, const int&) const;

 private:
  double cutsq_custom;

  typename AT::t_x_array_randomread x;

  typename AT::t_neighbors_2d_const d_neighbors_copy;
  typename AT::t_int_1d_const d_ilist_copy;
  typename AT::t_int_1d_const d_numneigh_copy;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d d_ilist;
  typename AT::t_int_1d d_numneigh;

  void trim_to_kokkos(class NeighList *);
  void trim_to_cpu(class NeighList *);
};

}

#endif
#endif

