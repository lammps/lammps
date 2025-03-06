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

#ifdef ANGLE_CLASS
// clang-format off
AngleStyle(hybrid/kk,AngleHybridKokkos);
AngleStyle(hybrid/kk/device,AngleHybridKokkos);
AngleStyle(hybrid/kk/host,AngleHybridKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_ANGLE_HYBRID_KOKKOS_H
#define LMP_ANGLE_HYBRID_KOKKOS_H

#include "angle_hybrid.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class AngleHybridKokkos : public AngleHybrid {
  friend class Force;

 public:
  AngleHybridKokkos(class LAMMPS *);
  ~AngleHybridKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double memory_usage() override;

 private:
  int maxangle_all;

  class NeighborKokkos *neighborKK;

  DAT::tdual_int_1d k_map;       // which style each angle type points to
  DAT::tdual_int_1d k_nanglelist; // # of angles in sub-style anglelists
  DAT::tdual_int_3d k_anglelist;  // anglelist for each sub-style

  void allocate() override;
  void deallocate() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
