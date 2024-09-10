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

#ifdef IMPROPER_CLASS
// clang-format off
ImproperStyle(hybrid/kk,ImproperHybridKokkos);
ImproperStyle(hybrid/kk/device,ImproperHybridKokkos);
ImproperStyle(hybrid/kk/host,ImproperHybridKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_IMPROPER_HYBRID_KOKKOS_H
#define LMP_IMPROPER_HYBRID_KOKKOS_H

#include "improper_hybrid.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class ImproperHybridKokkos : public ImproperHybrid {
  friend class Force;

 public:
  ImproperHybridKokkos(class LAMMPS *);
  ~ImproperHybridKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double memory_usage() override;

 private:
  int maximproper_all;

  class NeighborKokkos *neighborKK;

  DAT::tdual_int_1d k_map;       // which style each improper type points to
  DAT::tdual_int_1d k_nimproperlist; // # of impropers in sub-style improperlists
  DAT::tdual_int_3d k_improperlist;  // improperlist for each sub-style

  void allocate() override;
  void deallocate() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
