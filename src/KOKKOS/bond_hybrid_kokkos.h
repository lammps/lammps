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

#ifdef BOND_CLASS
// clang-format off
BondStyle(hybrid/kk,BondHybridKokkos);
BondStyle(hybrid/kk/device,BondHybridKokkos);
BondStyle(hybrid/kk/host,BondHybridKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_BOND_HYBRID_KOKKOS_H
#define LMP_BOND_HYBRID_KOKKOS_H

#include "bond_hybrid.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class BondHybridKokkos : public BondHybrid {
  friend class Force;

 public:
  BondHybridKokkos(class LAMMPS *);
  ~BondHybridKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double memory_usage() override;

 private:
  int maxbond_all;

  class NeighborKokkos *neighborKK;

  DAT::tdual_int_1d k_map;       // which style each bond type points to
  DAT::tdual_int_1d k_nbondlist; // # of bonds in sub-style bondlists
  DAT::tdual_int_3d k_bondlist;  // bondlist for each sub-style

  void allocate() override;
  void deallocate() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
