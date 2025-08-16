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

#ifdef DIHEDRAL_CLASS
// clang-format off
DihedralStyle(hybrid/kk,DihedralHybridKokkos);
DihedralStyle(hybrid/kk/device,DihedralHybridKokkos);
DihedralStyle(hybrid/kk/host,DihedralHybridKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_DIHEDRAL_HYBRID_KOKKOS_H
#define LMP_DIHEDRAL_HYBRID_KOKKOS_H

#include "dihedral_hybrid.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class DihedralHybridKokkos : public DihedralHybrid {
  friend class Force;

 public:
  DihedralHybridKokkos(class LAMMPS *);
  ~DihedralHybridKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double memory_usage() override;

 private:
  int maxdihedral_all;

  class NeighborKokkos *neighborKK;

  DAT::tdual_int_1d k_map;       // which style each dihedral type points to
  DAT::tdual_int_1d k_ndihedrallist; // # of dihedrals in sub-style dihedrallists
  DAT::tdual_int_3d k_dihedrallist;  // dihedrallist for each sub-style

  void allocate() override;
  void deallocate() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
