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

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(angle/kk,AtomVecAngleKokkos);
AtomStyle(angle/kk/device,AtomVecAngleKokkos);
AtomStyle(angle/kk/host,AtomVecAngleKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_ATOM_VEC_ANGLE_KOKKOS_H
#define LMP_ATOM_VEC_ANGLE_KOKKOS_H

#include "atom_vec_kokkos.h"
#include "atom_vec_angle.h"

namespace LAMMPS_NS {

class AtomVecAngleKokkos : public AtomVecKokkos, public AtomVecAngle {
 public:
  AtomVecAngleKokkos(class LAMMPS *);
  void init() override;

  void grow(int) override;
  void grow_pointers() override;
  void sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) override;
  void sync(ExecutionSpace space, uint64_t mask) override;
  void modified(ExecutionSpace space, uint64_t mask) override;
  void sync_pinned(ExecutionSpace space, uint64_t mask, int async_flag = 0) override;

 protected:
  tagint *molecule;
  tagint **special;
  tagint **bond_atom;
  tagint **angle_atom1,**angle_atom2,**angle_atom3;
};

}    // namespace LAMMPS_NS

#endif
#endif
