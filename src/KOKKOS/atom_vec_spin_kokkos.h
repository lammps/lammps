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
AtomStyle(spin/kk,AtomVecSpinKokkos);
AtomStyle(spin/kk/device,AtomVecSpinKokkos);
AtomStyle(spin/kk/host,AtomVecSpinKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_ATOM_VEC_SPIN_KOKKOS_H
#define LMP_ATOM_VEC_SPIN_KOKKOS_H

#include "atom_vec_kokkos.h"
#include "atom_vec_spin.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class AtomVecSpinKokkos : public AtomVecKokkos, public AtomVecSpin {
 public:
  AtomVecSpinKokkos(class LAMMPS *);
  void init() override;

  void grow(int) override;
  void grow_pointers() override;
  void force_clear(int, size_t) override;
  void sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) override;
  void sync(ExecutionSpace space, uint64_t mask) override;
  void modified(ExecutionSpace space, uint64_t mask) override;
  void sync_pinned(ExecutionSpace space, uint64_t mask, int async_flag = 0) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
