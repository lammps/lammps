/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale AtomicKokkos/Molecular Massively Parallel Simulator
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
AtomStyle(dpd/kk,AtomVecDPDKokkos);
AtomStyle(dpd/kk/device,AtomVecDPDKokkos);
AtomStyle(dpd/kk/host,AtomVecDPDKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_ATOM_VEC_DPD_KOKKOS_H
#define LMP_ATOM_VEC_DPD_KOKKOS_H

#include "atom_vec_kokkos.h"
#include "atom_vec_dpd.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class AtomVecDPDKokkos : public AtomVecKokkos, public AtomVecDPD {
 public:
  AtomVecDPDKokkos(class LAMMPS *);
  void init() override;

  void grow(int) override;
  void grow_pointers() override;
  void sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) override;
  void sync(ExecutionSpace space, uint64_t mask) override;
  void modified(ExecutionSpace space, uint64_t mask) override;
  void sync_pinned(ExecutionSpace space, uint64_t mask, int async_flag = 0) override;

  double *duChem;
};

}    // namespace LAMMPS_NS

#endif
#endif
