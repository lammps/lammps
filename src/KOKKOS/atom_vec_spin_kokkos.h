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
  void grow(int) override;
  void grow_pointers() override;
  void force_clear(int, size_t) override;
  void sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) override;
  int pack_border_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                         DAT::tdual_double_2d_lr buf,
                         int pbc_flag, int *pbc, ExecutionSpace space) override;
  void unpack_border_kokkos(const int &n, const int &nfirst,
                            const DAT::tdual_double_2d_lr &buf,
                            ExecutionSpace space) override;
  int pack_exchange_kokkos(const int &nsend,DAT::tdual_double_2d_lr &buf,
                           DAT::tdual_int_1d k_sendlist,
                           DAT::tdual_int_1d k_copylist,
                           ExecutionSpace space) override;
  int unpack_exchange_kokkos(DAT::tdual_double_2d_lr &k_buf, int nrecv,
                             int nlocal, int dim, double lo, double hi,
                             ExecutionSpace space,
                             DAT::tdual_int_1d &k_indices) override;

  void sync(ExecutionSpace space, unsigned int mask) override;
  void modified(ExecutionSpace space, unsigned int mask) override;
  void sync_pinned(ExecutionSpace space, unsigned int mask, int async_flag = 0) override;

 protected:
  DAT::t_tagint_1d d_tag;
  HAT::t_tagint_1d h_tag;

  DAT::t_int_1d d_type, d_mask;
  HAT::t_int_1d h_type, h_mask;

  DAT::t_imageint_1d d_image;
  HAT::t_imageint_1d h_image;

  DAT::t_kkfloat_1d_3_lr d_x;
  DAT::t_kkfloat_1d_3 d_v;
  DAT::t_kkacc_1d_3 d_f;

  DAT::t_kkfloat_1d_4 d_sp;
  DAT::t_kkacc_1d_3 d_fm;
  DAT::t_kkacc_1d_3 d_fm_long;

  HAT::t_kkfloat_1d_4 h_sp;
  HAT::t_kkacc_1d_3 h_fm;
  HAT::t_kkacc_1d_3 h_fm_long;
};

}    // namespace LAMMPS_NS

#endif
#endif
