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
AtomStyle(hybrid/kk,AtomVecHybridKokkos);
AtomStyle(hybrid/kk/device,AtomVecHybridKokkos);
AtomStyle(hybrid/kk/host,AtomVecHybridKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_ATOM_VEC_HYBRID_KOKKOS_H
#define LMP_ATOM_VEC_HYBRID_KOKKOS_H

#include "atom_vec_kokkos.h"
#include "atom_vec_hybrid.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class AtomVecHybridKokkos : public AtomVecKokkos, public AtomVecHybrid {
 public:
  AtomVecHybridKokkos(class LAMMPS *);

  void grow(int) override;

  int pack_comm_kokkos(const int &n, const DAT::tdual_int_2d &k_sendlist,
                       const int & iswap,
                       const DAT::tdual_xfloat_2d &buf,
                       const int &pbc_flag, const int pbc[]) override;
  void unpack_comm_kokkos(const int &n, const int &nfirst,
                          const DAT::tdual_xfloat_2d &buf) override;
  int pack_comm_self(const int &n, const DAT::tdual_int_2d &list,
                     const int & iswap, const int nfirst,
                     const int &pbc_flag, const int pbc[]) override;
  int pack_border_kokkos(int n, DAT::tdual_int_2d k_sendlist,
                         DAT::tdual_xfloat_2d buf,int iswap,
                         int pbc_flag, int *pbc, ExecutionSpace space) override;
  void unpack_border_kokkos(const int &n, const int &nfirst,
                            const DAT::tdual_xfloat_2d &buf,
                            ExecutionSpace space) override;
  int pack_exchange_kokkos(const int &nsend,DAT::tdual_xfloat_2d &buf,
                           DAT::tdual_int_1d k_sendlist,
                           DAT::tdual_int_1d k_copylist,
                           ExecutionSpace space) override;
  int unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf, int nrecv,
                             int nlocal, int dim, X_FLOAT lo, X_FLOAT hi,
                             ExecutionSpace space,
                             DAT::tdual_int_1d &k_indices) override;

  void sync(ExecutionSpace space, unsigned int mask) override;
  void modified(ExecutionSpace space, unsigned int mask) override;
  void sync_overlapping_device(ExecutionSpace space, unsigned int mask) override;

 private:
  DAT::t_tagint_1d d_tag;
  DAT::t_int_1d d_type, d_mask;
  HAT::t_tagint_1d h_tag;
  HAT::t_int_1d h_type, h_mask;

  DAT::t_imageint_1d d_image;
  HAT::t_imageint_1d h_image;

  DAT::t_x_array d_x;
  DAT::t_v_array d_v;
  DAT::t_f_array d_f;
  HAT::t_x_array h_x;
  HAT::t_v_array h_v;
  HAT::t_f_array h_f;

  DAT::t_v_array d_omega, d_angmom;
  HAT::t_v_array h_omega, h_angmom;
};

} // namespace LAMMPS_NS

#endif
#endif
