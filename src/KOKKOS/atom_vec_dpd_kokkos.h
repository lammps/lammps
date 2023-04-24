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

  void grow(int) override;
  void grow_pointers() override;
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

  double *duChem;

 protected:
  DAT::t_efloat_1d d_uCond, d_uMech, d_uChem, d_uCG, d_uCGnew,d_rho,d_dpdTheta,d_duChem;
  HAT::t_efloat_1d h_uCond, h_uMech, h_uChem, h_uCG, h_uCGnew,h_rho,h_dpdTheta,h_duChem;

  DAT::t_tagint_1d d_tag;
  HAT::t_tagint_1d h_tag;
  DAT::t_imageint_1d d_image;
  HAT::t_imageint_1d h_image;
  DAT::t_int_1d d_type, d_mask;
  HAT::t_int_1d h_type, h_mask;

  DAT::t_x_array d_x;
  DAT::t_v_array d_v;
  DAT::t_f_array d_f;
};

}    // namespace LAMMPS_NS

#endif
#endif
