// clang-format off
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

#ifndef LMP_ATOM_VEC_KOKKOS_H
#define LMP_ATOM_VEC_KOKKOS_H

#include "atom_vec.h"           //  IWYU pragma: export

#include "kokkos_type.h"
#include <type_traits>

#include <Kokkos_Sort.hpp>

namespace LAMMPS_NS {

class AtomVecKokkos : virtual public AtomVec {
 public:
  AtomVecKokkos(class LAMMPS *);
  ~AtomVecKokkos() override;

  using KeyViewType = DAT::t_kkfloat_1d_3_lr;
  using BinOp = BinOp3DLAMMPS<KeyViewType>;
  virtual void
    sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) = 0;

  virtual void sync(ExecutionSpace space, uint64_t mask) = 0;
  virtual void modified(ExecutionSpace space, uint64_t mask) = 0;
  virtual void sync_pinned(ExecutionSpace space, uint64_t mask, int async_flag = 0) = 0;

  int pack_comm_self(const int &n, const DAT::tdual_int_1d &list,
                     const int nfirst,
                     const int &pbc_flag, const int pbc[]);

  int pack_comm_self_fused(const int &n, const DAT::tdual_int_2d_lr &list,
                           const DAT::tdual_int_1d &sendnum_scan,
                           const DAT::tdual_int_1d &firstrecv,
                           const DAT::tdual_int_1d &pbc_flag,
                           const DAT::tdual_int_2d &pbc,
                           const DAT::tdual_int_1d &g2l);

  int pack_comm_kokkos(const int &n, const DAT::tdual_int_1d &list,
                       const DAT::tdual_double_2d_lr &buf,
                       const int &pbc_flag, const int pbc[]);

  void unpack_comm_kokkos(const int &n, const int &nfirst,
                          const DAT::tdual_double_2d_lr &buf);

  int pack_comm_vel_kokkos(const int &n, const DAT::tdual_int_1d &list,
                           const DAT::tdual_double_2d_lr &buf,
                           const int &pbc_flag, const int pbc[]);

  void unpack_comm_vel_kokkos(const int &n, const int &nfirst,
                              const DAT::tdual_double_2d_lr &buf);

  int pack_reverse_self(const int &n, const DAT::tdual_int_1d &list,
                        const int nfirst);

  int pack_reverse_kokkos(const int &n, const int &nfirst,
                          const DAT::tdual_double_2d_lr &buf);

  void unpack_reverse_kokkos(const int &n, const DAT::tdual_int_1d &list,
                             const DAT::tdual_double_2d_lr &buf);

  int pack_border_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                         DAT::tdual_double_2d_lr buf,
                         int pbc_flag, int *pbc, ExecutionSpace space);

  void unpack_border_kokkos(const int &n, const int &nfirst,
                            const DAT::tdual_double_2d_lr &buf,
                            ExecutionSpace space);

  int pack_border_vel_kokkos(int /*n*/, DAT::tdual_int_1d /*k_sendlist*/,
                             DAT::tdual_double_2d_lr /*buf*/,
                             int /*pbc_flag*/, int * /*pbc*/, ExecutionSpace /*space*/);

  void unpack_border_vel_kokkos(const int &/*n*/, const int & /*nfirst*/,
                                const DAT::tdual_double_2d_lr & /*buf*/,
                                ExecutionSpace /*space*/);

  int pack_exchange_kokkos(const int &nsend, DAT::tdual_double_2d_lr &buf,
                           DAT::tdual_int_1d k_sendlist,
                           DAT::tdual_int_1d k_copylist,
                           ExecutionSpace space);

  int unpack_exchange_kokkos(DAT::tdual_double_2d_lr &k_buf, int nrecv,
                             int nlocal, int dim, double lo, double hi,
                             ExecutionSpace space,
                             DAT::tdual_int_1d &k_indices);

  int size_exchange,size_exchange_default;

  uint64_t datamask_grow;
  uint64_t datamask_comm;
  uint64_t datamask_comm_vel;
  uint64_t datamask_reverse;
  uint64_t datamask_border;
  uint64_t datamask_border_vel;
  uint64_t datamask_exchange;

 protected:
  DAT::t_tagint_1d d_tag;
  DAT::t_int_1d d_type, d_mask;
  HAT::t_tagint_1d h_tag;
  HAT::t_int_1d h_type, h_mask;

  DAT::t_imageint_1d d_image;
  HAT::t_imageint_1d h_image;

  DAT::t_kkfloat_1d_3_lr d_x;
  DAT::t_kkfloat_1d_3 d_v;
  DAT::t_kkacc_1d_3 d_f;
  HAT::t_kkfloat_1d_3_lr h_x;
  HAT::t_kkfloat_1d_3 h_v;
  HAT::t_kkacc_1d_3 h_f;

  DAT::t_kkfloat_1d_3 d_omega, d_angmom;
  HAT::t_kkfloat_1d_3 h_omega, h_angmom;

  // FULL

  DAT::t_kkfloat_1d d_q;
  HAT::t_kkfloat_1d h_q;

  DAT::t_tagint_1d d_molecule;
  DAT::t_int_2d d_nspecial;
  DAT::t_tagint_2d d_special;
  DAT::t_int_1d d_num_bond;
  DAT::t_int_2d d_bond_type;
  DAT::t_tagint_2d d_bond_atom;

  HAT::t_tagint_1d h_molecule;
  HAT::t_int_2d h_nspecial;
  HAT::t_tagint_2d h_special;
  HAT::t_int_1d h_num_bond;
  HAT::t_int_2d h_bond_type;
  HAT::t_tagint_2d h_bond_atom;

  DAT::t_int_1d d_num_angle;
  DAT::t_int_2d d_angle_type;
  DAT::t_tagint_2d d_angle_atom1,d_angle_atom2,d_angle_atom3;

  HAT::t_int_1d h_num_angle;
  HAT::t_int_2d h_angle_type;
  HAT::t_tagint_2d h_angle_atom1,h_angle_atom2,h_angle_atom3;

  DAT::t_int_1d d_num_dihedral;
  DAT::t_int_2d d_dihedral_type;
  DAT::t_tagint_2d d_dihedral_atom1,d_dihedral_atom2,
    d_dihedral_atom3,d_dihedral_atom4;
  DAT::t_int_1d d_num_improper;
  DAT::t_int_2d d_improper_type;
  DAT::t_tagint_2d d_improper_atom1,d_improper_atom2,
    d_improper_atom3,d_improper_atom4;

  HAT::t_int_1d h_num_dihedral;
  HAT::t_int_2d h_dihedral_type;
  HAT::t_tagint_2d h_dihedral_atom1,h_dihedral_atom2,
    h_dihedral_atom3,h_dihedral_atom4;
  HAT::t_int_1d h_num_improper;
  HAT::t_int_2d h_improper_type;
  HAT::t_tagint_2d h_improper_atom1,h_improper_atom2,
    h_improper_atom3,h_improper_atom4;

  DAT::t_kkfloat_1d_4 d_mu;
  HAT::t_kkfloat_1d_4 h_mu;

  DAT::t_kkfloat_1d_4 d_sp;
  DAT::t_kkacc_1d_3 d_fm;
  DAT::t_kkacc_1d_3 d_fm_long;
  HAT::t_kkfloat_1d_4 h_sp;
  HAT::t_kkacc_1d_3 h_fm;
  HAT::t_kkacc_1d_3 h_fm_long;

  DAT::t_kkfloat_1d d_radius;
  HAT::t_kkfloat_1d h_radius;
  DAT::t_kkfloat_1d d_rmass;
  HAT::t_kkfloat_1d h_rmass;
  DAT::t_kkacc_1d_3 d_torque;
  HAT::t_kkacc_1d_3 h_torque;

  DAT::t_kkfloat_1d d_uCond, d_uMech, d_uChem, d_uCG, d_uCGnew,d_rho,d_dpdTheta,d_duChem;
  HAT::t_kkfloat_1d h_uCond, h_uMech, h_uChem, h_uCG, h_uCGnew,h_rho,h_dpdTheta,h_duChem;

  size_t buffer_size;
  void* buffer;

  DAT::tdual_int_1d k_count;

  uint64_t field2mask(std::string);
  int field2size(std::string);
  void set_atom_masks();
  void set_size_exchange();

 public:

  #ifdef LMP_KOKKOS_GPU
  template<class ViewType>
  void perform_pinned_copy(ViewType& src, unsigned int space, int async_flag = 0) {
    typedef Kokkos::View<typename ViewType::data_type,
                 typename ViewType::array_layout,
                 typename std::conditional<
                   std::is_same_v<typename ViewType::execution_space,LMPDeviceType>,
                   LMPPinnedHostType,typename ViewType::memory_space>::type,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged> > mirror_type;
    if (buffer_size == 0) {
       buffer_size = src.span()*sizeof(typename ViewType::value_type);
       buffer = Kokkos::kokkos_malloc<LMPPinnedHostType>(buffer_size);
    } else if (buffer_size < src.span()*sizeof(typename ViewType::value_type)) {
       buffer_size = src.span()*sizeof(typename ViewType::value_type);
       buffer = Kokkos::kokkos_realloc<LMPPinnedHostType>(buffer,buffer_size);
    }

    mirror_type tmp_view((typename ViewType::value_type*)buffer, src.view_device().layout());

    if (src.view_device().data() && space == Device) {
      Kokkos::deep_copy(LMPHostType(),tmp_view,src.view_host()),
      Kokkos::deep_copy(LMPHostType(),src.view_device(),tmp_view);
      src.clear_sync_state();
      if (!async_flag) Kokkos::fence(); // change to less agressive fence?
    } else if (src.view_host().data()) {
      Kokkos::deep_copy(LMPHostType(),tmp_view,src.view_device()),
      Kokkos::deep_copy(LMPHostType(),src.view_host(),tmp_view);
      src.clear_sync_state();
      if (!async_flag) Kokkos::fence(); // change to less agressive fence?
    }
  }
  #else
  template<class ViewType>
  void perform_pinned_copy(ViewType& src, unsigned int space, int /*async_flag*/ = 0) {
    if (space == Device)
      src.sync_device();
    else
      src.sync_host();
  }
  #endif

  #ifdef LMP_KOKKOS_GPU
  template<class TransformViewType>
  void perform_pinned_copy_transform(TransformViewType& src, unsigned int space, int async_flag = 0) {
    typedef typename TransformViewType::kk_view ViewType;
    typedef Kokkos::View<typename ViewType::data_type,
                 typename ViewType::array_layout,
                 typename std::conditional<
                   std::is_same_v<typename ViewType::execution_space,LMPDeviceType>,
                   LMPPinnedHostType,typename ViewType::memory_space>::type,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged> > mirror_type;
    if (buffer_size == 0) {
       buffer_size = src.view_device().span()*sizeof(typename ViewType::value_type);
       buffer = Kokkos::kokkos_malloc<LMPPinnedHostType>(buffer_size);
    } else if (buffer_size < src.view_device().span()*sizeof(typename ViewType::value_type)) {
       buffer_size = src.view_device().span()*sizeof(typename ViewType::value_type);
       buffer = Kokkos::kokkos_realloc<LMPPinnedHostType>(buffer,buffer_size);
    }

    if (space == Device)
      src.sync_device(buffer,async_flag);
    else if (space == Host)
      src.sync_host(buffer,async_flag);
    else if (space == HostKK)
      src.sync_hostkk(buffer,async_flag);
  }
  #else
  template<class TransformViewType>
  void perform_pinned_copy_transform(TransformViewType& src, unsigned int space, int /*async_flag*/ = 0) {
    if (space == Device)
      src.sync_device();
    else if (space == Host)
      src.sync_host();
    else if (space == HostKK)
      src.sync_hostkk();
  }
  #endif
};

}

#endif
