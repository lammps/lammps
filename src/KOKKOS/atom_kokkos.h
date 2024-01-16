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

#include "atom.h"               // IWYU pragma: export
#include "kokkos_type.h"

#include <Kokkos_Sort.hpp>

#ifndef LMP_ATOM_KOKKOS_H
#define LMP_ATOM_KOKKOS_H

namespace LAMMPS_NS {

class AtomKokkos : public Atom {
 public:
  bool sort_classic;
  int nprop_atom;
  class FixPropertyAtomKokkos **fix_prop_atom;

  DAT::tdual_tagint_1d k_tag;
  DAT::tdual_int_1d k_type, k_mask;
  DAT::tdual_imageint_1d k_image;
  DAT::tdual_x_array k_x;
  DAT::tdual_v_array k_v;
  DAT::tdual_f_array k_f;

  DAT::tdual_float_1d k_mass;

  DAT::tdual_float_1d k_q;
  DAT::tdual_float_1d k_radius;
  DAT::tdual_float_1d k_rmass;
  DAT::tdual_float_1d_4 k_mu;
  DAT::tdual_v_array k_omega;
  DAT::tdual_v_array k_angmom;
  DAT::tdual_f_array k_torque;
  DAT::tdual_tagint_1d k_molecule;
  DAT::tdual_int_2d k_nspecial;
  DAT::tdual_tagint_2d k_special;
  DAT::tdual_int_1d k_num_bond;
  DAT::tdual_int_2d k_bond_type;
  DAT::tdual_tagint_2d k_bond_atom;
  DAT::tdual_int_1d k_num_angle;
  DAT::tdual_int_2d k_angle_type;
  DAT::tdual_tagint_2d k_angle_atom1, k_angle_atom2, k_angle_atom3;
  DAT::tdual_int_1d k_num_dihedral;
  DAT::tdual_int_2d k_dihedral_type;
  DAT::tdual_tagint_2d k_dihedral_atom1, k_dihedral_atom2, k_dihedral_atom3, k_dihedral_atom4;
  DAT::tdual_int_1d k_num_improper;
  DAT::tdual_int_2d k_improper_type;
  DAT::tdual_tagint_2d k_improper_atom1, k_improper_atom2, k_improper_atom3, k_improper_atom4;

  DAT::tdual_float_2d k_dvector;

  // SPIN package

  DAT::tdual_float_1d_4 k_sp;
  DAT::tdual_f_array k_fm;
  DAT::tdual_f_array k_fm_long;

// DPD-REACT package
  DAT::tdual_efloat_1d k_uCond, k_uMech, k_uChem, k_uCG, k_uCGnew,
                       k_rho,k_dpdTheta,k_duChem;


  AtomKokkos(class LAMMPS *);
  ~AtomKokkos() override;

  void map_init(int check = 1) override;
  void map_clear() override;
  void map_set() override;
  void map_one(tagint, int) override;
  void map_delete() override;
  int map_find_hash(tagint) override;

  DAT::tdual_int_scalar k_error_flag;
  DAT::tdual_int_1d k_sametag;
  DAT::tdual_int_1d k_map_array;
  dual_hash_type k_map_hash;

  DAT::t_tagint_1d d_tag_sorted;
  DAT::t_int_1d d_i_sorted;

  typedef Kokkos::DualView<tagint[2], LMPDeviceType::array_layout, LMPDeviceType> tdual_tagint_2;
  typedef tdual_tagint_2::t_dev t_tagint_2;
  typedef tdual_tagint_2::t_host t_host_tagint_2;

  t_tagint_2 d_tag_min_max;
  t_host_tagint_2 h_tag_min_max;

  DAT::t_tagint_scalar d_tag_min,d_tag_max;
  HAT::t_tagint_scalar h_tag_min,h_tag_max;

  using MapKeyViewType = decltype(d_tag_sorted);
  using BinOpMap = Kokkos::BinOp1D<MapKeyViewType>;
  BinOpMap mapBinner;
  Kokkos::BinSort<MapKeyViewType, BinOpMap> mapSorter;

  class AtomVecKokkos* avecKK;

  // map lookup function inlined for efficiency
  // return -1 if no map defined

  inline int map(tagint global) override
  {
    if (map_style == 1) {
      k_map_array.sync_host();
      return map_array[global];
    } else if (map_style == 2)
      return map_find_hash(global);
    else
      return -1;
  };

  template<class DeviceType>
  KOKKOS_INLINE_FUNCTION
  static int map_kokkos(tagint global, int map_style, const DAT::tdual_int_1d &k_map_array, const dual_hash_type &k_map_hash)
  {
    if (map_style == 1)
      return k_map_array.view<DeviceType>()(global);
    else if (map_style == 2)
      return AtomKokkos::map_find_hash_kokkos<DeviceType>(global,k_map_hash);
    else
      return -1;
  }

  template<class DeviceType>
  KOKKOS_INLINE_FUNCTION
  static int map_find_hash_kokkos(tagint global, const dual_hash_type &k_map_hash)
  {
    int local = -1;
    auto& d_map_hash = k_map_hash.const_view<DeviceType>();
    auto index = d_map_hash.find(global);
    if (d_map_hash.valid_at(index))
      local = d_map_hash.value_at(index);
    return local;
  }

  void init() override;
  void update_property_atom();
  void allocate_type_arrays() override;
  void sync(const ExecutionSpace space, unsigned int mask);
  void modified(const ExecutionSpace space, unsigned int mask);
  void sync_overlapping_device(const ExecutionSpace space, unsigned int mask);
  void sort() override;
  virtual void grow(unsigned int mask);
  int add_custom(const char *, int, int, int border = 0) override;
  void remove_custom(int, int, int) override;
  virtual void deallocate_topology();
 private:
  void sort_device();
  class AtomVec *new_avec(const std::string &, int, int &) override;
};

template<class ViewType, class IndexView>
struct SortFunctor {
  typedef typename ViewType::device_type device_type;
  ViewType source;
  Kokkos::View<typename ViewType::non_const_data_type,typename ViewType::array_type,device_type> dest;
  IndexView index;
  SortFunctor(ViewType src, std::enable_if_t<ViewType::dynamic_rank==1,IndexView> ind):source(src),index(ind) {
    dest = Kokkos::View<typename ViewType::non_const_data_type,typename ViewType::array_type,device_type>("",src.extent(0));
  }
  SortFunctor(ViewType src, std::enable_if_t<ViewType::dynamic_rank==2,IndexView> ind):source(src),index(ind) {
    dest = Kokkos::View<typename ViewType::non_const_data_type,typename ViewType::array_type,device_type>("",src.extent(0),src.extent(1));
  }
  SortFunctor(ViewType src, std::enable_if_t<ViewType::dynamic_rank==3,IndexView> ind):source(src),index(ind) {
    dest = Kokkos::View<typename ViewType::non_const_data_type,typename ViewType::array_type,device_type>("",src.extent(0),src.extent(1),src.extent(2));
  }
  SortFunctor(ViewType src, std::enable_if_t<ViewType::dynamic_rank==4,IndexView> ind):source(src),index(ind) {
    dest = Kokkos::View<typename ViewType::non_const_data_type,typename ViewType::array_type,device_type>("",src.extent(0),src.extent(1),src.extent(2),src.extent(3));
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(const std::enable_if_t<ViewType::rank==1, int>& i) {
    dest(i) = source(index(i));
  }
  void operator()(const std::enable_if_t<ViewType::rank==2, int>& i) {
    for (int j=0; j < (int)source.extent(1); j++)
      dest(i,j) = source(index(i),j);
  }
  void operator()(const std::enable_if_t<ViewType::rank==3, int>& i) {
    for (int j=0; j < (int)source.extent(1); j++)
      for (int k=0; k < (int)source.extent(2); k++)
        dest(i,j,k) = source(index(i),j,k);
  }
  void operator()(const std::enable_if_t<ViewType::rank==4, int>& i) {
    for (int j=0; j < (int)source.extent(1); j++)
      for (int k=0; k < (int)source.extent(2); k++)
        for (int l=0; l < (int)source.extent(3); l++)
          dest(i,j,k,l) = source(index(i),j,k,l);
  }
};

}

#endif

