/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom.h"
#include "kokkos_type.h"

#ifndef LMP_ATOM_KOKKOS_H
#define LMP_ATOM_KOKKOS_H

namespace LAMMPS_NS {

class AtomKokkos : public Atom {
 public:
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


// USER-DPD package
  DAT::tdual_efloat_1d k_uCond, k_uMech, k_uChem, k_uCG, k_uCGnew,
                       k_rho,k_dpdTheta,k_duChem;


  AtomKokkos(class LAMMPS *);
  ~AtomKokkos();

  virtual void allocate_type_arrays();
  void sync(const ExecutionSpace space, unsigned int mask);
  void modified(const ExecutionSpace space, unsigned int mask);
  void sync_overlapping_device(const ExecutionSpace space, unsigned int mask);
  virtual void sort();
  virtual void grow(unsigned int mask);
  int add_custom(const char *, int);
  void remove_custom(int, int);
  virtual void deallocate_topology();
  void sync_modify(ExecutionSpace, unsigned int, unsigned int);
 private:
   class AtomVec *new_avec(const char *, int, int &);
};

template<class ViewType, class IndexView>
class SortFunctor {
  typedef typename ViewType::device_type device_type;
  ViewType source;
  Kokkos::View<typename ViewType::non_const_data_type,typename ViewType::array_type,device_type> dest;
  IndexView index;
  SortFunctor(ViewType src, typename Kokkos::Impl::enable_if<ViewType::dynamic_rank==1,IndexView>::type ind):source(src),index(ind){
    dest = Kokkos::View<typename ViewType::non_const_data_type,typename ViewType::array_type,device_type>("",src.extent(0));
  }
  SortFunctor(ViewType src, typename Kokkos::Impl::enable_if<ViewType::dynamic_rank==2,IndexView>::type ind):source(src),index(ind){
    dest = Kokkos::View<typename ViewType::non_const_data_type,typename ViewType::array_type,device_type>("",src.extent(0),src.extent(1));
  }
  SortFunctor(ViewType src, typename Kokkos::Impl::enable_if<ViewType::dynamic_rank==3,IndexView>::type ind):source(src),index(ind){
    dest = Kokkos::View<typename ViewType::non_const_data_type,typename ViewType::array_type,device_type>("",src.extent(0),src.extent(1),src.extent(2));
  }
  SortFunctor(ViewType src, typename Kokkos::Impl::enable_if<ViewType::dynamic_rank==4,IndexView>::type ind):source(src),index(ind){
    dest = Kokkos::View<typename ViewType::non_const_data_type,typename ViewType::array_type,device_type>("",src.extent(0),src.extent(1),src.extent(2),src.extent(3));
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::Impl::enable_if<ViewType::rank==1, int>::type& i) {
    dest(i) = source(index(i));
  }
  void operator()(const typename Kokkos::Impl::enable_if<ViewType::rank==2, int>::type& i) {
    for(int j=0;j<source.extent(1);j++)
      dest(i,j) = source(index(i),j);
  }
  void operator()(const typename Kokkos::Impl::enable_if<ViewType::rank==3, int>::type& i) {
    for(int j=0;j<source.extent(1);j++)
    for(int k=0;k<source.extent(2);k++)
      dest(i,j,k) = source(index(i),j,k);
  }
  void operator()(const typename Kokkos::Impl::enable_if<ViewType::rank==4, int>::type& i) {
    for(int j=0;j<source.extent(1);j++)
    for(int k=0;k<source.extent(2);k++)
    for(int l=0;l<source.extent(3);l++)
      dest(i,j,k,l) = source(index(i),j,k,l);
  }
};

}

#endif

/* ERROR/WARNING messages:

E: KOKKOS package requires a kokkos enabled atom_style

Self-explanatory.

*/
