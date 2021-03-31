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

#ifdef COMPUTE_CLASS

ComputeStyle(orientorder/atom/kk,ComputeOrientOrderAtomKokkos<LMPDeviceType>)
ComputeStyle(orientorder/atom/kk/device,ComputeOrientOrderAtomKokkos<LMPDeviceType>)
ComputeStyle(orientorder/atom/kk/host,ComputeOrientOrderAtomKokkos<LMPHostType>)

#else

#ifndef LMP_COMPUTE_ORIENTORDER_ATOM_KOKKOS_H
#define LMP_COMPUTE_ORIENTORDER_ATOM_KOKKOS_H

#include "compute_orientorder_atom.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagComputeOrientOrderAtomNeigh{};
struct TagComputeOrientOrderAtomSelect3{};
struct TagComputeOrientOrderAtomBOOP1{};
struct TagComputeOrientOrderAtomBOOP2{};

template<class DeviceType>
class ComputeOrientOrderAtomKokkos : public ComputeOrientOrderAtom {
 public:
  typedef Kokkos::View<int*, DeviceType> t_sna_1i;
  typedef Kokkos::View<double*, DeviceType> t_sna_1d;
  typedef Kokkos::View<double*, DeviceType, Kokkos::MemoryTraits<Kokkos::Atomic> > t_sna_1d_atomic;
  typedef Kokkos::View<int**, Kokkos::LayoutRight, DeviceType> t_sna_2i_lr;
  typedef Kokkos::View<int**, Kokkos::LayoutRight, DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> > t_sna_2i_lr_um;
  typedef Kokkos::View<int**, DeviceType> t_sna_2i;
  typedef Kokkos::View<double**, DeviceType> t_sna_2d;
  typedef Kokkos::View<double**, Kokkos::LayoutRight, DeviceType> t_sna_2d_lr;
  typedef Kokkos::DualView<double**, Kokkos::LayoutRight, DeviceType> tdual_sna_2d_lr;
  typedef Kokkos::View<double**, Kokkos::LayoutRight, DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> > t_sna_2d_lr_um;
  typedef Kokkos::View<double***, DeviceType> t_sna_3d;
  typedef Kokkos::View<double***, Kokkos::LayoutRight, DeviceType> t_sna_3d_lr;
  typedef Kokkos::View<double***, Kokkos::LayoutRight, DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> > t_sna_3d_lr_um;
  typedef Kokkos::View<double***[3], DeviceType> t_sna_4d;
  typedef Kokkos::View<double**[3], DeviceType> t_sna_3d3;
  typedef Kokkos::View<double*****, DeviceType> t_sna_5d;

  typedef Kokkos::View<SNAcomplex*, DeviceType> t_sna_1c;
  typedef Kokkos::View<SNAcomplex*, DeviceType, Kokkos::MemoryTraits<Kokkos::Atomic> > t_sna_1c_atomic;
  typedef Kokkos::View<SNAcomplex**, DeviceType> t_sna_2c;
  typedef Kokkos::View<SNAcomplex**, Kokkos::LayoutRight, DeviceType> t_sna_2c_lr;
  typedef Kokkos::View<SNAcomplex***, DeviceType> t_sna_3c;
  typedef Kokkos::View<SNAcomplex***[3], DeviceType> t_sna_4c;
  typedef Kokkos::View<SNAcomplex**[3], DeviceType> t_sna_3c3;
  typedef Kokkos::View<SNAcomplex*****, DeviceType> t_sna_5c;

  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef int value_type;

  ComputeOrientOrderAtomKokkos(class LAMMPS *, int, char **);
  ~ComputeOrientOrderAtomKokkos();
  void init();
  void compute_peratom();
  t_sna_1i d_qlist;

  template<class TagStyle>
  void check_team_size_for(int, int&, int);

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeOrientOrderAtomNeigh, const typename Kokkos::TeamPolicy<DeviceType, TagComputeOrientOrderAtomNeigh>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeOrientOrderAtomSelect3, const int& ii) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeOrientOrderAtomBOOP1, const typename Kokkos::TeamPolicy<DeviceType, TagComputeOrientOrderAtomBOOP1>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeOrientOrderAtomBOOP2, const int& ii) const;

  DAT::tdual_float_2d k_qnarray;
  typename AT::t_float_2d d_qnarray;

 private:
  int inum,chunk_size,chunk_offset;
  int host_flag;

  typename AT::t_x_array_randomread x;
  typename ArrayTypes<DeviceType>::t_int_1d mask;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  t_sna_1i d_ncount;
  t_sna_2d_lr d_distsq;
  t_sna_2i_lr d_nearest;
  t_sna_3d_lr d_rlist;

  t_sna_2d_lr_um d_distsq_um;
  t_sna_2i_lr_um d_nearest_um;
  t_sna_3d_lr_um d_rlist_um;

  t_sna_3c d_qnm;

  KOKKOS_INLINE_FUNCTION
  void select3(int, int, int) const;

  KOKKOS_INLINE_FUNCTION
  void calc_boop1(int, int, int) const;

  KOKKOS_INLINE_FUNCTION
  void calc_boop2(int, int) const;

  KOKKOS_INLINE_FUNCTION
  double polar_prefactor(int, int, double) const;

  KOKKOS_INLINE_FUNCTION
  double associated_legendre(int, int, double) const;

  void init_clebsch_gordan();
  t_sna_1d d_cglist;                     // Clebsch-Gordan coeffs
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
