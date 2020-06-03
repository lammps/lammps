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

#ifndef LMP_NEIGHBOR_KOKKOS_H
#define LMP_NEIGHBOR_KOKKOS_H

#include "neighbor.h"
#include "neigh_list_kokkos.h"
#include "neigh_bond_kokkos.h"
#include "kokkos_type.h"
#include <cmath>

namespace LAMMPS_NS {

template<ExecutionSpace Space>
struct TagNeighborCheckDistance{};

template<ExecutionSpace Space>
struct TagNeighborXhold{};

class NeighborKokkos : public Neighbor {
 public:
  typedef int value_type;

  NeighborKokkos(class LAMMPS *);
  ~NeighborKokkos();
  void init();
  void init_topology();
  void build_topology();

  template<ExecutionSpace Space>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighborCheckDistance<Space>, const int&, int&) const;

  template<ExecutionSpace Space>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighborXhold<Space>, const int&) const;

  DAT::tdual_float_2d k_cutneighsq;

  DAT::tdual_int_1d k_ex1_type,k_ex2_type;
  DAT::tdual_int_2d k_ex_type;
  DAT::tdual_int_1d k_ex1_group,k_ex2_group;
  DAT::tdual_int_1d k_ex1_bit,k_ex2_bit;
  DAT::tdual_int_1d k_ex_mol_group;
  DAT::tdual_int_1d k_ex_mol_bit;
  DAT::tdual_int_1d k_ex_mol_intra;

  NeighBondKokkos<Host> neighbond_host;
  NeighBondKokkos<Device> neighbond_device;

  DAT::tdual_int_2d k_bondlist;
  DAT::tdual_int_2d k_anglelist;
  DAT::tdual_int_2d k_dihedrallist;
  DAT::tdual_int_2d k_improperlist;

 private:

  DAT::tdual_float_1d_3 x;
  DAT::tdual_float_1d_3 xhold;

  KK_FLOAT deltasq;
  int device_flag;

  void init_cutneighsq_kokkos(int);
  void create_kokkos_list(int);
  void init_ex_type_kokkos(int);
  void init_ex_bit_kokkos();
  void init_ex_mol_bit_kokkos();
  void grow_ex_mol_intra_kokkos();
  virtual int check_distance();
  template<ExecutionSpace Space> int check_distance_kokkos();
  virtual void build(int);
  template<ExecutionSpace Space> void build_kokkos(int);
  void setup_bins_kokkos(int);
  void modify_ex_type_grow_kokkos();
  void modify_ex_group_grow_kokkos();
  void modify_mol_group_grow_kokkos();
  void modify_mol_intra_grow_kokkos();
  void set_binsize_kokkos();
};

}

#endif

/* ERROR/WARNING messages:

E: KOKKOS package only supports 'bin' neighbor lists

Self-explanatory.

E: Too many local+ghost atoms for neighbor list

The number of nlocal + nghost atoms on a processor
is limited by the size of a 32-bit integer with 2 bits
removed for masking 1-2, 1-3, 1-4 neighbors.

*/
