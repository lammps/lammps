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

#ifndef LMP_NEIGH_BOND_KOKKOS_H
#define LMP_NEIGH_BOND_KOKKOS_H

#include "neighbor.h"
#include "kokkos_type.h"
#include "domain_kokkos.h"
#include "pointers.h"
#include <Kokkos_UnorderedMap.hpp>

namespace LAMMPS_NS {

struct TagNeighBondBondAll{};
struct TagNeighBondBondPartial{};
struct TagNeighBondBondCheck{};
struct TagNeighBondAngleAll{};
struct TagNeighBondAnglePartial{};
struct TagNeighBondAngleCheck{};
struct TagNeighBondDihedralAll{};
struct TagNeighBondDihedralPartial{};
struct TagNeighBondDihedralCheck{};
struct TagNeighBondImproperAll{};
struct TagNeighBondImproperPartial{};

template<class DeviceType>
class NeighBondKokkos : protected Pointers  {
 public:
  typedef ArrayTypes<DeviceType> AT;
  typedef int value_type;

  NeighBondKokkos(class LAMMPS *);
  ~NeighBondKokkos() override = default;
  void init_topology_kk();
  void build_topology_kk();

  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighBondBondAll, const int&, int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighBondBondPartial, const int&, int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighBondBondCheck, const int&, int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighBondAngleAll, const int&, int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighBondAnglePartial, const int&, int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighBondAngleCheck, const int&, int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighBondDihedralAll, const int&, int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighBondDihedralPartial, const int&, int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighBondDihedralCheck, const int&, int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighBondImproperAll, const int&, int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighBondImproperPartial, const int&, int&) const;

  DAT::tdual_int_2d k_bondlist;
  DAT::tdual_int_2d k_anglelist;
  DAT::tdual_int_2d k_dihedrallist;
  DAT::tdual_int_2d k_improperlist;

  // KOKKOS host/device flag and data masks
  ExecutionSpace execution_space;
  unsigned int datamask_read,datamask_modify;

  int maxbond,maxangle,maxdihedral,maximproper;   // size of bond lists
  int me,nprocs;

 private:
  int map_style;
  DAT::tdual_int_1d k_sametag;
  typename AT::t_int_1d d_sametag;
  DAT::tdual_int_1d k_map_array;
  dual_hash_type k_map_hash;

  typename AT::t_int_2d v_bondlist;
  typename AT::t_int_2d v_anglelist;
  typename AT::t_int_2d v_dihedrallist;
  typename AT::t_int_2d v_improperlist;
  typename AT::t_int_2d list;

  typename AT::t_x_array_randomread x;
  typename AT::t_tagint_1d_randomread tag;

  typename AT::t_int_1d num_bond;
  typename AT::t_int_2d bond_type;
  typename AT::t_tagint_2d bond_atom;

  typename AT::t_int_1d num_angle;
  typename AT::t_int_2d angle_type;
  typename AT::t_tagint_2d angle_atom1,angle_atom2,angle_atom3;

  typename AT::t_int_1d num_dihedral;
  typename AT::t_int_2d dihedral_type;
  typename AT::t_tagint_2d dihedral_atom1,dihedral_atom2,
    dihedral_atom3,dihedral_atom4;

  typename AT::t_int_1d num_improper;
  typename AT::t_int_2d improper_type;
  typename AT::t_tagint_2d improper_atom1,improper_atom2,
    improper_atom3,improper_atom4;

  typename AT::t_int_1d d_scalars;
  HAT::t_int_1d h_scalars;
  typename AT::t_int_scalar d_nlist;
  HAT::t_int_scalar h_nlist;
  typename AT::t_int_scalar d_fail_flag;
  HAT::t_int_scalar h_fail_flag;

  KOKKOS_INLINE_FUNCTION
  int closest_image(const int, int) const;

  KOKKOS_INLINE_FUNCTION
  void minimum_image(X_FLOAT &dx, X_FLOAT &dy, X_FLOAT &dz) const;

  void update_class_variables();

  // topology build functions

  typedef void (NeighBondKokkos::*BondPtr)();   // ptrs to topology build functions

  BondPtr bond_build_kk;                 // ptr to bond list functions
  void bond_all();                    // bond list with all bonds
  void bond_template();               // bond list with templated bonds
  void bond_partial();                // exclude certain bonds
  void bond_check();

  BondPtr angle_build_kk;                // ptr to angle list functions
  void angle_all();                   // angle list with all angles
  void angle_template();              // angle list with templated bonds
  void angle_partial();               // exclude certain angles
  void angle_check();

  BondPtr dihedral_build_kk;             // ptr to dihedral list functions
  void dihedral_all();                // dihedral list with all dihedrals
  void dihedral_template();           // dihedral list with templated bonds
  void dihedral_partial();            // exclude certain dihedrals
  void dihedral_check(int, typename AT::t_int_2d list);

  BondPtr improper_build_kk;             // ptr to improper list functions
  void improper_all();                // improper list with all impropers
  void improper_template();           // improper list with templated bonds
  void improper_partial();            // exclude certain impropers

  int nlocal,newton_bond,lostbond,nmissing;

  int triclinic;
  int xperiodic,yperiodic,zperiodic;
  X_FLOAT xprd_half,yprd_half,zprd_half;
  X_FLOAT xprd,yprd,zprd;
  X_FLOAT xy,xz,yz;
};

}

#endif

