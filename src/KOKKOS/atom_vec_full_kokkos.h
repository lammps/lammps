/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(full/kk,AtomVecFullKokkos);
AtomStyle(full/kk/device,AtomVecFullKokkos);
AtomStyle(full/kk/host,AtomVecFullKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_ATOM_VEC_FULL_KOKKOS_H
#define LMP_ATOM_VEC_FULL_KOKKOS_H

#include "atom_vec_kokkos.h"

namespace LAMMPS_NS {

class AtomVecFullKokkos : public AtomVecKokkos {
 public:
  AtomVecFullKokkos(class LAMMPS *);

  void grow(int) override;
  void copy(int, int, int) override;
  int pack_border(int, int *, double *, int, int *) override;
  int pack_border_vel(int, int *, double *, int, int *) override;
  int pack_border_hybrid(int, int *, double *) override;
  void unpack_border(int, int, double *) override;
  void unpack_border_vel(int, int, double *) override;
  int unpack_border_hybrid(int, int, double *) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(double *) override;
  int size_restart() override;
  int pack_restart(int, double *) override;
  int unpack_restart(double *) override;
  void create_atom(int, double *) override;
  void data_atom(double *, imageint, const std::vector<std::string> &, std::string &) override;
  int data_atom_hybrid(int, const std::vector<std::string> &, int) override;
  void pack_data(double **) override;
  int pack_data_hybrid(int, double *) override;
  void write_data(FILE *, int, double **) override;
  int write_data_hybrid(FILE *, double *) override;
  double memory_usage() override;

  void grow_pointers() override;
  int pack_border_kokkos(int n, DAT::tdual_int_2d k_sendlist,
                         DAT::tdual_xfloat_2d buf,int iswap,
                         int pbc_flag, int *pbc, ExecutionSpace space) override;
  void unpack_border_kokkos(const int &n, const int &nfirst,
                            const DAT::tdual_xfloat_2d &buf,
                            ExecutionSpace space) override;
  int pack_exchange_kokkos(const int &nsend,DAT::tdual_xfloat_2d &buf,
                           DAT::tdual_int_1d k_sendlist,
                           DAT::tdual_int_1d k_copylist,
                           ExecutionSpace space, int dim,
                           X_FLOAT lo, X_FLOAT hi) override;
  int unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf, int nrecv,
                             int nlocal, int dim, X_FLOAT lo, X_FLOAT hi,
                             ExecutionSpace space) override;

  void sync(ExecutionSpace space, unsigned int mask) override;
  void modified(ExecutionSpace space, unsigned int mask) override;
  void sync_overlapping_device(ExecutionSpace space, unsigned int mask) override;

 protected:

  tagint *tag;
  int *type,*mask;
  imageint *image;
  double **x,**v,**f;

  double *q;

  tagint *molecule;
  int **nspecial;
  tagint **special;
  int *num_bond;
  int **bond_type;
  tagint **bond_atom;

  int *num_angle;
  int **angle_type;
  tagint **angle_atom1,**angle_atom2,**angle_atom3;

  int *num_dihedral;
  int **dihedral_type;
  tagint **dihedral_atom1,**dihedral_atom2,**dihedral_atom3,**dihedral_atom4;
  int *num_improper;
  int **improper_type;
  tagint **improper_atom1,**improper_atom2,**improper_atom3,**improper_atom4;

  DAT::t_tagint_1d d_tag;
  DAT::t_int_1d d_type, d_mask;
  HAT::t_tagint_1d h_tag;
  HAT::t_int_1d h_type, h_mask;

  DAT::t_imageint_1d d_image;
  HAT::t_imageint_1d h_image;

  DAT::t_x_array d_x;
  DAT::t_v_array d_v;
  DAT::t_f_array d_f;

  DAT::t_float_1d d_q;
  HAT::t_float_1d h_q;

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

  HAT::tdual_int_1d k_count;

};

}

#endif
#endif

