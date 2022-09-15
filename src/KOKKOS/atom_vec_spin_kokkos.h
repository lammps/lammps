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
AtomStyle(spin/kk,AtomVecSpinKokkos);
AtomStyle(spin/kk/device,AtomVecSpinKokkos);
AtomStyle(spin/kk/host,AtomVecSpinKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_ATOM_VEC_SPIN_KOKKOS_H
#define LMP_ATOM_VEC_SPIN_KOKKOS_H

#include "atom_vec_kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class AtomVecSpinKokkos : public AtomVecKokkos {
 public:
  AtomVecSpinKokkos(class LAMMPS *);
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

  // clear magnetic and mechanic forces

  void force_clear(int, size_t) override;

  void grow_pointers() override;
  // input lists to be checked
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
  double **x,**v,**f;           // lattice quantities

                                // spin quantities
  double **sp;                  // sp[i][0-2] direction of the spin i
                                // sp[i][3] atomic magnetic moment of the spin i
  double **fm;                  // fm[i][0-2] direction of magnetic precession
  double **fm_long;             // storage of long-range spin prec. components

  DAT::t_tagint_1d d_tag;
  HAT::t_tagint_1d h_tag;

  DAT::t_int_1d d_type, d_mask;
  HAT::t_int_1d h_type, h_mask;

  DAT::t_imageint_1d d_image;
  HAT::t_imageint_1d h_image;

  DAT::t_x_array d_x;
  DAT::t_v_array d_v;
  DAT::t_f_array d_f;

  DAT::t_sp_array d_sp;
  DAT::t_fm_array d_fm;
  DAT::t_fm_long_array d_fm_long;

  HAT::t_sp_array h_sp;
  HAT::t_fm_array h_fm;
  HAT::t_fm_long_array h_fm_long;

  DAT::tdual_int_1d k_count;
};

}

#endif
#endif

