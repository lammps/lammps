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

#ifdef ATOM_CLASS

AtomStyle(sphere/kk,AtomVecSphereKokkos)
AtomStyle(sphere/kk/device,AtomVecSphereKokkos)
AtomStyle(sphere/kk/host,AtomVecSphereKokkos)

#else

#ifndef LMP_ATOM_VEC_SPHERE_KOKKOS_H
#define LMP_ATOM_VEC_SPHERE_KOKKOS_H

#include "atom_vec_kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class AtomVecSphereKokkos : public AtomVecKokkos {
 public:
  AtomVecSphereKokkos(class LAMMPS *);
  ~AtomVecSphereKokkos() {}
  void init();
  void grow(int);
  void grow_pointers();
  void copy(int, int, int);
  int pack_comm(int, int *, double *, int, int *);
  int pack_comm_vel(int, int *, double *, int, int *);
  int pack_comm_hybrid(int, int *, double *);
  void unpack_comm(int, int, double *);
  void unpack_comm_vel(int, int, double *);
  int unpack_comm_hybrid(int, int, double *);
  int pack_reverse(int, int, double *);
  int pack_reverse_hybrid(int, int, double *);
  void unpack_reverse(int, int *, double *);
  int unpack_reverse_hybrid(int, int *, double *);
  int pack_border(int, int *, double *, int, int *);
  int pack_border_vel(int, int *, double *, int, int *);
  int pack_border_hybrid(int, int *, double *);
  void unpack_border(int, int, double *);
  void unpack_border_vel(int, int, double *);
  int unpack_border_hybrid(int, int, double *);
  int pack_exchange(int, double *);
  int unpack_exchange(double *);
  int size_restart();
  int pack_restart(int, double *);
  int unpack_restart(double *);
  void create_atom(int, double *);
  void data_atom(double *, imageint, char **);
  int data_atom_hybrid(int, char **);
  void data_vel(int, char **);
  int data_vel_hybrid(int, char **);
  void pack_data(double **);
  int pack_data_hybrid(int, double *);
  void write_data(FILE *, int, double **);
  int write_data_hybrid(FILE *, double *);
  void pack_vel(double **);
  int pack_vel_hybrid(int, double *);
  void write_vel(FILE *, int, double **);
  int write_vel_hybrid(FILE *, double *);
  bigint memory_usage();

  int pack_comm_kokkos(const int &n, const DAT::tdual_int_2d &k_sendlist,
                       const int & iswap,
                       const DAT::tdual_xfloat_2d &buf,
                       const int &pbc_flag, const int pbc[]);
  void unpack_comm_kokkos(const int &n, const int &nfirst,
                          const DAT::tdual_xfloat_2d &buf);
  int pack_comm_vel_kokkos(const int &n, const DAT::tdual_int_2d &k_sendlist,
                           const int & iswap,
                           const DAT::tdual_xfloat_2d &buf,
                           const int &pbc_flag, const int pbc[]);
  void unpack_comm_vel_kokkos(const int &n, const int &nfirst,
                              const DAT::tdual_xfloat_2d &buf);
  int pack_comm_self(const int &n, const DAT::tdual_int_2d &list,
                     const int & iswap, const int nfirst,
                     const int &pbc_flag, const int pbc[]);
  int pack_border_kokkos(int n, DAT::tdual_int_2d k_sendlist,
                         DAT::tdual_xfloat_2d buf,int iswap,
                         int pbc_flag, int *pbc, ExecutionSpace space);
  void unpack_border_kokkos(const int &n, const int &nfirst,
                            const DAT::tdual_xfloat_2d &buf,
                            ExecutionSpace space);
  int pack_border_vel_kokkos(int n, DAT::tdual_int_2d k_sendlist,
                             DAT::tdual_xfloat_2d buf,int iswap,
                             int pbc_flag, int *pbc, ExecutionSpace space);
  void unpack_border_vel_kokkos(const int &n, const int &nfirst,
                                const DAT::tdual_xfloat_2d &buf,
                                ExecutionSpace space);
  int pack_exchange_kokkos(const int &nsend,DAT::tdual_xfloat_2d &buf,
                           DAT::tdual_int_1d k_sendlist,
                           DAT::tdual_int_1d k_copylist,
                           ExecutionSpace space, int dim,
                           X_FLOAT lo, X_FLOAT hi);
  int unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf, int nrecv,
                             int nlocal, int dim, X_FLOAT lo, X_FLOAT hi,
                             ExecutionSpace space);

  void sync(ExecutionSpace space, unsigned int mask);
  void modified(ExecutionSpace space, unsigned int mask);
  void sync_overlapping_device(ExecutionSpace space, unsigned int mask);

 private:
  tagint *tag;
  int *type,*mask;
  imageint *image;
  double **x,**v,**f;
  double *radius,*rmass;
  double **omega,**torque;
  int radvary;

  DAT::t_tagint_1d d_tag;
  HAT::t_tagint_1d h_tag;
  DAT::t_imageint_1d d_image;
  HAT::t_imageint_1d h_image;
  DAT::t_int_1d d_type, d_mask;
  HAT::t_int_1d h_type, h_mask;

  DAT::t_x_array d_x;
  DAT::t_v_array d_v;
  DAT::t_f_array d_f;
  DAT::t_float_1d d_radius;
  HAT::t_float_1d h_radius;
  DAT::t_float_1d d_rmass;
  HAT::t_float_1d h_rmass;
  DAT::t_v_array d_omega;
  HAT::t_v_array h_omega;
  DAT::t_f_array d_torque;
  HAT::t_f_array h_torque;

  DAT::tdual_int_1d k_count;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Per-processor system is too big

The number of owned atoms plus ghost atoms on a single
processor must fit in 32-bit integer.

E: Invalid atom type in Atoms section of data file

Atom types must range from 1 to specified # of types.

E: Invalid radius in Atoms section of data file

Radius must be >= 0.0.

E: Invalid density in Atoms section of data file

Density value cannot be <= 0.0.

*/
