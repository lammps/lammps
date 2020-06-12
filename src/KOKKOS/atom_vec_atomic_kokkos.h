/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale AtomicKokkos/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS

AtomStyle(atomic/kk,AtomVecAtomicKokkos)
AtomStyle(atomic/kk/device,AtomVecAtomicKokkos)
AtomStyle(atomic/kk/host,AtomVecAtomicKokkos)

#else

#ifndef LMP_ATOM_VEC_ATOMIC_KOKKOS_H
#define LMP_ATOM_VEC_ATOMIC_KOKKOS_H

#include "atom_vec_kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class AtomVecAtomicKokkos : public AtomVecKokkos {
 public:
  AtomVecAtomicKokkos(class LAMMPS *);
  virtual ~AtomVecAtomicKokkos() {}
  void grow(int);
  void copy(int, int, int);
  int pack_border(int, int *, double *, int, int *);
  int pack_border_vel(int, int *, double *, int, int *);
  void unpack_border(int, int, double *);
  void unpack_border_vel(int, int, double *);
  int pack_exchange(int, double *);
  int unpack_exchange(double *);
  int size_restart();
  int pack_restart(int, double *);
  int unpack_restart(double *);
  void create_atom(int, double *);
  void data_atom(double *, tagint, char **);
  void pack_data(double **);
  void write_data(FILE *, int, double **);
  bigint memory_usage();

  void grow_pointers();
  int pack_border_kokkos(int n, DAT::tdual_int_2d k_sendlist,
                         DAT::tdual_double_2d_lr buf,int iswap,
                         int pbc_flag, int *pbc, ExecutionSpace space);
  void unpack_border_kokkos(const int &n, const int &nfirst,
                            const DAT::tdual_double_2d_lr &buf,
                            ExecutionSpace space);
  int pack_exchange_kokkos(const int &nsend,DAT::tdual_double_2d_lr &buf,
                           DAT::tdual_int_1d k_sendlist,
                           DAT::tdual_int_1d k_copylist,
                           ExecutionSpace space, int dim,
                           KK_FLOAT lo, KK_FLOAT hi);
  int unpack_exchange_kokkos(DAT::tdual_double_2d_lr &k_buf, int nrecv,
                             int nlocal, int dim, KK_FLOAT lo, KK_FLOAT hi,
                             ExecutionSpace space);

  void sync(ExecutionSpace space, unsigned int mask);
  void modified(ExecutionSpace space, unsigned int mask);
  void sync_overlapping_device(ExecutionSpace space, unsigned int mask);
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

*/
