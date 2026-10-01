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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(hybrid/scaled/kk,PairHybridScaledKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_HYBRID_SCALED_KOKKOS_H
#define LMP_PAIR_HYBRID_SCALED_KOKKOS_H

#include "pair_hybrid_kokkos.h"

#include "kokkos_base.h"

namespace LAMMPS_NS {

class PairHybridScaledKokkos : public PairHybridKokkos, public KokkosBase {
 public:
  PairHybridScaledKokkos(class LAMMPS *);
  ~PairHybridScaledKokkos() override;
  void compute(int, int) override;
  void init_style() override;
  void settings(int, char **) override;
  void coeff(int, char **) override;

  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void born_matrix(int, int, int, int, double, double, double, double &, double &) override;

  void init_svector() override;
  void copy_svector(int, int) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  // nvcc's extended __host__ __device__ lambda extension requires the
  // enclosing function to have public access

  int pack_forward_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_double_1d &, int, int *) override;
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_double_1d &) override;

  template <class DeviceType> void save_forces(int);
  template <class DeviceType> void clear_forces(int);
  template <class DeviceType> void accumulate_forces(int, double, int);
  template <class DeviceType> void restore_forces(int);

protected:
  double *scaleval;
  int *scaleidx;
  std::vector<std::string> scalevars;
  int nmaxfsum;
  int nmaxatomscale;
  int *atomvar;         // indices of atom-style variables
  double *atomscale;    // scratch vector for evaluating atom-style variables

  void update_atomscale(int);

  // accumulators for the scaled sum over the sub-style forces and torques
  // and the per-atom scale factors, all in the accumulation memory space

  DAT::tdual_kkacc_1d_3 k_fsum, k_tsum;
  DAT::tdual_kkacc_1d k_atomscale;

  // memory space the scaled sum is carried out in, see init_style()

  ExecutionSpace accum_space;

  DAT::t_kkfloat_1d_3_lr_randomread x;
  DAT::t_kkacc_1d_3 f;
  friend void pair_virial_fdotr_compute<PairHybridScaledKokkos>(PairHybridScaledKokkos*);
};

}    // namespace LAMMPS_NS

#endif
#endif


