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

/* ----------------------------------------------------------------------
   Contributing author: Matt Bettencourt (NVIDIA)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(mliap/kk,PairMLIAPKokkos<LMPDeviceType>);
PairStyle(mliap/kk/device,PairMLIAPKokkos<LMPDeviceType>);
PairStyle(mliap/kk/host,PairMLIAPKokkos<LMPHostType>);
// clang-format off
#else

#ifndef LMP_PAIR_MLIAP_KOKKOS_H
#define LMP_PAIR_MLIAP_KOKKOS_H

#include "pair_mliap.h"
#include "pair_kokkos.h"
#include "kokkos_type.h"
#include "kokkos_base.h"
#include "comm.h"
#include <variant>

namespace LAMMPS_NS {

template<class DeviceType>
class PairMLIAPKokkos : public PairMLIAP, public KokkosBase  {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  PairMLIAPKokkos(class LAMMPS*);
  ~PairMLIAPKokkos() override;
  void settings(int narg, char ** arg) override;
  void init_style() override;

  void compute(int, int) override;
  void e_tally(MLIAPData* data);

  void allocate() override;

  void coeff(int narg, char **arg) override;

  //Outward facing functions to be invoked by the ML layer via MLIAPDataKokkosDevice
  template <typename CommType>
  int forward_comm(CommType* copy_from, CommType* copy_to, const int vec_len);
  template <typename CommType>
  int reverse_comm(CommType* copy_from, CommType* copy_to, const int vec_len);


  //Virtual functions overwriting the kokkos virtual functions in Pair. The
  //virtual functions call the templated versions using the work-around cache
  //variables below
  int pack_forward_comm_kokkos(int nv, DAT::tdual_int_1d iii, DAT::tdual_xfloat_1d &fill,
                               int int2, int *intp) override;
  template <typename CommType>
  int pack_forward_comm_kokkos(int nv, DAT::tdual_int_1d iii, DAT::tdual_xfloat_1d &fill,
                               int int2, int *intp, CommType *copy_to);
  int pack_forward_comm(int, int *, double *, int, int *) override;
  template <typename CommType>
  int pack_forward_comm(int, int *, double *, int, int *, CommType *copy_to) ;

  void unpack_forward_comm_kokkos(int nv, int first_up, DAT::tdual_xfloat_1d& fill) override;
  template <typename CommType>
  void unpack_forward_comm_kokkos(int nv, int first_up, DAT::tdual_xfloat_1d& fill, CommType *copy_to);

  void unpack_forward_comm(int, int, double *) override;
  template <typename CommType>
  void unpack_forward_comm(int, int, double *, CommType *copy_to);

  int pack_reverse_comm_kokkos(int, int, DAT::tdual_xfloat_1d&) override;
  template <typename CommType>
  int pack_reverse_comm_kokkos(int, int, DAT::tdual_xfloat_1d&, CommType *copy_to);

  int pack_reverse_comm(int, int, double *) override;
  template <typename CommType>
  int pack_reverse_comm(int, int, double *,CommType *copy_to);

  void unpack_reverse_comm_kokkos(int, DAT::tdual_int_1d idx,DAT::tdual_xfloat_1d& buf) override;
  template <typename CommType>
  void unpack_reverse_comm_kokkos(int, DAT::tdual_int_1d idx,DAT::tdual_xfloat_1d& buf, CommType *copy_to);
  void unpack_reverse_comm(int, int *, double *)override;
  template <typename CommType>
  void unpack_reverse_comm(int, int *, double *,CommType *copy_to);


  // Variables to store to/from pointers and types between
  // forward_comm/reverse_comm and virutal functions. Ugly work-around which
  // could be partially be removed by templating LAMMPS comms on float type
  std::variant<float*, double*> copy_to, copy_from;
  enum class COMM_TYPE {FLOAT=0, DOUBLE, UNSET} comm_type;
  int vec_len;

  typename AT::t_x_array_randomread x;
  typename AT::t_x_array_randomread v;
  typename AT::t_f_array f;
  DAT::tdual_int_1d k_map;
  DAT::tdual_double_2d k_cutsq;
  DAT::tdual_int_2d k_setflag;
  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_double_2d k_vatom;


  friend void pair_virial_fdotr_compute<PairMLIAPKokkos>(PairMLIAPKokkos*);
};

}
#endif
#endif
