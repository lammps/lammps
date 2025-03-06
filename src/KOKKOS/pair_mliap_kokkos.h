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

namespace LAMMPS_NS {

template<class DeviceType>
class PairMLIAPKokkos : public PairMLIAP {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  PairMLIAPKokkos(class LAMMPS*);
  ~PairMLIAPKokkos();
  void settings(int narg, char ** arg);
  void init_style();

  void compute(int, int);
  void e_tally(MLIAPData* data);

  void allocate();

  void coeff(int narg, char **arg);
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
