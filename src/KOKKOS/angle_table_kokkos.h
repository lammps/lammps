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

#ifdef ANGLE_CLASS
// clang-format off
AngleStyle(table/kk,AngleTableKokkos<LMPDeviceType>);
AngleStyle(table/kk/device,AngleTableKokkos<LMPDeviceType>);
AngleStyle(table/kk/host,AngleTableKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_ANGLE_TABLE_KOKKOS_H
#define LMP_ANGLE_TABLE_KOKKOS_H

#include "angle_table.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagAngleTableCompute{};

template<class DeviceType>
class AngleTableKokkos : public AngleTable {

 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  AngleTableKokkos(class LAMMPS *);
  ~AngleTableKokkos() override;
  void compute(int, int) override;
  void init_style() override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void uf_lookup_kk(const int &type, const KK_FLOAT &x, KK_FLOAT &u, KK_FLOAT &mdu) const;

  template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagAngleTableCompute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagAngleTableCompute<NEWTON_BOND,EVFLAG>, const int&) const;

  //template<int NEWTON_BOND>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int i, const int j, const int k,
                     KK_FLOAT &eangle, KK_FLOAT *f1, KK_FLOAT *f3,
                     const KK_FLOAT &delx1, const KK_FLOAT &dely1, const KK_FLOAT &delz1,
                     const KK_FLOAT &delx2, const KK_FLOAT &dely2, const KK_FLOAT &delz2) const;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;

 protected:

  class NeighborKokkos *neighborKK;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_2d_lr anglelist;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int nlocal,newton_bond;
  int eflag,vflag;

  // device copies of the tabulated data.  the per-table scalars and the
  // tablength-long arrays are packed into views indexed by table number,
  // so the kernel only needs tabindex[type] to reach them.  the arrays
  // carry one extra element so that the spline branch can read itable+1
  // at the last bin, where its weight is exactly zero

  DAT::tdual_int_1d k_tabindex;
  DAT::tdual_kkfloat_1d k_invdelta, k_deltasq6;
  DAT::tdual_kkfloat_2d k_ang, k_e, k_de, k_f, k_df, k_e2, k_f2;

  typename AT::t_int_1d d_tabindex;
  typename AT::t_kkfloat_1d d_invdelta, d_deltasq6;
  typename AT::t_kkfloat_2d d_ang, d_e, d_de, d_f, d_df, d_e2, d_f2;

  void setup_tables();
};

}

#endif
#endif

