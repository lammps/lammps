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

#ifdef DIHEDRAL_CLASS
// clang-format off
DihedralStyle(table/kk,DihedralTableKokkos<LMPDeviceType>);
DihedralStyle(table/kk/device,DihedralTableKokkos<LMPDeviceType>);
DihedralStyle(table/kk/host,DihedralTableKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_DIHEDRAL_TABLE_KOKKOS_H
#define LMP_DIHEDRAL_TABLE_KOKKOS_H

#include "dihedral_table.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagDihedralTableCompute{};

template<class DeviceType>
class DihedralTableKokkos : public DihedralTable {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  DihedralTableKokkos(class LAMMPS *);
  ~DihedralTableKokkos() override;
  void compute(int, int) override;
  void init_style() override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void uf_lookup_kk(const int &type, const KK_FLOAT &x, KK_FLOAT &u, KK_FLOAT &mdu) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void minimum_image(KK_FLOAT &dx, KK_FLOAT &dy, KK_FLOAT &dz) const;

  template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDihedralTableCompute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDihedralTableCompute<NEWTON_BOND,EVFLAG>, const int&) const;

  //template<int NEWTON_BOND>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int i1, const int i2, const int i3, const int i4,
                          KK_FLOAT &edihedral, KK_FLOAT *f1, KK_FLOAT *f3, KK_FLOAT *f4,
                          const KK_FLOAT &vb1x, const KK_FLOAT &vb1y, const KK_FLOAT &vb1z,
                          const KK_FLOAT &vb2x, const KK_FLOAT &vb2y, const KK_FLOAT &vb2z,
                          const KK_FLOAT &vb3x, const KK_FLOAT &vb3y, const KK_FLOAT &vb3z) const;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;

 protected:

  class NeighborKokkos *neighborKK;
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_2d_lr dihedrallist;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int nlocal,newton_bond;
  int eflag,vflag;

  // domain data for the device-side minimum_image

  int triclinic;
  int xperiodic,yperiodic,zperiodic;
  KK_FLOAT xprd,yprd,zprd;
  KK_FLOAT xprd_half,yprd_half,zprd_half;
  KK_FLOAT xy,xz,yz;

  // device copies of the tabulated data, indexed by table number.  the
  // dihedral tables are cyclic, so no extra element is needed: the kernel
  // wraps itable+1 back to 0 exactly as the CPU style does

  DAT::tdual_int_1d k_tabindex, k_f_unspecified;
  DAT::tdual_kkfloat_1d k_delta, k_invdelta, k_deltasq6;
  DAT::tdual_kkfloat_2d k_e, k_de, k_f_tab, k_df, k_e2, k_f2;

  typename AT::t_int_1d d_tabindex, d_f_unspecified;
  typename AT::t_kkfloat_1d d_delta, d_invdelta, d_deltasq6;
  typename AT::t_kkfloat_2d d_e, d_de, d_f_tab, d_df, d_e2, d_f2;

  void setup_tables();
};

}

#endif
#endif

