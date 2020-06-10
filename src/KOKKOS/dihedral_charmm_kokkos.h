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

#ifdef DIHEDRAL_CLASS

DihedralStyle(charmm/kk,DihedralCharmmKokkos<Device>)
DihedralStyle(charmm/kk/device,DihedralCharmmKokkos<Device>)
DihedralStyle(charmm/kk/host,DihedralCharmmKokkos<Host>)

#else

#ifndef LMP_DIHEDRAL_CHARMM_KOKKOS_H
#define LMP_DIHEDRAL_CHARMM_KOKKOS_H

#include "dihedral_charmm.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct s_EVM_FLOAT {
  KK_FLOAT evdwl;
  KK_FLOAT ecoul;
  KK_FLOAT emol;
  KK_FLOAT v[6];
  KK_FLOAT vp[6];
  KOKKOS_INLINE_FUNCTION
  s_EVM_FLOAT() {
          evdwl = 0;
          ecoul = 0;
          emol = 0;
          v[0] = 0; v[1] = 0; v[2] = 0;
          v[3] = 0; v[4] = 0; v[5] = 0;
          vp[0] = 0; vp[1] = 0; vp[2] = 0;
          vp[3] = 0; vp[4] = 0; vp[5] = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const s_EVM_FLOAT &rhs) {
    evdwl += rhs.evdwl;
    ecoul += rhs.ecoul;
    emol += rhs.emol;
    v[0] += rhs.v[0];
    v[1] += rhs.v[1];
    v[2] += rhs.v[2];
    v[3] += rhs.v[3];
    v[4] += rhs.v[4];
    v[5] += rhs.v[5];
    vp[0] += rhs.vp[0];
    vp[1] += rhs.vp[1];
    vp[2] += rhs.vp[2];
    vp[3] += rhs.vp[3];
    vp[4] += rhs.vp[4];
    vp[5] += rhs.vp[5];
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const volatile s_EVM_FLOAT &rhs) volatile {
    evdwl += rhs.evdwl;
    ecoul += rhs.ecoul;
    emol += rhs.emol;
    v[0] += rhs.v[0];
    v[1] += rhs.v[1];
    v[2] += rhs.v[2];
    v[3] += rhs.v[3];
    v[4] += rhs.v[4];
    v[5] += rhs.v[5];
    vp[0] += rhs.vp[0];
    vp[1] += rhs.vp[1];
    vp[2] += rhs.vp[2];
    vp[3] += rhs.vp[3];
    vp[4] += rhs.vp[4];
    vp[5] += rhs.vp[5];
  }
};
typedef struct s_EVM_FLOAT EVM_FLOAT;

template<int NEWTON_BOND, int EVFLAG>
struct TagDihedralCharmmCompute{};

template<ExecutionSpace Space>
class DihedralCharmmKokkos : public DihedralCharmm {
 public:
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef EVM_FLOAT value_type;
  typedef ArrayTypes<Space> AT;

  DihedralCharmmKokkos(class LAMMPS *);
  virtual ~DihedralCharmmKokkos();
  void compute(int, int);
  void coeff(int, char **);
  void init_style();
  void read_restart(FILE *);

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDihedralCharmmCompute<NEWTON_BOND,EVFLAG>, const int&, EVM_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDihedralCharmmCompute<NEWTON_BOND,EVFLAG>, const int&) const;

  //template<int NEWTON_BOND>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EVM_FLOAT &evm, const int i1, const int i2, const int i3, const int i4,
                          KK_FLOAT &edihedral, KK_FLOAT *f1, KK_FLOAT *f3, KK_FLOAT *f4,
                          const KK_FLOAT &vb1x, const KK_FLOAT &vb1y, const KK_FLOAT &vb1z,
                          const KK_FLOAT &vb2x, const KK_FLOAT &vb2y, const KK_FLOAT &vb2z,
                          const KK_FLOAT &vb3x, const KK_FLOAT &vb3y, const KK_FLOAT &vb3z) const;

  KOKKOS_INLINE_FUNCTION
  void ev_tally(EVM_FLOAT &evm, const int i, const int j,
        const KK_FLOAT &evdwl, const KK_FLOAT &ecoul, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                const KK_FLOAT &dely, const KK_FLOAT &delz) const;

 protected:

  class NeighborKokkos *neighborKK;

  typename AT::t_float_1d_3_lr_randomread x;
  typename AT::t_int_1d_randomread atomtype;
  typename AT::t_float_1d_randomread q;
  typename AT::t_float_1d_3 f;
  typename AT::t_int_2d dihedrallist;

  typedef typename KKDevice<DeviceType>::value KKDeviceType;
  DAT::tdual_float_1d k_eatom;
  DAT::tdual_float_1d_6 k_vatom;
  Kokkos::View<typename AT::t_float_1d::data_type,typename AT::t_float_1d::array_layout,KKDeviceType,Kokkos::MemoryTraits<Kokkos::Atomic> > d_eatom;
  Kokkos::View<typename AT::t_float_1d_6::data_type,typename AT::t_float_1d_6::array_layout,KKDeviceType,Kokkos::MemoryTraits<Kokkos::Atomic> > d_vatom;

  DAT::tdual_float_1d k_eatom_pair;
  DAT::tdual_float_1d_6 k_vatom_pair;
  Kokkos::View<typename AT::t_float_1d::data_type,typename AT::t_float_1d::array_layout,KKDeviceType,Kokkos::MemoryTraits<Kokkos::Atomic> > d_eatom_pair;
  Kokkos::View<typename AT::t_float_1d_6::data_type,typename AT::t_float_1d_6::array_layout,KKDeviceType,Kokkos::MemoryTraits<Kokkos::Atomic> > d_vatom_pair;

  int nlocal,newton_bond;
  int eflag,vflag;
  KK_FLOAT qqrd2e;

  Kokkos::DualView<int,DeviceType> k_warning_flag;
  typename Kokkos::DualView<int,DeviceType>::t_dev d_warning_flag;
  typename Kokkos::DualView<int,DeviceType>::t_host h_warning_flag;

  typename AT::t_float_2d d_lj14_1;
  typename AT::t_float_2d d_lj14_2;
  typename AT::t_float_2d d_lj14_3;
  typename AT::t_float_2d d_lj14_4;

  typename AT::t_float_1d d_k;
  typename AT::t_float_1d d_multiplicity;
  typename AT::t_float_1d d_shift;
  typename AT::t_float_1d d_sin_shift;
  typename AT::t_float_1d d_cos_shift;
  typename AT::t_float_1d d_weight;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: Dihedral problem

Conformation of the 4 listed dihedral atoms is extreme; you may want
to check your simulation geometry.

*/
