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

DihedralStyle(charmm/kk,DihedralCharmmKokkos<LMPDeviceType>)
DihedralStyle(charmm/kk/device,DihedralCharmmKokkos<LMPDeviceType>)
DihedralStyle(charmm/kk/host,DihedralCharmmKokkos<LMPHostType>)

#else

#ifndef LMP_DIHEDRAL_CHARMM_KOKKOS_H
#define LMP_DIHEDRAL_CHARMM_KOKKOS_H

#include "dihedral_charmm.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct s_EVM_FLOAT {
  E_FLOAT evdwl;
  E_FLOAT ecoul;
  E_FLOAT emol;
  F_FLOAT v[6];
  F_FLOAT vp[6];
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

template<class DeviceType>
class DihedralCharmmKokkos : public DihedralCharmm {
 public:
  typedef DeviceType device_type;
  typedef EVM_FLOAT value_type;
  typedef ArrayTypes<DeviceType> AT;

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
                          F_FLOAT &edihedral, F_FLOAT *f1, F_FLOAT *f3, F_FLOAT *f4,
                          const F_FLOAT &vb1x, const F_FLOAT &vb1y, const F_FLOAT &vb1z,
                          const F_FLOAT &vb2x, const F_FLOAT &vb2y, const F_FLOAT &vb2z,
                          const F_FLOAT &vb3x, const F_FLOAT &vb3y, const F_FLOAT &vb3z) const;

  KOKKOS_INLINE_FUNCTION
  void ev_tally(EVM_FLOAT &evm, const int i, const int j,
        const F_FLOAT &evdwl, const F_FLOAT &ecoul, const F_FLOAT &fpair, const F_FLOAT &delx,
                const F_FLOAT &dely, const F_FLOAT &delz) const;

 protected:

  class NeighborKokkos *neighborKK;

  typename AT::t_x_array_randomread x;
  typename AT::t_int_1d_randomread atomtype;
  typename AT::t_ffloat_1d_randomread q;
  typename AT::t_f_array f;
  typename AT::t_int_2d dihedrallist;

  typedef typename KKDevice<DeviceType>::value KKDeviceType;
  Kokkos::DualView<E_FLOAT*,Kokkos::LayoutRight,KKDeviceType> k_eatom;
  Kokkos::DualView<F_FLOAT*[6],Kokkos::LayoutRight,KKDeviceType> k_vatom;
  Kokkos::View<E_FLOAT*,Kokkos::LayoutRight,KKDeviceType,Kokkos::MemoryTraits<Kokkos::Atomic> > d_eatom;
  Kokkos::View<F_FLOAT*[6],Kokkos::LayoutRight,KKDeviceType,Kokkos::MemoryTraits<Kokkos::Atomic> > d_vatom;

  Kokkos::DualView<E_FLOAT*,Kokkos::LayoutRight,KKDeviceType> k_eatom_pair;
  Kokkos::DualView<F_FLOAT*[6],Kokkos::LayoutRight,KKDeviceType> k_vatom_pair;
  Kokkos::View<E_FLOAT*,Kokkos::LayoutRight,KKDeviceType,Kokkos::MemoryTraits<Kokkos::Atomic> > d_eatom_pair;
  Kokkos::View<F_FLOAT*[6],Kokkos::LayoutRight,KKDeviceType,Kokkos::MemoryTraits<Kokkos::Atomic> > d_vatom_pair;

  int nlocal,newton_bond;
  int eflag,vflag;
  double qqrd2e;

  Kokkos::DualView<int,DeviceType> k_warning_flag;
  typename Kokkos::DualView<int,DeviceType>::t_dev d_warning_flag;
  typename Kokkos::DualView<int,DeviceType>::t_host h_warning_flag;

  typename AT::t_ffloat_2d d_lj14_1;
  typename AT::t_ffloat_2d d_lj14_2;
  typename AT::t_ffloat_2d d_lj14_3;
  typename AT::t_ffloat_2d d_lj14_4;

  typename AT::t_ffloat_1d d_k;
  typename AT::t_ffloat_1d d_multiplicity;
  typename AT::t_ffloat_1d d_shift;
  typename AT::t_ffloat_1d d_sin_shift;
  typename AT::t_ffloat_1d d_cos_shift;
  typename AT::t_ffloat_1d d_weight;

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
