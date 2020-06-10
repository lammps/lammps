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

#ifdef PAIR_CLASS

PairStyle(tersoff/zbl/kk,PairTersoffZBLKokkos<Device>)
PairStyle(tersoff/zbl/kk/device,PairTersoffZBLKokkos<Device>)
PairStyle(tersoff/zbl/kk/host,PairTersoffZBLKokkos<Host>)

#else

#ifndef LMP_PAIR_TERSOFF_ZBL_KOKKOS_H
#define LMP_PAIR_TERSOFF_ZBL_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_tersoff_zbl.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<int NEIGHFLAG, int EVFLAG>
struct TagPairTersoffZBLComputeHalf{};

template<int NEIGHFLAG, int EVFLAG>
struct TagPairTersoffZBLComputeFullA{};

template<int NEIGHFLAG, int EVFLAG>
struct TagPairTersoffZBLComputeFullB{};

struct TagPairTersoffZBLComputeShortNeigh{};

template<ExecutionSpace Space>
class PairTersoffZBLKokkos : public PairTersoffZBL {
 public:
  enum {EnabledNeighFlags=FULL};
  enum {COUL_FLAG=0};
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;
  typedef EV_FLOAT value_type;

  PairTersoffZBLKokkos(class LAMMPS *);
  virtual ~PairTersoffZBLKokkos();
  virtual void compute(int, int);
  void init_style();

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffZBLComputeHalf<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffZBLComputeHalf<NEIGHFLAG,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffZBLComputeFullA<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffZBLComputeFullA<NEIGHFLAG,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffZBLComputeFullB<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffZBLComputeFullB<NEIGHFLAG,EVFLAG>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffZBLComputeShortNeigh, const int&) const;
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT ters_fc_k(const int &i, const int &j, const int &k, const KK_FLOAT &r) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT ters_dfc(const int &i, const int &j, const int &k, const KK_FLOAT &r) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT ters_fa_k(const int &i, const int &j, const int &k, const KK_FLOAT &r) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT ters_dfa(const int &i, const int &j, const int &k, const KK_FLOAT &r) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT ters_bij_k(const int &i, const int &j, const int &k, const KK_FLOAT &bo) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT ters_dbij(const int &i, const int &j, const int &k, const KK_FLOAT &bo) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT bondorder(const int &i, const int &j, const int &k,
        const KK_FLOAT &rij, const KK_FLOAT &dx1, const KK_FLOAT &dy1, const KK_FLOAT &dz1,
        const KK_FLOAT &rik, const KK_FLOAT &dx2, const KK_FLOAT &dy2, const KK_FLOAT &dz2) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT ters_gijk(const int &i, const int &j, const int &k, const KK_FLOAT &cos) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT ters_dgijk(const int &i, const int &j, const int &k, const KK_FLOAT &cos) const;

  KOKKOS_INLINE_FUNCTION
  void ters_dthb(const int &i, const int &j, const int &k, const KK_FLOAT &prefactor,
        const KK_FLOAT &rij, const KK_FLOAT &dx1, const KK_FLOAT &dy1, const KK_FLOAT &dz1,
        const KK_FLOAT &rik, const KK_FLOAT &dx2, const KK_FLOAT &dy2, const KK_FLOAT &dz2,
        KK_FLOAT *fi, KK_FLOAT *fj, KK_FLOAT *fk) const;

  KOKKOS_INLINE_FUNCTION
  void ters_dthbj(const int &i, const int &j, const int &k, const KK_FLOAT &prefactor,
        const KK_FLOAT &rij, const KK_FLOAT &dx1, const KK_FLOAT &dy1, const KK_FLOAT &dz1,
        const KK_FLOAT &rik, const KK_FLOAT &dx2, const KK_FLOAT &dy2, const KK_FLOAT &dz2,
        KK_FLOAT *fj, KK_FLOAT *fk) const;

  KOKKOS_INLINE_FUNCTION
  void ters_dthbk(const int &i, const int &j, const int &k, const KK_FLOAT &prefactor,
        const KK_FLOAT &rij, const KK_FLOAT &dx1, const KK_FLOAT &dy1, const KK_FLOAT &dz1,
        const KK_FLOAT &rik, const KK_FLOAT &dx2, const KK_FLOAT &dy2, const KK_FLOAT &dz2,
        KK_FLOAT *fk) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT vec3_dot(const KK_FLOAT x[3], const KK_FLOAT y[3]) const {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
  }

  KOKKOS_INLINE_FUNCTION
  void vec3_add(const KK_FLOAT x[3], const KK_FLOAT y[3], KK_FLOAT * const z) const {
    z[0] = x[0]+y[0]; z[1] = x[1]+y[1]; z[2] = x[2]+y[2];
  }

  KOKKOS_INLINE_FUNCTION
  void vec3_scale(const KK_FLOAT k, const KK_FLOAT x[3], KK_FLOAT y[3]) const {
    y[0] = k*x[0]; y[1] = k*x[1]; y[2] = k*x[2];
  }

  KOKKOS_INLINE_FUNCTION
  void vec3_scaleadd(const KK_FLOAT k, const KK_FLOAT x[3], const KK_FLOAT y[3], KK_FLOAT * const z) const {
    z[0] = k*x[0]+y[0]; z[1] = k*x[1]+y[1]; z[2] = k*x[2]+y[2];
  }

  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const;

  struct params_ters{
    KOKKOS_INLINE_FUNCTION
    params_ters(){powerm=0;gamma=0;lam3=0;c=0;d=0;h=0;powern=0;beta=0;lam2=0;bigb=0;
          bigr=0;bigd=0;lam1=0;biga=0;cutsq=0;c1=0;c2=0;c3=0;c4=0;Z_i=0;Z_j=0;ZBLcut=0;ZBLexpscale=0;};
    KOKKOS_INLINE_FUNCTION
    params_ters(int i){powerm=0;gamma=0;lam3=0;c=0;d=0;h=0;powern=0;beta=0;lam2=0;bigb=0;
          bigr=0;bigd=0;lam1=0;biga=0;cutsq=0;c1=0;c2=0;c3=0;c4=0;Z_i=0;Z_j=0;ZBLcut=0;ZBLexpscale=0;};
    KK_FLOAT powerm, gamma, lam3, c, d, h, powern, beta, lam2, bigb, bigr,
      bigd, lam1, biga, cutsq, c1, c2, c3, c4, Z_i, Z_j, ZBLcut, ZBLexpscale;
  };

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                  const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void v_tally3(EV_FLOAT &ev, const int &i, const int &j, const int &k,
    KK_FLOAT *fj, KK_FLOAT *fk, KK_FLOAT *drij, KK_FLOAT *drik) const;

  KOKKOS_INLINE_FUNCTION
  void v_tally3_atom(EV_FLOAT &ev, const int &i, const int &j, const int &k,
    KK_FLOAT *fj, KK_FLOAT *fk, KK_FLOAT *drji, KK_FLOAT *drjk) const;

  void allocate();
  void setup_params();

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT fermi_k(const int &i, const int &j, const int &k, const KK_FLOAT &r) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT fermi_d_k(const int &i, const int &j, const int &k, const KK_FLOAT &r) const;

 protected:
  void cleanup_copy();

  typedef Kokkos::DualView<int***,DeviceType> tdual_int_3d;
  Kokkos::DualView<params_ters***,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_ters***,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um paramskk;
  // hardwired to space for 12 atom types
  //params_ters m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  int inum;
  typename AT::t_float_1d_3_lr_randomread x;
  typename AT::t_float_1d_3 f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_tagint_1d tag;

  DAT::tdual_float_1d k_eatom;
  DAT::tdual_float_1d_6 k_vatom;
  typename AT::t_float_1d d_eatom;
  typename AT::t_float_1d_6 d_vatom;

  int need_dup;
  Kokkos::Experimental::ScatterView<typename AT::t_float_1d_3::data_type, typename AT::t_float_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_f;
  Kokkos::Experimental::ScatterView<typename AT::t_float_1d::data_type, typename AT::t_float_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_eatom;
  Kokkos::Experimental::ScatterView<typename AT::t_float_1d_6::data_type, typename AT::t_float_1d_6::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_vatom;
  Kokkos::Experimental::ScatterView<typename AT::t_float_1d_3::data_type, typename AT::t_float_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_f;
  Kokkos::Experimental::ScatterView<typename AT::t_float_1d::data_type, typename AT::t_float_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_eatom;
  Kokkos::Experimental::ScatterView<typename AT::t_float_1d_6::data_type, typename AT::t_float_1d_6::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_vatom;


  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;
  //NeighListKokkos<Space> k_list;

  int neighflag,newton_pair;
  int nlocal,nall,eflag,vflag;

  Kokkos::View<int**,DeviceType> d_neighbors_short;
  Kokkos::View<int*,DeviceType> d_numneigh_short;

  // ZBL
  KK_FLOAT global_a_0;                // Bohr radius for Coulomb repulsion
  KK_FLOAT global_epsilon_0;        // permittivity of vacuum for Coulomb repulsion
  KK_FLOAT global_e;                // proton charge (negative of electron charge)

  friend void pair_virial_fdotr_compute<Space,PairTersoffZBLKokkos>(PairTersoffZBLKokkos*);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Pair tersoff/zbl/kk requires metal or real units

This is a current restriction of this pair potential.

E: Cannot (yet) use full neighbor list style with tersoff/zbl/kk

Self-explanatory.

E: Cannot use chosen neighbor list style with tersoff/zbl/kk

Self-explanatory.

*/
