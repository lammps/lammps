/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(granular/kk,PairGranularKokkos<LMPDeviceType>);
PairStyle(granular/kk/device,PairGranularKokkos<LMPDeviceType>);
PairStyle(granular/kk/host,PairGranularKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_PAIR_GRANULAR_KOKKOS_H
#define LMP_PAIR_GRANULAR_KOKKOS_H

#include "pair_granular.h"
#include "pair_kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template <class DeviceType> class FixNeighHistoryKokkos;

template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int HISTORYUPDATE>
struct TagPairGranularCompute {};

template <class DeviceType>
class PairGranularKokkos : public PairGranular {
 public:
  using device_type = DeviceType;
  using AT = ArrayTypes<DeviceType>;
  using value_type = EV_FLOAT;

  PairGranularKokkos(class LAMMPS *);
  ~PairGranularKokkos() override;
  void compute(int, int) override;
  void init_style() override;
  double init_one(int, int) override;

  template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int HISTORYUPDATE>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairGranularCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,HISTORYUPDATE>,
                  const int, EV_FLOAT &) const;
  template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int HISTORYUPDATE>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairGranularCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,HISTORYUPDATE>,
                  const int) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &, int, int, KK_FLOAT, KK_FLOAT, KK_FLOAT,
                    KK_FLOAT, KK_FLOAT, KK_FLOAT) const;

 protected:
  enum {N_HOOKE=1, N_HERTZ=2};
  enum {D_VELOCITY=1, D_MASS_VELOCITY=2, D_VISCOELASTIC=3, D_TSUJI=4};
  enum {T_LINEAR_NOHISTORY=1, T_LINEAR_HISTORY=2, T_LINEAR_CLASSIC=3,
        T_MINDLIN_CLASSIC=4, T_MINDLIN=5};
  enum {MI_NORMAL=0, MI_DAMPING, MI_TANGENTIAL, MI_HISTORY_INDEX,
        MI_LIMIT_DAMPING, MI_MINDLIN_FORCE, MI_MINDLIN_RESCALE, MI_COUNT};
  enum {MF_NORMAL_K=0, MF_DAMP, MF_TANGENTIAL_K, MF_XT, MF_MU, MF_COUNT};

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3_randomread v, omega;
  typename AT::t_kkacc_1d_3 f, torque;
  typename AT::t_int_1d_randomread type, mask;
  typename AT::t_kkfloat_1d_randomread rmass, radius;
  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist, d_numneigh;
  typename AT::t_int_2d d_firsttouch;
  typename AT::t_kkfloat_2d d_firsthistory;

  Kokkos::View<int**,Kokkos::LayoutRight,DeviceType> d_model_i;
  Kokkos::View<KK_FLOAT**,Kokkos::LayoutRight,DeviceType> d_model_f;
  Kokkos::View<int**,Kokkos::LayoutRight,DeviceType> d_type_model;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  FixNeighHistoryKokkos<DeviceType> *fix_historyKK;
  int neighflag, newton_pair, nlocal, nall, eflag, vflag;
  int use_history_kk, size_history_kk;
  KK_FLOAT dt_kk;

  void create_kokkos_history();
  void pack_models();

  KOKKOS_INLINE_FUNCTION
  int sbmask(const int &j) const { return j >> SBBITS & 3; }

  friend void pair_virial_fdotr_compute<PairGranularKokkos>(PairGranularKokkos*);
};

} // namespace LAMMPS_NS

#endif
#endif
