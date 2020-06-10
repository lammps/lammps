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

PairStyle(snap/kk,PairSNAPKokkos<Device>)
PairStyle(snap/kk/device,PairSNAPKokkos<Device>)
PairStyle(snap/kk/host,PairSNAPKokkos<Host>)

#else

#ifndef LMP_PAIR_SNAP_KOKKOS_H
#define LMP_PAIR_SNAP_KOKKOS_H

#include "pair_snap.h"
#include "kokkos_type.h"
#include "neigh_list_kokkos.h"
#include "sna_kokkos.h"
#include "pair_kokkos.h"

namespace LAMMPS_NS {

template<int NEIGHFLAG, int EVFLAG>
struct TagPairSNAPComputeForce{};

struct TagPairSNAPBeta{};
struct TagPairSNAPComputeNeigh{};
struct TagPairSNAPPreUi{};
struct TagPairSNAPComputeUi{};
struct TagPairSNAPComputeUiCPU{};
struct TagPairSNAPComputeZi{};
struct TagPairSNAPComputeBi{};
struct TagPairSNAPZeroYi{};
struct TagPairSNAPComputeYi{};
struct TagPairSNAPComputeFusedDeidrj{};
struct TagPairSNAPComputeDuidrjCPU{};
struct TagPairSNAPComputeDeidrjCPU{};

template<ExecutionSpace Space>
class PairSNAPKokkos : public PairSNAP {
public:
  enum {EnabledNeighFlags=FULL|HALF|HALFTHREAD};
  enum {COUL_FLAG=0};
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;
  typedef EV_FLOAT value_type;

  PairSNAPKokkos(class LAMMPS *);
  ~PairSNAPKokkos();

  void coeff(int, char**);
  void init_style();
  double init_one(int, int);
  void compute(int, int);
  double memory_usage();

  template<class TagStyle>
  void check_team_size_for(int, int&, int);

  template<class TagStyle>
  void check_team_size_reduce(int, int&, int);

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG>,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG> >::member_type& team) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG>,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG> >::member_type& team, EV_FLOAT&) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeNeigh,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeNeigh>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPPreUi,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPPreUi>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeUi,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeUi>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeUiCPU,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeUiCPU>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeZi,const int& ii) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeBi,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeBi>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPZeroYi,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPZeroYi>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeYi,const int& ii) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeFusedDeidrj,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeFusedDeidrj>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeDuidrjCPU,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeDuidrjCPU>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeDeidrjCPU,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeDeidrjCPU>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPBeta,const int& ii) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void v_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &fx, const KK_FLOAT &fy, const KK_FLOAT &fz,
      const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const;

protected:
  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  DAT::tdual_float_1d k_eatom;
  DAT::tdual_float_1d_6 k_vatom;
  typename AT::t_float_1d d_eatom;
  typename AT::t_float_1d_6 d_vatom;

  typedef Kokkos::View<KK_FLOAT**> t_bvec;
  t_bvec bvec;
  typedef Kokkos::View<KK_FLOAT***> t_dbvec;
  t_dbvec dbvec;
  SNAKokkos<Space> snaKK;

  int inum,max_neighs,chunk_size,chunk_offset;
  int host_flag;

  int eflag,vflag;

  void allocate();
  //void read_files(char *, char *);
  /*template<ExecutionSpace Space>
inline int equal(KK_FLOAT* x,KK_FLOAT* y);
  template<ExecutionSpace Space>
inline KK_FLOAT dist2(KK_FLOAT* x,KK_FLOAT* y);
  KK_FLOAT extra_cutoff();
  void load_balance();
  void set_sna_to_shared(int snaid,int i);
  void build_per_atom_arrays();*/

  int neighflag;

  Kokkos::View<int*, DeviceType> ilistmast;
  Kokkos::View<int*, DeviceType> ghostilist;
  Kokkos::View<int*, DeviceType> ghostnumneigh;
  Kokkos::View<int*, DeviceType> ghostneighs;
  Kokkos::View<int*, DeviceType> ghostfirstneigh;

  Kokkos::View<int**, Kokkos::LayoutRight, DeviceType> i_pairs;
  Kokkos::View<int***, Kokkos::LayoutRight, DeviceType> i_rij;
  Kokkos::View<int**, Kokkos::LayoutRight, DeviceType> i_inside;
  Kokkos::View<KK_FLOAT**, Kokkos::LayoutRight, DeviceType> i_wj;
  Kokkos::View<KK_FLOAT***, Kokkos::LayoutRight, DeviceType>i_rcutij;
  Kokkos::View<int*, DeviceType> i_ninside;
  Kokkos::View<KK_FLOAT****, Kokkos::LayoutRight, DeviceType> i_uarraytot_r, i_uarraytot_i;
  Kokkos::View<KK_FLOAT******, Kokkos::LayoutRight, DeviceType> i_zarray_r, i_zarray_i;

  Kokkos::View<KK_FLOAT*, DeviceType> d_radelem;              // element radii
  Kokkos::View<KK_FLOAT*, DeviceType> d_wjelem;               // elements weights
  Kokkos::View<KK_FLOAT**, Kokkos::LayoutRight, DeviceType> d_coeffelem;           // element bispectrum coefficients
  Kokkos::View<int*, DeviceType> d_map;                    // mapping from atom types to elements
  Kokkos::View<int*, DeviceType> d_ninside;                // ninside for all atoms in list
  Kokkos::View<KK_FLOAT**, DeviceType> d_beta;                // betas for all atoms in list
  Kokkos::View<KK_FLOAT**, DeviceType> d_bispectrum;          // bispectrum components for all atoms in list

  typedef Kokkos::DualView<KK_FLOAT**, DeviceType> tdual_fparams;
  tdual_fparams k_cutsq;
  typedef Kokkos::View<const KK_FLOAT**, DeviceType,
      Kokkos::MemoryTraits<Kokkos::RandomAccess> > t_fparams_rnd;
  t_fparams_rnd rnd_cutsq;

  typename AT::t_float_1d_3_lr_randomread x;
  typename AT::t_float_1d_3 f;
  typename AT::t_int_1d_randomread type;

  int need_dup;
  Kokkos::Experimental::ScatterView<typename AT::t_float_1d_3::data_type, typename AT::t_float_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_f;
  Kokkos::Experimental::ScatterView<typename AT::t_float_1d_6::data_type, typename AT::t_float_1d_6::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_vatom;
  Kokkos::Experimental::ScatterView<typename AT::t_float_1d_3::data_type, typename AT::t_float_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_f;
  Kokkos::Experimental::ScatterView<typename AT::t_float_1d_6::data_type, typename AT::t_float_1d_6::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_vatom;

  friend void pair_virial_fdotr_compute<Space,PairSNAPKokkos>(PairSNAPKokkos*);

};

}

#endif
#endif
