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

#ifndef LMP_FIX_RIGID_NH_SMALL_KOKKOS_H
#define LMP_FIX_RIGID_NH_SMALL_KOKKOS_H

#include "fix_rigid_nh_small.h"
#include "kokkos_type.h"
#include "kokkos_few.h"

namespace LAMMPS_NS {

struct TagRigidNHInitialIntegrate {};
struct TagRigidNHFinalIntegrate {};
struct TagRigidNHAccumKE {};

template<int TRICLINIC, int NEIGHFLAG, int EVFLAG>
struct TagRigidNHSetXV {};

template<int TRICLINIC, int NEIGHFLAG, int EVFLAG>
struct TagRigidNHSetV {};

struct TagRigidNHComputeForcesTorquesZero {};
struct TagRigidNHComputeForcesTorques {};
struct TagRigidNHImageShift {};
struct TagRigidNHEnforce2d {};

template<class DeviceType>
class FixRigidNHSmallKokkos : public FixRigidNHSmall {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  FixRigidNHSmallKokkos(class LAMMPS *, int, char **);
  ~FixRigidNHSmallKokkos() override;

  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void pre_neighbor() override;
  void setup_pre_neighbor() override;
  void post_force(int) override;

  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  void set_molecule(int, tagint, int, double *, double *, double *) override;

  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  void reset_dt() override;
  void zero_momentum() override;
  void zero_rotation() override;

  double compute_scalar() override;

  void grow_body() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidNHInitialIntegrate, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidNHFinalIntegrate, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidNHAccumKE, const int&, double&, double&) const;

  template<int TRICLINIC, int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidNHSetXV<TRICLINIC,NEIGHFLAG,EVFLAG>, const int&) const;

  template<int TRICLINIC, int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidNHSetXV<TRICLINIC,NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT &) const;

  template<int TRICLINIC, int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidNHSetV<TRICLINIC,NEIGHFLAG,EVFLAG>, const int&) const;

  template<int TRICLINIC, int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidNHSetV<TRICLINIC,NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT &) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidNHComputeForcesTorquesZero, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidNHComputeForcesTorques, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidNHImageShift, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidNHEnforce2d, const int&) const;

 protected:
  class AtomKokkos *atomKK;
  class DomainKokkos *domainKK;
  ExecutionSpace execution_space;

  Kokkos::DualView<Body*, Kokkos::LayoutRight, DeviceType> k_body;
  typename Kokkos::DualView<Body*, Kokkos::LayoutRight, DeviceType>::t_dev d_body;

  DAT::tdual_int_1d k_bodyown;
  DAT::tdual_tagint_1d k_bodytag;
  DAT::tdual_int_1d k_atom2body;
  DAT::tdual_imageint_1d k_xcmimage;
  DAT::tdual_kkfloat_2d k_displace;
  DAT::tdual_int_1d k_eflags;

  typename AT::t_int_1d d_bodyown;
  typename AT::t_tagint_1d d_bodytag;
  typename AT::t_int_1d d_atom2body;
  typename AT::t_imageint_1d d_xcmimage;
  typename AT::t_kkfloat_2d d_displace;
  typename AT::t_int_1d d_eflags;

  typename AT::t_kkfloat_1d_3_lr d_x;
  typename AT::t_kkfloat_1d_3 d_v;
  typename AT::t_kkacc_1d_3 d_f;
  typename AT::t_kkfloat_1d d_rmass;
  typename AT::t_kkfloat_1d d_mass;
  typename AT::t_int_1d d_type;
  typename AT::t_int_1d d_mask;
  typename AT::t_imageint_1d d_image;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  Few<double,3> d_prd;
  Few<double,6> d_h;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView =
      KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView =
      KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<KK_ACC_FLOAT *[6], typename AT::t_kkacc_1d_6::array_layout> dup_vatom;
  NonDupScatterView<KK_ACC_FLOAT *[6], typename AT::t_kkacc_1d_6::array_layout> ndup_vatom;

  KK_FLOAT d_scale_t[3], d_scale_v[3], d_scale_r;
  KK_FLOAT d_dtf2;

  void remap();
  void compute_forces_and_torques_kokkos();
  void enforce2d_kokkos();
  void image_shift_kokkos();
  void grow_body_kokkos();

  void sync_fix_data_device();

  template<int TRICLINIC, int EVFLAG>
  void set_xv_kokkos();

  template<int TRICLINIC, int EVFLAG>
  void set_v_kokkos();
};

}    // namespace LAMMPS_NS

#endif
