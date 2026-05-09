/* -*- c++ -*- ----------------------------------------------------------
   Shared Kokkos rigid body storage and legacy transforms for
   fix_rigid_small_kokkos and fix_rigid_nh_small_kokkos.
------------------------------------------------------------------------- */

#ifndef LMP_FIX_RIGID_BASE_KOKKOS_H
#define LMP_FIX_RIGID_BASE_KOKKOS_H

#include "fix_rigid.h"
#include "fix_rigid_nh.h"
#include "fix_rigid_small.h"
#include "fix_rigid_nh_small.h"

#include "atom_kokkos.h"
#include "kokkos.h"
#include "kokkos_base.h"
#include "kokkos_few.h"
#include "kokkos_type.h"

// langflag
#include "Kokkos_Random.hpp"
#include "rand_pool_wrap_kokkos.h"

namespace LAMMPS_NS {

using Kokkos::fma;

template<class DeviceType, class FixRigidBase>
class FixRigidBaseKokkos : public FixRigidBase, KokkosBase {
 public:

  FixRigidBaseKokkos(class LAMMPS *, int, char **);
  ~FixRigidBaseKokkos();

protected:

  using Pointers::atom;
  using Pointers::atomKK;
  using Pointers::comm;
  using Pointers::domain;
  using Pointers::error;
  using Pointers::force;
  using Pointers::group;
  using Pointers::lmp;
  using Pointers::memory;
  using Pointers::memoryKK;
  using Pointers::neighbor;
  using Pointers::update;
  using Pointers::world;

  // should be in Pointers but im not willing to die on that hill
  class DomainKokkos *domainKK;

  using Fix::style;
  using Fix::kokkosable;
  using Fix::copymode;
  using Fix::execution_space;
  using Fix::datamask_read;
  using Fix::datamask_modify;
  using Fix::forward_comm_device;
  using Fix::reverse_comm_device;
  using Fix::exchange_comm_device;
  using Fix::sort_device;

  using FixRigidBase::commflag;
  using FixRigidBase::eflags;
  using FixRigidBase::extended;

  using FixRigidBase::body;
  using FixRigidBase::atom2body;
  using FixRigidBase::xcmimage;

  using FixRigidBase::triclinic;

  using FixRigidBase::t_start;
  using FixRigidBase::t_stop;
  using FixRigidBase::langextra;
  using FixRigidBase::t_period;
  using FixRigidBase::nlocal_body;
  // nbody_total defined below as a method (not in FixRigidSmall)
  using FixRigidBase::dtq;
  using FixRigidBase::allremap;
  using FixRigidBase::p_flag;
  using FixRigidBase::itensor;
  using FixRigidBase::inpfile;
  using FixRigidBase::readfile;

  using FixRigidBase::nmol;

  using FixRigidBase::onemols;

  using FixRigidBase::p_start;
  using FixRigidBase::p_stop;
  using FixRigidBase::p_freq;
  using FixRigidBase::displace;
  using FixRigidBase::counts;

  // FixRigidSmall protected members
  using FixRigidBase::bodyown;
  using FixRigidBase::bodytag;
  using FixRigidBase::nmax_body;
  using FixRigidBase::nghost_body;
  using FixRigidBase::seed;
  using FixRigidBase::maxextent;
  using FixRigidBase::langflag;
  using FixRigidBase::id_gravity;
  using FixRigidBase::gvec;
  using FixRigidBase::maxlang;
  using FixRigidBase::orientflag;
  using FixRigidBase::dorientflag;
  using FixRigidBase::orient;
  using FixRigidBase::dorient;
  using FixRigidBase::avec_ellipsoid;
  using FixRigidBase::avec_line;
  using FixRigidBase::avec_tri;
  using FixRigidBase::dilate_group_bit;
  using FixRigidBase::setupflag;
  using FixRigidBase::earlyflag;
  using FixRigidBase::nlinear;
  using FixRigidBase::dtv;
  using FixRigidBase::dtf;
  using FixRigidBase::nbody;
  using FixRigidBase::tstat_flag;
  using FixRigidBase::pstat_flag;
  using FixRigidBase::t_chain;
  using FixRigidBase::t_iter;
  using FixRigidBase::t_order;
  using FixRigidBase::pstyle;
  using FixRigidBase::p_chain;

  // Fix protected members
  using Fix::evflag;
  using Fix::vflag_global;
  using Fix::vflag_atom;
  using Fix::virial;
  using Fix::vatom;
  using Fix::maxvatom;
  using Fix::v_init;
  using Fix::dynamic;

  // nbody_total helper (nlocal_body + nghost_body; not a member of FixRigidSmall)
  int nbody_total() const { return FixRigidBase::nlocal_body + FixRigidBase::nghost_body; }


  // fix methods
  void post_constructor() override;
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void pre_neighbor() override;
  void final_integrate() override;
  void zero_momentum() override;
  void zero_rotation() override;
  double compute_scalar() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  void set_molecule(int, tagint, int, double *, double *, double *) override;
  bigint dof(int) override;
  void deform(int) override;

  // fix rigid protected methods
  void setup_bodies_static() override;
  void setup_bodies_dynamic() override;
  void compute_forces_and_torques() override;
  void enforce2d() override;
  void compute_dof() override;
  void remap() override;
  void grow_body() override;
  void grow_body(int);
  void image_shift() override;
  void reset_atom2body() override;

  template<bool EVFLAG>
  void set_xv();

  template<bool EVFLAG>
  void set_v();

  // Curiously Repeating Template Pattern (CRTP)
  static constexpr bool is_nh = std::is_base_of_v<FixRigidNHSmall, FixRigidBase>;

  KOKKOS_INLINE_FUNCTION
  FixRigidNHSmall* nh() {
    if constexpr (is_nh) return static_cast<FixRigidNHSmall*>(this);
    else return nullptr;
  }

  int nlocal() { return atomKK->nlocal; }

  // kokkos views
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  DAT::tdual_int_1d k_atom2body, k_bodyown, k_eflags;
  DAT::tdual_int_2d k_counts;
  DAT::tdual_tagint_1d k_bodytag;
  DAT::tdual_imageint_1d k_xcmimage;
  DAT::ttransform_kkfloat_2d k_displace, k_itensor;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView =
      KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView =
      KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<KK_ACC_FLOAT *[6], typename AT::t_kkacc_1d_6::array_layout> dup_vatom;
  NonDupScatterView<KK_ACC_FLOAT *[6], typename AT::t_kkacc_1d_6::array_layout> ndup_vatom;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  void modify_host();
  void sync_host();

  // LANGFLAG

#ifndef LMP_KOKKOS_DEBUG_RNG
    Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
    typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;
#else
    RandPoolWrap rand_pool;
    typedef RandWrap rand_type;
#endif

  TransformView<KK_FLOAT**, double**, Kokkos::LayoutRight, DeviceType> k_langextra;
  void apply_langevin_thermostat() override;

  // HOST COMM

  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  // KOKKOS BASE

  int pack_forward_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_double_1d &,
                               int, int *) override;

  void unpack_forward_comm_kokkos(int, int, DAT::tdual_double_1d &) override;

  int pack_reverse_comm_kokkos(int, int, DAT::tdual_double_1d &) override;

  void unpack_reverse_comm_kokkos(int, DAT::tdual_int_1d,
                                          DAT::tdual_double_1d &) override;

  int pack_exchange_kokkos(const int &, DAT::tdual_double_2d_lr &,
                           DAT::tdual_int_1d, DAT::tdual_int_1d, ExecutionSpace) override;

  void unpack_exchange_kokkos(DAT::tdual_double_2d_lr &, DAT::tdual_int_1d &,
                              int, int, int, ExecutionSpace) override;

  void sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &) override;

  // BODY KOKKOS

  using Body = FixRigidSmall::Body;

  struct BodyKokkos {
    int natoms, ilocal;
    KK_FLOAT mass, xcm[3], xgc[3], vcm[3], fcm[3], torque[3], quat[4], inertia[3];
    KK_FLOAT ex_space[3], ey_space[3], ez_space[3];
    KK_FLOAT xgc_body[3], angmom[3], omega[3], conjqm[4];
    imageint image;

    KOKKOS_INLINE_FUNCTION
    BodyKokkos() = default;

    KOKKOS_INLINE_FUNCTION
    BodyKokkos(const Body &b) {
    natoms = b.natoms;
    ilocal = b.ilocal;
    mass   = static_cast<KK_FLOAT>(b.mass);
    xcm[0] = static_cast<KK_FLOAT>(b.xcm[0]);
    xcm[1] = static_cast<KK_FLOAT>(b.xcm[1]);
    xcm[2] = static_cast<KK_FLOAT>(b.xcm[2]);
    xgc[0] = static_cast<KK_FLOAT>(b.xgc[0]);
    xgc[1] = static_cast<KK_FLOAT>(b.xgc[1]);
    xgc[2] = static_cast<KK_FLOAT>(b.xgc[2]);
    vcm[0] = static_cast<KK_FLOAT>(b.vcm[0]);
    vcm[1] = static_cast<KK_FLOAT>(b.vcm[1]);
    vcm[2] = static_cast<KK_FLOAT>(b.vcm[2]);
    fcm[0] = static_cast<KK_FLOAT>(b.fcm[0]);
    fcm[1] = static_cast<KK_FLOAT>(b.fcm[1]);
    fcm[2] = static_cast<KK_FLOAT>(b.fcm[2]);
    torque[0] = static_cast<KK_FLOAT>(b.torque[0]);
    torque[1] = static_cast<KK_FLOAT>(b.torque[1]);
    torque[2] = static_cast<KK_FLOAT>(b.torque[2]);
    quat[0] = static_cast<KK_FLOAT>(b.quat[0]);
    quat[1] = static_cast<KK_FLOAT>(b.quat[1]);
    quat[2] = static_cast<KK_FLOAT>(b.quat[2]);
    quat[3] = static_cast<KK_FLOAT>(b.quat[3]);
    inertia[0] = static_cast<KK_FLOAT>(b.inertia[0]);
    inertia[1] = static_cast<KK_FLOAT>(b.inertia[1]);
    inertia[2] = static_cast<KK_FLOAT>(b.inertia[2]);
    ex_space[0] = static_cast<KK_FLOAT>(b.ex_space[0]);
    ex_space[1] = static_cast<KK_FLOAT>(b.ex_space[1]);
    ex_space[2] = static_cast<KK_FLOAT>(b.ex_space[2]);
    ey_space[0] = static_cast<KK_FLOAT>(b.ey_space[0]);
    ey_space[1] = static_cast<KK_FLOAT>(b.ey_space[1]);
    ey_space[2] = static_cast<KK_FLOAT>(b.ey_space[2]);
    ez_space[0] = static_cast<KK_FLOAT>(b.ez_space[0]);
    ez_space[1] = static_cast<KK_FLOAT>(b.ez_space[1]);
    ez_space[2] = static_cast<KK_FLOAT>(b.ez_space[2]);
    xgc_body[0] = static_cast<KK_FLOAT>(b.xgc_body[0]);
    xgc_body[1] = static_cast<KK_FLOAT>(b.xgc_body[1]);
    xgc_body[2] = static_cast<KK_FLOAT>(b.xgc_body[2]);
    angmom[0] = static_cast<KK_FLOAT>(b.angmom[0]);
    angmom[1] = static_cast<KK_FLOAT>(b.angmom[1]);
    angmom[2] = static_cast<KK_FLOAT>(b.angmom[2]);
    omega[0] = static_cast<KK_FLOAT>(b.omega[0]);
    omega[1] = static_cast<KK_FLOAT>(b.omega[1]);
    omega[2] = static_cast<KK_FLOAT>(b.omega[2]);
    conjqm[0] = static_cast<KK_FLOAT>(b.conjqm[0]);
    conjqm[1] = static_cast<KK_FLOAT>(b.conjqm[1]);
    conjqm[2] = static_cast<KK_FLOAT>(b.conjqm[2]);
    conjqm[3] = static_cast<KK_FLOAT>(b.conjqm[3]);
    image = b.image;
  }

  KOKKOS_INLINE_FUNCTION
  operator Body() const {
    Body b;
    b.natoms = natoms;
    b.ilocal = ilocal;
    b.mass   = static_cast<double>(mass);
    b.xcm[0] = static_cast<double>(xcm[0]);
    b.xcm[1] = static_cast<double>(xcm[1]);
    b.xcm[2] = static_cast<double>(xcm[2]);
    b.xgc[0] = static_cast<double>(xgc[0]);
    b.xgc[1] = static_cast<double>(xgc[1]);
    b.xgc[2] = static_cast<double>(xgc[2]);
    b.vcm[0] = static_cast<double>(vcm[0]);
    b.vcm[1] = static_cast<double>(vcm[1]);
    b.vcm[2] = static_cast<double>(vcm[2]);
    b.fcm[0] = static_cast<double>(fcm[0]);
    b.fcm[1] = static_cast<double>(fcm[1]);
    b.fcm[2] = static_cast<double>(fcm[2]);
    b.torque[0] = static_cast<double>(torque[0]);
    b.torque[1] = static_cast<double>(torque[1]);
    b.torque[2] = static_cast<double>(torque[2]);
    b.quat[0] = static_cast<double>(quat[0]);
    b.quat[1] = static_cast<double>(quat[1]);
    b.quat[2] = static_cast<double>(quat[2]);
    b.quat[3] = static_cast<double>(quat[3]);
    b.inertia[0] = static_cast<double>(inertia[0]);
    b.inertia[1] = static_cast<double>(inertia[1]);
    b.inertia[2] = static_cast<double>(inertia[2]);
    b.ex_space[0] = static_cast<double>(ex_space[0]);
    b.ex_space[1] = static_cast<double>(ex_space[1]);
    b.ex_space[2] = static_cast<double>(ex_space[2]);
    b.ey_space[0] = static_cast<double>(ey_space[0]);
    b.ey_space[1] = static_cast<double>(ey_space[1]);
    b.ey_space[2] = static_cast<double>(ey_space[2]);
    b.ez_space[0] = static_cast<double>(ez_space[0]);
    b.ez_space[1] = static_cast<double>(ez_space[1]);
    b.ez_space[2] = static_cast<double>(ez_space[2]);
    b.xgc_body[0] = static_cast<double>(xgc_body[0]);
    b.xgc_body[1] = static_cast<double>(xgc_body[1]);
    b.xgc_body[2] = static_cast<double>(xgc_body[2]);
    b.angmom[0] = static_cast<double>(angmom[0]);
    b.angmom[1] = static_cast<double>(angmom[1]);
    b.angmom[2] = static_cast<double>(angmom[2]);
    b.omega[0] = static_cast<double>(omega[0]);
    b.omega[1] = static_cast<double>(omega[1]);
    b.omega[2] = static_cast<double>(omega[2]);
    b.conjqm[0] = static_cast<double>(conjqm[0]);
    b.conjqm[1] = static_cast<double>(conjqm[1]);
    b.conjqm[2] = static_cast<double>(conjqm[2]);
    b.conjqm[3] = static_cast<double>(conjqm[3]);
    b.image = image;
    return b;
  }
  }; // struct BodyKokkos

  TransformView<BodyKokkos*, Body*, Kokkos::LayoutRight, DeviceType> k_body;

  template<typename To, typename From>
  struct Transform {
    static constexpr bool is_identity = std::is_same_v<To, From>;

    KOKKOS_INLINE_FUNCTION
    static To transform(const From &x) {
      if constexpr (std::is_same_v<To, From>)
        return x;
      else if constexpr (std::is_same_v<To, BodyKokkos> && std::is_same_v<From, Body>)
        return BodyKokkos(x);
      else if constexpr (std::is_same_v<To, Body> && std::is_same_v<From, BodyKokkos>)
        return static_cast<Body>(x);
    }
  };

}; // class FixRigidBaseKokkos


// helper for printf debugging
static const std::string commflag_string(int commflag_) { return std::vector<std::string>({
    "FULL_BODY", "INITIAL", "FINAL", "FORCE_TORQUE", "VCM_ANGMOM", "XCM_MASS", "ITENSOR", "DOF"
    })[commflag_];
  }

}    // namespace LAMMPS_NS

#endif // !LMP_FIX_RIGID_BASE_KOKKOS_H
