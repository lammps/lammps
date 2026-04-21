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

#include "kokkos.h"
#include "kokkos_base.h"
#include "kokkos_few.h"
#include "kokkos_type.h"

// langflag
#include "Kokkos_Random.hpp"
#include "rand_pool_wrap_kokkos.h"

namespace LAMMPS_NS {

using Kokkos::fma;

struct TagRigidResetAtom2Body {};

template<int TRICLINIC, int NEIGHFLAG, int EVFLAG>
struct TagRigidSetXV {};

template<int TRICLINIC, int NEIGHFLAG, int EVFLAG>
struct TagRigidSetV {};

template<class DeviceType, class FixRigidBase>
class FixRigidBaseKokkos : public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  FixRigidBaseKokkos(Atom*, Domain*);
  ~FixRigidBaseKokkos();

  // templated tagged operators

  template<bool TRICLINIC, bool NEIGHFLAG, bool EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidSetXV<TRICLINIC,NEIGHFLAG,EVFLAG>, const int&) const;

  template<bool TRICLINIC, bool NEIGHFLAG, bool EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidSetXV<TRICLINIC,NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT&) const;

  template<bool TRICLINIC, bool NEIGHFLAG, bool EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidSetV<TRICLINIC,NEIGHFLAG,EVFLAG>, const int&) const;

  template<bool TRICLINIC, bool NEIGHFLAG, bool EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidSetV<TRICLINIC,NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidResetAtom2Body, const int &i) const;

protected:

  // fix methods
  void setup_base(int);
  void setup_pre_neighbor_base();
  void initial_integrate_base(int);
  void pre_neighbor_base();
  void final_integrate_base();
  void zero_momentum_base();
  void zero_rotation_base();
  double compute_scalar_base();
  void deform_base(int);

/*
  virtual void post_constructor() {}
  virtual void init() {}
  void init_flags();
  virtual void init_list(int, class NeighList *) {}
  virtual void setup(int) {}
  virtual void setup_pre_exchange() {}
  virtual void setup_pre_neighbor() {}
  virtual void setup_post_neighbor() {}
  virtual void setup_pre_force(int) {}
  virtual void setup_pre_reverse(int, int) {}
  virtual void min_setup(int) {}
  virtual void initial_integrate(int) {}
  virtual void post_integrate() {}
  virtual void pre_exchange() {}
  virtual void pre_neighbor() {}
  virtual void post_neighbor() {}
  virtual void pre_force(int) {}
  virtual void pre_reverse(int, int) {}
  virtual void post_force(int) {}
  virtual void final_integrate() {}
  virtual void fused_integrate(int) {}
  virtual void end_of_step() {}
  virtual void post_run() {}
  virtual void write_restart(FILE *) {}
  virtual void write_restart_file(const char *) {}
  virtual void restart(char *) {}

  virtual void grow_arrays(int) {}
  virtual void copy_arrays(int, int, int) {}
  virtual void set_arrays(int) {}
  virtual void update_arrays(int, int) {}
  virtual void set_molecule(int, tagint, int, double *, double *, double *);
  virtual void clear_bonus() {}

  virtual int pack_border(int, int *, double *) { return 0; }
  virtual int unpack_border(int, int, double *) { return 0; }
  virtual int pack_exchange(int, double *) { return 0; }
  virtual int unpack_exchange(int, double *) { return 0; }
  virtual int pack_restart(int, double *) { return 0; }
  virtual void unpack_restart(int, int) {}
  virtual int size_restart(int) { return 0; }
  virtual int maxsize_restart() { return 0; }

  virtual void setup_pre_force_respa(int, int) {}
  virtual void initial_integrate_respa(int, int, int) {}
  virtual void post_integrate_respa(int, int) {}
  virtual void pre_force_respa(int, int, int) {}
  virtual void post_force_respa(int, int, int) {}
  virtual void final_integrate_respa(int, int) {}

  virtual void min_pre_exchange() {}
  virtual void min_pre_neighbor() {}
  virtual void min_post_neighbor() {}
  virtual void min_pre_force(int) {}
  virtual void min_pre_reverse(int, int) {}
  virtual void min_post_force(int) {}

  virtual double min_energy(double *) { return 0.0; }
  virtual void min_store() {}
  virtual void min_clearstore() {}
  virtual void min_pushstore() {}
  virtual void min_popstore() {}
  virtual int min_reset_ref() { return 0; }
  virtual void min_step(double, double *) {}
  virtual double max_alpha(double *) { return 0.0; }
  virtual int min_dof() { return 0; }
  virtual int pack_forward_comm(int, int *, double *, int, int *) { return 0; }
  virtual void unpack_forward_comm(int, int, double *) {}
  virtual int pack_reverse_comm_size(int, int) { return 0; }
  virtual int pack_reverse_comm(int, int, double *) { return 0; }
  virtual void unpack_reverse_comm(int, int *, double *) {}
  virtual void reset_grid() {};
  virtual void pack_forward_grid(int, void *, int, int *) {};
  virtual void unpack_forward_grid(int, void *, int, int *) {};
  virtual void pack_reverse_grid(int, void *, int, int *) {};
  virtual void unpack_reverse_grid(int, void *, int, int *) {};
  virtual void pack_remap_grid(int, void *, int, int *) {};
  virtual void unpack_remap_grid(int, void *, int, int *) {};
  virtual int unpack_read_grid(int, char *) { return 0; };
  virtual void pack_write_grid(int, void *) {};
  virtual void unpack_write_grid(int, void *, int *) {};
  virtual int get_grid_by_name(const std::string &, int &) { return -1; };
  virtual void *get_grid_by_index(int) { return nullptr; };
  virtual int get_griddata_by_name(int, const std::string &, int &) { return -1; };
  virtual void *get_griddata_by_index(int) { return nullptr; };
  virtual double compute_scalar() { return 0.0; }
  virtual bigint dof(int) { return 0; }
  virtual void deform(int) {}
  virtual void reset_target(double) {}
  virtual void reset_dt() {}
  virtual void read_data_header(char *) {}
  virtual void read_data_section(char *, int, char *, tagint) {}
  virtual bigint read_data_skip_lines(char *) { return 0; }
  virtual void write_data_header(FILE *, int) {}
  virtual void write_data_section_size(int, int &, int &) {}
  virtual void write_data_section_pack(int, double **) {}
  virtual void write_data_section_keyword(int, FILE *) {}
  virtual void write_data_section(int, FILE *, int, double **, int) {}
  virtual void rebuild_special() {}
  virtual int image(int *&, double **&) { return 0; }
  virtual int modify_param(int, char **) { return 0; }
  virtual void *extract(const char *, int &) { return nullptr; }
*/


  // protected methods
  void compute_forces_and_torques_base();
  void enforce2d_base();
  void remap_base();
  void grow_arrays_base(int);
  void grow_body_base(int);
  void image_shift_base();
  void reset_atom2body_base();

  template<bool TRICLINIC, bool EVFLAG>
  void set_xv_base();

  template<bool TRICLINIC, bool EVFLAG>
  void set_v_base();

  // CRTP accessors — every concrete KK class inherits from
  // both FixRigidBase and FixRigidBaseKokkos<DeviceType, FixRigidBase>
  FixRigidBase* base() { return static_cast<FixRigidBase*>(this); }
  const FixRigidBase* base() const { return static_cast<const FixRigidBase*>(this); }

  static constexpr bool is_nh    = std::is_base_of_v<FixRigidNH,    FixRigidBase>;
  static constexpr bool is_small = std::is_base_of_v<FixRigidSmall, FixRigidBase>;

  int nbody_total() { return base()->nlocal_body + base()->nghost_body; }

  class AtomKokkos *atomKK;
  class DomainKokkos *domainKK;
  ExecutionSpace execution_space;

  DAT::tdual_int_1d k_atom2body, k_bodyown, k_eflags;
  DAT::tdual_tagint_1d k_bodytag;
  DAT::tdual_imageint_1d k_xcmimage;
  
  TransformView<KK_FLOAT**, double**, Kokkos::LayoutRight, DeviceType> k_displace;
  typename AT::t_kkfloat_2d d_displace;

  int map_style;
  DAT::tdual_int_1d k_map_array;
  dual_hash_type k_map_hash;

  int comm_me;
  bigint ntimestep;

  typename AT::t_kkfloat_1d_3_lr d_x;
  typename AT::t_kkfloat_1d_3 d_v;
  typename AT::t_kkacc_1d_3 d_f;
  typename AT::t_kkfloat_1d d_rmass, d_mass;
  typename AT::t_int_1d d_type, d_mask, d_atom2body, d_bodyown, d_eflags;
  typename AT::t_tagint_1d d_tag, d_bodytag;
  typename AT::t_imageint_1d d_image, d_xcmimage;

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



  void modify_host_base();
  void modify_device_base();
  void sync_host_base();
  void sync_device_base();

  // LANGFLAG

#ifndef LMP_KOKKOS_DEBUG_RNG
    Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
    typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;
#else
    RandPoolWrap rand_pool;
    typedef RandWrap rand_type;
#endif

  TransformView<KK_FLOAT**, double**, Kokkos::LayoutRight, DeviceType> k_langextra;
  void apply_langevin_thermostat_base();

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
  Kokkos::View <BodyKokkos*,        Kokkos::LayoutRight, DeviceType> d_body;

  template<typename To, typename From>
  struct Transform {
    static constexpr bool is_identity = std::is_same_v<To, From>;
    KOKKOS_INLINE_FUNCTION
    static To transform(const From &x) { return x; }
  };

  template<>
  struct Transform<BodyKokkos, Body> {
    static constexpr bool is_identity = false;
    KOKKOS_INLINE_FUNCTION
    static BodyKokkos transform(const Body &b) {
      return BodyKokkos(b);
    }
  };

  template<>
  struct Transform<Body, BodyKokkos> {
    static constexpr bool is_identity = false;
    KOKKOS_INLINE_FUNCTION
    static Body transform(const BodyKokkos &bk) {
      return static_cast<Body>(bk);
    }
  };

}; // class FixRigidBaseKokkos

}    // namespace LAMMPS_NS

#endif // !LMP_FIX_RIGID_BASE_KOKKOS_H
