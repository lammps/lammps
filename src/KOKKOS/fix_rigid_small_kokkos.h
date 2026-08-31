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

#ifdef FIX_CLASS
// clang-format off
FixStyle(rigid/small/kk,FixRigidSmallKokkos<LMPDeviceType>);
FixStyle(rigid/small/kk/device,FixRigidSmallKokkos<LMPDeviceType>);
FixStyle(rigid/small/kk/host,FixRigidSmallKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_FIX_RIGID_SMALL_KOKKOS_H
#define LMP_FIX_RIGID_SMALL_KOKKOS_H

#include "Kokkos_Random.hpp"
#include "atom_masks.h"
#include "atom_vec_ellipsoid_kokkos.h"
#include "comm_brick_kokkos.h"
#include "fix_rigid_small.h"
#include "kokkos_base.h"
#ifdef LMP_KOKKOS_DEBUG_RNG
#include "rand_pool_wrap_kokkos.h"
#endif

struct TagInitialIntegrate {};
template <int SETXFLAG> struct TagSetXV {};
struct TagUpdateXGC {};

namespace LAMMPS_NS {

// Kokkos port of fix rigid/small.
//
// Design notes:
//  - Body setup (create_bodies / rendezvous_body) is intentionally kept on the
//    host: it runs only at setup, uses STL maps and an MPI rendezvous callback,
//    and has no per-timestep cost, so a device reimplementation would add
//    complexity with no runtime benefit.  The tied DualViews make the host
//    setup results directly visible to the device kernels.
//  - Atom exchange (migration) runs on whichever side CommBrickKokkos/AtomKokkos
//    select: the device path (pack/unpack_exchange_kokkos) when comm and sort
//    are both on device (GPU default, or "-pk kokkos comm device sort device"),
//    and the host FixRigidSmall::pack/unpack_exchange path otherwise (CPU
//    default), with the tied DualViews flushed to the host in pre_exchange() and
//    re-synced in pre_neighbor().  The exchange and the sort must run on the same
//    side -- a mixed configuration would let a host sort permute the per-atom
//    arrays out from under a device exchange -- so setup_pre_neighbor()/
//    setup_device_push() tie sort_device to the exchange side: if another fix
//    forces a host exchange (e.g. fix property/atom, which has no device
//    exchange), the fix follows it onto the host sort.  Forward/reverse comm run
//    on the device during the run.
template <class DeviceType> class FixRigidSmallKokkos : public FixRigidSmall, public KokkosBase {

 public:
  typedef EV_FLOAT value_type;
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixRigidSmallKokkos(class LAMMPS *, int, char **);
  ~FixRigidSmallKokkos();
  int setmask() override;
  void pre_exchange() override;

  // Verlet::setup() runs comm->exchange() on every run, and on the second and
  // later runs the device copies are the current ones, so the host exchange
  // needs the same flush the exchange during a run gets from pre_exchange().
  // Modify::setup_pre_exchange() is the only hook that runs before it.
  void setup_pre_exchange() override { pre_exchange(); }
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void post_force(int) override;
  void final_integrate() override;
  void write_restart_file(const char *) override;
  void pre_neighbor() override;

  void grow_arrays(int) override;
  void grow_body() override;
  void set_molecule(int, tagint, int, double *, double *, double *) override;
  void resample_momenta(int, int, class RanPark *, double) override;

  int pack_exchange_kokkos(const int &nsend, DAT::tdual_double_2d_lr &buf,
                           DAT::tdual_int_1d k_sendlist, DAT::tdual_int_1d k_copylist,
                           ExecutionSpace space) override;

  void unpack_exchange_kokkos(DAT::tdual_double_2d_lr &k_buf, DAT::tdual_int_1d &indices, int nrecv,
                              int, int, ExecutionSpace space) override;
  int pack_forward_comm_kokkos(int n, DAT::tdual_int_1d k_sendlist, DAT::tdual_double_1d &k_buf,
                               int pbc_flag, int *pbc) override;
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_double_1d &) override;

  int pack_reverse_comm_kokkos(int, int, DAT::tdual_double_1d &) override;
  void unpack_reverse_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_double_1d &) override;
  // reverse comm handled by host,
  // only happens when body and bodyown
  // are already on host

  void setup_pre_neighbor() override;
  bigint dof(int) override;
  void deform(int) override;
  void enforce2d() override;
  void zero_momentum() override;
  void zero_rotation() override;
  void *extract(const char *, int &) override;
  double extract_ke() override;
  double extract_erotational() override;
  double compute_scalar() override;
  void reset_atom2body() override;
  void image_shift() override;
  void sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) override;
  void sync_host_for_sort() override;
  void modified_host_for_sort() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagInitialIntegrate, const int) const;

  template <int SETXFLAG>
  KOKKOS_INLINE_FUNCTION void operator()(TagSetXV<SETXFLAG>, const int, EV_FLOAT &ev) const;
  KOKKOS_INLINE_FUNCTION
  void operator()(TagUpdateXGC, const int) const;

  void compute_forces_and_torques_kokkos();

  // nvcc's extended __host__ __device__ lambda extension requires the
  // enclosing function to have public access
  void setup_device_push();
  void apply_langevin_thermostat_kokkos();

 protected:
  void set_xv_kokkos(int);

  using ImageIntView1D = typename AT::t_imageint_1d;
  using TagIntView1D = typename AT::t_tagint_1d;
  using IntView1D = typename AT::t_int_1d;
  using View2D = typename AT::t_double_2d_lr;

  using Range1D = Kokkos::RangePolicy<DeviceType>;

  // Everything the host FixRigidSmall code reads when it rebuilds the bodies or
  // runs its setup: setup_bodies_static() reads radius/rmass/mass/type/image/
  // mask and the extended orientation arrays, compute_forces_and_torques() reads
  // f and (for extended bodies) torque, and set_xv/set_v read v and rmass.  The
  // per-step datamask_read is deliberately narrower than this, so the host-path
  // entry points have to ask for the difference explicitly.  Masks for arrays a
  // given atom style does not have are ignored by the atom-vec sync.
  static constexpr uint64_t host_rebuild_datamask =
      X_MASK | V_MASK | F_MASK | VIRIAL_MASK | TYPE_MASK | TAG_MASK | MASK_MASK |
      IMAGE_MASK | RMASS_MASK | RADIUS_MASK | TORQUE_MASK | OMEGA_MASK |
      ANGMOM_MASK | MU_MASK;

  void resize_body_views();
  void copy_body_host();
  void copy_body_device();
  void refresh_atom_views();
  KOKKOS_INLINE_FUNCTION
  void v_tally(EV_FLOAT &, int, double[6], double[3], double[3], double[3]) const;
  KOKKOS_INLINE_FUNCTION
  void v_tally(EV_FLOAT &, int, double[6]) const;

  // per-atom DualViews, tied to the FixRigidSmall host pointers via grow_kokkos
  DAT::tdual_int_1d k_bodyown;
  DAT::tdual_tagint_1d k_bodytag;
  DAT::tdual_int_1d k_atom2body;
  DAT::tdual_imageint_1d k_xcmimage;
  DAT::tdual_double_2d_lr k_displace, k_vatom, k_langextra;

  IntView1D d_bodyown;
  TagIntView1D d_bodytag;
  IntView1D d_atom2body;
  ImageIntView1D d_xcmimage;
  View2D d_displace, d_vatom, d_langextra;

  // extended-particle per-atom data (tied DualViews), only allocated when the
  // body contains finite-size particles (sphere/ellipsoid/line/tri/dipole).
  // eflags: per-atom bit flags; orient: orientation rel. to body (ellipsoid/tri/
  // quat, orientflag cols); dorient: dipole orientation rel. to body (3 cols).
  DAT::tdual_int_1d k_eflags;
  DAT::tdual_double_2d_lr k_orient, k_dorient;
  IntView1D d_eflags;
  View2D d_orient, d_dorient;

  // atom-style views written per-step for extended particles (set_xv/set_v):
  // omega (sphere/line), angmom (ellipsoid/tri), mu (dipole); torque is read
  // for the extended force accumulation.
  typename AT::t_kkfloat_1d_3 d_omega, d_angmom;
  typename AT::t_kkfloat_1d_4 d_mu;
  typename AT::t_kkacc_1d_3 d_torque;

  // ellipsoid bonus (shape read, quat written) reached on the device via the
  // ellipsoid atom-vec; accessed only through dynamic_cast + inline templates so
  // the KOKKOS-package fix still links when the ASPHERE package is absent.
  typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d d_bonus;
  IntView1D d_ellipsoid;

  // extended size added per body atom to the device-exchange payload, and the
  // atomKK datamask of per-atom views the extended set_xv/set_v kernels modify
  int extended_per_atom = 0;
  uint64_t extended_datamask = 0;

  // 1 once grow_kokkos owns the base per-atom pointers
  bool tied_initialized = false;

  // set while the host body rebuild runs in setup_pre_neighbor() on the 2nd and
  // later runs, so the pre_neighbor()/reset_atom2body()/image_shift() overrides
  // take their host fallback even though setupflag is already 1 by then
  bool setup_host_rebuild = false;

  // scratch reused across calls: allocating a device View (and, for the counter,
  // a fenced readback) inside every pack_exchange_kokkos / reset_atom2body call
  // costs an allocation per dimension per reneighbor
  Kokkos::View<int, DeviceType> d_exchange_count;
  IntView1D d_from_indices;

  // Per-swap ghost-body bookkeeping, rebuilt every reneighbor by the
  // BODY_SENDLIST comm and read by the INITIAL/FINAL/FULL_BODY packers on every
  // timestep.  Flat vectors scanned linearly rather than std::map: nswap is
  // small, and this keeps the per-step lookups out of a red-black tree (and off
  // its node allocations).  The send side is still keyed by the k_sendlist row
  // pointer because the packer signature carries no swap index; both vectors are
  // cleared wherever the swaps are rebuilt, so a stale key cannot survive.
  struct BodySend {
    int *key;
    int n;
    IntView1D list;
  };
  struct BodyRecv {
    int first;
    int n;
    int start;
    bool n_set;   // n and start are registered independently, as the two maps
  };              // this replaced were
  std::vector<BodySend> body_sends;
  std::vector<BodyRecv> body_recvs;

  BodySend &body_send(int *key) {
    for (auto &e : body_sends)
      if (e.key == key) return e;
    // reuse a retired slot: its list View is already allocated, and the whole
    // point of caching these is to avoid a device allocation per swap per
    // reneighbor
    for (auto &e : body_sends)
      if (e.key == nullptr) {
        e.key = key;
        e.n = 0;
        return e;
      }
    body_sends.push_back(BodySend{key, 0, IntView1D()});
    return body_sends.back();
  }
  // retire the keys at a reneighbor without freeing the cached list Views
  void retire_body_sends() {
    for (auto &e : body_sends) {
      e.key = nullptr;
      e.n = 0;
    }
  }
  BodyRecv &body_recv(int first) {
    for (auto &e : body_recvs)
      if (e.first == first) return e;
    body_recvs.push_back(BodyRecv{first, 0, 0, false});
    return body_recvs.back();
  }
  bool has_body_recv(int first) const {
    for (const auto &e : body_recvs)
      if (e.first == first) return e.n_set;
    return false;
  }

  IntView1D d_sendlist;
  typename AT::t_double_1d_um d_buf;
  int first;

  CommBrickKokkos *commKK;

  // body DualView (struct array); not tied to base `body` (different allocator),
  // bridged by copy_body_host()/copy_body_device()
  Kokkos::DualView<Body *, DeviceType> k_body;
  typename Kokkos::DualView<Body *, DeviceType>::t_dev d_body;

  double xbox, ybox, zbox, xprd, yprd, zprd, xy, xz, yz;
  typename AT::t_kkfloat_1d_3_lr d_x;
  typename AT::t_kkfloat_1d_3 d_v;
  typename AT::t_kkacc_1d_3 d_f;

  typename AT::t_kkfloat_1d d_rmass, d_mass;
  IntView1D d_type;

  // RNG pool for the Langevin thermostat.  By default this is the on-device
  // Kokkos XorShift pool (fast, but its stream differs from the host RanMars,
  // so Kokkos Langevin runs are statistically -- not bit-for-bit -- comparable
  // to the non-Kokkos style).  Building with -DLMP_KOKKOS_DEBUG_RNG swaps in a
  // wrapper around the host RanMars (the same instance the base class uses for
  // thread 0), so a Serial single-thread Kokkos run reproduces the non-Kokkos
  // Langevin trajectory bit-for-bit -- used for validation only.
#ifndef LMP_KOKKOS_DEBUG_RNG
  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;
#else
  RandPoolWrap rand_pool;
  typedef RandWrap rand_type;
#endif
};

KOKKOS_INLINE_FUNCTION
void copy_body(FixRigidSmall::Body *dest, FixRigidSmall::Body *src)
{
  memcpy(dest, src, sizeof(FixRigidSmall::Body));
}

}    // namespace LAMMPS_NS

#endif
#endif
