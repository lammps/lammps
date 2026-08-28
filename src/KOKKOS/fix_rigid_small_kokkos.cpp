// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_rigid_small_kokkos.h"

#include "atom_kokkos.h"
#include "atom_vec_kokkos.h"
#include "atom_masks.h"
#include "kokkos.h"
#include "modify.h"
#include "molecule.h"
#include "math_const.h"
#include "math_eigen.h"
#include "math_extra.h"
#include "rigid_const.h"
#include "memory.h"
#include "memory_kokkos.h"
#include "update.h"
#include "force.h"
#include "random_mars.h"
#include "domain.h"
#include "error.h"
#include "math_extra_kokkos.h"
#include "kokkos_few.h"
#include "domain_kokkos.h"

#include <cmath>
#include <cstring>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace RigidConst;

#define RVOUS 1   // 0 for irregular, 1 for all2all

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidSmallKokkos<DeviceType>::FixRigidSmallKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRigidSmall(lmp, narg, arg),
#ifdef LMP_KOKKOS_DEBUG_RNG
  rand_pool(12345 + comm->me, lmp)
#else
  rand_pool(12345 + comm->me)
#endif
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  commKK = (CommKokkos *) comm;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | V_MASK | VIRIAL_MASK | TYPE_MASK | TAG_MASK;
  datamask_modify = X_MASK | V_MASK | VIRIAL_MASK;

  // Variable-stride device exchange (pack/unpack_exchange_kokkos) writes a
  // per-sent-atom header (1 double, the offset of that atom's data) followed by
  // the atom's variable-length payload.  The largest payload is a body owner:
  // bodytag + xcmimage + displace[3] + vatom[6] + own-flag (12) + Body
  // (bodysize).  maxexchange must bound header + payload so CommKokkos sizes the
  // send buffer (nsend*maxexchange) large enough.  (The host path uses the base
  // class's own smaller variable packing; over-sizing the buffer is harmless.)
  maxexchange = 1 + 12 + bodysize;

  grow_arrays(atom->nmax);

  // For the device instantiation, declare device-exchange support immediately so
  // CommKokkos::exchange() does not permanently switch to exchange_comm_legacy=true
  // on the first (setup-time) call, and declare sort_device before the first
  // AtomKokkos::sort() (which Verlet::setup() runs before modify->setup() reaches
  // setup_device_push()).  init() then pairs the two -- see the reasoning there --
  // and may move both back to the host before either is first used.
  //
  // Note the base constructor has already run create_bodies() by this point, so
  // the per-atom body bookkeeping is real from here on: the setup-time exchange
  // is NOT a no-op for body data, and the host arrays create_bodies() wrote are
  // the authoritative copy until they are first pushed to the device.
  //
  // forward/reverse comm flags are dispatched per-call by CommKokkos and so can
  // safely be set later in setup_device_push().
  if (std::is_same<DeviceType,LMPDeviceType>::value) {
    exchange_comm_device = 1;
    sort_device = 1;
  }
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidSmallKokkos<DeviceType>::~FixRigidSmallKokkos()
{
  if (copymode) return;

  // No atomKK->sync() here.  Modify::delete_fix() deletes the fix before it
  // compacts modify->fix[], and AtomKokkos::sync() walks that list for the
  // property/atom fixes -- so a sync from this destructor reads a fix that has
  // already been freed, which AddressSanitizer reports as a heap-use-after-free
  // at the end of any run whose input has a fix property/atom.  Nothing reads
  // the atom arrays after this point, and no other KOKKOS fix syncs from its
  // destructor.

  // free the Kokkos-owned buffers and null the base pointers so the
  // FixRigidSmall destructor's memory->destroy() calls become no-ops
  // (these per-atom arrays are owned by grow_kokkos, not memory->grow)
  memoryKK->destroy_kokkos(k_bodyown, bodyown);
  memoryKK->destroy_kokkos(k_bodytag, bodytag);
  memoryKK->destroy_kokkos(k_atom2body, atom2body);
  memoryKK->destroy_kokkos(k_xcmimage, xcmimage);
  memoryKK->destroy_kokkos(k_displace, displace);
  memoryKK->destroy_kokkos(k_vatom, vatom);
  if (extended) {
    memoryKK->destroy_kokkos(k_eflags, eflags);
    if (orientflag) memoryKK->destroy_kokkos(k_orient, orient);
    if (dorientflag) memoryKK->destroy_kokkos(k_dorient, dorient);
  }

#ifdef LMP_KOKKOS_DEBUG_RNG
  if (langflag) rand_pool.destroy();
#endif
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::init()
{
  FixRigidSmall::init();
  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot yet use respa with Kokkos");

  // Keep the per-reneighbor atom exchange and the atom sort on the same side.
  // They permute the same per-atom bookkeeping: a device sort would reorder the
  // arrays a host exchange just rebuilt, and a host sort permutes the base
  // class' arrays while our tied DualViews stay device-modified, so
  // pre_neighbor()'s device branch then discards the permutation -- silently
  // corrupting the body<->atom mapping.
  //
  // This has to be decided here.  Verlet::setup() runs comm->exchange() and
  // atom->sort() before it ever reaches setup_pre_neighbor(), so a choice made
  // there (or later, in setup_device_push()) comes too late for the setup-time
  // pair.
  //
  // AtomKokkos::sort() drops to the legacy host sort -- permanently, for the
  // rest of the run -- when the atom style carries bonus data
  // (ellipsoid/line/tri/body), for atom_style hybrid, or when any fix that grows
  // per-atom arrays does not itself support the device sort.  None of those
  // consult exchange_comm_legacy, so a device exchange can otherwise end up
  // paired with a host sort with no warning; extended-particle rigid bodies hit
  // the bonus-data case routinely.  Follow the sort onto the host when that is
  // going to happen, and equally follow a forced host exchange onto the host
  // sort.  exchange_comm_on_host ("comm host") also routes through the
  // exchange_device<LMPHostType> Kokkos path, which leaves atom-vec arrays
  // host-modified; if the device sort then runs on a stale snapshot the
  // resulting atom ordering differs from what a true device exchange would give,
  // producing a trajectory that deviates from the correct physics.
  // The all-host combination is a supported path: the tied DualViews are
  // flushed in pre_exchange() and pushed back in pre_neighbor().
  if (std::is_same<DeviceType,LMPDeviceType>::value) {
    bool sort_on_host = lmp->kokkos->sort_legacy || atomKK->hybrid_flag ||
      atom->ellipsoid_flag || atom->line_flag || atom->tri_flag || atom->body_flag;
    for (int i = 0; i < atom->nextra_grow; i++)
      if (!modify->fix[atom->extra_grow[i]]->sort_device) sort_on_host = true;

    if (sort_on_host || lmp->kokkos->exchange_comm_legacy ||
        lmp->kokkos->exchange_comm_on_host) {
      exchange_comm_device = 0;
      sort_device = 0;
    }
  }

  // 2d is not supported: the base setup()/final_integrate() call enforce2d()
  // virtually, which dispatches to the device kk version operating on d_body
  // before it is populated from the host, leaving the host body's 2d
  // components un-zeroed for the host set_v.  Error rather than run silently
  // wrong until the 2d path is validated.
  if (domain->dimension == 2)
    error->all(FLERR,"fix rigid/small/kk does not yet support 2d systems");

  // set up the Langevin RNG pool
  if (langflag && random) {
#ifdef LMP_KOKKOS_DEBUG_RNG
    // validation build: wrap the base-class host RanMars (random_thr[0] = random)
    // so a Serial single-thread run draws the identical stream, in the identical
    // per-body order, as the non-Kokkos apply_langevin_thermostat().  Do not draw
    // from random here -- that would desync the stream from the non-Kokkos style.
    rand_pool.init(random, 12345 + comm->me);
#else
    // production build: seed the on-device XorShift pool from the host RNG (which
    // carries the user's langevin seed) so results are reproducible per input.
    int s = (int)(random->uniform() * 900000000.0) + 1;
    rand_pool = Kokkos::Random_XorShift64_Pool<DeviceType>(s + comm->me);
#endif
  }
}

/* ----------------------------------------------------------------------
   add PRE_EXCHANGE so the tied DualViews can be flushed to the host before
   the (host) atom migration in Comm::exchange
------------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::setmask()
{
  return FixRigidSmall::setmask() | PRE_EXCHANGE;
}

/* ----------------------------------------------------------------------
   flush body and per-atom state to the host before host atom exchange
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::pre_exchange()
{
  if (!setupflag) return;

  // Device exchange path: pack_exchange_kokkos reads from the device DualViews
  // which are always authoritative during the run; no host flush needed.
  // Only skip the flush when CommKokkos is actually using the device exchange
  // path (exchange_comm_legacy==false).  For OpenMP/Serial builds the global
  // default is exchange_comm_legacy=1 even though we set exchange_comm_device=1,
  // so CommBrick::exchange() is called and the base-class pack_exchange() reads
  // from the host arrays -- we must not skip the flush in that case.
  if (exchange_comm_device && !commKK->exchange_comm_legacy) return;

  // Host exchange path: flush device→host so the base-class
  // pack_exchange/copy_arrays/unpack_exchange see valid data.
  copy_body_host();
  k_bodytag.sync_host();
  k_bodyown.sync_host();
  k_atom2body.sync_host();
  k_xcmimage.sync_host();
  k_displace.sync_host();
  k_vatom.sync_host();
  // the extended-particle arrays migrate with the atoms too, and pre_neighbor()
  // marks them host-modified after the host exchange, so they have to be flushed
  // here as well or that modify_host() collides with an outstanding device claim
  if (extended) {
    k_eflags.sync_host();
    if (orientflag) k_orient.sync_host();
    if (dorientflag) k_dorient.sync_host();
  }
}

/* ----------------------------------------------------------------------
   setup static/dynamic properties of rigid bodies, using current atom info.
   if reinitflag is not set, do the initialization only once, b/c properties
   may not be re-computable especially if overlapping particles or bodies
   are inserted from mol template.
     do not do dynamic init if read body properties from inpfile. this
   is b/c the inpfile defines the static and dynamic properties and may not
   be computable if contain overlapping particles setup_bodies_static()
   reads inpfile itself.
     cannot do this until now, b/c requires comm->setup() to have setup stencil
   invoke pre_neighbor() to ensure body xcmimage flags are reset
     needed if Verlet::setup::pbc() has remapped/migrated atoms for 2nd run
     setup_bodies_static() invokes pre_neighbor itself
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::setup_pre_neighbor()
{
  // setup_bodies_static() rebuilds the bodies from the unwrapped atom
  // positions, so it also reads atom->image and atom->mask (and rmass for
  // finite-size particles), which are not in this fix's per-step datamask.
  atomKK->sync(Host, host_rebuild_datamask);

  // On the 2nd and later runs the host rebuild below (setup_bodies_static/
  // _dynamic) also re-derives the extended-particle state from the per-atom
  // orientation data -- omega/angmom/mu, covered by extended_datamask above, and
  // for ellipsoids the bonus quaternion.  set_xv_kokkos() wrote all of it on the
  // device, so pull the bonus array back too; the atom-style arrays are handled
  // by the sync above.
  if (setupflag && extended && atom->ellipsoid_flag) {
    auto *avecEllipKK = dynamic_cast<AtomVecEllipsoidKokkos *>(atom->style_match("ellipsoid"));
    if (avecEllipKK) avecEllipKK->k_bonus.sync_host();
  }

  // Re-assert the exchange/sort pairing chosen in init() (see the reasoning
  // there).  This only catches a side that changed after init() -- e.g.
  // CommKokkos deciding on a host exchange once it has seen every fix.  It
  // cannot protect the setup-time sort, which Verlet::setup() has already run by
  // the time we get here; that is why the decision is made in init().
  if (std::is_same<DeviceType,LMPDeviceType>::value) {
    if (lmp->kokkos->exchange_comm_legacy || lmp->kokkos->exchange_comm_on_host ||
        lmp->kokkos->sort_legacy) {
      exchange_comm_device = 0;
      sort_device = 0;
    }
  }

  // Retire the setup-time device claim, once, at its origin.  create_bodies()
  // ran in the base constructor and wrote the per-atom body bookkeeping through
  // the base class' raw pointers -- i.e. into the host mirrors of the tied
  // DualViews, without touching the modify flags.  Nothing has pushed that to
  // the device yet, so the device copies hold uninitialized data; the setup-time
  // exchange nevertheless ran pack/unpack_exchange_kokkos over them and marked
  // the views device-modified.  The host copies are therefore the authoritative
  // ones here, and that device claim is not.  Left in place it makes every later
  // reader defend itself: a sync_host() would pull the uninitialized copy over
  // the host arrays, and a modify_host() would trip the DualView
  // concurrent-modification guard.  Clearing it here -- after the exchange,
  // before the host build below and before the sort, dof() and
  // setup_device_push() that follow -- lets all of them use a plain
  // modify_host().  Only for the first run: once setupflag is set the device
  // copies carry real data and their claim must be honoured (the flush below
  // syncs it down instead).  No-op when host space == device space.
  if (!setupflag) {
    k_bodyown.clear_sync_state();
    k_bodytag.clear_sync_state();
    k_atom2body.clear_sync_state();
    k_xcmimage.clear_sync_state();
    k_displace.clear_sync_state();
    k_vatom.clear_sync_state();
    if (extended) {
      k_eflags.clear_sync_state();
      if (orientflag) k_orient.clear_sync_state();
      if (dorientflag) k_dorient.clear_sync_state();
    }
  }

  // On the 2nd and later runs reinitflag re-runs the host body build
  // (setup_bodies_static/dynamic), which reads bodytag/bodyown/atom2body/
  // xcmimage and rewrites body[].  After a preceding run that state is live on
  // the device, so flush it down first or the host rebuilds from stale
  // bookkeeping.  Note the base sets setupflag = 1 before returning, so latch it
  // here: on the first run the device views do not exist yet and must not be
  // touched.
  const int rebuild_on_host = setupflag;
  if (rebuild_on_host) {
    copy_body_host();
    k_bodytag.sync_host();
    k_bodyown.sync_host();
    k_atom2body.sync_host();
    k_xcmimage.sync_host();
    k_displace.sync_host();
    k_vatom.sync_host();     // setup_device_push() marks this host-modified below
    if (extended) {
      k_eflags.sync_host();
      if (orientflag) k_orient.sync_host();
      if (dorientflag) k_dorient.sync_host();
    }
  }

  // The host build also drives pre_neighbor/reset_atom2body/image_shift and
  // non-FORCE_TORQUE forward/reverse comms.  Run 1's setup_device_push() left
  // *_comm_device=1 and setupflag=1, which would route all of that to the device
  // packers (FORCE_TORQUE/body-sendlist only) and abort.  Restore the run-1
  // conditions around the base call: force the host fallback in the overrides
  // (setup_host_rebuild) and the host comm path; setup_device_push() re-enables
  // device comm for the run.
  const int saved_forward_comm_device = forward_comm_device;
  const int saved_reverse_comm_device = reverse_comm_device;
  forward_comm_device = 0;
  reverse_comm_device = 0;
  setup_host_rebuild = true;
  FixRigidSmall::setup_pre_neighbor();
  setup_host_rebuild = false;
  forward_comm_device = saved_forward_comm_device;
  reverse_comm_device = saved_reverse_comm_device;

  // the host rebuild rewrote body[] and the per-atom body arrays; push them
  // back to the device (setup() -> setup_device_push() does the rest)
  if (rebuild_on_host) {
    const int nbody_all = nlocal_body + nghost_body;
    if (k_body.view_device().extent_int(0) < nbody_all) {
      k_body.sync_host();
      k_body.resize(nbody_all > nmax_body ? nbody_all : nmax_body);
      // refill the host mirror resize() re-created empty, and clear the device
      // claim it left, before the host rebuild below writes body[] (see grow_body)
      k_body.sync_host();
    }
    copy_body_device();
    k_bodytag.modify_host();
    k_bodyown.modify_host();
    k_atom2body.modify_host();
    k_xcmimage.modify_host();
    k_displace.modify_host();
    k_vatom.modify_host();
    k_bodytag.template sync<DeviceType>();
    k_bodyown.template sync<DeviceType>();
    k_atom2body.template sync<DeviceType>();
    k_xcmimage.template sync<DeviceType>();
    k_displace.template sync<DeviceType>();
    k_vatom.template sync<DeviceType>();
    refresh_atom_views();
  }

  atomKK->modified(Host, datamask_modify);
  atomKK->sync(execution_space, datamask_read);

  // The host body arrays just (re)written above are authoritative until setup()
  // pushes them to the device: mark them so.  ModifyKokkos::setup() runs
  // compute temp -> dof() (which sync_host()s them) BEFORE the fix setup() loop,
  // and without this the sync would pull the device copy over the good host
  // arrays and segfault.  A plain modify_host() suffices here because no device
  // claim is outstanding: on the first run it was retired at the top of this
  // function, and on later runs the flush above synced it down (which clears it).
  k_bodyown.modify_host();
  k_bodytag.modify_host();
  k_atom2body.modify_host();
  k_xcmimage.modify_host();
  k_displace.modify_host();
  if (extended) {
    k_eflags.modify_host();
    if (orientflag) k_orient.modify_host();
    if (dorientflag) k_dorient.modify_host();
  }
}

/* ----------------------------------------------------------------------
   compute initial fcm and torque on bodies, also initial virial
   reset all particle velocities to be consistent with vcm and omega
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::setup(int vflag)
{
  atomKK->sync(Host, host_rebuild_datamask);

  // FixRigidSmall::setup() is host code throughout: compute_forces_and_torques()
  // does a commflag=FORCE_TORQUE reverse comm and then reads the host body[],
  // and a commflag=FINAL forward comm updates the ghost bodies that the
  // following set_v() reads -- again from the host body[].  On the 2nd and later
  // runs the device comm flags are still set from run 1's setup_device_push(),
  // which would route both comms to the device packers: they would accumulate
  // into d_body and update the ghost bodies on the device while the host body[]
  // that setup() actually reads stayed stale, giving wrong setup forces and
  // velocities for atoms in bodies that span processors.  Clear the flags for
  // the duration of the base call, exactly as dof() and setup_pre_neighbor() do.
  const int saved_forward_comm_device = forward_comm_device;
  const int saved_reverse_comm_device = reverse_comm_device;
  forward_comm_device = 0;
  reverse_comm_device = 0;

  FixRigidSmall::setup(vflag);

  forward_comm_device = saved_forward_comm_device;
  reverse_comm_device = saved_reverse_comm_device;

  atomKK->modified(Host, datamask_modify);

  // The base class did its work through the host atom arrays, so the modified()
  // above marks the host side of atom:x/v.  Reconcile to the fix's
  // execution_space (matching setup_pre_neighbor) before returning, so the
  // host data is pushed to the device and the host-modified flag cleared.
  // Otherwise ModifyKokkos::setup() then marks atom:x modified(Device, ...)
  // while the host is still flagged, which aborts: on a GPU build with a
  // host/device DualView, and on a mixed/single-precision build with the
  // legacy/Kokkos host TransformViews.  No-op when host space == device space
  // and no transform.
  atomKK->sync(execution_space, datamask_read);

  setup_device_push();
}

/* ----------------------------------------------------------------------
   push the host body/per-atom state populated by the (host) setup routines
   to the device, then enable device comm/sort for the run and seed the
   ghost-body sendlists.  Shared by the base setup() and the Nose-Hoover
   subclass setup() (which inserts its own scalar setup before this call).
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::setup_device_push()
{
  // FixRigidSmall::setup() populated the host per-atom arrays, which are the
  // host mirrors of the tied DualViews -> mark host-modified and push to device.
  // setup_pre_neighbor() always runs earlier in the same setup sequence and
  // leaves no device claim outstanding (on the first run it retires the one the
  // setup-time exchange made, and on later runs it syncs the device copy down),
  // so a plain modify_host() cannot trip the concurrent-modification guard.
  k_bodyown.modify_host();
  k_bodytag.modify_host();
  k_atom2body.modify_host();
  k_xcmimage.modify_host();
  k_displace.modify_host();
  k_vatom.modify_host();
  k_bodyown.template sync<DeviceType>();
  k_bodytag.template sync<DeviceType>();
  k_atom2body.template sync<DeviceType>();
  k_xcmimage.template sync<DeviceType>();
  k_displace.template sync<DeviceType>();
  k_vatom.template sync<DeviceType>();
  refresh_atom_views();

  // extended-particle arrays: setup_bodies_static() filled eflags/orient/dorient
  // on the host; push them to the device (they then migrate in the device
  // exchange) and size the exchange payload to carry them.
  extended_per_atom = 0;
  if (extended) {
    // The device set_xv/set_v kernels handle finite spheres, ellipsoids and
    // point dipoles.  Line and triangle particles (and a per-atom-quaternion
    // atom_style) are not handled, so error rather than silently skipping their
    // per-particle orientation update.  eflags was filled on the host by
    // setup_bodies_static().
    int has_line = 0, has_tri = 0;
    for (int i = 0; i < atom->nlocal; i++) {
      if (eflags[i] & LINE) has_line = 1;
      if (eflags[i] & TRIANGLE) has_tri = 1;
    }
    int loc[2] = {has_line, has_tri}, glob[2];
    MPI_Allreduce(loc, glob, 2, MPI_INT, MPI_MAX, world);
    if (glob[0])
      error->all(FLERR, "fix rigid/small/kk does not yet support line particles "
                 "in rigid bodies");
    if (glob[1])
      error->all(FLERR, "fix rigid/small/kk does not yet support triangle "
                 "particles in rigid bodies");
    if (atom->quat_flag)
      error->all(FLERR, "fix rigid/small/kk does not yet support a per-atom "
                 "quaternion (atom_style storing quat) in rigid bodies");
    if (atom->ellipsoid_flag &&
        dynamic_cast<AtomVecEllipsoidKokkos *>(atom->style_match("ellipsoid")) == nullptr)
      error->all(FLERR, "fix rigid/small/kk requires the Kokkos ellipsoid atom "
                 "style; run with the '-sf kk' suffix");

    k_eflags.modify_host();
    k_eflags.template sync<DeviceType>();
    d_eflags = k_eflags.template view<DeviceType>();
    extended_per_atom = 1;                      // eflags
    if (orientflag) {
      k_orient.modify_host();
      k_orient.template sync<DeviceType>();
      d_orient = k_orient.template view<DeviceType>();
      extended_per_atom += orientflag;          // orient cols
    }
    if (dorientflag) {
      k_dorient.modify_host();
      k_dorient.template sync<DeviceType>();
      d_dorient = k_dorient.template view<DeviceType>();
      extended_per_atom += 3;                    // dorient cols
    }
    // each migrating body atom now also carries extended_per_atom doubles
    maxexchange = 1 + 12 + extended_per_atom + bodysize;

    // per-step set_xv/set_v writes these atom-style views for extended particles
    extended_datamask = 0;
    if (atom->omega_flag)  extended_datamask |= OMEGA_MASK;
    if (atom->angmom_flag) extended_datamask |= ANGMOM_MASK;
    if (atom->mu_flag)     extended_datamask |= MU_MASK;

    // The host setup (FixRigidSmall::setup -> set_xv/set_v) already wrote the
    // per-atom omega/angmom/mu on the host; push them to the device so the first
    // device kernel does not see a stale host claim and trip the concurrent-
    // modification guard.
    atomKK->modified(Host, extended_datamask);
    atomKK->sync(execution_space, extended_datamask);
  }

  // size the body DualView and push host body[] to device
  k_body.sync_host();
  k_body.resize(nmax_body);
  copy_body_device();

  // size the per-body Langevin buffer to match (grow_body does this during the
  // run, but it is not necessarily called before the first post_force)
  if (langflag) {
    k_langextra.sync_host();
    k_langextra.resize(nmax_body, 6);
    k_langextra.sync_host();          // refill the mirror resize() re-created (see grow_body)
    k_langextra.modify_host();
    k_langextra.template sync<DeviceType>();
    d_langextra = k_langextra.template view<DeviceType>();
  }

  // Before this, the non-Kokkos (host) setup routines populated the host
  // arrays; from now on the run uses device communication and sorting.
  // Only the device instantiation enables these flags; CommKokkos / AtomKokkos
  // dispatch per-fix on them, so no global "comm device"/"sort device"
  // command-line override is needed and other fixes may use their own setting.
  // (The rigid/small/host instantiation keeps host comm/sort, execution_space
  // == Host, which CommKokkos already forces.)
  // exchange_comm_device is already 1 for the device instantiation (set in
  // the constructor) so CommKokkos doesn't fall back to the legacy path.
  if (std::is_same<DeviceType,LMPDeviceType>::value) {
    forward_comm_device = 1;
    reverse_comm_device = 1;

    // Re-assert the exchange/sort pairing decided in init() (see there): keep
    // both on the same side for the run, following whichever one has been forced
    // onto the host.  forward/reverse comm above are dispatched per call and are
    // unaffected by this.
    if (lmp->kokkos->exchange_comm_legacy || lmp->kokkos->exchange_comm_on_host ||
        lmp->kokkos->sort_legacy) {
      exchange_comm_device = 0;
      sort_device = 0;
    }
  }

  // The setup-time atom sort permutes the per-atom arrays (bodyown/bodytag) to
  // follow the atoms, but the body array is not permuted, so each body's ilocal
  // (its owner's local atom index) is left pointing at the pre-sort ordering.
  // d_bodyown stays correct, but the device exchange body-compaction relies on
  // body.ilocal to find a moved body's owner, so rebuild ilocal from the
  // authoritative per-atom bodyown here.  pack/unpack/sort_kokkos keep it
  // consistent thereafter.
  {
    auto d_body_l = k_body.view_device();
    auto d_bodyown_l = this->d_bodyown;
    int nlocal = atom->nlocal;
    Kokkos::parallel_for(
      "fix rigid/small rebuild body ilocal",
      Range1D(0, nlocal),
      KOKKOS_LAMBDA(const int i){
        if (d_bodyown_l(i) >= 0) d_body_l(d_bodyown_l(i)).ilocal = i;
      }
    );
  }

  // pre_neighbor isn't called again until necessary during the run,
  // need to set up the body sendlist and send the ghost bodies now
  nghost_body = 0;
  body_recvs.clear();
  retire_body_sends();

  commflag = BODY_SENDLIST;
  commKK->forward_comm_device<DeviceType>(this, 1);
  commflag = FULL_BODY;
  commKK->forward_comm_device<DeviceType>(this, bodysize);
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::pre_neighbor(){
  // setup_host_rebuild: setup_bodies_static() calls pre_neighbor() (virtually)
  // in the middle of re-deriving the bodies on the host, to remap each xcm into
  // the box and reset the atom xcmimage flags -- and then unwraps the atom
  // coordinates with those host xcmimage values to build the inertia tensor.
  // The device path below updates only d_xcmimage, so taking it here would
  // leave the host flags stale and the unwrap would place the atoms of any
  // body straddling a periodic boundary a full box length away, producing a
  // garbage inertia tensor ("Bad principal moments") on the second and later
  // runs.  Use the host path while the host setup owns the data.
  if (!setupflag || setup_host_rebuild) {
    FixRigidSmall::pre_neighbor();
    return;
  }
  Kokkos::Profiling::pushRegion("rigid/small pre_neighbor");

  nghost_body = 0;

  // True when CommKokkos is actually executing the device exchange path.
  // For OpenMP/Serial, exchange_comm_legacy defaults to 1 so CommBrick::exchange()
  // runs host pack/unpack even though exchange_comm_device=1; for GPU builds
  // with exchange_comm_legacy=0 the device path is used.
  bool using_device_exchange = exchange_comm_device && !commKK->exchange_comm_legacy;

  if (!using_device_exchange) {
    // Host exchange path: Comm::exchange migrated atoms and bodies on the host;
    // push the updated state back to the device.
    k_bodytag.modify_host();
    k_bodyown.modify_host();
    k_atom2body.modify_host();
    k_xcmimage.modify_host();
    k_displace.modify_host();
    k_vatom.modify_host();
    k_bodytag.template sync<DeviceType>();
    k_bodyown.template sync<DeviceType>();
    k_atom2body.template sync<DeviceType>();
    k_xcmimage.template sync<DeviceType>();
    k_displace.template sync<DeviceType>();
    k_vatom.template sync<DeviceType>();
    if (extended) {                 // base-class host exchange updated eflags/...
      k_eflags.modify_host();
      k_eflags.template sync<DeviceType>();
      if (orientflag) { k_orient.modify_host(); k_orient.template sync<DeviceType>(); }
      if (dorientflag) { k_dorient.modify_host(); k_dorient.template sync<DeviceType>(); }
      // The legacy exchange/sort/borders left the atom-style extended-particle
      // arrays (omega, angmom, mu) host-modified.  Retire that host claim now by
      // pushing them to the device.  On a build with strict host/device memory
      // separation (e.g. HIP without UVM), a subsequent modify_device() call for
      // any of these arrays would trip the concurrent-modification guard if the
      // host-modified flag is still set at that point.
      if (extended_datamask) atomKK->sync(execution_space, extended_datamask);
    }
    refresh_atom_views();
    copy_body_device();  // push host body[] → d_body
  } else {
    // Device exchange path: pack/unpack_exchange_kokkos already updated all
    // per-atom DualViews and d_body on the device.  Do NOT push stale host
    // data to the device — just mark the device side as authoritative and
    // refresh the d_* view aliases.
    //
    // Clear any outstanding host claim first.  comm->borders() runs between the
    // exchange and here and can grow the per-atom arrays, and grow_arrays() ends
    // by marking these host-modified; marking the device on top of that leaves
    // both flags set and Kokkos::abort()s on a build where the two memory spaces
    // differ.  Dropping the host claim is correct on this path: the device
    // copies are the ones the exchange just wrote.
    k_body.clear_sync_state();
    k_bodytag.clear_sync_state();
    k_bodyown.clear_sync_state();
    k_atom2body.clear_sync_state();
    k_xcmimage.clear_sync_state();
    k_displace.clear_sync_state();
    k_vatom.clear_sync_state();
    k_body.modify_device();
    k_bodytag.template modify<DeviceType>();
    k_bodyown.template modify<DeviceType>();
    k_atom2body.template modify<DeviceType>();
    k_xcmimage.template modify<DeviceType>();
    k_displace.template modify<DeviceType>();
    k_vatom.template modify<DeviceType>();
    if (extended) {
      k_eflags.clear_sync_state();
      k_eflags.template modify<DeviceType>();
      if (orientflag) { k_orient.clear_sync_state(); k_orient.template modify<DeviceType>(); }
      if (dorientflag) { k_dorient.clear_sync_state(); k_dorient.template modify<DeviceType>(); }
    }
    refresh_atom_views();
    // d_body is already the correct device view; no copy needed
    d_body = k_body.view_device();
  }

  int triclinic = domain->triclinic;
  int xperiodic = domain->xperiodic;
  int yperiodic = domain->yperiodic;
  int zperiodic = domain->zperiodic;
  Few<double,6> h(domain->h);
  Few<double,6> h_inv(domain->h_inv);
  Few<double,3> boxlo(domain->boxlo);
  Few<double,3> lo, hi, period;

  if (triclinic == 0) {
    lo = Few<double,3>(domain->boxlo);
    hi = Few<double,3>(domain->boxhi);
    period = Few<double,3>(domain->prd);
  } else {
    lo = Few<double,3>(domain->boxlo_lamda);
    hi = Few<double,3>(domain->boxhi_lamda);
    period = Few<double,3>(domain->prd_lamda);
  }

  auto d_body = this->d_body;

  Kokkos::parallel_for(
    "fix rigid/small remap",
    Range1D(0, nlocal_body),
    KOKKOS_LAMBDA(const int ibody) {
      Body &b = d_body(ibody);
      Few<double,3> x(b.xcm);
      Few<double,3> coord;
      if (triclinic == 0) {
        coord = Few<double,3>(x);
      } else {
        coord = DomainKokkos::x2lamda(boxlo, h_inv, x);
      }
      // remap the center of mass back into the box and step the body image
      // flags accordingly.  Both results have to be stored: remap() returns the
      // remapped coordinate and updates image through its reference argument.
      // Dropping either lets xcm integrate unbounded, so xcmimage grows until
      // the IMGMASK field overflows and set_xv() loses precision to two large
      // cancelling terms.
      coord = DomainKokkos::remap(lo, hi, period, xperiodic, yperiodic, zperiodic, coord, b.image);
      if (triclinic) {
        x = DomainKokkos::lamda2x(boxlo, h, coord);
        for (int i = 0; i < 3; i++) b.xcm[i] = x[i];
      } else {
        for (int i = 0; i < 3; i++) b.xcm[i] = coord[i];
      }
    }
  );

  body_recvs.clear();
  retire_body_sends();

  commflag = BODY_SENDLIST;
  commKK->forward_comm_device<DeviceType>(this, 1);
  commflag = FULL_BODY;
  commKK->forward_comm_device<DeviceType>(this, bodysize);
  reset_atom2body();

  image_shift();
  Kokkos::Profiling::popRegion();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::initial_integrate(int vflag)
{
  Kokkos::Profiling::pushRegion("rigid/small initial integrate");
  copymode = 1;
  Kokkos::parallel_for("fix rigid/small/kk initial_integrate",
    Kokkos::RangePolicy<DeviceType, TagInitialIntegrate>(0, nlocal_body), *this
  );
  copymode = 0;

  // virial setup before call to set_xv

  v_init(vflag);

  // forward communicate updated info of all bodies

  commflag = INITIAL;
  commKK->forward_comm_device<DeviceType>(this,29);

  // set coords/orient and velocity/rotation of atoms in rigid bodies

  set_xv_kokkos(1);
  Kokkos::Profiling::popRegion();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixRigidSmallKokkos<DeviceType>::operator()(TagInitialIntegrate, const int ibody) const {
  Body &b = d_body(ibody);

  // update vcm by 1/2 step

  double dtfm = dtf / b.mass;
  b.vcm[0] += dtfm * b.fcm[0];
  b.vcm[1] += dtfm * b.fcm[1];
  b.vcm[2] += dtfm * b.fcm[2];

  // update xcm by full step

  b.xcm[0] += dtv * b.vcm[0];
  b.xcm[1] += dtv * b.vcm[1];
  b.xcm[2] += dtv * b.vcm[2];

  // update angular momentum by 1/2 step

  b.angmom[0] += dtf * b.torque[0];
  b.angmom[1] += dtf * b.torque[1];
  b.angmom[2] += dtf * b.torque[2];

  // compute omega at 1/2 step from angmom at 1/2 step and current q
  // update quaternion a full step via Richardson iteration
  // returns new normalized quaternion, also updated omega at 1/2 step
  // update ex,ey,ez to reflect new quaternion

  // rotational quantities computed in KK_FLOAT; the quaternion stays double
  KK_FLOAT angmom[3] = {(KK_FLOAT)b.angmom[0],(KK_FLOAT)b.angmom[1],(KK_FLOAT)b.angmom[2]};
  KK_FLOAT ex[3] = {(KK_FLOAT)b.ex_space[0],(KK_FLOAT)b.ex_space[1],(KK_FLOAT)b.ex_space[2]};
  KK_FLOAT ey[3] = {(KK_FLOAT)b.ey_space[0],(KK_FLOAT)b.ey_space[1],(KK_FLOAT)b.ey_space[2]};
  KK_FLOAT ez[3] = {(KK_FLOAT)b.ez_space[0],(KK_FLOAT)b.ez_space[1],(KK_FLOAT)b.ez_space[2]};
  KK_FLOAT inertia[3] = {(KK_FLOAT)b.inertia[0],(KK_FLOAT)b.inertia[1],(KK_FLOAT)b.inertia[2]};
  KK_FLOAT omega[3];
  MathExtraKokkos::angmom_to_omega(angmom,ex,ey,ez,inertia,omega);
  MathExtraKokkos::richardson(b.quat,angmom,omega,inertia,dtq);
  MathExtraKokkos::q_to_exyz(b.quat,ex,ey,ez);
  b.ex_space[0]=ex[0]; b.ex_space[1]=ex[1]; b.ex_space[2]=ex[2];
  b.ey_space[0]=ey[0]; b.ey_space[1]=ey[1]; b.ey_space[2]=ey[2];
  b.ez_space[0]=ez[0]; b.ez_space[1]=ez[1]; b.ez_space[2]=ez[2];
  b.omega[0]=omega[0]; b.omega[1]=omega[1]; b.omega[2]=omega[2];
}

/* ----------------------------------------------------------------------
   apply Langevin thermostat to all 6 DOF of rigid bodies I own
   unlike fix langevin, this stores extra force in extra arrays,
     which are added in when a new fcm/torque are calculated
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::apply_langevin_thermostat_kokkos()
{
  // Langevin forces/torques are computed per-body on the device.  k_langextra
  // is sized to nmax_body by grow_body(), so it always covers nlocal_body.

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  const double t_target = t_start + delta * (t_stop - t_start);
  const double tsqrt = sqrt(t_target);
  const double boltz = force->boltz;
  const double dt = update->dt;
  const double mvv2e = force->mvv2e;
  const double ftm2v = force->ftm2v;
  const double tp = t_period;

  auto d_body = this->d_body;
  auto l_d_langextra = d_langextra;
  auto l_rand_pool = rand_pool;

  Kokkos::parallel_for("fix rigid/small langevin",
    Range1D(0, nlocal_body),
    KOKKOS_LAMBDA(const int ibody) {
      rand_type rand_gen = l_rand_pool.get_state();
      Body &b = d_body(ibody);

      double gamma1 = -b.mass / tp / ftm2v;
      double gamma2 = sqrt(b.mass) * tsqrt * sqrt(24.0*boltz/tp/dt/mvv2e) / ftm2v;
      l_d_langextra(ibody,0) = gamma1*b.vcm[0] + gamma2*(rand_gen.drand()-0.5);
      l_d_langextra(ibody,1) = gamma1*b.vcm[1] + gamma2*(rand_gen.drand()-0.5);
      l_d_langextra(ibody,2) = gamma1*b.vcm[2] + gamma2*(rand_gen.drand()-0.5);

      gamma1 = -1.0 / tp / ftm2v;
      gamma2 = tsqrt * sqrt(24.0*boltz/tp/dt/mvv2e) / ftm2v;

      // convert omega from space frame to body frame, compute body-frame torque
      // (rotational quantities in KK_FLOAT)

      KK_FLOAT ex[3] = {(KK_FLOAT)b.ex_space[0],(KK_FLOAT)b.ex_space[1],(KK_FLOAT)b.ex_space[2]};
      KK_FLOAT ey[3] = {(KK_FLOAT)b.ey_space[0],(KK_FLOAT)b.ey_space[1],(KK_FLOAT)b.ey_space[2]};
      KK_FLOAT ez[3] = {(KK_FLOAT)b.ez_space[0],(KK_FLOAT)b.ez_space[1],(KK_FLOAT)b.ez_space[2]};
      KK_FLOAT omega[3] = {(KK_FLOAT)b.omega[0],(KK_FLOAT)b.omega[1],(KK_FLOAT)b.omega[2]};
      KK_FLOAT wbody[3], tbody[3], lang[3];
      MathExtraKokkos::transpose_matvec(ex,ey,ez,omega,wbody);
      tbody[0] = b.inertia[0]*gamma1*wbody[0] + sqrt(b.inertia[0])*gamma2*(rand_gen.drand()-0.5);
      tbody[1] = b.inertia[1]*gamma1*wbody[1] + sqrt(b.inertia[1])*gamma2*(rand_gen.drand()-0.5);
      tbody[2] = b.inertia[2]*gamma1*wbody[2] + sqrt(b.inertia[2])*gamma2*(rand_gen.drand()-0.5);

      // convert langevin torque from body frame back to space frame

      MathExtraKokkos::matvec(ex,ey,ez,tbody,lang);
      l_d_langextra(ibody,3) = lang[0];
      l_d_langextra(ibody,4) = lang[1];
      l_d_langextra(ibody,5) = lang[2];

      l_rand_pool.free_state(rand_gen);
    });
  k_langextra.modify_device();
}

/* ----------------------------------------------------------------------
   called from FixEnforce post_force() for 2d problems
   zero all body values that should be zero for 2d model
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::enforce2d()
{
  // NOTE: the langextra columns zeroed below follow FixRigid::enforce2d();
  // FixRigidSmall::enforce2d() deliberately leaves langextra alone.  The
  // divergence is currently unreachable -- init() rejects 2d systems -- but it
  // has to be reconciled with the CPU style before that restriction is lifted,
  // or the Langevin buffer will differ between the two.

  auto d_body = this->d_body;
  auto langflag = this->langflag && (langextra!=nullptr);
  auto d_langextra = this->d_langextra;

  Kokkos::parallel_for(
      "fix rigid/small zero z components",
      Range1D(0, nlocal_body),
      KOKKOS_LAMBDA (const int ibody) {
        Body &b = d_body(ibody);
        b.xcm[2] = 0.0;
        b.vcm[2] = 0.0;
        b.fcm[2] = 0.0;
        b.xgc[2] = 0.0;
        b.torque[0] = 0.0;
        b.torque[1] = 0.0;
        b.angmom[0] = 0.0;
        b.angmom[1] = 0.0;
        b.omega[0] = 0.0;
        b.omega[1] = 0.0;
        if(langflag) {
          d_langextra(ibody,2) = 0.0;
          d_langextra(ibody,3) = 0.0;
          d_langextra(ibody,4) = 0.0;
        }
      }
  );
  if (langflag) k_langextra.modify_device();
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::post_force(int /*vflag*/)
{
  Kokkos::Profiling::pushRegion("rigid/small post force");
  if (langflag) apply_langevin_thermostat_kokkos();
  if (earlyflag) compute_forces_and_torques_kokkos();
  Kokkos::Profiling::popRegion();
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::compute_forces_and_torques_kokkos()
{
  // sum over atoms to get force and torque on rigid body

  Kokkos::Profiling::pushRegion("rigix/small forces and torques");
  atomKK->sync(execution_space, X_MASK | F_MASK);
  auto d_x = atomKK->k_x.view<DeviceType>();
  auto d_f = atomKK->k_f.view<DeviceType>();
  int nlocal = atom->nlocal;

  auto d_body = this->d_body;

  // extended particles that carry their own torque (sphere/ellipsoid/line/tri)
  // contribute it to the body torque; fetch atom torque + eflags on device
  const bool ext_torque = extended;
  if (ext_torque) {
    atomKK->sync(execution_space, TORQUE_MASK);
    d_torque = atomKK->k_torque.template view<DeviceType>();
    d_eflags = k_eflags.template view<DeviceType>();
  }
  auto d_torque_l = d_torque;
  auto d_eflags_l = d_eflags;

  Kokkos::parallel_for(
    "fix rigid/small zero fcm&tcm",
    Range1D(0,nlocal_body+nghost_body),
    KOKKOS_LAMBDA (const int ibody) {
      Body &b = d_body(ibody);
      b.fcm[0] = b.fcm[1] = b.fcm[2] = 0.0;
      b.torque[0] = b.torque[1] = b.torque[2] = 0.0;
    }
  );

  auto d_xcmimage = this->d_xcmimage;
  auto d_atom2body = this->d_atom2body;

  auto prd = Few<double,3>(domain->prd);
  auto h = Few<double,6>(domain->h);
  auto triclinic = domain->triclinic;

  Kokkos::parallel_for(
    "fix rigid/small add tcm&fcm",
    Range1D(0,nlocal),
    KOKKOS_LAMBDA (const int i) {
      if (d_atom2body(i) < 0) return;
      Body &b = d_body(d_atom2body(i));

      Kokkos::atomic_add(&b.fcm[0], d_f(i,0));
      Kokkos::atomic_add(&b.fcm[1], d_f(i,1));
      Kokkos::atomic_add(&b.fcm[2], d_f(i,2));

      Few<double,3> x_i;
      x_i[0] = d_x(i,0);
      x_i[1] = d_x(i,1);
      x_i[2] = d_x(i,2);
      auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,d_xcmimage(i));
      double dx = unwrap[0] - b.xcm[0];
      double dy = unwrap[1] - b.xcm[1];
      double dz = unwrap[2] - b.xcm[2];

      Kokkos::atomic_add(&b.torque[0], dy*d_f(i,2) - dz*d_f(i,1));
      Kokkos::atomic_add(&b.torque[1], dz*d_f(i,0) - dx*d_f(i,2));
      Kokkos::atomic_add(&b.torque[2], dx*d_f(i,1) - dy*d_f(i,0));

      // extended particle's own torque (e.g. granular/dipole) adds to the body
      if (ext_torque && (d_eflags_l(i) & RigidConst::TORQUE)) {
        Kokkos::atomic_add(&b.torque[0], d_torque_l(i,0));
        Kokkos::atomic_add(&b.torque[1], d_torque_l(i,1));
        Kokkos::atomic_add(&b.torque[2], d_torque_l(i,2));
      }
    }
  );


  // reverse communicate fcm, torque of all bodies

  Kokkos::Profiling::pushRegion("reverse communicate");
  commflag = FORCE_TORQUE;
  commKK->reverse_comm_device<DeviceType>(this,6);
  Kokkos::Profiling::popRegion();

  // include Langevin thermostat forces and torques (computed on the device)

  if (langflag) {
    Kokkos::Profiling::pushRegion("rigid/small langevin");
    k_langextra.template sync<DeviceType>();
    auto l_d_langextra = d_langextra;
    Kokkos::parallel_for(
      "fix rigid/small add langevin",
      Range1D(0, nlocal_body),
      KOKKOS_LAMBDA (const int ibody) {
        Body &b = d_body(ibody);
        b.fcm[0] += l_d_langextra(ibody,0);
        b.fcm[1] += l_d_langextra(ibody,1);
        b.fcm[2] += l_d_langextra(ibody,2);
        b.torque[0] += l_d_langextra(ibody,3);
        b.torque[1] += l_d_langextra(ibody,4);
        b.torque[2] += l_d_langextra(ibody,5);
      }
    );
    Kokkos::Profiling::popRegion();
  }

  // add gravity force to COM of each body (matches FixRigidSmall::post_force).
  // gvec points into the gravity fix and is refreshed on the host each step, so
  // read its three components on the host and apply them in a device kernel.
  if (id_gravity) {
    const double gx = gvec[0], gy = gvec[1], gz = gvec[2];
    Kokkos::parallel_for(
      "fix rigid/small gravity",
      Range1D(0, nlocal_body),
      KOKKOS_LAMBDA(const int ibody) {
        Body &b = d_body(ibody);
        b.fcm[0] += gx*b.mass;
        b.fcm[1] += gy*b.mass;
        b.fcm[2] += gz*b.mass;
      }
    );
  }

  Kokkos::Profiling::popRegion();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::final_integrate()
{
  Kokkos::Profiling::pushRegion("rigid/small final integrate");

  if (!earlyflag) compute_forces_and_torques_kokkos();

  // update vcm and angmom, recompute omega

  auto nlocal_body = this->nlocal_body;
  auto d_body = this->d_body;
  auto dtf = this->dtf;

  Kokkos::parallel_for(
    "fix rigid/small update vcm+angmom+omega",
    Range1D(0, nlocal_body),
    KOKKOS_LAMBDA (const int ibody) {
      Body &b = d_body(ibody);

      // update vcm by 1/2 step

      double dtfm = dtf / b.mass;
      b.vcm[0] += dtfm * b.fcm[0];
      b.vcm[1] += dtfm * b.fcm[1];
      b.vcm[2] += dtfm * b.fcm[2];

      // update angular momentum by 1/2 step

      b.angmom[0] += dtf * b.torque[0];
      b.angmom[1] += dtf * b.torque[1];
      b.angmom[2] += dtf * b.torque[2];

      KK_FLOAT angmom[3] = {(KK_FLOAT)b.angmom[0],(KK_FLOAT)b.angmom[1],(KK_FLOAT)b.angmom[2]};
      KK_FLOAT ex[3] = {(KK_FLOAT)b.ex_space[0],(KK_FLOAT)b.ex_space[1],(KK_FLOAT)b.ex_space[2]};
      KK_FLOAT ey[3] = {(KK_FLOAT)b.ey_space[0],(KK_FLOAT)b.ey_space[1],(KK_FLOAT)b.ey_space[2]};
      KK_FLOAT ez[3] = {(KK_FLOAT)b.ez_space[0],(KK_FLOAT)b.ez_space[1],(KK_FLOAT)b.ez_space[2]};
      KK_FLOAT inertia[3] = {(KK_FLOAT)b.inertia[0],(KK_FLOAT)b.inertia[1],(KK_FLOAT)b.inertia[2]};
      KK_FLOAT omega[3];
      MathExtraKokkos::angmom_to_omega(angmom,ex,ey,ez,inertia,omega);
      b.omega[0]=omega[0]; b.omega[1]=omega[1]; b.omega[2]=omega[2];
    }
  );

  // forward communicate updated info of all bodies

  commflag = FINAL;
  commKK->forward_comm_device<DeviceType>(this,10);

  // set velocity/rotation of atoms in rigid bodies
  // virial is already setup from initial_integrate

  set_xv_kokkos(0);
  Kokkos::Profiling::popRegion();
}

/* ----------------------------------------------------------------------
   count # of DOF removed by rigid bodies for atoms in igroup
   return total count of DOF
------------------------------------------------------------------------- */

template<class DeviceType>
bigint FixRigidSmallKokkos<DeviceType>::dof(int tgroup)
{
  // FixRigidSmall::dof() returns 0 before setup; otherwise it reads the per-atom
  // atom2body/bodyown/bodytag (and eflags for extended particles) and atom->mask
  // on the host, and does a DOF reverse_comm.  On the device path a reneighbor
  // may have left those host mirrors stale, so flush them when they are actually
  // read (setupflag true); sync_host() is a no-op when the host copy is current.
  // It also reads body[ibody].inertia to detect linear bodies, so the body
  // array has to come down from the device as well.
  if (setupflag) {
    copy_body_host();
    k_atom2body.sync_host();
    k_bodyown.sync_host();
    k_bodytag.sync_host();
    if (extended) k_eflags.sync_host();
    atomKK->sync(Host, MASK_MASK);
  }

  // Route the DOF reverse_comm through the host comm path: with
  // reverse_comm_device set for the run, CommKokkos would send it to
  // pack_reverse_comm_kokkos, which only handles FORCE_TORQUE.  Clearing the
  // flag around the base call makes CommKokkos use CommBrick::reverse_comm and
  // the base host packer (which handles commflag == DOF).  dof() is an
  // infrequent diagnostic query (compute temp), so the host path is fine.
  const int saved_reverse_comm_device = reverse_comm_device;
  reverse_comm_device = 0;
  bigint ndof = FixRigidSmall::dof(tgroup);
  reverse_comm_device = saved_reverse_comm_device;
  return ndof;
}

/* ----------------------------------------------------------------------
   adjust xcm of each rigid body due to box deformation
   called by various fixes that change box size/shape
   flag = 0/1 means map from box to lamda coords or vice versa
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::deform(int flag)
{
  // Barostatted runs call remap() twice per step and each remap() calls this
  // once with flag=0 and once with flag=1, so a host round-trip here would move
  // the whole body array off and back onto the device eight times per timestep
  // to apply a pure affine transform.  Do it in place on the device.
  if (!setupflag) {
    copy_body_host();
    if (flag == 0)
      for (int ibody = 0; ibody < nlocal_body; ibody++)
        domain->x2lamda(body[ibody].xcm,body[ibody].xcm);
    else
      for (int ibody = 0; ibody < nlocal_body; ibody++)
        domain->lamda2x(body[ibody].xcm,body[ibody].xcm);
    copy_body_device();
    return;
  }

  Few<double,3> boxlo;
  boxlo[0] = domain->boxlo[0]; boxlo[1] = domain->boxlo[1]; boxlo[2] = domain->boxlo[2];
  Few<double,6> h, h_inv;
  for (int i = 0; i < 6; i++) { h[i] = domain->h[i]; h_inv[i] = domain->h_inv[i]; }

  auto d_body = this->d_body;
  Kokkos::parallel_for(
    "fix rigid/small deform",
    Range1D(0, nlocal_body),
    KOKKOS_LAMBDA(const int ibody) {
      Body &b = d_body(ibody);
      Few<double,3> x(b.xcm);
      Few<double,3> out = (flag == 0) ? DomainKokkos::x2lamda(boxlo, h_inv, x)
                                      : DomainKokkos::lamda2x(boxlo, h, x);
      for (int i = 0; i < 3; i++) b.xcm[i] = out[i];
    }
  );
}

/* ----------------------------------------------------------------------
   set space-frame coords and velocity of each atom in each rigid body
   set orientation and rotation of extended particles
   x = Q displace + Xcm, mapped back to periodic box
   v = Vcm + (W cross (x - Xcm))
   setxflag = whether to update positions
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::set_xv_kokkos(int setxflag)
{
  Kokkos::Profiling::pushRegion("rigid/small set_xv");

  // The per-atom centroid virial (cvatom) is filled on the host by the base
  // Fix::v_tally and has no tied device buffer, so the device set_xv kernel
  // cannot accumulate it without dereferencing a host pointer on the device.
  // Refuse rather than crash on GPU / silently corrupt.  The global and the
  // ordinary per-atom virial (d_vatom) are fully supported.
  if (cvflag_atom)
    error->all(FLERR,"fix rigid/small/kk does not support the per-atom centroid "
               "stress (e.g. compute centroid/stress/atom); use the non-Kokkos style");

  this->xprd = domain->xprd;
  this->yprd = domain->yprd;
  this->zprd = domain->zprd;
  this->xy = domain->xy;
  this->xz = domain->xz;
  this->yz = domain->yz;

  atomKK->sync(execution_space, datamask_read);
  if (extended) atomKK->sync(execution_space, extended_datamask);
  d_x = atomKK->k_x.view<DeviceType>();
  d_v = atomKK->k_v.view<DeviceType>();
  d_f = atomKK->k_f.view<DeviceType>();
  d_rmass = atomKK->k_rmass.view<DeviceType>();
  d_mass = atomKK->k_mass.view<DeviceType>();
  d_type = atomKK->k_type.view<DeviceType>();
  int nlocal = atom->nlocal;

  // extended particles: the kernel sets each finite-size particle's
  // omega/angmom (and dipole mu) from the body's rotational state
  AtomVecEllipsoidKokkos *avecEllipKK = nullptr;
  if (extended) {
    if (atom->omega_flag)  d_omega  = atomKK->k_omega.template view<DeviceType>();
    if (atom->angmom_flag) d_angmom = atomKK->k_angmom.template view<DeviceType>();
    if (atom->mu_flag)     d_mu     = atomKK->k_mu.template view<DeviceType>();
    if (atom->ellipsoid_flag) {
      avecEllipKK = dynamic_cast<AtomVecEllipsoidKokkos *>(atom->style_match("ellipsoid"));
      if (avecEllipKK) {
        avecEllipKK->k_bonus.template sync<DeviceType>();
        d_bonus = avecEllipKK->k_bonus.template view<DeviceType>();
        d_ellipsoid = atomKK->k_ellipsoid.template view<DeviceType>();
        d_rmass = atomKK->k_rmass.template view<DeviceType>();
      }
    }
    d_eflags = k_eflags.template view<DeviceType>();
  }

  // The kernel also reads atom-style arrays that datamask_read does not cover:
  // rmass (constraint-force virial and ellipsoid inertia), the per-atom
  // ellipsoid index, and mu[3] (the dipole magnitude).  Sync them to the device
  // or they are read stale -- or never initialised -- on a GPU (host==device
  // hides this on CPU builds).
  int extra_read = 0;
  if (atom->rmass) extra_read |= RMASS_MASK;
  if (extended) {
    if (atom->mu_flag) extra_read |= MU_MASK;
    if (atom->ellipsoid_flag && avecEllipKK) extra_read |= ELLIPSOID_MASK;
  }
  if (extra_read) atomKK->sync(execution_space, extra_read);
  // per-type mass is not a datamask-tracked per-atom array; sync the DualView
  // directly (no-op once it is device-current, as it is set only at setup)
  if (!atom->rmass) atomKK->k_mass.template sync<DeviceType>();

  EV_FLOAT ev;
  // The constraint virial is tallied in two halves per timestep: set_xv() from
  // initial_integrate() and set_v() from final_integrate().  Only zero the
  // accumulator at the start of the step (setxflag, i.e. the set_xv call), or
  // the set_v pass would discard the first half.
  if (vflag_atom && setxflag) {
    Kokkos::deep_copy(d_vatom, 0.0);
    k_vatom.modify_device();
  }
  copymode = 1;
  if (setxflag) {
    Kokkos::parallel_reduce("fix rigid/small setxv", Kokkos::RangePolicy<DeviceType, TagSetXV<1>>(0, nlocal), *this, ev);
    Kokkos::parallel_for("fix rigid/small update xgc", Kokkos::RangePolicy<DeviceType, TagUpdateXGC>(0, nlocal_body+nghost_body), *this);
  } else {
    Kokkos::parallel_reduce("fix rigid/small setv", Kokkos::RangePolicy<DeviceType, TagSetXV<0>>(0, nlocal), *this, ev);
  }
  copymode = 0;
  if(evflag){
    if(vflag_global){
      for(int i = 0; i < 6; i++){
        virial[i] += ev.v[i];
      }
    }
    if(vflag_atom){
      // d_vatom aliases the host vatom[] buffer (tied DualView); flush to host
      k_vatom.modify_device();
      k_vatom.sync_host();
    }
  }
  atomKK->modified(execution_space, datamask_modify);
  if (extended) {
    atomKK->modified(execution_space, extended_datamask);
    if (avecEllipKK) avecEllipKK->k_bonus.template modify<DeviceType>();
  }
  Kokkos::Profiling::popRegion();
}

template<class DeviceType>
template<int SETXFLAG>
KOKKOS_INLINE_FUNCTION
void FixRigidSmallKokkos<DeviceType>::operator()(TagSetXV<SETXFLAG>, const int i, EV_FLOAT &ev) const{
  if (d_atom2body(i) < 0) return;
  Body &b = d_body(d_atom2body(i));

  double xbox = (d_xcmimage(i) & IMGMASK) - IMGMAX;
  double ybox = (d_xcmimage(i) >> IMGBITS & IMGMASK) - IMGMAX;
  double zbox = (d_xcmimage(i) >> IMG2BITS) - IMGMAX;

  double x0, x1, x2, v0, v1, v2;

  // save old positions and velocities for virial

  if (evflag) {
    if (triclinic == 0) {
      x0 = d_x(i,0) + xbox*xprd;
      x1 = d_x(i,1) + ybox*yprd;
      x2 = d_x(i,2) + zbox*zprd;
    } else {
      x0 = d_x(i,0) + xbox*xprd + ybox*xy + zbox*xz;
      x1 = d_x(i,1) + ybox*yprd + zbox*yz;
      x2 = d_x(i,2) + zbox*zprd;
    }
    v0 = d_v(i,0);
    v1 = d_v(i,1);
    v2 = d_v(i,2);
  }

  // x = displacement from center-of-mass, based on body orientation
  // v = vcm + omega around center-of-mass

  KK_FLOAT ex[3] = {(KK_FLOAT)b.ex_space[0],(KK_FLOAT)b.ex_space[1],(KK_FLOAT)b.ex_space[2]};
  KK_FLOAT ey[3] = {(KK_FLOAT)b.ey_space[0],(KK_FLOAT)b.ey_space[1],(KK_FLOAT)b.ey_space[2]};
  KK_FLOAT ez[3] = {(KK_FLOAT)b.ez_space[0],(KK_FLOAT)b.ez_space[1],(KK_FLOAT)b.ez_space[2]};
  KK_FLOAT displace[3] = {(KK_FLOAT)d_displace(i,0),(KK_FLOAT)d_displace(i,1),(KK_FLOAT)d_displace(i,2)};
  KK_FLOAT delta[3];
  MathExtraKokkos::matvec(ex,ey,ez,displace,delta);

  d_v(i,0) = b.omega[1]*delta[2] - b.omega[2]*delta[1] + b.vcm[0];
  d_v(i,1) = b.omega[2]*delta[0] - b.omega[0]*delta[2] + b.vcm[1];
  d_v(i,2) = b.omega[0]*delta[1] - b.omega[1]*delta[0] + b.vcm[2];

  // add center of mass to displacement
  // map back into periodic box via xbox,ybox,zbox
  // for triclinic, add in box tilt factors as well

  if constexpr(SETXFLAG) {
    if (triclinic == 0) {
      d_x(i,0) = delta[0] + b.xcm[0] - xbox*xprd;
      d_x(i,1) = delta[1] + b.xcm[1] - ybox*yprd;
      d_x(i,2) = delta[2] + b.xcm[2] - zbox*zprd;
    } else {
      d_x(i,0) = delta[0] + b.xcm[0] - xbox*xprd - ybox*xy - zbox*xz;
      d_x(i,1) = delta[1] + b.xcm[1] - ybox*yprd - zbox*yz;
      d_x(i,2) = delta[2] + b.xcm[2] - zbox*zprd;
    }
  }

  // virial = unwrapped coords dotted into body constraint force
  // body constraint force = implied force due to v change minus f external
  // assume f does not include forces internal to body
  // 1/2 factor b/c final_integrate contributes other half
  // assume per-atom contribution is due to constraint force on that atom

  if (evflag) {
    double massone;
    if (d_rmass.data()) massone = d_rmass(i);
    else massone = d_mass(d_type(i));
    double fc0 = massone*(d_v(i,0) - v0)/dtf - d_f(i,0);
    double fc1 = massone*(d_v(i,1) - v1)/dtf - d_f(i,1);
    double fc2 = massone*(d_v(i,2) - v2)/dtf - d_f(i,2);

    double vr[6];
    vr[0] = 0.5*x0*fc0;
    vr[1] = 0.5*x1*fc1;
    vr[2] = 0.5*x2*fc2;
    vr[3] = 0.5*x0*fc1;
    vr[4] = 0.5*x0*fc2;
    vr[5] = 0.5*x1*fc2;

    double rlist[3] = {x0, x1, x2};
    double flist[3] = {0.5*fc0, 0.5*fc1, 0.5*fc2};
    v_tally(ev,i,vr,rlist,flist,b.xgc);
  }

  // extended particles: set per-particle rotational state from the body
  // (mirrors FixRigidSmall::set_xv / set_v extended block)
  if (extended) {
    const int ef = d_eflags(i);
    if (ef & RigidConst::SPHERE) {
      d_omega(i,0) = b.omega[0];
      d_omega(i,1) = b.omega[1];
      d_omega(i,2) = b.omega[2];
    }
    // ellipsoid: compose body orientation with the particle's body-frame
    // orientation to get its space-frame quaternion, and set its angmom
    if (ef & RigidConst::ELLIPSOID) {
      const int e = d_ellipsoid(i);
      double *shape = d_bonus(e).shape;
      double *quatatom = d_bonus(e).quat;
      double orient_i[4] = {d_orient(i,0), d_orient(i,1), d_orient(i,2), d_orient(i,3)};
      MathExtraKokkos::quatquat(b.quat, orient_i, quatatom);
      MathExtraKokkos::qnormalize(quatatom);
      const KK_FLOAT rm = d_rmass(i);
      KK_FLOAT ione[3];
      ione[0] = RigidConst::EINERTIA*rm * (shape[1]*shape[1] + shape[2]*shape[2]);
      ione[1] = RigidConst::EINERTIA*rm * (shape[0]*shape[0] + shape[2]*shape[2]);
      ione[2] = RigidConst::EINERTIA*rm * (shape[0]*shape[0] + shape[1]*shape[1]);
      KK_FLOAT exone[3], eyone[3], ezone[3], am[3];
      KK_FLOAT omegae[3] = {(KK_FLOAT)b.omega[0],(KK_FLOAT)b.omega[1],(KK_FLOAT)b.omega[2]};
      MathExtraKokkos::q_to_exyz(quatatom, exone, eyone, ezone);
      MathExtraKokkos::omega_to_angmom(omegae, exone, eyone, ezone, ione, am);
      d_angmom(i,0) = am[0];
      d_angmom(i,1) = am[1];
      d_angmom(i,2) = am[2];
    }
    // point dipole: rotate the body-frame dipole orientation into space frame
    if (ef & RigidConst::DIPOLE) {
      KK_FLOAT p[3][3];
      MathExtraKokkos::quat_to_mat(b.quat, p);
      KK_FLOAT dor[3] = {(KK_FLOAT)d_dorient(i,0), (KK_FLOAT)d_dorient(i,1), (KK_FLOAT)d_dorient(i,2)};
      KK_FLOAT muvec[3];
      MathExtraKokkos::matvec(p, dor, muvec);
      MathExtraKokkos::snormalize3(d_mu(i,3), muvec, muvec);
      d_mu(i,0) = muvec[0];
      d_mu(i,1) = muvec[1];
      d_mu(i,2) = muvec[2];
    }
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixRigidSmallKokkos<DeviceType>::operator()(TagUpdateXGC, const int ibody) const {
  Body &b = d_body(ibody);
  KK_FLOAT ex[3] = {(KK_FLOAT)b.ex_space[0],(KK_FLOAT)b.ex_space[1],(KK_FLOAT)b.ex_space[2]};
  KK_FLOAT ey[3] = {(KK_FLOAT)b.ey_space[0],(KK_FLOAT)b.ey_space[1],(KK_FLOAT)b.ey_space[2]};
  KK_FLOAT ez[3] = {(KK_FLOAT)b.ez_space[0],(KK_FLOAT)b.ez_space[1],(KK_FLOAT)b.ez_space[2]};
  KK_FLOAT xgc_body[3] = {(KK_FLOAT)b.xgc_body[0],(KK_FLOAT)b.xgc_body[1],(KK_FLOAT)b.xgc_body[2]};
  KK_FLOAT xgc[3];
  MathExtraKokkos::matvec(ex,ey,ez,xgc_body,xgc);
  b.xgc[0] = xgc[0] + b.xcm[0];
  b.xgc[1] = xgc[1] + b.xcm[1];
  b.xgc[2] = xgc[2] + b.xcm[2];
}


/* ----------------------------------------------------------------------
   resample body momenta (fix hmc); the base implementation is host code that
   rewrites body[], does a host FINAL forward comm and calls the host set_v()
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::resample_momenta(int groupbit_in, int mom_flag,
                                                       class RanPark *random_in, double KT)
{
  if (!setupflag) {
    FixRigidSmall::resample_momenta(groupbit_in, mom_flag, random_in, KT);
    return;
  }

  // bring the bodies and the atom arrays the base reads down to the host, and
  // put the comm back on the host for the duration, exactly as setup() and
  // dof() do -- otherwise the FINAL forward comm goes to the device packers
  // while the host body[] this reads and writes stays stale, and the resampled
  // momenta never reach the device at all
  atomKK->sync(Host, host_rebuild_datamask);
  copy_body_host();

  const int saved_forward_comm_device = forward_comm_device;
  const int saved_reverse_comm_device = reverse_comm_device;
  forward_comm_device = 0;
  reverse_comm_device = 0;
  setup_host_rebuild = true;

  FixRigidSmall::resample_momenta(groupbit_in, mom_flag, random_in, KT);

  setup_host_rebuild = false;
  forward_comm_device = saved_forward_comm_device;
  reverse_comm_device = saved_reverse_comm_device;

  // push the new body momenta and the per-atom velocities set_v() wrote back
  copy_body_device();
  atomKK->modified(Host, V_MASK);
  atomKK->sync(execution_space, V_MASK);
}

/* ----------------------------------------------------------------------
   write out restart info for mass, COM, inertia tensor to file
   identical format to inpfile option, so info can be read in when restarting
   each proc contributes info for rigid bodies it owns
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::write_restart_file(const char *file)
{
  // k_body is not sized until setup_device_push(); the base method returns
  // early before setup, so guard the flush the same way (cf. compute_scalar)
  if (setupflag) copy_body_host();
  FixRigidSmall::write_restart_file(file);
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::refresh_atom_views()
{
  d_bodyown   = k_bodyown.template view<DeviceType>();
  d_bodytag   = k_bodytag.template view<DeviceType>();
  d_atom2body = k_atom2body.template view<DeviceType>();
  d_xcmimage  = k_xcmimage.template view<DeviceType>();
  d_displace  = k_displace.template view<DeviceType>();
  d_vatom     = k_vatom.template view<DeviceType>();
  if (extended) {
    d_eflags  = k_eflags.template view<DeviceType>();
    if (orientflag) d_orient  = k_orient.template view<DeviceType>();
    if (dorientflag) d_dorient = k_dorient.template view<DeviceType>();
  }
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::grow_arrays(int nmax)
{
  Kokkos::Profiling::pushRegion("rigid/small grow arrays");
  int prev_size = k_bodyown.view_device().extent(0);

  // snapshot of existing base-pointer data to preserve across the one-time
  // hand-off from memory->grow ownership to grow_kokkos ownership (the base
  // constructor's create_bodies() has already populated bodytag/bodyown)
  std::vector<int> save_bodyown, save_atom2body;
  std::vector<tagint> save_bodytag;
  std::vector<imageint> save_xcmimage;
  std::vector<double> save_displace;
  int nsave = 0;

  // The grow below keeps the device copy of each array and refills the host
  // mirror from it (see the sync_host() calls further down).  That is right
  // during a run, when the device copy is the current one -- and wrong while
  // the host setup still owns this bookkeeping.  create_bodies() fills bodytag,
  // bodyown, atom2body, xcmimage and displace on the host through the base
  // class' raw pointers without touching the modify flags, and nothing pushes
  // them to the device until setup_device_push().  setup_bodies_static() calls
  // grow_arrays() in the middle of that -- but only for bodies of finite-size
  // particles, which is why this was invisible for point particles -- and the
  // refill then replaced all five arrays with the empty device allocation:
  // every bodytag went to zero, no atom was left assigned to a body, and every
  // body came out with zero mass, zero inertia and, through a zero degree-of-
  // freedom count, a nan temperature.  Carry the host side across the grow
  // whenever the host side is the one holding the values.
  const bool host_holds_bookkeeping = (!setupflag || setup_host_rebuild);

  if (!tied_initialized) {
    // the base ctor's virtual grow_arrays() allocated these as plain
    // memory->grow pointers before the Kokkos views existed; free them so the
    // grow_kokkos() calls below take the clean create_kokkos branch and own the
    // buffers (the base host pointers then alias the DualView host mirrors).
    // Preserve the already-populated contents first.
    if (bodyown != nullptr) {
      nsave = atom->nmax;
      save_bodyown.assign(bodyown, bodyown+nsave);
      save_bodytag.assign(bodytag, bodytag+nsave);
      save_atom2body.assign(atom2body, atom2body+nsave);
      save_xcmimage.assign(xcmimage, xcmimage+nsave);
      save_displace.resize((size_t)nsave*3);
      for (int i = 0; i < nsave; i++)
        for (int k = 0; k < 3; k++) save_displace[(size_t)i*3+k] = displace[i][k];
    }
    memory->destroy(bodyown);   bodyown = nullptr;
    memory->destroy(bodytag);   bodytag = nullptr;
    memory->destroy(atom2body); atom2body = nullptr;
    memory->destroy(xcmimage);  xcmimage = nullptr;
    memory->destroy(displace);  displace = nullptr;
    memory->destroy(vatom);     vatom = nullptr;
    maxvatom = 0;
    tied_initialized = true;
    prev_size = nsave;
  } else if (host_holds_bookkeeping && (bodyown != nullptr)) {
    nsave = MIN(prev_size, nmax);
    save_bodyown.assign(bodyown, bodyown+nsave);
    save_bodytag.assign(bodytag, bodytag+nsave);
    save_atom2body.assign(atom2body, atom2body+nsave);
    save_xcmimage.assign(xcmimage, xcmimage+nsave);
    save_displace.resize((size_t)nsave*3);
    for (int i = 0; i < nsave; i++)
      for (int k = 0; k < 3; k++) save_displace[(size_t)i*3+k] = displace[i][k];
  }

  // extended-particle per-atom arrays: tied DualViews so the host setup writes
  // (setup_bodies_static fills eflags/orient/dorient on the host) alias the
  // device data, and so they migrate in the device exchange.  Allocated lazily
  // the first time setup_bodies_static() sets extended (it calls grow_arrays).
  if (extended) {
    memoryKK->grow_kokkos(k_eflags, eflags, nmax, "rigid/small:eflags");
    if (orientflag) memoryKK->grow_kokkos(k_orient, orient, nmax, orientflag, "rigid/small:orient");
    if (dorientflag) memoryKK->grow_kokkos(k_dorient, dorient, nmax, 3, "rigid/small:dorient");
    d_eflags = k_eflags.template view<DeviceType>();
    if (orientflag) d_orient = k_orient.template view<DeviceType>();
    if (dorientflag) d_dorient = k_dorient.template view<DeviceType>();
  }

  memoryKK->grow_kokkos(k_bodyown, bodyown, nmax, "rigid/small:bodyown");
  memoryKK->grow_kokkos(k_bodytag, bodytag, nmax, "rigid/small:bodytag");
  memoryKK->grow_kokkos(k_atom2body, atom2body, nmax, "rigid/small:atom2body");
  memoryKK->grow_kokkos(k_xcmimage, xcmimage, nmax, "rigid/small:xcmimage");
  memoryKK->grow_kokkos(k_displace, displace, nmax, 3, "rigid/small:displace");
  if (nmax > maxvatom) {
    maxvatom = atom->nmax;
    memoryKK->grow_kokkos(k_vatom, vatom, maxvatom, 6, "fix:vatom");
  }

  refresh_atom_views();

  // grow_kokkos() resizes the DualView.  DualView::resize() resizes the device
  // view in place, replaces the host mirror with a fresh empty allocation, and
  // marks the device modified -- so after a grow the data is on the device and
  // the mirror is empty.  Every array grown above therefore needs sync_host():
  // it refills the mirror and clears the claim, which is also what keeps the
  // next modify_host() -- from here, setup_device_push(), pre_neighbor() or
  // sort_kokkos() -- from seeing both flags set and aborting.
  k_bodyown.sync_host();
  k_atom2body.sync_host();
  k_bodytag.sync_host();
  k_xcmimage.sync_host();
  k_displace.sync_host();
  k_vatom.sync_host();
  if (extended) {
    k_eflags.sync_host();
    if (orientflag) k_orient.sync_host();
    if (dorientflag) k_dorient.sync_host();
  }

  // restore the preserved contents into the new host mirrors
  for (int i = 0; i < nsave; i++) {
    bodyown[i] = save_bodyown[i];
    bodytag[i] = save_bodytag[i];
    atom2body[i] = save_atom2body[i];
    xcmimage[i] = save_xcmimage[i];
    for (int k = 0; k < 3; k++) displace[i][k] = save_displace[(size_t)i*3+k];
  }

  // Initialize new slots' sentinels on the HOST only (bodyown/atom2body = -1,
  // bodytag = 0). These are read by create_bodies/reset_atom2body during the
  // non-Kokkos setup; the device copies are synced from the host by setup()
  // (and atom2body is recomputed on device by reset_atom2body during a run).
  for (int i = prev_size; i < nmax; i++) {
    bodyown[i] = -1;
    atom2body[i] = -1;
    bodytag[i] = 0;
  }
  k_bodyown.modify_host();
  k_atom2body.modify_host();
  k_bodytag.modify_host();
  k_xcmimage.modify_host();
  k_displace.modify_host();

  Kokkos::Profiling::popRegion();
}

/* ----------------------------------------------------------------------
   initialize a molecule inserted by another fix, e.g. deposit or pour
   called when molecule is created
   nlocalprev = # of atoms on this proc before molecule inserted
   tagprev = atom ID previous to new atoms in the molecule
   xgeom = geometric center of new molecule
   vcm = COM velocity of new molecule
   quat = rotation of new molecule (around geometric center)
          relative to template in Molecule class
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::set_molecule(int nlocalprev, tagint tagprev, int imol,
                                 double *xgeom, double *vcm, double *quat)
{
  // Molecule insertion (fix deposit/pour) is done on the host: copy the current
  // bodies down, let the base class add the new body, then push the grown body
  // list and the new atoms' per-atom data back to the device.  Insertion is a
  // rare, host-side event, so the round-trip is not a per-step cost.

  // copy current bodies to host because they will later overwrite device
  copy_body_host();

  // Flush the per-atom arrays device->host first.  On the device-exchange path
  // these DualViews are modify-device (last written by unpack_exchange_kokkos),
  // so the host mirrors are stale; the base set_molecule only appends the new
  // atoms' rows, leaving the existing rows untouched.  Without this sync the
  // modify_host/sync<DeviceType> below would push the stale host array back over
  // the live device state, corrupting body bookkeeping for every existing atom.
  k_bodytag.sync_host();
  k_bodyown.sync_host();
  k_displace.sync_host();
  k_xcmimage.sync_host();
  if (extended) {
    k_eflags.sync_host();
    if (orientflag) k_orient.sync_host();
    if (dorientflag) k_dorient.sync_host();
  }

  // add new bodies on host
  FixRigidSmall::set_molecule(nlocalprev, tagprev, imol, xgeom, vcm, quat);
  // the base call already grew the capacity through the virtual grow_body() if
  // the body list was full; only match the device views to it here
  resize_body_views();
  copy_body_device();

  // base set_molecule wrote new atoms' bodytag/bodyown/displace/xcmimage (and,
  // for extended particles, eflags/orient/dorient) on the host; push back
  k_bodytag.modify_host();
  k_bodyown.modify_host();
  k_displace.modify_host();
  k_xcmimage.modify_host();
  k_bodytag.template sync<DeviceType>();
  k_bodyown.template sync<DeviceType>();
  k_displace.template sync<DeviceType>();
  k_xcmimage.template sync<DeviceType>();
  if (extended) {
    k_eflags.modify_host();
    k_eflags.template sync<DeviceType>();
    if (orientflag) { k_orient.modify_host(); k_orient.template sync<DeviceType>(); }
    if (dorientflag) { k_dorient.modify_host(); k_dorient.template sync<DeviceType>(); }
  }
  refresh_atom_views();
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::pack_exchange_kokkos (
    const int &nsend,DAT::tdual_double_2d_lr &k_buf,
    DAT::tdual_int_1d k_sendlist,
    DAT::tdual_int_1d k_copylist,
    ExecutionSpace space) {

  Kokkos::Profiling::pushRegion("rigid/small pack exchange");

  k_buf.sync<DeviceType>();
  k_copylist.sync<DeviceType>();
  k_sendlist.sync<DeviceType>();

  // bring the tied per-atom DualViews up to date on the device
  k_bodytag.template sync<DeviceType>();
  k_xcmimage.template sync<DeviceType>();
  k_displace.template sync<DeviceType>();
  k_vatom.template sync<DeviceType>();
  k_atom2body.template sync<DeviceType>();
  k_bodyown.template sync<DeviceType>();
  refresh_atom_views();
  d_body = k_body.view_device();

  auto d_buf = typename ArrayTypes<DeviceType>::t_double_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));
  auto d_copylist = k_copylist.view<DeviceType>();
  auto d_sendlist = k_sendlist.view<DeviceType>();

  auto d_vatom = this->d_vatom;
  auto d_bodytag = this->d_bodytag;
  auto d_xcmimage = this->d_xcmimage;
  auto d_displace = this->d_displace;
  auto d_atom2body = this->d_atom2body;
  auto d_bodyown = this->d_bodyown;
  auto d_body = this->d_body;
  const int bodysize = this->bodysize;
  const int nsend_local = nsend;

  // extended-particle per-atom data carried alongside each migrating body atom
  const int ext = extended_per_atom;            // 0 if not extended
  const int oflag = orientflag;
  const int dflag = dorientflag;
  auto d_eflags = this->d_eflags;
  auto d_orient = this->d_orient;
  auto d_dorient = this->d_dorient;
  if (ext) {
    k_eflags.template sync<DeviceType>();
    if (oflag) k_orient.template sync<DeviceType>();
    if (dflag) k_dorient.template sync<DeviceType>();
    d_eflags = k_eflags.template view<DeviceType>();
    if (oflag) d_orient = k_orient.template view<DeviceType>();
    if (dflag) d_dorient = k_dorient.template view<DeviceType>();
  }

  // Variable-stride packing (modeled on fix shake/kk): the first nsend doubles
  // form a header table (d_buf(mysend) = absolute offset of that atom's payload),
  // followed by the densely packed payloads.  A non-body atom costs 1 double, a
  // body member 12, and a body owner 12 + bodysize -- versus the old fixed
  // maxexchange stride for every atom.  d_count holds the total doubles written.
  if (!d_exchange_count.data())
    d_exchange_count = Kokkos::View<int, DeviceType>("rigid/small:exchange_count");
  auto d_count = d_exchange_count;
  Kokkos::deep_copy(d_count, 0);

  // count bodies whose owner is being sent (their slots are freed below); read
  // d_bodyown before the pack scan rewrites it via the copylist compaction
  int n_deleted_bodies = 0;
  Kokkos::parallel_reduce(
    "fix rigid/small count deleted bodies",
    Range1D(0, nsend),
    KOKKOS_LAMBDA(const int isend, int &n) {
      if (d_bodyown(d_sendlist(isend)) >= 0) n++;
    },
    n_deleted_bodies
  );

  Kokkos::parallel_scan(
    "fix rigid/small pack exchange",
    Range1D(0, nsend),
    KOKKOS_LAMBDA(const int isend, int &offset, const bool is_final) {
      const int i = d_sendlist(isend);
      const tagint btag = d_bodytag(i);
      const int bown = d_bodyown(i);

      int size;
      if (btag == 0) size = 1;
      else if (bown < 0) size = 12 + ext;
      else size = 12 + ext + bodysize;

      if (is_final) {
        int m = nsend_local + offset;
        d_buf(isend) = m;                       // header: offset of this payload
        d_buf(m++) = d_ubuf(btag).d;
        if (btag) {
          d_buf(m++) = d_ubuf(d_xcmimage(i)).d;
          d_buf(m++) = d_displace(i,0);
          d_buf(m++) = d_displace(i,1);
          d_buf(m++) = d_displace(i,2);
          for (int k = 0; k < 6; k++) d_buf(m++) = d_vatom(i,k);
          if (ext) {                            // extended-particle data
            d_buf(m++) = d_ubuf(d_eflags(i)).d;
            for (int k = 0; k < oflag; k++) d_buf(m++) = d_orient(i,k);
            if (dflag) { d_buf(m++) = d_dorient(i,0); d_buf(m++) = d_dorient(i,1); d_buf(m++) = d_dorient(i,2); }
          }
          if (bown < 0) {
            d_buf(m++) = 0;
          } else {
            d_buf(m++) = 1;
            memcpy(&d_buf(m),&d_body(bown),sizeof(Body));
            m += bodysize;
            d_body(bown).natoms = -1;           // mark owner's body deleted
          }
        }
        if (isend == nsend_local-1) d_count() = m;

        // copylist compaction: backfill freed slot i with kept atom j
        const int j = d_copylist(isend);
        if (j >= 0) {
          d_bodytag(i) = d_bodytag(j);
          d_xcmimage(i) = d_xcmimage(j);
          d_bodyown(i) = d_bodyown(j);
          for (int k = 0; k < 3; k++) d_displace(i,k) = d_displace(j,k);
          for (int k = 0; k < 6; k++) d_vatom(i,k) = d_vatom(j,k);
          if (ext) {
            d_eflags(i) = d_eflags(j);
            for (int k = 0; k < oflag; k++) d_orient(i,k) = d_orient(j,k);
            if (dflag) { d_dorient(i,0)=d_dorient(j,0); d_dorient(i,1)=d_dorient(j,1); d_dorient(i,2)=d_dorient(j,2); }
          }
          if (d_bodyown(i) >= 0) d_body(d_bodyown(i)).ilocal = i;
          d_bodyown(j) = -1;
          d_bodytag(j) = 0;   // 0 == "not in a body"; -1 is not a valid sentinel
          d_atom2body(j) = -1;
        }
      }
      offset += size;
    }
  );

  // Need to pack remaining bodies densely
  int new_nlocal_body = nlocal_body - n_deleted_bodies;
  if (d_from_indices.extent_int(0) < n_deleted_bodies)
    d_from_indices = IntView1D("rigid/small:from_idx", n_deleted_bodies);
  auto from_indices = d_from_indices;

  // count bodies that need to be moved
  // and do cumulative sum to determine
  // their relative position
  int n_to_move = 0;
  Kokkos::parallel_scan(
    "fix rigid/small count bodies to move",
    Range1D(new_nlocal_body, nlocal_body),
    KOKKOS_LAMBDA(const int ibody, int &count, const bool is_final){
      if (d_body(ibody).natoms > 0) {
        if (is_final) {
          from_indices(count) = ibody;
        }
        count++;
      }
    },
    n_to_move
  );


  // count open slots, fill in corresponding body
  // need copylist for updating bodyown
  Kokkos::parallel_scan(
    "fix rigid/small create body copylist",
    Range1D(0, new_nlocal_body),
    KOKKOS_LAMBDA(const int i, int &count, const bool is_final){
      if (d_body(i).natoms<0) {
        if (is_final) {
          copy_body(&d_body(i), &d_body(from_indices(count)));
          d_bodyown(d_body(i).ilocal) = i;
        }
        count++;
      }
    }
  );

  nlocal_body = new_nlocal_body;

  k_bodytag.template modify<DeviceType>();
  k_xcmimage.template modify<DeviceType>();
  k_displace.template modify<DeviceType>();
  k_vatom.template modify<DeviceType>();
  k_atom2body.template modify<DeviceType>();
  k_bodyown.template modify<DeviceType>();
  if (ext) {
    k_eflags.template modify<DeviceType>();
    if (oflag) k_orient.template modify<DeviceType>();
    if (dflag) k_dorient.template modify<DeviceType>();
  }

  k_buf.template modify<DeviceType>();
  if (space == HostKK) k_buf.sync_host();
  else k_buf.sync_device();

  int total = 0;
  Kokkos::deep_copy(total, d_count);

  Kokkos::Profiling::popRegion();

  return total;
}


/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::unpack_exchange_kokkos(DAT::tdual_double_2d_lr &k_buf,
                              DAT::tdual_int_1d &k_indices,int nrecv,
                              int nrecv1, int nrecv1extra,
                              ExecutionSpace space)
{
  Kokkos::Profiling::pushRegion("rigid/small unpack exchange");
  k_buf.sync<DeviceType>();
  k_indices.sync<DeviceType>();

  k_bodytag.template sync<DeviceType>();
  k_xcmimage.template sync<DeviceType>();
  k_displace.template sync<DeviceType>();
  k_vatom.template sync<DeviceType>();
  k_bodyown.template sync<DeviceType>();
  refresh_atom_views();

  const int ext = extended_per_atom;
  const int oflag = orientflag;
  const int dflag = dorientflag;
  auto d_eflags = this->d_eflags;
  auto d_orient = this->d_orient;
  auto d_dorient = this->d_dorient;
  if (ext) {
    k_eflags.template sync<DeviceType>();
    if (oflag) k_orient.template sync<DeviceType>();
    if (dflag) k_dorient.template sync<DeviceType>();
    d_eflags = k_eflags.template view<DeviceType>();
    if (oflag) d_orient = k_orient.template view<DeviceType>();
    if (dflag) d_dorient = k_dorient.template view<DeviceType>();
  }

  auto d_buf = typename ArrayTypes<DeviceType>::t_double_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));
  auto d_indices = k_indices.view<DeviceType>();

  auto d_bodytag = this->d_bodytag;
  auto d_xcmimage = this->d_xcmimage;
  auto d_displace = this->d_displace;
  auto d_vatom = this->d_vatom;
  auto d_bodyown = this->d_bodyown;

  // Variable-stride layout written by pack_exchange_kokkos: the buffer received
  // from a single neighbor holds a header table (one offset per atom) followed
  // by densely packed payloads.  When a swap exchanges with two neighbors
  // (procgrid > 2) the receive buffer concatenates two such blocks: atoms
  // [0,nrecv1) come from block 1 (header at d_buf(irecv)), atoms [nrecv1,nrecv)
  // come from block 2 which starts at offset nrecv1extra, so their header lives
  // at d_buf(nrecv1extra + irecv - nrecv1) and offsets are relative to
  // nrecv1extra.  This mirrors fix shake/kk's nrecv1/nextrarecv1 indirection.
  const int nrecv1_local = nrecv1;
  const int nrecv1extra_local = nrecv1extra;

  int n_new_body = 0;
  Kokkos::parallel_reduce(
    "fix rigid/small count incoming bodies",
    Range1D(0, nrecv),
    KOKKOS_LAMBDA(const int irecv, int &count){
      int i = d_indices(irecv);
      if (i < 0) return;

      int m = (int) d_buf(irecv);
      if (irecv >= nrecv1_local)
        m = nrecv1extra_local + (int) d_buf(nrecv1extra_local + irecv - nrecv1_local);

      d_bodyown(i) = -1; // set later
      d_bodytag(i) = (tagint) d_ubuf(d_buf(m++)).i;
      if (d_bodytag(i)) {
        d_xcmimage(i) = (imageint) d_ubuf(d_buf(m++)).i;
        for(int k = 0; k < 3; k++){
          d_displace(i,k) = d_buf(m++);
        }
        for(int k = 0; k < 6; k++){
          d_vatom(i,k) = d_buf(m++);
        }
        if (ext) {
          d_eflags(i) = (int) d_ubuf(d_buf(m++)).i;
          for (int k = 0; k < oflag; k++) d_orient(i,k) = d_buf(m++);
          if (dflag) { d_dorient(i,0)=d_buf(m++); d_dorient(i,1)=d_buf(m++); d_dorient(i,2)=d_buf(m++); }
        }
        if (d_buf(m++) > 0) count++;
      }
    },
    n_new_body
  );

  while (nlocal_body + n_new_body > nmax_body) {
    grow_body();
  }
  auto d_body = this->d_body;
  auto nlocal_body = this->nlocal_body;

  Kokkos::parallel_scan(
    "fix rigid/small copy incoming bodies",
    Range1D(0, nrecv),
    KOKKOS_LAMBDA(const int irecv, int &count, const bool is_final){
      int i = d_indices(irecv);
      if(i < 0) return;

      if(d_bodytag(i) <= 0) return;

      int m = (int) d_buf(irecv);
      if (irecv >= nrecv1_local)
        m = nrecv1extra_local + (int) d_buf(nrecv1extra_local + irecv - nrecv1_local);
      // payload: bodytag(1) xcmimage(1) displace(3) vatom(6) extended(ext) -> flag
      m += 11 + ext;

      // if owning body
      if (d_buf(m) > 0) {
        if (is_final) {
          int body_idx = nlocal_body + count;
          d_bodyown(i) = body_idx;
          memcpy(&d_body(body_idx),&d_buf(m+1),sizeof(Body));
          d_body(body_idx).ilocal = i;
        }
        count++;
      }
    }
  );


  this->nlocal_body += n_new_body;

  k_bodytag.template modify<DeviceType>();
  k_xcmimage.template modify<DeviceType>();
  k_displace.template modify<DeviceType>();
  k_vatom.template modify<DeviceType>();
  k_bodyown.template modify<DeviceType>();
  if (ext) {
    k_eflags.template modify<DeviceType>();
    if (oflag) k_orient.template modify<DeviceType>();
    if (dflag) k_dorient.template modify<DeviceType>();
  }

  Kokkos::Profiling::popRegion();
}

/* ----------------------------------------------------------------------
   only pack body info if own or ghost atom owns the body
   for FULL_BODY, send 0/1 flag with every atom
------------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::pack_forward_comm_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                                        DAT::tdual_double_1d &k_buf,
                                                        int pbc_flag, int* pbc)

{
  Kokkos::Profiling::pushRegion("rigid/small pack forward");

  auto d_sendlist = k_sendlist.view<DeviceType>();
  auto d_buf = k_buf.view<DeviceType>();
  int *iswap = d_sendlist.data(); // used as key by which to store body_sendlist etc.

  int m = 0;

  auto d_body = this->d_body;
  auto d_bodyown = this->d_bodyown;

  if (commflag == INITIAL) {
    const auto &send = body_send(iswap);
    int n_body = send.n;
    auto d_body_sendlist = send.list;
    Kokkos::parallel_for("fix rigid/small pack forward comm initial",
      Range1D(0, n_body),
      KOKKOS_LAMBDA (const int ibodysend) {
        int ibody = d_body_sendlist(ibodysend);
        Body &b = d_body(ibody);
        int m = 29*ibodysend;
        d_buf(m++) = b.xcm[0];
        d_buf(m++) = b.xcm[1];
        d_buf(m++) = b.xcm[2];
        d_buf(m++) = b.xgc[0];
        d_buf(m++) = b.xgc[1];
        d_buf(m++) = b.xgc[2];
        d_buf(m++) = b.vcm[0];
        d_buf(m++) = b.vcm[1];
        d_buf(m++) = b.vcm[2];
        d_buf(m++) = b.quat[0];
        d_buf(m++) = b.quat[1];
        d_buf(m++) = b.quat[2];
        d_buf(m++) = b.quat[3];
        d_buf(m++) = b.omega[0];
        d_buf(m++) = b.omega[1];
        d_buf(m++) = b.omega[2];
        d_buf(m++) = b.ex_space[0];
        d_buf(m++) = b.ex_space[1];
        d_buf(m++) = b.ex_space[2];
        d_buf(m++) = b.ey_space[0];
        d_buf(m++) = b.ey_space[1];
        d_buf(m++) = b.ey_space[2];
        d_buf(m++) = b.ez_space[0];
        d_buf(m++) = b.ez_space[1];
        d_buf(m++) = b.ez_space[2];
        d_buf(m++) = b.conjqm[0];
        d_buf(m++) = b.conjqm[1];
        d_buf(m++) = b.conjqm[2];
        d_buf(m++) = b.conjqm[3];
      }
    );
    Kokkos::Profiling::popRegion();
    return 29*n_body;
  }
  if (commflag == FINAL) {
    const auto &send = body_send(iswap);
    int n_body = send.n;
    auto d_body_sendlist = send.list;
    Kokkos::parallel_for("fix rigid/small pack forward comm final",
      Range1D(0, n_body),
      KOKKOS_LAMBDA (const int ibodysend) {
        int ibody = d_body_sendlist(ibodysend);
        Body &b = d_body(ibody);
        int m = 10*ibodysend;
        d_buf(m++) = b.vcm[0];
        d_buf(m++) = b.vcm[1];
        d_buf(m++) = b.vcm[2];
        d_buf(m++) = b.omega[0];
        d_buf(m++) = b.omega[1];
        d_buf(m++) = b.omega[2];
        d_buf(m++) = b.conjqm[0];
        d_buf(m++) = b.conjqm[1];
        d_buf(m++) = b.conjqm[2];
        d_buf(m++) = b.conjqm[3];
      }
    );
    Kokkos::Profiling::popRegion();
    return 10*n_body;
  } else if (commflag == FULL_BODY) {
    auto bodysize = this->bodysize;
    const auto &send = body_send(iswap);
    int n_body = send.n;
    auto d_body_sendlist = send.list;
    Kokkos::parallel_for(
      "fix rigid/small full body pack forward",
      Range1D(0, n_body),
      KOKKOS_LAMBDA(const int ibodysend) {
        int ibody = d_body_sendlist(ibodysend);
        int m = (bodysize)*ibodysend;

        memcpy(&d_buf(m), &d_body(ibody), sizeof(Body));
      }
    );
    Kokkos::Profiling::popRegion();
    return bodysize*n_body;
  } else if (commflag == BODY_SENDLIST) {
    int n_sent = 0;
    Kokkos::parallel_reduce(
      "fix rigid/small count bodies sent",
      Range1D(0, n),
      KOKKOS_LAMBDA(const int isend, int &tmp) {
        int i = d_sendlist(isend);
        if (d_bodyown(i) >= 0) {
          tmp++;
          d_buf(isend) = 1;
        } else {
          d_buf(isend) = 0;
        }
      },
      n_sent
    );
    auto &send = body_send(iswap);
    send.n = n_sent;

    if (send.list.extent_int(0) < n_sent) send.list = IntView1D("body sendlist", n_sent);
    auto d_body_sendlist = send.list;

    Kokkos::parallel_scan(
      "fix rigid/small create body sendlist",
      Range1D(0, n),
      KOKKOS_LAMBDA(const int isend, int &count, const bool is_final) {
        const int i = d_sendlist(isend);
        if (d_bodyown(i) >= 0) {
          if (is_final) {
            d_body_sendlist(count) = d_bodyown(i);
          }
          count++;
        }
      }
    );
    Kokkos::Profiling::popRegion();
    return n;
  }

  Kokkos::Profiling::popRegion();
  return m;
}

/* ----------------------------------------------------------------------
   only ghost atoms are looped over
   for FULL_BODY, store a new ghost body if this atom owns it
   for other commflag values, only unpack body info if atom owns it
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::unpack_forward_comm_kokkos(int n, int first, DAT::tdual_double_1d &k_buf)
{
  if (n==0) return;
  Kokkos::Profiling::pushRegion("rigid/small unpack forward");
  this->first = first;
  auto d_buf = k_buf.view<DeviceType>();
  auto d_bodyown = this->d_bodyown;
  auto d_body = this->d_body;

  if (commflag == INITIAL) {
    const auto &recv = body_recv(first);
    int n_incoming_bodies = recv.n;
    int start_body = recv.start;
    Kokkos::parallel_for("fix rigid/small unpack forward comm initial",
      Range1D(0, n_incoming_bodies),
      KOKKOS_LAMBDA(const int ibodyrecv) {
        int ibody = ibodyrecv + start_body;
        int m = 29*ibodyrecv;
        Body &b = d_body(ibody);
        b.xcm[0] = d_buf(m++);
        b.xcm[1] = d_buf(m++);
        b.xcm[2] = d_buf(m++);
        b.xgc[0] = d_buf(m++);
        b.xgc[1] = d_buf(m++);
        b.xgc[2] = d_buf(m++);
        b.vcm[0] = d_buf(m++);
        b.vcm[1] = d_buf(m++);
        b.vcm[2] = d_buf(m++);
        b.quat[0] = d_buf(m++);
        b.quat[1] = d_buf(m++);
        b.quat[2] = d_buf(m++);
        b.quat[3] = d_buf(m++);
        b.omega[0] = d_buf(m++);
        b.omega[1] = d_buf(m++);
        b.omega[2] = d_buf(m++);
        b.ex_space[0] = d_buf(m++);
        b.ex_space[1] = d_buf(m++);
        b.ex_space[2] = d_buf(m++);
        b.ey_space[0] = d_buf(m++);
        b.ey_space[1] = d_buf(m++);
        b.ey_space[2] = d_buf(m++);
        b.ez_space[0] = d_buf(m++);
        b.ez_space[1] = d_buf(m++);
        b.ez_space[2] = d_buf(m++);
        b.conjqm[0] = d_buf(m++);
        b.conjqm[1] = d_buf(m++);
        b.conjqm[2] = d_buf(m++);
        b.conjqm[3] = d_buf(m++);
      }
    );
  } else if (commflag == FINAL) {
    const auto &recv = body_recv(first);
    int n_incoming_bodies = recv.n;
    int start_body = recv.start;
    Kokkos::parallel_for("fix rigid/small/kk unpack forward comm final",
      Range1D(0, n_incoming_bodies),
      KOKKOS_LAMBDA(const int ibodyrecv) {
        int ibody = ibodyrecv + start_body;
        int m = 10*ibodyrecv;
        Body &b = d_body(ibody);
        b.vcm[0] = d_buf(m++);
        b.vcm[1] = d_buf(m++);
        b.vcm[2] = d_buf(m++);
        b.omega[0] = d_buf(m++);
        b.omega[1] = d_buf(m++);
        b.omega[2] = d_buf(m++);
        b.conjqm[0] = d_buf(m++);
        b.conjqm[1] = d_buf(m++);
        b.conjqm[2] = d_buf(m++);
        b.conjqm[3] = d_buf(m++);
      }
    );
  } else if (commflag == FULL_BODY) {
    Kokkos::Profiling::pushRegion("unpack forward full body");
    auto bodysize = this->bodysize;
    const auto &recv = body_recv(first);
    int n_incoming_bodies = recv.n;
    int start_body = recv.start;

    Kokkos::parallel_for(
      "fix rigid/small pack incoming ghost bodies",
      Range1D(0, n_incoming_bodies),
      KOKKOS_LAMBDA(const int ibodyrecv) {
        int m = ibodyrecv*(bodysize);
        memcpy(&d_body(ibodyrecv+start_body), &d_buf(m), sizeof(Body));
      }
    );
    Kokkos::Profiling::popRegion();
  } else if (commflag == BODY_SENDLIST) {
    Kokkos::Profiling::pushRegion("unpack forward body sendlist");
    int n_curr_bodies = this->nlocal_body + this->nghost_body;
    body_recv(first).start = n_curr_bodies;
    int n_incoming_bodies = 0;
    Kokkos::parallel_scan(
      "fix rigid/small count incoming bodies",
      Range1D(0, n),
      KOKKOS_LAMBDA(const int irecv, int &count, const bool is_final) {
        int i = irecv+first;
        int m = irecv;
        d_bodyown(i) = d_buf(m);
        if (d_bodyown(i)) {
          if (is_final) d_bodyown(i) = n_curr_bodies + count;
          count++;
        } else {
          d_bodyown(i) = -1;
        }
      },
      n_incoming_bodies
    );
    if (has_body_recv(first))
      error->one(FLERR, "Internal error in fix rigid/small/kk: receive buffer offset {} "
                 "is already registered", first);
    auto &recv_entry = body_recv(first);
    recv_entry.n = n_incoming_bodies;
    recv_entry.n_set = true;
    while (n_curr_bodies+n_incoming_bodies > nmax_body) {
      grow_body();
    }

    this->nghost_body += n_incoming_bodies;
    k_bodyown.template modify<DeviceType>();
    Kokkos::Profiling::popRegion();
  }
  Kokkos::Profiling::popRegion();
}

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::pack_reverse_comm_kokkos(int n, int first, DAT::tdual_double_1d &k_buf)
{
  if (commflag != FORCE_TORQUE) {
    error->all(FLERR, "attempting invalid reverse comm on device");
  }

  if (n==0) return 0;

  auto d_buf = k_buf.view<DeviceType>();
  auto d_bodyown = this->d_bodyown;
  auto d_body = this->d_body;

  const auto &recv = body_recv(first);
  int n_body = recv.n;
  int start_body = recv.start;

  Kokkos::parallel_for(
    "fix rigid/small pack reverse comm",
    Range1D(0, n_body),
    KOKKOS_LAMBDA(const int ibodysend) {
      int ibody = ibodysend + start_body;
      int m = ibodysend*6;
      Body &b = d_body(ibody);

      d_buf(m++) = b.fcm[0];
      d_buf(m++) = b.fcm[1];
      d_buf(m++) = b.fcm[2];
      d_buf(m++) = b.torque[0];
      d_buf(m++) = b.torque[1];
      d_buf(m++) = b.torque[2];
    }
  );

  return 6*n_body;
}


template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::unpack_reverse_comm_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                          DAT::tdual_double_1d &k_buf)
{
  if (commflag != FORCE_TORQUE) {
    error->all(FLERR, "attempting invalid reverse comm on device");
  }

  auto d_buf = k_buf.view<DeviceType>();
  auto d_bodyown = this->d_bodyown;
  auto d_body = this->d_body;
  auto d_sendlist = k_sendlist.view<DeviceType>();
  int *iswap = d_sendlist.data();

  const auto &send = body_send(iswap);
  int n_body = send.n;
  auto d_body_sendlist = send.list;
  Kokkos::parallel_for(
    "fix rigid/small unpack reverse comm",
    Range1D(0, n_body),
    KOKKOS_LAMBDA (const int ibodyrecv) {
      int ibody = d_body_sendlist(ibodyrecv);
      Body &b = d_body(ibody);
      int m = 6*ibodyrecv;

      b.fcm[0] += d_buf(m++);
      b.fcm[1] += d_buf(m++);
      b.fcm[2] += d_buf(m++);
      b.torque[0] += d_buf(m++);
      b.torque[1] += d_buf(m++);
      b.torque[2] += d_buf(m++);
    }
  );
}

/* ----------------------------------------------------------------------
   grow body data structure
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::resize_body_views()
{
  // Match the device views to nmax_body without changing the capacity itself.
  //
  // DualView::resize() takes its resize_on_device branch (both modify flags are
  // 0 after the sync above it), which resizes the device view -- preserving its
  // contents -- and then replaces the host mirror with a fresh, empty
  // allocation, marking the device modified.  That device claim is real: the
  // data is on the device and the mirror is empty.  sync_host() refills the
  // mirror and clears the claim; clear_sync_state() would instead leave the
  // empty mirror standing as authoritative and the next modify_host()/
  // sync<DeviceType>() would push zeros back over the live device data.
  if (k_body.view_device().extent_int(0) != nmax_body) {
    k_body.sync_host();
    k_body.resize(nmax_body);
    k_body.sync_host();
  }
  d_body = k_body.view_device();

  if (langflag && k_langextra.view_device().extent_int(0) != nmax_body) {
    k_langextra.sync_host();
    k_langextra.resize(nmax_body,6);
    k_langextra.sync_host();
    d_langextra = k_langextra.template view<DeviceType>();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::grow_body()
{
  Kokkos::Profiling::pushRegion("rigid/small grow body");

  // Always grow the capacity: callers loop on `while (... > nmax_body)
  // grow_body();`, so nmax_body has to advance on every call.  (Guarding this
  // on the device extent instead made set_molecule() -- which the base class
  // already calls grow_body() from when the body list is full -- add another
  // DELTA_BODY slots on top for every single molecule insertion.)
  FixRigidSmall::grow_body();
  resize_body_views();

  Kokkos::Profiling::popRegion();
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::reset_atom2body()
{
  if (!setupflag || setup_host_rebuild) {
    // called during setup_bodies
    FixRigidSmall::reset_atom2body();
    return;
  }
  Kokkos::Profiling::pushRegion("rigid/small reset_atom2body");

  int nlocal = atom->nlocal;

  auto d_atom2body = this->d_atom2body;
  auto d_bodytag = this->d_bodytag;
  auto d_bodyown = this->d_bodyown;
  auto atomKK = this->atomKK;

  auto map_style = atom->map_style;
  decltype(atomKK->k_map_array) k_map_array;
  decltype(atomKK->k_map_hash) k_map_hash;
  if (map_style == Atom::MAP_ARRAY) {
    k_map_array = atomKK->k_map_array;
    k_map_array.template sync<DeviceType>();
  } else if (map_style == Atom::MAP_HASH) {
    k_map_hash = atomKK->k_map_hash;
    k_map_hash.template sync<DeviceType>();
  }


  // count atoms whose body owner is missing on this proc (iowner < 0); the host
  // FixRigidSmall::reset_atom2body() treats this as a fatal error.  Doing it on
  // the device avoids the out-of-bounds d_bodyown(-1) read the previous code
  // performed, and lets us reproduce the host error afterwards.
  // reduce the missing-owner count rather than atomically incrementing one
  // device address from every atom, which serialises the whole kernel
  int nmissing = 0;
  Kokkos::parallel_reduce(
    "fix rigid/small reset atom2body",
    Range1D(0, nlocal),
    KOKKOS_LAMBDA(const int i, int &missing){
      d_atom2body(i) = -1;
      if (d_bodytag(i)) {
        int iowner = AtomKokkos::map_kokkos<DeviceType>(d_bodytag(i),map_style,k_map_array,k_map_hash);
        if (iowner >= 0) d_atom2body(i) = d_bodyown(iowner);
        else missing++;
      }
    },
    nmissing
  );
  k_atom2body.template modify<DeviceType>();
  if (nmissing)
    error->one(FLERR, "Rigid body atoms missing on proc {} at step {} "
               "(fix rigid/small/kk)", comm->me, update->ntimestep);
  Kokkos::Profiling::popRegion();
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::image_shift()
{
  if (!setupflag || setup_host_rebuild) {
    // called during setup_bodies
    FixRigidSmall::image_shift();
    return;
  }
  Kokkos::Profiling::pushRegion("rigid/small image shift");

  // the kernel reads atom->image; sync it to the device first (on the host
  // exchange path image was last updated on the host, so the device view is
  // stale -- invisible on CPU where host and device share the buffer)
  atomKK->sync(execution_space, IMAGE_MASK);
  ImageIntView1D d_image = atomKK->k_image.view<DeviceType>();
  int nlocal = atom->nlocal;
  auto d_body = this->d_body;
  auto d_xcmimage = this->d_xcmimage;
  auto d_atom2body = this->d_atom2body;

  Kokkos::parallel_for(
    "fix rigid/small image shift",
    Range1D(0, nlocal),
    KOKKOS_LAMBDA(const int i) {
      if (d_atom2body(i) < 0) return;
      Body &b = d_body(d_atom2body(i));

      imageint tdim,bdim,xdim[3];
      tdim = d_image(i) & IMGMASK;
      bdim = b.image & IMGMASK;
      xdim[0] = IMGMAX + tdim - bdim;
      tdim = (d_image(i) >> IMGBITS) & IMGMASK;
      bdim = (b.image >> IMGBITS) & IMGMASK;
      xdim[1] = IMGMAX + tdim - bdim;
      tdim = d_image(i) >> IMG2BITS;
      bdim = b.image >> IMG2BITS;
      xdim[2] = IMGMAX + tdim - bdim;

      d_xcmimage(i) = (xdim[2] << IMG2BITS) | (xdim[1] << IMGBITS) | xdim[0];
    }
  );
  k_xcmimage.template modify<DeviceType>();
  Kokkos::Profiling::popRegion();
}

/* ----------------------------------------------------------------------
   The legacy sort permutes these arrays through copy_arrays() on the host,
   where sort_kokkos() below permutes them on the device.  Bring them to the
   host for it, and claim the host side once it has run, or the device copies
   keep the pre-sort ordering: atom2body then resolves an atom to the wrong
   body and set_xv() places it at that body's centre.
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::sync_host_for_sort()
{
  k_bodytag.sync_host();
  k_bodyown.sync_host();
  k_atom2body.sync_host();
  k_xcmimage.sync_host();
  k_displace.sync_host();
  k_vatom.sync_host();
  if (extended) {
    k_eflags.sync_host();
    if (orientflag) k_orient.sync_host();
    if (dorientflag) k_dorient.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::modified_host_for_sort()
{
  k_bodytag.modify_host();
  k_bodyown.modify_host();
  k_atom2body.modify_host();
  k_xcmimage.modify_host();
  k_displace.modify_host();
  k_vatom.modify_host();
  if (extended) {
    k_eflags.modify_host();
    if (orientflag) k_orient.modify_host();
    if (dorientflag) k_dorient.modify_host();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  Kokkos::Profiling::pushRegion("rigid/small sort");

  // At setup (setupflag still 0) this sort fires after create_bodies() wrote the
  // host body arrays, so those host copies are the ones that must be sorted:
  // make them authoritative, and the sync<DeviceType>() below then carries them
  // to the device to be permuted alongside the atoms the sort reorders.
  // Otherwise the host bookkeeping no longer matches the atom order:
  // reset_atom2body leaves bodyown[map(bodytag)] < 0 and set_v() skips those
  // atoms, leaving their velocities unprojected.
  //
  // The device claim from the setup-time exchange is still outstanding here and
  // has to be retired first: Verlet::setup() runs comm->exchange() and
  // atom->sort() back to back, both *before* the first setup_pre_neighbor(), so
  // the equivalent clear there has not run yet at this point.  Nothing is lost
  // by dropping it -- nothing had pushed the host arrays create_bodies() wrote
  // to the device, so what the exchange moved derives from an uninitialized
  // copy.  No-op when host space == device space.
  if (!setupflag) {
    k_bodytag.clear_sync_state();   k_bodytag.modify_host();
    k_bodyown.clear_sync_state();   k_bodyown.modify_host();
    k_atom2body.clear_sync_state(); k_atom2body.modify_host();
    k_xcmimage.clear_sync_state();  k_xcmimage.modify_host();
    k_displace.clear_sync_state();  k_displace.modify_host();
    k_vatom.clear_sync_state();     k_vatom.modify_host();
    if (extended) {
      k_eflags.clear_sync_state();  k_eflags.modify_host();
      if (orientflag)  { k_orient.clear_sync_state();  k_orient.modify_host(); }
      if (dorientflag) { k_dorient.clear_sync_state(); k_dorient.modify_host(); }
    }
  }

  // sort the device side of each tied DualView in place
  auto space = LMPDeviceType();
  k_bodytag.template sync<LMPDeviceType>();
  k_bodyown.template sync<LMPDeviceType>();
  k_xcmimage.template sync<LMPDeviceType>();
  k_displace.template sync<LMPDeviceType>();
  k_vatom.template sync<LMPDeviceType>();
  k_atom2body.template sync<LMPDeviceType>();

  Sorter.sort(space, k_bodytag.template view<LMPDeviceType>());
  Sorter.sort(space, k_bodyown.template view<LMPDeviceType>());
  Sorter.sort(space, k_xcmimage.template view<LMPDeviceType>());
  Sorter.sort(space, k_displace.template view<LMPDeviceType>());
  Sorter.sort(space, k_vatom.template view<LMPDeviceType>());
  Sorter.sort(space, k_atom2body.template view<LMPDeviceType>());

  k_bodytag.modify_device();
  k_bodyown.modify_device();
  k_xcmimage.modify_device();
  k_displace.modify_device();
  k_vatom.modify_device();
  k_atom2body.modify_device();

  // the extended per-atom arrays are indexed by local atom in set_xv, so they
  // must follow the same permutation (an extended system that ever reaches the
  // device sort path -- this is currently gated to the host sort, but keep the
  // sort self-consistent regardless)
  if (extended) {
    k_eflags.template sync<LMPDeviceType>();
    Sorter.sort(space, k_eflags.template view<LMPDeviceType>());
    k_eflags.modify_device();
    if (orientflag) {
      k_orient.template sync<LMPDeviceType>();
      Sorter.sort(space, k_orient.template view<LMPDeviceType>());
      k_orient.modify_device();
    }
    if (dorientflag) {
      k_dorient.template sync<LMPDeviceType>();
      Sorter.sort(space, k_dorient.template view<LMPDeviceType>());
      k_dorient.modify_device();
    }
  }

  // At setup the host setup path (setup_bodies_static -> reset_atom2body / set_v)
  // reads these right after, so pull the freshly sorted device copies back to the
  // host -- they now follow the same permutation as the reordered atoms.  No
  // body.ilocal re-link yet: bodies are not fully wired until setup() completes.
  // During a run setupflag is 1 and the device copies stay authoritative (the
  // host is synced on demand).
  if (!setupflag) {
    k_bodytag.sync_host();
    k_bodyown.sync_host();
    k_atom2body.sync_host();
    k_xcmimage.sync_host();
    k_displace.sync_host();
    k_vatom.sync_host();
    if (extended) {
      k_eflags.sync_host();
      if (orientflag)  k_orient.sync_host();
      if (dorientflag) k_dorient.sync_host();
    }
    Kokkos::Profiling::popRegion();
    return;
  }

  // refresh d_body from k_body: grow_body() during a preceding exchange may have
  // reallocated the body buffer, leaving this->d_body stale.  The body.ilocal
  // back-pointers must be written into the live buffer that pack_exchange reads.
  this->d_body = k_body.view_device();
  auto d_body = this->d_body;
  auto d_bodyown = this->d_bodyown;
  int nlocal = atom->nlocal;

  Kokkos::parallel_for(
    "fix rigid/small update body.ilocal after sort",
    Range1D(0, nlocal),
    KOKKOS_LAMBDA(const int i){
      if (d_bodyown(i) < 0) return;
      d_body(d_bodyown(i)).ilocal = i;
    }
  );
  Kokkos::Profiling::popRegion();
}


/* ----------------------------------------------------------------------
   zero linear momentum of each rigid body
   set Vcm to 0.0, then reset velocities of particles via set_v()
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::zero_momentum()
{
  if (!setupflag) {
    FixRigidSmall::zero_momentum();
    return;
  }
  auto d_body = this->d_body;
  auto nlocal_body = this->nlocal_body;
  auto nghost_body = this->nghost_body;
  Kokkos::parallel_for(
      "fix rigid/small zero momentum",
      Range1D(0, nlocal_body+nghost_body),
      KOKKOS_LAMBDA(const int ibody){
        double *vcm = d_body(ibody).vcm;
        vcm[0] = vcm[1] = vcm[2] = 0.0;
      }
  );

  // forward communicate vcm to all ghost copies

  commflag = FINAL;
  commKK->forward_comm_device<DeviceType>(this,10);

  // set velocity of atoms in rigid bodues

  evflag = 0;
  set_xv_kokkos(0);
}

/* ----------------------------------------------------------------------
   zero angular momentum of each rigid body
   set angmom/omega to 0.0, then reset velocities of particles via set_v()
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::zero_rotation()
{
  if (!setupflag) {
    FixRigidSmall::zero_rotation();
    return;
  }
  auto d_body = this->d_body;
  Kokkos::parallel_for(
    "fix rigid/small zero rotation",
    Range1D(0, nlocal_body+nghost_body),
    KOKKOS_LAMBDA(const int ibody){
      double *angmom = d_body(ibody).angmom;
      angmom[0] = angmom[1] = angmom[2] = 0.0;
      double *omega = d_body(ibody).omega;
      omega[0] = omega[1] = omega[2] = 0.0;
    }
  );

  // forward communicate of omega to all ghost copies

  commflag = FINAL;
  commKK->forward_comm_device<DeviceType>(this,10);

  // set velocity of atoms in rigid bodues

  evflag = 0;
  set_xv_kokkos(0);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void *FixRigidSmallKokkos<DeviceType>::extract(const char *str, int &dim)
{
  dim = 0;

  if (strcmp(str,"body") == 0) {
    if (!setupflag) return nullptr;
    dim = 1;
    k_atom2body.sync_host();
    return atom2body;
  }

  if (strcmp(str,"onemol") == 0) {
    dim = 0;
    return onemols;
  }

  // return vector of rigid body masses, for owned+ghost bodies
  // used by granular pair styles, indexed by atom2body

  if (strcmp(str,"masstotal") == 0) {
    if (!setupflag) return nullptr;
    dim = 1;

    copy_body_host();
    if (nmax_mass < nmax_body) {
      memory->destroy(mass_body);
      nmax_mass = nmax_body;
      memory->create(mass_body,nmax_mass,"rigid:mass_body");
    }

    int n = nlocal_body + nghost_body;
    for (int i = 0; i < n; i++)
      mass_body[i] = body[i].mass;

    return mass_body;
  }

  return nullptr;
}

/* ----------------------------------------------------------------------
   return translational KE for all rigid bodies
   KE = 1/2 M Vcm^2
   sum local body results across procs
------------------------------------------------------------------------- */

template<class DeviceType>
double FixRigidSmallKokkos<DeviceType>::extract_ke()
{
  if (!setupflag) {
    return FixRigidSmall::extract_ke();
  }

  auto d_body = this->d_body;

  double ke = 0.0;
  Kokkos::parallel_reduce(
    "fix rigid/small ke",
    Range1D(0, nlocal_body),
    KOKKOS_LAMBDA(const int i, double &ke){
      double *vcm = d_body(i).vcm;
      ke += d_body(i).mass * (vcm[0]*vcm[0] + vcm[1]*vcm[1] + vcm[2]*vcm[2]);
    },
    ke
  );

  double keall;
  MPI_Allreduce(&ke,&keall,1,MPI_DOUBLE,MPI_SUM,world);

  return 0.5*keall;
}

/* ----------------------------------------------------------------------
   return rotational KE for all rigid bodies
   Erotational = 1/2 I wbody^2
------------------------------------------------------------------------- */

template<class DeviceType>
double FixRigidSmallKokkos<DeviceType>::extract_erotational()
{
  if (!setupflag) {
    return FixRigidSmall::extract_erotational();
  }

  double erotate = 0.0;
  auto d_body = this->d_body;
  auto nlocal_body = this->nlocal_body;
  Kokkos::parallel_reduce(
    "fix rigid/small erotational",
    Range1D(0, nlocal_body),
    KOKKOS_LAMBDA(const int i, double &erotate){
      KK_FLOAT wbody[3],rot[3][3];
      double *inertia;

      // for Iw^2 rotational term, need wbody = angular velocity in body frame
      // not omega = angular velocity in space frame

      inertia = d_body(i).inertia;
      KK_FLOAT angmom[3] = {(KK_FLOAT)d_body(i).angmom[0],(KK_FLOAT)d_body(i).angmom[1],(KK_FLOAT)d_body(i).angmom[2]};
      MathExtraKokkos::quat_to_mat(d_body(i).quat,rot);
      MathExtraKokkos::transpose_matvec(rot,angmom,wbody);
      if (inertia[0] == 0.0) wbody[0] = 0.0;
      else wbody[0] /= inertia[0];
      if (inertia[1] == 0.0) wbody[1] = 0.0;
      else wbody[1] /= inertia[1];
      if (inertia[2] == 0.0) wbody[2] = 0.0;
      else wbody[2] /= inertia[2];

      erotate += inertia[0]*wbody[0]*wbody[0] + inertia[1]*wbody[1]*wbody[1] +
        inertia[2]*wbody[2]*wbody[2];
    },
    erotate
  );

  double erotateall;
  MPI_Allreduce(&erotate,&erotateall,1,MPI_DOUBLE,MPI_SUM,world);

  return 0.5*erotateall;
}

/* ----------------------------------------------------------------------
   return temperature of collection of rigid bodies
   non-active DOF are removed by fflag/tflag and in tfactor
------------------------------------------------------------------------- */

template<class DeviceType>
double FixRigidSmallKokkos<DeviceType>::compute_scalar()
{
  if (!setupflag) {
    return FixRigidSmall::compute_scalar();
  }

  double t = 0.0;
  auto d_body = this->d_body;
  auto nlocal_body = this->nlocal_body;

  Kokkos::parallel_reduce(
    "fix rigid/small compute scalar",
    Range1D(0, nlocal_body),
    KOKKOS_LAMBDA(const int i, double &t) {
      KK_FLOAT wbody[3],rot[3][3];
      double *vcm,*inertia;

      vcm = d_body(i).vcm;
      t += d_body(i).mass * (vcm[0]*vcm[0] + vcm[1]*vcm[1] + vcm[2]*vcm[2]);

      // for Iw^2 rotational term, need wbody = angular velocity in body frame
      // not omega = angular velocity in space frame

      inertia = d_body(i).inertia;
      KK_FLOAT angmom[3] = {(KK_FLOAT)d_body(i).angmom[0],(KK_FLOAT)d_body(i).angmom[1],(KK_FLOAT)d_body(i).angmom[2]};
      MathExtraKokkos::quat_to_mat(d_body(i).quat,rot);
      MathExtraKokkos::transpose_matvec(rot,angmom,wbody);
      if (inertia[0] == 0.0) wbody[0] = 0.0;
      else wbody[0] /= inertia[0];
      if (inertia[1] == 0.0) wbody[1] = 0.0;
      else wbody[1] /= inertia[1];
      if (inertia[2] == 0.0) wbody[2] = 0.0;
      else wbody[2] /= inertia[2];

      t += inertia[0]*wbody[0]*wbody[0] + inertia[1]*wbody[1]*wbody[1] +
        inertia[2]*wbody[2]*wbody[2];
    },
    t
  );

  double tall;
  MPI_Allreduce(&t,&tall,1,MPI_DOUBLE,MPI_SUM,world);

  double tfactor = force->mvv2e / ((6.0*nbody - nlinear) * force->boltz);
  tall *= tfactor;
  return tall;
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::copy_body_host(){
  Kokkos::Profiling::pushRegion("rigid/small copy body host");
  // k_body is not sized until setup_device_push().  Callers can run before that
  // even with setupflag already set -- setup_pre_neighbor() sets it, and
  // ModifyKokkos::setup() then runs compute temp -> dof() before the fix setup()
  // loop is reached.  Until the device array exists the host body[] that
  // create_bodies() filled is the authoritative copy, so there is nothing to
  // bring down.
  if (k_body.view_device().extent_int(0) < nlocal_body + nghost_body) {
    Kokkos::Profiling::popRegion();
    return;
  }
  // d_body is written directly by device kernels (integrate/comm) without
  // updating DualView modify flags, so copy explicitly device -> host
  Kokkos::deep_copy(k_body.view_host(), k_body.view_device());
  for(int ibody = 0; ibody < nlocal_body + nghost_body; ibody++){
    copy_body(&body[ibody], &k_body.view_host()(ibody));
  }
  Kokkos::Profiling::popRegion();
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::copy_body_device(){
  Kokkos::Profiling::pushRegion("rigid/small copy body device");
  // Size the device array to what we are about to write.  nlocal_body +
  // nghost_body can exceed nmax_body (the base set_molecule() only grows when
  // nlocal_body == nmax_body, which misses the ghost bodies), and k_body does
  // not exist at all before setup_device_push(); either way writing
  // k_body.view_host()(ibody) below would run past the extent.
  const int nbody_all = nlocal_body + nghost_body;
  if (k_body.view_host().extent_int(0) < nbody_all) {
    k_body.sync_host();
    k_body.resize(nbody_all);
    k_body.sync_host();
    d_body = k_body.view_device();
  }
  for(int ibody = 0; ibody < nlocal_body + nghost_body; ibody++){
    copy_body(&k_body.view_host()(ibody), &body[ibody]);
  }
  Kokkos::deep_copy(k_body.view_device(), k_body.view_host());
  d_body = k_body.view_device();
  Kokkos::Profiling::popRegion();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixRigidSmallKokkos<DeviceType>::v_tally(EV_FLOAT &ev, int i, double vtot[6], double r[3], double f[3], double center[3]) const{
  // the centroid (r, f, center) arguments are only needed for the per-atom
  // centroid virial, which is unsupported on the device and rejected up front
  // in set_xv_kokkos(); here we only tally the global / ordinary per-atom virial
  (void) r; (void) f; (void) center;
  v_tally(ev, i, vtot);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixRigidSmallKokkos<DeviceType>::v_tally(EV_FLOAT &ev, int i, double vtot[6]) const{
  if (vflag_global) {
    ev.v[0] += vtot[0];
    ev.v[1] += vtot[1];
    ev.v[2] += vtot[2];
    ev.v[3] += vtot[3];
    ev.v[4] += vtot[4];
    ev.v[5] += vtot[5];
  }

  if (vflag_atom) {
    d_vatom(i,0) += vtot[0];
    d_vatom(i,1) += vtot[1];
    d_vatom(i,2) += vtot[2];
    d_vatom(i,3) += vtot[3];
    d_vatom(i,4) += vtot[4];
    d_vatom(i,5) += vtot[5];
  }
}

namespace LAMMPS_NS {
template class FixRigidSmallKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixRigidSmallKokkos<LMPHostType>;
#endif
}

