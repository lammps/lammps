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

/* ----------------------------------------------------------------------
   Contributing author: Mitch Murphy (alphataubio@gmail.com)
------------------------------------------------------------------------- */

#include "fix_rigid_small_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "comm.h"
#include "domain_kokkos.h"
#include "error.h"
#include "force.h"
#include "math_extra_kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "rigid_const.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace RigidConst;

static constexpr int DELTA_BODY = 10000;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidSmallKokkos<DeviceType>::FixRigidSmallKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRigidSmall(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  // re-allocate per-atom data as Kokkos DualViews
  // base constructor already allocated with memory->grow via grow_arrays()

  int nmax = atom->nmax;

  memoryKK->destroy_kokkos(k_bodyown, bodyown);
  memoryKK->destroy_kokkos(k_bodytag, bodytag);
  memoryKK->destroy_kokkos(k_atom2body, atom2body);
  memoryKK->destroy_kokkos(k_xcmimage, xcmimage);
  memoryKK->destroy_kokkos(k_displace, displace);

  memoryKK->create_kokkos(k_bodyown, bodyown, nmax, "rigid/small:bodyown");
  memoryKK->create_kokkos(k_bodytag, bodytag, nmax, "rigid/small:bodytag");
  memoryKK->create_kokkos(k_atom2body, atom2body, nmax, "rigid/small:atom2body");
  memoryKK->create_kokkos(k_xcmimage, xcmimage, nmax, "rigid/small:xcmimage");
  memoryKK->create_kokkos(k_displace, displace, nmax, 3, "rigid/small:displace");

  d_bodyown = k_bodyown.template view<DeviceType>();
  d_bodytag = k_bodytag.template view<DeviceType>();
  d_atom2body = k_atom2body.template view<DeviceType>();
  d_xcmimage = k_xcmimage.template view<DeviceType>();
  d_displace = k_displace.template view<DeviceType>();

  if (extended) {
    memoryKK->destroy_kokkos(k_eflags, eflags);
    memoryKK->create_kokkos(k_eflags, eflags, nmax, "rigid/small:eflags");
    d_eflags = k_eflags.template view<DeviceType>();
  }

  // body DualView
  k_body = Kokkos::DualView<Body*, Kokkos::LayoutRight, DeviceType>("rigid/small:body", nmax_body);
  if (nlocal_body > 0)
    Kokkos::deep_copy(Kokkos::subview(k_body.h_view, Kokkos::make_pair(0, nlocal_body)),
                      Kokkos::View<Body*, Kokkos::LayoutRight, Kokkos::HostSpace>(body, nlocal_body));
  memory->sfree(body);
  body = k_body.h_view.data();
  d_body = k_body.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidSmallKokkos<DeviceType>::~FixRigidSmallKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_bodyown, bodyown);
  memoryKK->destroy_kokkos(k_bodytag, bodytag);
  memoryKK->destroy_kokkos(k_atom2body, atom2body);
  memoryKK->destroy_kokkos(k_xcmimage, xcmimage);
  memoryKK->destroy_kokkos(k_displace, displace);
  if (extended)
    memoryKK->destroy_kokkos(k_eflags, eflags);

  // null base class pointers so base destructor won't double-free
  bodyown = nullptr;
  bodytag = nullptr;
  atom2body = nullptr;
  xcmimage = nullptr;
  displace = nullptr;
  eflags = nullptr;

  // body managed by k_body DualView
  body = nullptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::init()
{
  FixRigidSmall::init();

  atomKK->k_mass.modify_host();
  atomKK->k_mass.template sync<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::setup(int vflag)
{
  atomKK->sync(Host, ALL_MASK);
  FixRigidSmall::setup(vflag);
  atomKK->modified(Host, X_MASK | V_MASK);

  k_body.modify_host();
  sync_fix_data_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::initial_integrate(int vflag)
{
  // sync atom data to device
  atomKK->sync(execution_space, X_MASK | V_MASK | F_MASK | MASK_MASK |
               TYPE_MASK | RMASS_MASK | IMAGE_MASK);

  d_x = atomKK->k_x.template view<DeviceType>();
  d_v = atomKK->k_v.template view<DeviceType>();
  d_f = atomKK->k_f.template view<DeviceType>();
  d_rmass = atomKK->k_rmass.template view<DeviceType>();
  d_mass = atomKK->k_mass.template view<DeviceType>();
  d_type = atomKK->k_type.template view<DeviceType>();
  d_mask = atomKK->k_mask.template view<DeviceType>();
  d_image = atomKK->k_image.template view<DeviceType>();

  // sync fix data to device
  k_body.template sync<DeviceType>();
  k_atom2body.template sync<DeviceType>();
  k_displace.template sync<DeviceType>();
  k_xcmimage.template sync<DeviceType>();
  d_body = k_body.template view<DeviceType>();
  d_atom2body = k_atom2body.template view<DeviceType>();
  d_displace = k_displace.template view<DeviceType>();
  d_xcmimage = k_xcmimage.template view<DeviceType>();

  // body integration loop on device
  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,
    TagRigidSmallInitialIntegrate>(0, nlocal_body), *this);
  copymode = 0;

  k_body.template modify<DeviceType>();
  k_body.sync_host();

  // virial setup
  v_init(vflag);

  // forward communicate body info
  commflag = INITIAL;
  comm->forward_comm(this, 29);

  k_body.modify_host();
  k_body.template sync<DeviceType>();
  d_body = k_body.template view<DeviceType>();

  // set coords/velocity of atoms from body state -- device kernel
  d_prd = Few<double,3>(domain->prd);
  d_h = Few<double,6>(domain->h);

  if (vflag_atom) {
    if (atom->nmax > (int) d_vatom.extent(0)) {
      memoryKK->destroy_kokkos(k_vatom, vatom);
      memoryKK->create_kokkos(k_vatom, vatom, atom->nmax, 6, "fix:vatom");
      d_vatom = k_vatom.template view<DeviceType>();
    } else {
      k_vatom.template sync<DeviceType>();
    }
    Kokkos::deep_copy(d_vatom, 0.0);
  }

  if (triclinic) {
    if (evflag)
      set_xv_kokkos<1,1>();
    else
      set_xv_kokkos<1,0>();
  } else {
    if (evflag)
      set_xv_kokkos<0,1>();
    else
      set_xv_kokkos<0,0>();
  }

  // extended particles: host fallback
  if (extended) {
    atomKK->sync(Host, ALL_MASK);
    k_body.sync_host();
    // extended particle code runs on host via base class set_xv internals
    // handled inside set_xv_kokkos after device kernel
  }

  atomKK->modified(execution_space, X_MASK | V_MASK);

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixRigidSmallKokkos<DeviceType>::operator()(TagRigidSmallInitialIntegrate, const int &ibody) const
{
  Body b = d_body[ibody];

  KK_FLOAT dtfm = dtf / b.mass;
  b.vcm[0] += dtfm * b.fcm[0];
  b.vcm[1] += dtfm * b.fcm[1];
  b.vcm[2] += dtfm * b.fcm[2];

  b.xcm[0] += dtv * b.vcm[0];
  b.xcm[1] += dtv * b.vcm[1];
  b.xcm[2] += dtv * b.vcm[2];

  b.angmom[0] += dtf * b.torque[0];
  b.angmom[1] += dtf * b.torque[1];
  b.angmom[2] += dtf * b.torque[2];

  KK_FLOAT omega_k[3];
  MathExtraKokkos::angmom_to_omega(b.angmom, b.ex_space, b.ey_space,
                                   b.ez_space, b.inertia, omega_k);
  b.omega[0] = omega_k[0];
  b.omega[1] = omega_k[1];
  b.omega[2] = omega_k[2];

  MathExtraKokkos::richardson(b.quat, b.angmom, b.omega, b.inertia, dtq);
  MathExtraKokkos::q_to_exyz(b.quat, b.ex_space, b.ey_space, b.ez_space);

  d_body[ibody] = b;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::post_force(int /*vflag*/)
{
  if (langflag) {
    atomKK->sync(Host, ALL_MASK);
    k_body.sync_host();
    apply_langevin_thermostat();
    k_body.modify_host();
  }
  if (earlyflag) {
    compute_forces_and_torques_kokkos();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::final_integrate()
{
  atomKK->sync(execution_space, X_MASK | V_MASK | F_MASK | MASK_MASK |
               TYPE_MASK | RMASS_MASK | IMAGE_MASK);

  d_x = atomKK->k_x.template view<DeviceType>();
  d_v = atomKK->k_v.template view<DeviceType>();
  d_f = atomKK->k_f.template view<DeviceType>();
  d_rmass = atomKK->k_rmass.template view<DeviceType>();
  d_mass = atomKK->k_mass.template view<DeviceType>();
  d_type = atomKK->k_type.template view<DeviceType>();
  d_mask = atomKK->k_mask.template view<DeviceType>();
  d_image = atomKK->k_image.template view<DeviceType>();

  if (!earlyflag) compute_forces_and_torques_kokkos();
  if (domain->dimension == 2) enforce2d_kokkos();

  // body integration loop on device
  k_body.template sync<DeviceType>();
  d_body = k_body.template view<DeviceType>();

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,
    TagRigidSmallFinalIntegrate>(0, nlocal_body), *this);
  copymode = 0;

  k_body.template modify<DeviceType>();
  k_body.sync_host();

  commflag = FINAL;
  comm->forward_comm(this, 10);

  k_body.modify_host();
  k_body.template sync<DeviceType>();
  d_body = k_body.template view<DeviceType>();

  // sync fix data
  k_atom2body.template sync<DeviceType>();
  k_displace.template sync<DeviceType>();
  k_xcmimage.template sync<DeviceType>();
  d_atom2body = k_atom2body.template view<DeviceType>();
  d_displace = k_displace.template view<DeviceType>();
  d_xcmimage = k_xcmimage.template view<DeviceType>();

  d_prd = Few<double,3>(domain->prd);
  d_h = Few<double,6>(domain->h);

  if (vflag_atom) {
    k_vatom.template sync<DeviceType>();
  }

  if (triclinic) {
    if (evflag)
      set_v_kokkos<1,1>();
    else
      set_v_kokkos<1,0>();
  } else {
    if (evflag)
      set_v_kokkos<0,1>();
    else
      set_v_kokkos<0,0>();
  }

  // extended: host fallback
  if (extended) {
    atomKK->sync(Host, ALL_MASK);
    k_body.sync_host();
  }

  atomKK->modified(execution_space, V_MASK);

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixRigidSmallKokkos<DeviceType>::operator()(TagRigidSmallFinalIntegrate, const int &ibody) const
{
  Body b = d_body[ibody];

  KK_FLOAT dtfm = dtf / b.mass;
  b.vcm[0] += dtfm * b.fcm[0];
  b.vcm[1] += dtfm * b.fcm[1];
  b.vcm[2] += dtfm * b.fcm[2];

  b.angmom[0] += dtf * b.torque[0];
  b.angmom[1] += dtf * b.torque[1];
  b.angmom[2] += dtf * b.torque[2];

  KK_FLOAT omega_k[3];
  MathExtraKokkos::angmom_to_omega(b.angmom, b.ex_space, b.ey_space,
                                   b.ez_space, b.inertia, omega_k);
  b.omega[0] = omega_k[0];
  b.omega[1] = omega_k[1];
  b.omega[2] = omega_k[2];

  d_body[ibody] = b;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int TRICLINIC, int EVFLAG>
void FixRigidSmallKokkos<DeviceType>::set_xv_kokkos()
{
  int nlocal = atomKK->nlocal;

  double result[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  copymode = 1;
  if (EVFLAG) {
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,
      TagRigidSmallSetXV<TRICLINIC,1>>(0, nlocal), *this, result);
  } else {
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,
      TagRigidSmallSetXV<TRICLINIC,0>>(0, nlocal), *this, result);
  }
  copymode = 0;

  if (EVFLAG && vflag_global) {
    virial[0] += result[0];
    virial[1] += result[1];
    virial[2] += result[2];
    virial[3] += result[3];
    virial[4] += result[4];
    virial[5] += result[5];
  }

  // update geometric center of bodies
  int nbody_total = nlocal_body + nghost_body;
  for (int ibody = 0; ibody < nbody_total; ibody++) {
    Body *b = &body[ibody];
    double xgc_tmp[3];
    MathExtra::matvec(b->ex_space, b->ey_space, b->ez_space,
                      b->xgc_body, xgc_tmp);
    b->xgc[0] = xgc_tmp[0] + b->xcm[0];
    b->xgc[1] = xgc_tmp[1] + b->xcm[1];
    b->xgc[2] = xgc_tmp[2] + b->xcm[2];
  }
  k_body.modify_host();

  // extended particles: host fallback
  if (extended) {
    atomKK->sync(Host, ALL_MASK);
    k_body.sync_host();

    double theta_body, theta;
    double *shape, *quatatom, *inertiaatom;
    double ione[3], exone[3], eyone[3], ezone[3], p[3][3];

    AtomVecEllipsoid::Bonus *ebonus;
    if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
    AtomVecLine::Bonus *lbonus;
    if (avec_line) lbonus = avec_line->bonus;
    AtomVecTri::Bonus *tbonus;
    if (avec_tri) tbonus = avec_tri->bonus;
    double **omega_atom = atom->omega;
    double **angmom_atom = atom->angmom;
    double **mu = atom->mu;
    double *rmass_h = atom->rmass;
    int *ellipsoid = atom->ellipsoid;
    int *line = atom->line;
    int *tri = atom->tri;
    int nlocal_h = atom->nlocal;

    for (int i = 0; i < nlocal_h; i++) {
      if (atom2body[i] < 0) continue;
      Body *b = &body[atom2body[i]];

      if (eflags[i] & SPHERE) {
        omega_atom[i][0] = b->omega[0];
        omega_atom[i][1] = b->omega[1];
        omega_atom[i][2] = b->omega[2];
      } else if (eflags[i] & ELLIPSOID) {
        shape = ebonus[ellipsoid[i]].shape;
        quatatom = ebonus[ellipsoid[i]].quat;
        MathExtra::quatquat(b->quat, orient[i], quatatom);
        MathExtra::qnormalize(quatatom);
        ione[0] = EINERTIA*rmass_h[i] * (shape[1]*shape[1] + shape[2]*shape[2]);
        ione[1] = EINERTIA*rmass_h[i] * (shape[0]*shape[0] + shape[2]*shape[2]);
        ione[2] = EINERTIA*rmass_h[i] * (shape[0]*shape[0] + shape[1]*shape[1]);
        MathExtra::q_to_exyz(quatatom, exone, eyone, ezone);
        MathExtra::omega_to_angmom(b->omega, exone, eyone, ezone, ione, angmom_atom[i]);
      } else if (eflags[i] & LINE) {
        if (b->quat[3] >= 0.0) theta_body = 2.0*acos(b->quat[0]);
        else theta_body = -2.0*acos(b->quat[0]);
        theta = orient[i][0] + theta_body;
        while (theta <= -MY_PI) theta += MY_2PI;
        while (theta > MY_PI) theta -= MY_2PI;
        lbonus[line[i]].theta = theta;
        omega_atom[i][0] = b->omega[0];
        omega_atom[i][1] = b->omega[1];
        omega_atom[i][2] = b->omega[2];
      } else if (eflags[i] & TRIANGLE) {
        inertiaatom = tbonus[tri[i]].inertia;
        quatatom = tbonus[tri[i]].quat;
        MathExtra::quatquat(b->quat, orient[i], quatatom);
        MathExtra::qnormalize(quatatom);
        MathExtra::q_to_exyz(quatatom, exone, eyone, ezone);
        MathExtra::omega_to_angmom(b->omega, exone, eyone, ezone,
                                   inertiaatom, angmom_atom[i]);
      }
      if (atom->quat_flag) {
        quatatom = atom->quat[i];
        MathExtra::quatquat(b->quat, orient[i], quatatom);
        MathExtra::qnormalize(quatatom);
      }
      if (eflags[i] & DIPOLE) {
        MathExtra::quat_to_mat(b->quat, p);
        MathExtra::matvec(p, dorient[i], mu[i]);
        MathExtra::snormalize3(mu[i][3], mu[i], mu[i]);
      }
    }
    atomKK->modified(Host, ALL_MASK);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int TRICLINIC, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void FixRigidSmallKokkos<DeviceType>::operator()(TagRigidSmallSetXV<TRICLINIC,EVFLAG>,
                                                  const int &i, value_type result) const
{
  const int ibody = d_atom2body(i);
  if (ibody < 0) return;

  const Body &b = d_body[ibody];

  const int xbox = (d_xcmimage(i) & IMGMASK) - IMGMAX;
  const int ybox = (d_xcmimage(i) >> IMGBITS & IMGMASK) - IMGMAX;
  const int zbox = (d_xcmimage(i) >> IMG2BITS) - IMGMAX;

  const double deltax = xbox*d_prd[0] + (TRICLINIC ? ybox*d_h[5] + zbox*d_h[4] : 0.0);
  const double deltay = ybox*d_prd[1] + (TRICLINIC ? zbox*d_h[3] : 0.0);
  const double deltaz = zbox*d_prd[2];

  double x0, x1, x2, vx, vy, vz;
  if (EVFLAG) {
    x0 = d_x(i,0) + deltax;
    x1 = d_x(i,1) + deltay;
    x2 = d_x(i,2) + deltaz;
    vx = d_v(i,0);
    vy = d_v(i,1);
    vz = d_v(i,2);
  }

  KK_FLOAT disp[3] = {(KK_FLOAT)d_displace(i,0), (KK_FLOAT)d_displace(i,1), (KK_FLOAT)d_displace(i,2)};
  KK_FLOAT xnew[3];
  MathExtraKokkos::matvec(b.ex_space, b.ey_space, b.ez_space, disp, xnew);

  d_v(i,0) = b.omega[1]*xnew[2] - b.omega[2]*xnew[1] + b.vcm[0];
  d_v(i,1) = b.omega[2]*xnew[0] - b.omega[0]*xnew[2] + b.vcm[1];
  d_v(i,2) = b.omega[0]*xnew[1] - b.omega[1]*xnew[0] + b.vcm[2];

  d_x(i,0) = xnew[0] + b.xcm[0] - deltax;
  d_x(i,1) = xnew[1] + b.xcm[1] - deltay;
  d_x(i,2) = xnew[2] + b.xcm[2] - deltaz;

  if (EVFLAG) {
    double massone;
    if (d_rmass.data()) massone = d_rmass(i);
    else massone = d_mass(d_type(i));

    const double fc0 = 0.5*(massone*(d_v(i,0) - vx)/dtf - d_f(i,0));
    const double fc1 = 0.5*(massone*(d_v(i,1) - vy)/dtf - d_f(i,1));
    const double fc2 = 0.5*(massone*(d_v(i,2) - vz)/dtf - d_f(i,2));

    if (vflag_global) {
      result[0] += x0*fc0;
      result[1] += x1*fc1;
      result[2] += x2*fc2;
      result[3] += x0*fc1;
      result[4] += x0*fc2;
      result[5] += x1*fc2;
    }

    if (vflag_atom) {
      Kokkos::atomic_add(&d_vatom(i,0), x0*fc0);
      Kokkos::atomic_add(&d_vatom(i,1), x1*fc1);
      Kokkos::atomic_add(&d_vatom(i,2), x2*fc2);
      Kokkos::atomic_add(&d_vatom(i,3), x0*fc1);
      Kokkos::atomic_add(&d_vatom(i,4), x0*fc2);
      Kokkos::atomic_add(&d_vatom(i,5), x1*fc2);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int TRICLINIC, int EVFLAG>
void FixRigidSmallKokkos<DeviceType>::set_v_kokkos()
{
  int nlocal = atomKK->nlocal;

  double result[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  copymode = 1;
  if (EVFLAG) {
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,
      TagRigidSmallSetV<TRICLINIC,1>>(0, nlocal), *this, result);
  } else {
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,
      TagRigidSmallSetV<TRICLINIC,0>>(0, nlocal), *this, result);
  }
  copymode = 0;

  if (EVFLAG && vflag_global) {
    virial[0] += result[0];
    virial[1] += result[1];
    virial[2] += result[2];
    virial[3] += result[3];
    virial[4] += result[4];
    virial[5] += result[5];
  }

  // extended particles: host fallback
  if (extended) {
    atomKK->sync(Host, ALL_MASK);
    k_body.sync_host();

    double ione[3], exone[3], eyone[3], ezone[3];
    double *shape, *quatatom, *inertiaatom;

    AtomVecEllipsoid::Bonus *ebonus;
    if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
    AtomVecTri::Bonus *tbonus;
    if (avec_tri) tbonus = avec_tri->bonus;
    double **omega_atom = atom->omega;
    double **angmom_atom = atom->angmom;
    double *rmass_h = atom->rmass;
    int *ellipsoid = atom->ellipsoid;
    int *tri = atom->tri;
    int nlocal_h = atom->nlocal;

    for (int i = 0; i < nlocal_h; i++) {
      if (atom2body[i] < 0) continue;
      Body *b = &body[atom2body[i]];

      if (eflags[i] & SPHERE) {
        omega_atom[i][0] = b->omega[0];
        omega_atom[i][1] = b->omega[1];
        omega_atom[i][2] = b->omega[2];
      } else if (eflags[i] & ELLIPSOID) {
        shape = ebonus[ellipsoid[i]].shape;
        quatatom = ebonus[ellipsoid[i]].quat;
        ione[0] = EINERTIA*rmass_h[i] * (shape[1]*shape[1] + shape[2]*shape[2]);
        ione[1] = EINERTIA*rmass_h[i] * (shape[0]*shape[0] + shape[2]*shape[2]);
        ione[2] = EINERTIA*rmass_h[i] * (shape[0]*shape[0] + shape[1]*shape[1]);
        MathExtra::q_to_exyz(quatatom, exone, eyone, ezone);
        MathExtra::omega_to_angmom(b->omega, exone, eyone, ezone, ione, angmom_atom[i]);
      } else if (eflags[i] & LINE) {
        omega_atom[i][0] = b->omega[0];
        omega_atom[i][1] = b->omega[1];
        omega_atom[i][2] = b->omega[2];
      } else if (eflags[i] & TRIANGLE) {
        inertiaatom = tbonus[tri[i]].inertia;
        quatatom = tbonus[tri[i]].quat;
        MathExtra::q_to_exyz(quatatom, exone, eyone, ezone);
        MathExtra::omega_to_angmom(b->omega, exone, eyone, ezone,
                                   inertiaatom, angmom_atom[i]);
      }
    }
    atomKK->modified(Host, ALL_MASK);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int TRICLINIC, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void FixRigidSmallKokkos<DeviceType>::operator()(TagRigidSmallSetV<TRICLINIC,EVFLAG>,
                                                  const int &i, value_type result) const
{
  const int ibody = d_atom2body(i);
  if (ibody < 0) return;

  const Body &b = d_body[ibody];

  KK_FLOAT disp[3] = {(KK_FLOAT)d_displace(i,0), (KK_FLOAT)d_displace(i,1), (KK_FLOAT)d_displace(i,2)};
  KK_FLOAT delta[3];
  MathExtraKokkos::matvec(b.ex_space, b.ey_space, b.ez_space, disp, delta);

  double vx, vy, vz;
  if (EVFLAG) {
    vx = d_v(i,0);
    vy = d_v(i,1);
    vz = d_v(i,2);
  }

  d_v(i,0) = b.omega[1]*delta[2] - b.omega[2]*delta[1] + b.vcm[0];
  d_v(i,1) = b.omega[2]*delta[0] - b.omega[0]*delta[2] + b.vcm[1];
  d_v(i,2) = b.omega[0]*delta[1] - b.omega[1]*delta[0] + b.vcm[2];

  if (EVFLAG) {
    double massone;
    if (d_rmass.data()) massone = d_rmass(i);
    else massone = d_mass(d_type(i));

    const double fc0 = 0.5*(massone*(d_v(i,0) - vx)/dtf - d_f(i,0));
    const double fc1 = 0.5*(massone*(d_v(i,1) - vy)/dtf - d_f(i,1));
    const double fc2 = 0.5*(massone*(d_v(i,2) - vz)/dtf - d_f(i,2));

    const int xbox = (d_xcmimage(i) & IMGMASK) - IMGMAX;
    const int ybox = (d_xcmimage(i) >> IMGBITS & IMGMASK) - IMGMAX;
    const int zbox = (d_xcmimage(i) >> IMG2BITS) - IMGMAX;

    const double x0 = d_x(i,0) + xbox*d_prd[0] + (TRICLINIC ? ybox*d_h[5] + zbox*d_h[4] : 0.0);
    const double x1 = d_x(i,1) + ybox*d_prd[1] + (TRICLINIC ? zbox*d_h[3] : 0.0);
    const double x2 = d_x(i,2) + zbox*d_prd[2];

    if (vflag_global) {
      result[0] += x0*fc0;
      result[1] += x1*fc1;
      result[2] += x2*fc2;
      result[3] += x0*fc1;
      result[4] += x0*fc2;
      result[5] += x1*fc2;
    }

    if (vflag_atom) {
      Kokkos::atomic_add(&d_vatom(i,0), x0*fc0);
      Kokkos::atomic_add(&d_vatom(i,1), x1*fc1);
      Kokkos::atomic_add(&d_vatom(i,2), x2*fc2);
      Kokkos::atomic_add(&d_vatom(i,3), x0*fc1);
      Kokkos::atomic_add(&d_vatom(i,4), x0*fc2);
      Kokkos::atomic_add(&d_vatom(i,5), x1*fc2);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::compute_forces_and_torques_kokkos()
{
  atomKK->sync(execution_space, X_MASK | F_MASK | IMAGE_MASK);
  d_x = atomKK->k_x.template view<DeviceType>();
  d_f = atomKK->k_f.template view<DeviceType>();
  d_image = atomKK->k_image.template view<DeviceType>();

  k_body.template sync<DeviceType>();
  k_atom2body.template sync<DeviceType>();
  k_xcmimage.template sync<DeviceType>();
  d_body = k_body.template view<DeviceType>();
  d_atom2body = k_atom2body.template view<DeviceType>();
  d_xcmimage = k_xcmimage.template view<DeviceType>();

  d_prd = Few<double,3>(domain->prd);
  d_h = Few<double,6>(domain->h);

  int nbody_total = nlocal_body + nghost_body;
  int nlocal = atomKK->nlocal;

  // zero body forces/torques
  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,
    TagRigidSmallComputeForcesTorquesZero>(0, nbody_total), *this);

  // accumulate from atoms
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,
    TagRigidSmallComputeForcesTorques>(0, nlocal), *this);
  copymode = 0;

  // extended particles torque: host
  if (extended) {
    k_body.template modify<DeviceType>();
    k_body.sync_host();
    double **torque_one = atom->torque;
    atomKK->sync(Host, TORQUE_MASK);
    int nlocal_h = atom->nlocal;
    for (int i = 0; i < nlocal_h; i++) {
      if (atom2body[i] < 0) continue;
      if (eflags[i] & TORQUE) {
        double *tcm = body[atom2body[i]].torque;
        tcm[0] += torque_one[i][0];
        tcm[1] += torque_one[i][1];
        tcm[2] += torque_one[i][2];
      }
    }
    k_body.modify_host();
  } else {
    k_body.template modify<DeviceType>();
    k_body.sync_host();
  }

  // reverse comm
  commflag = FORCE_TORQUE;
  comm->reverse_comm(this, 6);
  k_body.modify_host();

  // Langevin + gravity on host (body-level, small count)
  if (langflag) {
    for (int ibody = 0; ibody < nlocal_body; ibody++) {
      double *fcm = body[ibody].fcm;
      fcm[0] += langextra[ibody][0];
      fcm[1] += langextra[ibody][1];
      fcm[2] += langextra[ibody][2];
      double *tcm = body[ibody].torque;
      tcm[0] += langextra[ibody][3];
      tcm[1] += langextra[ibody][4];
      tcm[2] += langextra[ibody][5];
    }
  }

  if (id_gravity) {
    for (int ibody = 0; ibody < nlocal_body; ibody++) {
      double gmass = body[ibody].mass;
      double *fcm = body[ibody].fcm;
      fcm[0] += gvec[0]*gmass;
      fcm[1] += gvec[1]*gmass;
      fcm[2] += gvec[2]*gmass;
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixRigidSmallKokkos<DeviceType>::operator()(TagRigidSmallComputeForcesTorquesZero,
                                                  const int &ibody) const
{
  Body b = d_body[ibody];
  b.fcm[0] = b.fcm[1] = b.fcm[2] = 0.0;
  b.torque[0] = b.torque[1] = b.torque[2] = 0.0;
  d_body[ibody] = b;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixRigidSmallKokkos<DeviceType>::operator()(TagRigidSmallComputeForcesTorques,
                                                  const int &i) const
{
  const int ibody = d_atom2body(i);
  if (ibody < 0) return;

  // accumulate force
  Kokkos::atomic_add(&d_body[ibody].fcm[0], d_f(i,0));
  Kokkos::atomic_add(&d_body[ibody].fcm[1], d_f(i,1));
  Kokkos::atomic_add(&d_body[ibody].fcm[2], d_f(i,2));

  // unwrap position and compute torque
  Few<double,3> x_i;
  x_i[0] = d_x(i,0); x_i[1] = d_x(i,1); x_i[2] = d_x(i,2);
  Few<double,3> unwrap = DomainKokkos::unmap(d_prd, d_h, triclinic, x_i, d_xcmimage(i));

  const double dx = unwrap[0] - d_body[ibody].xcm[0];
  const double dy = unwrap[1] - d_body[ibody].xcm[1];
  const double dz = unwrap[2] - d_body[ibody].xcm[2];

  Kokkos::atomic_add(&d_body[ibody].torque[0], dy*d_f(i,2) - dz*d_f(i,1));
  Kokkos::atomic_add(&d_body[ibody].torque[1], dz*d_f(i,0) - dx*d_f(i,2));
  Kokkos::atomic_add(&d_body[ibody].torque[2], dx*d_f(i,1) - dy*d_f(i,0));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::enforce2d_kokkos()
{
  k_body.template sync<DeviceType>();
  d_body = k_body.template view<DeviceType>();

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,
    TagRigidSmallEnforce2d>(0, nlocal_body), *this);
  copymode = 0;

  k_body.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixRigidSmallKokkos<DeviceType>::operator()(TagRigidSmallEnforce2d,
                                                  const int &ibody) const
{
  Body b = d_body[ibody];
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
  d_body[ibody] = b;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::image_shift_kokkos()
{
  atomKK->sync(execution_space, IMAGE_MASK);
  d_image = atomKK->k_image.template view<DeviceType>();

  k_body.template sync<DeviceType>();
  k_atom2body.template sync<DeviceType>();
  d_body = k_body.template view<DeviceType>();
  d_atom2body = k_atom2body.template view<DeviceType>();

  int nlocal = atomKK->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,
    TagRigidSmallImageShift>(0, nlocal), *this);
  copymode = 0;

  k_xcmimage.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixRigidSmallKokkos<DeviceType>::operator()(TagRigidSmallImageShift,
                                                  const int &i) const
{
  if (d_atom2body(i) < 0) return;
  const Body &b = d_body[d_atom2body(i)];

  imageint tdim, bdim;
  imageint xdim[3];

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

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::setup_pre_neighbor()
{
  atomKK->sync(Host, ALL_MASK);
  k_body.sync_host();
  k_bodyown.sync_host();
  k_bodytag.sync_host();
  k_atom2body.sync_host();
  k_xcmimage.sync_host();

  FixRigidSmall::setup_pre_neighbor();

  k_body.modify_host();
  k_bodyown.modify_host();
  k_bodytag.modify_host();
  k_atom2body.modify_host();
  k_xcmimage.modify_host();
  atomKK->modified(Host, X_MASK | IMAGE_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::pre_neighbor()
{
  k_body.sync_host();
  for (int ibody = 0; ibody < nlocal_body; ibody++) {
    Body *b = &body[ibody];
    domain->remap(b->xcm, b->image);
  }
  k_body.modify_host();

  nghost_body = 0;
  commflag = FULL_BODY;
  comm->forward_comm(this);

  k_body.modify_host();
  k_bodyown.sync_host();
  k_bodytag.sync_host();
  k_atom2body.sync_host();

  reset_atom2body();

  k_atom2body.modify_host();

  image_shift_kokkos();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::grow_arrays(int nmax)
{
  k_bodyown.template sync<DeviceType>();
  k_bodytag.template sync<DeviceType>();
  k_atom2body.template sync<DeviceType>();
  k_xcmimage.template sync<DeviceType>();
  k_displace.template sync<DeviceType>();

  memoryKK->grow_kokkos(k_bodyown, bodyown, nmax, "rigid/small:bodyown");
  memoryKK->grow_kokkos(k_bodytag, bodytag, nmax, "rigid/small:bodytag");
  memoryKK->grow_kokkos(k_atom2body, atom2body, nmax, "rigid/small:atom2body");
  memoryKK->grow_kokkos(k_xcmimage, xcmimage, nmax, "rigid/small:xcmimage");
  memoryKK->grow_kokkos(k_displace, displace, nmax, 3, "rigid/small:displace");

  d_bodyown = k_bodyown.template view<DeviceType>();
  d_bodytag = k_bodytag.template view<DeviceType>();
  d_atom2body = k_atom2body.template view<DeviceType>();
  d_xcmimage = k_xcmimage.template view<DeviceType>();
  d_displace = k_displace.template view<DeviceType>();

  if (extended) {
    k_eflags.template sync<DeviceType>();
    memoryKK->grow_kokkos(k_eflags, eflags, nmax, "rigid/small:eflags");
    d_eflags = k_eflags.template view<DeviceType>();
    if (orientflag) memory->grow(orient, nmax, orientflag, "rigid/small:orient");
    if (dorientflag) memory->grow(dorient, nmax, 3, "rigid/small:dorient");
  }

  if (nmax > maxvatom) {
    maxvatom = nmax;
    memoryKK->destroy_kokkos(k_vatom, vatom);
    memoryKK->create_kokkos(k_vatom, vatom, maxvatom, 6, "fix:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::grow_body_kokkos()
{
  nmax_body += DELTA_BODY;
  k_body.resize(nmax_body);
  body = k_body.h_view.data();
  d_body = k_body.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::copy_arrays(int i, int j, int delflag)
{
  k_bodyown.sync_host();
  k_bodytag.sync_host();
  k_atom2body.sync_host();
  k_xcmimage.sync_host();
  k_displace.sync_host();
  k_body.sync_host();
  if (extended) k_eflags.sync_host();

  FixRigidSmall::copy_arrays(i, j, delflag);

  k_bodyown.modify_host();
  k_bodytag.modify_host();
  k_atom2body.modify_host();
  k_xcmimage.modify_host();
  k_displace.modify_host();
  k_body.modify_host();
  if (extended) k_eflags.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::set_arrays(int i)
{
  k_bodyown.sync_host();
  k_bodytag.sync_host();
  k_atom2body.sync_host();
  k_xcmimage.sync_host();
  k_displace.sync_host();

  FixRigidSmall::set_arrays(i);

  k_bodyown.modify_host();
  k_bodytag.modify_host();
  k_atom2body.modify_host();
  k_xcmimage.modify_host();
  k_displace.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::set_molecule(int nlocalprev, tagint tagprev,
                                                    int imol, double *xgeom,
                                                    double *vcm, double *quat)
{
  k_bodyown.sync_host();
  k_bodytag.sync_host();
  k_displace.sync_host();
  k_body.sync_host();
  if (extended) k_eflags.sync_host();

  FixRigidSmall::set_molecule(nlocalprev, tagprev, imol, xgeom, vcm, quat);

  k_bodyown.modify_host();
  k_bodytag.modify_host();
  k_displace.modify_host();
  k_body.modify_host();
  if (extended) k_eflags.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::pack_exchange(int i, double *buf)
{
  k_bodytag.sync_host();
  k_xcmimage.sync_host();
  k_displace.sync_host();
  k_bodyown.sync_host();
  k_body.sync_host();
  if (extended) k_eflags.sync_host();

  return FixRigidSmall::pack_exchange(i, buf);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::unpack_exchange(int nlocal, double *buf)
{
  int result = FixRigidSmall::unpack_exchange(nlocal, buf);

  k_bodytag.modify_host();
  k_xcmimage.modify_host();
  k_displace.modify_host();
  k_bodyown.modify_host();
  k_body.modify_host();
  if (extended) k_eflags.modify_host();

  return result;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::pack_forward_comm(int n, int *list,
                                                        double *buf, int pbc_flag, int *pbc)
{
  k_bodyown.sync_host();
  k_body.sync_host();
  return FixRigidSmall::pack_forward_comm(n, list, buf, pbc_flag, pbc);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  FixRigidSmall::unpack_forward_comm(n, first, buf);
  k_body.modify_host();
  k_bodyown.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  k_body.sync_host();
  return FixRigidSmall::pack_reverse_comm(n, first, buf);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  FixRigidSmall::unpack_reverse_comm(n, list, buf);
  k_body.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::reset_dt()
{
  FixRigidSmall::reset_dt();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::zero_momentum()
{
  k_body.sync_host();
  FixRigidSmall::zero_momentum();
  k_body.modify_host();
  atomKK->modified(Host, V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::zero_rotation()
{
  k_body.sync_host();
  FixRigidSmall::zero_rotation();
  k_body.modify_host();
  atomKK->modified(Host, V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::sync_body_device()
{
  k_body.template sync<DeviceType>();
  d_body = k_body.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::sync_body_host()
{
  k_body.sync_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::sync_fix_data_device()
{
  k_bodyown.template sync<DeviceType>();
  k_bodytag.template sync<DeviceType>();
  k_atom2body.template sync<DeviceType>();
  k_xcmimage.template sync<DeviceType>();
  k_displace.template sync<DeviceType>();
  k_body.template sync<DeviceType>();

  d_bodyown = k_bodyown.template view<DeviceType>();
  d_bodytag = k_bodytag.template view<DeviceType>();
  d_atom2body = k_atom2body.template view<DeviceType>();
  d_xcmimage = k_xcmimage.template view<DeviceType>();
  d_displace = k_displace.template view<DeviceType>();
  d_body = k_body.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixRigidSmallKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixRigidSmallKokkos<LMPHostType>;
#endif
}
