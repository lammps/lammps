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
   Contributing author: Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#include "fix_rigid_nh_small_kokkos.h"

#include "atom_kokkos.h"
#include "kokkos.h"
#include "atom_masks.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "comm.h"
#include "compute.h"
#include "domain_kokkos.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "math_extra.h"
#include "math_extra_kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "rigid_const.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <type_traits>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace MathExtra;
using namespace RigidConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidNHSmallKokkos<DeviceType>::FixRigidNHSmallKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRigidNHSmall(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  domainKK = (DomainKokkos *) domain;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
  int nmax = atom->nmax;
  int nlocal = atom->nlocal;
  // save bodytag and bodyown filled by the base constructor's create_bodies()
  int *old_bodyown = bodyown;   bodyown = nullptr;
  tagint *old_bodytag = bodytag; bodytag = nullptr;
  memoryKK->destroy_kokkos(k_bodyown, bodyown);
  memoryKK->destroy_kokkos(k_bodytag, bodytag);
  memoryKK->destroy_kokkos(k_atom2body, atom2body);
  memoryKK->destroy_kokkos(k_xcmimage, xcmimage);
  {
    double **old_displace = displace;
    std::vector<double> displace_backup((bigint) nmax * 3);
    for (int i = 0; i < nmax; i++)
      for (int j = 0; j < 3; j++) displace_backup[(bigint) i * 3 + j] = old_displace[i][j];
    memory->destroy(displace);
    k_displace =
        TransformView<KK_FLOAT **, double **, Kokkos::LayoutRight, DeviceType>("rigid/small:displace", nmax, 3);
    double *dh = const_cast<double *>(k_displace.view_host().data());
    memcpy(dh, displace_backup.data(), displace_backup.size() * sizeof(double));
    bigint nbytes = ((bigint) sizeof(double *)) * nmax;
    displace = (double **) memory->smalloc(nbytes, "rigid/small:displace");
    for (int i = 0; i < nmax; i++) displace[i] = &k_displace.view_host()(i, 0);
    k_displace.modify_host();
    k_displace.sync_device();
  }
  memoryKK->create_kokkos(k_bodyown, bodyown, nmax, "rigid/small:bodyown");
  memoryKK->create_kokkos(k_bodytag, bodytag, nmax, "rigid/small:bodytag");
  memoryKK->create_kokkos(k_atom2body, atom2body, nmax, "rigid/small:atom2body");
  memoryKK->create_kokkos(k_xcmimage, xcmimage, nmax, "rigid/small:xcmimage");
  if (nlocal > 0) {
    memcpy(bodyown, old_bodyown, nlocal * sizeof(int));
    memcpy(bodytag, old_bodytag, nlocal * sizeof(tagint));
    k_bodyown.modify_host();
    k_bodytag.modify_host();
  }
  memory->sfree(old_bodyown);
  memory->sfree(old_bodytag);
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
  k_body = TransformView<BodyKokkos*, Body*, Kokkos::LayoutRight, DeviceType>("rigid/small:body", nmax_body);
  if (nmax_body > 0 && body != nullptr) {
    memcpy(k_body.view_host().data(), body, (bigint) nmax_body * sizeof(Body));
    memory->sfree(body);
    body = k_body.view_host().data();
    k_body.modify_host();
    k_body.sync_device();
  }
  d_body = k_body.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidNHSmallKokkos<DeviceType>::~FixRigidNHSmallKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_bodyown, bodyown);
  memoryKK->destroy_kokkos(k_bodytag, bodytag);
  memoryKK->destroy_kokkos(k_atom2body, atom2body);
  memoryKK->destroy_kokkos(k_xcmimage, xcmimage);
  if (displace) {
    memory->sfree(displace);
    displace = nullptr;
  }
  memoryKK->destroy_kokkos(k_displace);
  if (extended) memoryKK->destroy_kokkos(k_eflags, eflags);
  body = nullptr;
  bodyown = nullptr;
  bodytag = nullptr;
  atom2body = nullptr;
  xcmimage = nullptr;
  eflags = nullptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::init()
{
  FixRigidNHSmall::init();
  atomKK->k_mass.modify_host();
  atomKK->k_mass.template sync<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::setup(int vflag)
{
  atomKK->sync(Host, ALL_MASK);
  // Do not sync_host() here: device-side body buffers can be uninitialized before
  // FixRigidSmall::setup() and would overwrite valid host rigid-body data.

  FixRigidSmall::setup(vflag);
  // Match FixRigidSmallKokkos::setup after CPU rigid/small setup
  atomKK->modified(Host, X_MASK | V_MASK);
  k_body.modify_host();
  k_body.sync_device();

  compute_dof();

  double mbody[3];
  akin_t = akin_r = 0.0;
  for (int ibody = 0; ibody < nlocal_body; ibody++) {
    Body *b = &body[ibody];
    MathExtra::transpose_matvec(b->ex_space, b->ey_space, b->ez_space,
                                b->angmom, mbody);
    MathExtra::quatvec(b->quat, mbody, b->conjqm);
    b->conjqm[0] *= 2.0;
    b->conjqm[1] *= 2.0;
    b->conjqm[2] *= 2.0;
    b->conjqm[3] *= 2.0;

    if (tstat_flag || pstat_flag) {
      akin_t += b->mass*(b->vcm[0]*b->vcm[0] + b->vcm[1]*b->vcm[1] +
        b->vcm[2]*b->vcm[2]);
      akin_r += b->angmom[0]*b->omega[0] + b->angmom[1]*b->omega[1] +
        b->angmom[2]*b->omega[2];
    }
  }

  if (tstat_flag || pstat_flag) {
    double ke[2], keall[2];
    ke[0] = akin_t;
    ke[1] = akin_r;
    MPI_Allreduce(ke, keall, 2, MPI_DOUBLE, MPI_SUM, world);
    akin_t = keall[0];
    akin_r = keall[1];
  }

  if (tstat_flag) compute_temp_target();
  else if (pstat_flag) {
    t0 = temperature->compute_scalar();
    if (t0 == 0.0) {
      if (strcmp(update->unit_style, "lj") == 0) t0 = 1.0;
      else t0 = 300.0;
    }
    t_target = t0;
  }

  if (pstat_flag) {
    compute_press_target();
    if (pstyle == ISO) {
      temperature->compute_scalar();
      pressure->compute_scalar();
    } else {
      temperature->compute_vector();
      pressure->compute_vector();
    }
    couple();
    pressure->addstep(update->ntimestep+1);
  }

  double kt, t_mass, tb_mass;
  kt = boltz * t_target;

  if (tstat_flag) {
    t_mass = kt / (t_freq*t_freq);
    q_t[0] = nf_t * t_mass;
    q_r[0] = nf_r * t_mass;
    for (int i = 1; i < t_chain; i++)
      q_t[i] = q_r[i] = t_mass;
    for (int i = 1; i < t_chain; i++) {
      f_eta_t[i] = (q_t[i-1] * eta_dot_t[i-1] * eta_dot_t[i-1] - kt)/q_t[i];
      f_eta_r[i] = (q_r[i-1] * eta_dot_r[i-1] * eta_dot_r[i-1] - kt)/q_r[i];
    }
  }

  int dimension = domain->dimension;

  if (pstat_flag) {
    for (int i = 0; i < 3; i++)
      if (p_flag[i]) {
        epsilon_mass[i] = (g_f + dimension) * kt / (p_freq[i]*p_freq[i]);
        epsilon[i] = log(vol0)/dimension;
      }

    tb_mass = kt / (p_freq_max * p_freq_max);
    q_b[0] = dimension * dimension * tb_mass;
    for (int i = 1; i < p_chain; i++) {
      q_b[i] = tb_mass;
      f_eta_b[i] = (q_b[i] * eta_dot_b[i-1] * eta_dot_b[i-1] - kt)/q_b[i];
    }
  }

  if (tstat_flag || pstat_flag) {
    for (int i = 0; i < t_order; i++) {
      wdti1[i] = w[i] * dtv / t_iter;
      wdti2[i] = wdti1[i] / 2.0;
      wdti4[i] = wdti1[i] / 4.0;
    }
  }

  if (pstat_flag) {
    compute_press_target();
    nh_epsilon_dot();
  }

  k_body.modify_host();
  atomKK->modified(Host, X_MASK | V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::initial_integrate(int vflag)
{
  double tmp;

  d_scale_t[0] = d_scale_t[1] = d_scale_t[2] = 1.0;
  d_scale_v[0] = d_scale_v[1] = d_scale_v[2] = 1.0;
  d_scale_r = 1.0;
  d_dtf2 = dtf * 2.0;

  if (tstat_flag) {
    tmp = exp(-dtq * eta_dot_t[0]);
    d_scale_t[0] = d_scale_t[1] = d_scale_t[2] = tmp;
    tmp = exp(-dtq * eta_dot_r[0]);
    d_scale_r = tmp;
  }

  if (pstat_flag) {
    d_scale_t[0] *= exp(-dtq * (epsilon_dot[0] + mtk_term2));
    d_scale_t[1] *= exp(-dtq * (epsilon_dot[1] + mtk_term2));
    d_scale_t[2] *= exp(-dtq * (epsilon_dot[2] + mtk_term2));
    d_scale_r *= exp(-dtq * (pdim * mtk_term2));

    tmp = dtq * epsilon_dot[0];
    d_scale_v[0] = dtv * exp(tmp) * maclaurin_series(tmp);
    tmp = dtq * epsilon_dot[1];
    d_scale_v[1] = dtv * exp(tmp) * maclaurin_series(tmp);
    tmp = dtq * epsilon_dot[2];
    d_scale_v[2] = dtv * exp(tmp) * maclaurin_series(tmp);
  }

  // sync data to device
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

  k_body.template sync<DeviceType>();
  k_atom2body.template sync<DeviceType>();
  k_displace.template sync<DeviceType>();
  k_xcmimage.template sync<DeviceType>();
  d_body = k_body.template view<DeviceType>();
  d_atom2body = k_atom2body.template view<DeviceType>();
  d_displace = k_displace.template view<DeviceType>();
  d_xcmimage = k_xcmimage.template view<DeviceType>();

  // NH body integration on device
  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,
    TagRigidNHInitialIntegrate>(0, nlocal_body), *this);
  copymode = 0;

  k_body.template modify<DeviceType>();
  k_body.sync_host();

  // forward communicate
  commflag = INITIAL;
  comm->forward_comm(this, 29);
  k_body.sync_device();

  // accumulate kinetic energies
  if (tstat_flag || pstat_flag) {
    copymode = 1;
    auto l_body = d_body;
    KK_ACC_FLOAT ke[2], keall[2];
    Kokkos::parallel_reduce(nlocal_body,
      KOKKOS_LAMBDA(const int &ibody, KK_ACC_FLOAT &l_akin_t, KK_ACC_FLOAT &l_akin_r) {
        BodyKokkos &bk = l_body[ibody];
        l_akin_t += bk.mass*(bk.vcm[0]*bk.vcm[0] + bk.vcm[1]*bk.vcm[1] +
                    bk.vcm[2]*bk.vcm[2]);
        l_akin_r += bk.angmom[0]*bk.omega[0] + bk.angmom[1]*bk.omega[1] +
                    bk.angmom[2]*bk.omega[2];
      }, ke[0], ke[1]);
    copymode = 0;
    MPI_Allreduce(ke, keall, 2, MPI_KK_ACC_FLOAT, MPI_SUM, world);
    akin_t = keall[0];
    akin_r = keall[1];
  }

  if (tstat_flag) {
    compute_temp_target();
    if (dynamic) compute_dof();
    nhc_temp_integrate();
  }

  if (pstat_flag) nhc_press_integrate();

  // virial setup
  v_init(vflag);

  // remap
  if (pstat_flag) remap();

  // set_xv on device
  k_body.template sync<DeviceType>();
  d_body = k_body.template view<DeviceType>();

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
    if (evflag) set_xv_kokkos<1,1>();
    else set_xv_kokkos<1,0>();
  } else {
    if (evflag) set_xv_kokkos<0,1>();
    else set_xv_kokkos<0,0>();
  }

  if (extended) {
    atomKK->sync(Host, ALL_MASK);
    k_body.sync_host();
  }

  // remap again + kspace
  if (pstat_flag) {
    remap();
    if (kspace_flag) {
      atomKK->sync(Host, X_MASK | V_MASK);
      force->kspace->setup();
    }
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
void FixRigidNHSmallKokkos<DeviceType>::operator()(TagRigidNHInitialIntegrate,
                                                    const int &ibody) const
{
  BodyKokkos &bk = d_body[ibody];
  KK_FLOAT mbody[3], tbody[3], fquat[4];

  // step 1.1 - update vcm by 1/2 step
  const KK_FLOAT dtfm = dtf / bk.mass;
  bk.vcm[0] += dtfm * bk.fcm[0];
  bk.vcm[1] += dtfm * bk.fcm[1];
  bk.vcm[2] += dtfm * bk.fcm[2];

  if (tstat_flag || pstat_flag) {
    bk.vcm[0] *= d_scale_t[0];
    bk.vcm[1] *= d_scale_t[1];
    bk.vcm[2] *= d_scale_t[2];
  }

  // step 1.2 - update xcm by full step
  if (!pstat_flag) {
    bk.xcm[0] += dtv * bk.vcm[0];
    bk.xcm[1] += dtv * bk.vcm[1];
    bk.xcm[2] += dtv * bk.vcm[2];
  } else {
    bk.xcm[0] += d_scale_v[0] * bk.vcm[0];
    bk.xcm[1] += d_scale_v[1] * bk.vcm[1];
    bk.xcm[2] += d_scale_v[2] * bk.vcm[2];
  }

  // step 1.3 - apply torque to quaternion momentum
  MathExtraKokkos::transpose_matvec(bk.ex_space, bk.ey_space, bk.ez_space,
                                    bk.torque, tbody);
  MathExtraKokkos::quatvec(bk.quat, tbody, fquat);

  bk.conjqm[0] += d_dtf2 * fquat[0];
  bk.conjqm[1] += d_dtf2 * fquat[1];
  bk.conjqm[2] += d_dtf2 * fquat[2];
  bk.conjqm[3] += d_dtf2 * fquat[3];

  if (tstat_flag || pstat_flag) {
    bk.conjqm[0] *= d_scale_r;
    bk.conjqm[1] *= d_scale_r;
    bk.conjqm[2] *= d_scale_r;
    bk.conjqm[3] *= d_scale_r;
  }

  // step 1.4 to 1.13 - no_squish_rotate
  MathExtraKokkos::no_squish_rotate(3, bk.conjqm, bk.quat, bk.inertia, dtq);
  MathExtraKokkos::no_squish_rotate(2, bk.conjqm, bk.quat, bk.inertia, dtq);
  MathExtraKokkos::no_squish_rotate(1, bk.conjqm, bk.quat, bk.inertia, dtv);
  MathExtraKokkos::no_squish_rotate(2, bk.conjqm, bk.quat, bk.inertia, dtq);
  MathExtraKokkos::no_squish_rotate(3, bk.conjqm, bk.quat, bk.inertia, dtq);

  // update exyz_space, transform p back to angmom, update omega
  MathExtraKokkos::q_to_exyz(bk.quat, bk.ex_space, bk.ey_space, bk.ez_space);
  MathExtraKokkos::invquatvec(bk.quat, bk.conjqm, mbody);
  MathExtraKokkos::matvec(bk.ex_space, bk.ey_space, bk.ez_space, mbody, bk.angmom);

  bk.angmom[0] *= 0.5;
  bk.angmom[1] *= 0.5;
  bk.angmom[2] *= 0.5;

  MathExtraKokkos::angmom_to_omega(bk.angmom, bk.ex_space, bk.ey_space,
                                   bk.ez_space, bk.inertia, bk.omega);

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::final_integrate()
{
  double tmp;
  KK_FLOAT scale_t[3], scale_r_local;
  KK_FLOAT dtf2_local = dtf * 2.0;

  scale_t[0] = scale_t[1] = scale_t[2] = 1.0;
  scale_r_local = 1.0;

  if (tstat_flag) {
    tmp = exp(-1.0 * dtq * eta_dot_t[0]);
    scale_t[0] = scale_t[1] = scale_t[2] = tmp;
    scale_r_local = exp(-1.0 * dtq * eta_dot_r[0]);
  }

  if (pstat_flag) {
    scale_t[0] *= exp(-dtq * (epsilon_dot[0] + mtk_term2));
    scale_t[1] *= exp(-dtq * (epsilon_dot[1] + mtk_term2));
    scale_t[2] *= exp(-dtq * (epsilon_dot[2] + mtk_term2));
    scale_r_local *= exp(-dtq * (pdim * mtk_term2));
  }

  d_scale_t[0] = scale_t[0]; d_scale_t[1] = scale_t[1]; d_scale_t[2] = scale_t[2];
  d_scale_r = scale_r_local;
  d_dtf2 = dtf2_local;

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

  if (!earlyflag) compute_forces_and_torques();
  if (domain->dimension == 2) enforce2d_kokkos();

  k_body.template sync<DeviceType>();
  d_body = k_body.template view<DeviceType>();

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,
    TagRigidNHFinalIntegrate>(0, nlocal_body), *this);
  copymode = 0;

  k_body.template modify<DeviceType>();
  k_body.sync_host();

  commflag = FINAL;
  comm->forward_comm(this, 10);

  k_body.modify_host();

  // accumulate kinetic energies for pstat
  if (pstat_flag) {
    akin_t = akin_r = 0.0;
    for (int ibody = 0; ibody < nlocal_body; ibody++) {
      Body *b = &body[ibody];
      akin_t += b->mass*(b->vcm[0]*b->vcm[0] + b->vcm[1]*b->vcm[1] +
        b->vcm[2]*b->vcm[2]);
      akin_r += b->angmom[0]*b->omega[0] + b->angmom[1]*b->omega[1] +
        b->angmom[2]*b->omega[2];
    }
    double ke[2], keall[2];
    ke[0] = akin_t;
    ke[1] = akin_r;
    MPI_Allreduce(ke, keall, 2, MPI_DOUBLE, MPI_SUM, world);
    akin_t = keall[0];
    akin_r = keall[1];
  }

  // set_v on device
  k_body.template sync<DeviceType>();
  d_body = k_body.template view<DeviceType>();
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
    if (evflag) set_v_kokkos<1,1>();
    else set_v_kokkos<1,0>();
  } else {
    if (evflag) set_v_kokkos<0,1>();
    else set_v_kokkos<0,0>();
  }

  if (extended) {
    atomKK->sync(Host, ALL_MASK);
    k_body.sync_host();
  }

  atomKK->modified(execution_space, V_MASK);

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  // compute temperature
  if (tcomputeflag) {
    atomKK->sync(temperature->execution_space, temperature->datamask_read);
    t_current = temperature->compute_scalar();
    atomKK->modified(temperature->execution_space, temperature->datamask_modify);
    atomKK->sync(execution_space, temperature->datamask_modify);
  }

  // pressure
  if (pstat_flag) {
    if (pstyle == ISO) {
      atomKK->sync(temperature->execution_space, temperature->datamask_read);
      temperature->compute_scalar();
      atomKK->modified(temperature->execution_space, temperature->datamask_modify);
      pressure->compute_scalar();
    } else {
      atomKK->sync(temperature->execution_space, temperature->datamask_read);
      temperature->compute_vector();
      atomKK->modified(temperature->execution_space, temperature->datamask_modify);
      pressure->compute_vector();
    }
    couple();
    pressure->addstep(update->ntimestep+1);
    compute_press_target();
    nh_epsilon_dot();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixRigidNHSmallKokkos<DeviceType>::operator()(TagRigidNHFinalIntegrate,
                                                    const int &ibody) const
{
  BodyKokkos &bk = d_body[ibody];
  KK_FLOAT mbody[3], tbody[3], fquat[4];

  const KK_FLOAT dtfm = dtf / bk.mass;
  if (tstat_flag || pstat_flag) {
    bk.vcm[0] *= d_scale_t[0];
    bk.vcm[1] *= d_scale_t[1];
    bk.vcm[2] *= d_scale_t[2];
  }

  bk.vcm[0] += dtfm * bk.fcm[0];
  bk.vcm[1] += dtfm * bk.fcm[1];
  bk.vcm[2] += dtfm * bk.fcm[2];

  MathExtraKokkos::transpose_matvec(bk.ex_space, bk.ey_space,
                                    bk.ez_space, bk.torque, tbody);
  MathExtraKokkos::quatvec(bk.quat, tbody, fquat);

  if (tstat_flag || pstat_flag) {
    bk.conjqm[0] = d_scale_r * bk.conjqm[0] + d_dtf2 * fquat[0];
    bk.conjqm[1] = d_scale_r * bk.conjqm[1] + d_dtf2 * fquat[1];
    bk.conjqm[2] = d_scale_r * bk.conjqm[2] + d_dtf2 * fquat[2];
    bk.conjqm[3] = d_scale_r * bk.conjqm[3] + d_dtf2 * fquat[3];
  } else {
    bk.conjqm[0] += d_dtf2 * fquat[0];
    bk.conjqm[1] += d_dtf2 * fquat[1];
    bk.conjqm[2] += d_dtf2 * fquat[2];
    bk.conjqm[3] += d_dtf2 * fquat[3];
  }

  MathExtraKokkos::invquatvec(bk.quat, bk.conjqm, mbody);
  MathExtraKokkos::matvec(bk.ex_space, bk.ey_space, bk.ez_space, mbody, bk.angmom);

  bk.angmom[0] *= 0.5;
  bk.angmom[1] *= 0.5;
  bk.angmom[2] *= 0.5;

  MathExtraKokkos::angmom_to_omega(bk.angmom, bk.ex_space, bk.ey_space,
                                   bk.ez_space, bk.inertia, bk.omega);

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixRigidNHSmallKokkos<DeviceType>::operator()(TagRigidNHAccumKE,
                                                    const int &ibody,
                                                    double &akin_t_local,
                                                    double &akin_r_local) const
{
  const BodyKokkos &bk = d_body[ibody];
  akin_t_local += bk.mass*(bk.vcm[0]*bk.vcm[0] + bk.vcm[1]*bk.vcm[1] +
    bk.vcm[2]*bk.vcm[2]);
  akin_r_local += bk.angmom[0]*bk.omega[0] + bk.angmom[1]*bk.omega[1] +
    bk.angmom[2]*bk.omega[2];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::remap()
{
  int i;
  double oldlo, oldhi, ctr, expfac;

  int nlocal = atom->nlocal;

  for (i = 0; i < 3; i++) epsilon[i] += dtq * epsilon_dot[i];

  // Use Kokkos domain remapping on device coordinates (same pattern as fix nh/kk).
  // CPU domain->x2lamda on atom->x bypasses the atom coord TransformView and
  // breaks sync with Kokkos device buffers.
  if (allremap)
    domainKK->x2lamda(nlocal);
  else
    domainKK->x2lamda(nlocal, dilate_group_bit);

  for (auto &ifix : rfix) ifix->deform(0);

  for (i = 0; i < 3; i++) {
    if (p_flag[i]) {
      oldlo = domain->boxlo[i];
      oldhi = domain->boxhi[i];
      ctr = 0.5 * (oldlo + oldhi);
      expfac = exp(dtq * epsilon_dot[i]);
      domain->boxlo[i] = (oldlo-ctr)*expfac + ctr;
      domain->boxhi[i] = (oldhi-ctr)*expfac + ctr;
    }
  }

  domain->set_global_box();
  domain->set_local_box();

  if (allremap)
    domainKK->lamda2x(nlocal);
  else
    domainKK->lamda2x(nlocal, dilate_group_bit);

  for (auto &ifix : rfix) ifix->deform(1);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::post_force(int /*vflag*/)
{
  if (langflag) {
    atomKK->sync(Host, ALL_MASK);
    k_body.sync_host();
    apply_langevin_thermostat();
    k_body.modify_host();
  }
  if (earlyflag) compute_forces_and_torques();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::compute_forces_and_torques()
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

  d_prd = Few<KK_FLOAT,3>(domain->prd);
  d_h = Few<KK_FLOAT,6>(domain->h);

  int nbody_total = nlocal_body + nghost_body;
  int nlocal = atomKK->nlocal;

  copymode = 1;
  auto l_body = d_body;
  Kokkos::parallel_for(nbody_total, KOKKOS_LAMBDA(const int &ibody) {
    BodyKokkos &bk = l_body[ibody];
    bk.fcm[0] = bk.fcm[1] = bk.fcm[2] = 0.0;
    bk.torque[0] = bk.torque[1] = bk.torque[2] = 0.0;
  });
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,
    TagRigidNHComputeForcesTorques>(0, nlocal), *this);
  copymode = 0;

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

  commflag = FORCE_TORQUE;
  comm->reverse_comm(this, 6);
  k_body.modify_host();

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
void FixRigidNHSmallKokkos<DeviceType>::operator()(TagRigidNHComputeForcesTorques,
                                                    const int &i) const
{
  const int ibody = d_atom2body(i);
  if (ibody < 0) return;
  BodyKokkos &bk = d_body[ibody];
  Kokkos::atomic_add(&bk.fcm[0], d_f(i,0));
  Kokkos::atomic_add(&bk.fcm[1], d_f(i,1));
  Kokkos::atomic_add(&bk.fcm[2], d_f(i,2));

  Few<KK_FLOAT,3> x_i;
  x_i[0] = d_x(i,0); x_i[1] = d_x(i,1); x_i[2] = d_x(i,2);
  Few<KK_FLOAT,3> unwrap = DomainKokkos::unmap(d_prd, d_h, triclinic, x_i, d_xcmimage(i));

  const KK_FLOAT dx = unwrap[0] - bk.xcm[0];
  const KK_FLOAT dy = unwrap[1] - bk.xcm[1];
  const KK_FLOAT dz = unwrap[2] - bk.xcm[2];

  Kokkos::atomic_add(&bk.torque[0], dy*d_f(i,2) - dz*d_f(i,1));
  Kokkos::atomic_add(&bk.torque[1], dz*d_f(i,0) - dx*d_f(i,2));
  Kokkos::atomic_add(&bk.torque[2], dx*d_f(i,1) - dy*d_f(i,0));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::enforce2d_kokkos()
{
  k_body.template sync<DeviceType>();
  d_body = k_body.template view<DeviceType>();
  copymode = 1;
  auto l_body = d_body;
  Kokkos::parallel_for(nlocal_body, KOKKOS_LAMBDA(const int &ibody) {
    BodyKokkos &bk = l_body[ibody];
    bk.xcm[2] = bk.vcm[2] = bk.fcm[2] = bk.xgc[2] = 0.0;
    bk.torque[0] = bk.torque[1] = 0.0;
    bk.angmom[0] = bk.angmom[1] = 0.0;
    bk.omega[0] = bk.omega[1] = 0.0;
  });
  copymode = 0;
  k_body.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::image_shift_kokkos()
{
  atomKK->sync(execution_space, IMAGE_MASK);
  k_body.template sync<DeviceType>();
  k_atom2body.template sync<DeviceType>();
  auto l_image = atomKK->k_image.template view<DeviceType>();
  auto l_xcmimage = d_xcmimage;
  auto l_body = d_body;
  auto l_atom2body = d_atom2body;
  copymode = 1;
  Kokkos::parallel_for(atomKK->nlocal, KOKKOS_LAMBDA(const int &i) {
    if (l_atom2body(i) < 0) return;
    const BodyKokkos &bk = l_body[l_atom2body(i)];
    imageint tdim = l_image(i) & IMGMASK;
    imageint bdim = bk.image & IMGMASK;
    const imageint xdim0 = IMGMAX + tdim - bdim;
    tdim = (l_image(i) >> IMGBITS) & IMGMASK;
    bdim = (bk.image >> IMGBITS) & IMGMASK;
    const imageint xdim1 = IMGMAX + tdim - bdim;
    tdim = l_image(i) >> IMG2BITS;
    bdim = bk.image >> IMG2BITS;
    const imageint xdim2 = IMGMAX + tdim - bdim;
    l_xcmimage(i) = (xdim2 << IMG2BITS) | (xdim1 << IMGBITS) | xdim0;
  });
  copymode = 0;
  k_xcmimage.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int TRICLINIC, int EVFLAG>
void FixRigidNHSmallKokkos<DeviceType>::set_xv_kokkos()
{
  int nlocal = atomKK->nlocal;

  copymode = 1;
  if constexpr (!EVFLAG) {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,
      TagRigidNHSetXV<TRICLINIC,HALF,0>>(0, nlocal), *this);
  } else {
    int neighflag = lmp->kokkos->neighflag;
    if (neighflag == FULL) {
      neighflag =
        (lmp->kokkos->nthreads > 1 || lmp->kokkos->ngpus > 0) ? HALFTHREAD : HALF;
    }

    const bool need_dup =
      (neighflag != HALF) &&
      std::is_same_v<NeedDup_v<HALFTHREAD, DeviceType>,
                     Kokkos::Experimental::ScatterDuplicated>;

    if (vflag_atom) {
      if (need_dup)
        dup_vatom = Kokkos::Experimental::create_scatter_view<
          Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
      else
        ndup_vatom = Kokkos::Experimental::create_scatter_view<
          Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
    }

    EV_FLOAT ev{};
    if (neighflag == HALF) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,
        TagRigidNHSetXV<TRICLINIC,HALF,1>>(0, nlocal), *this, ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,
        TagRigidNHSetXV<TRICLINIC,HALFTHREAD,1>>(0, nlocal), *this, ev);
    }
    if (vflag_global) {
      virial[0] += static_cast<double>(ev.v[0]);
      virial[1] += static_cast<double>(ev.v[1]);
      virial[2] += static_cast<double>(ev.v[2]);
      virial[3] += static_cast<double>(ev.v[3]);
      virial[4] += static_cast<double>(ev.v[4]);
      virial[5] += static_cast<double>(ev.v[5]);
    }
    if (vflag_atom && need_dup) Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    if (vflag_atom && need_dup) dup_vatom = {};
  }
  copymode = 0;

  // update geometric center
  k_body.sync_host();
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
}

template<class DeviceType>
template<int TRICLINIC, int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void FixRigidNHSmallKokkos<DeviceType>::operator()(TagRigidNHSetXV<TRICLINIC,NEIGHFLAG,EVFLAG>,
                                                    const int &i) const
{
  EV_FLOAT ev;
  this->template operator()(TagRigidNHSetXV<TRICLINIC,NEIGHFLAG,EVFLAG>(), i, ev);
}

template<class DeviceType>
template<int TRICLINIC, int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void FixRigidNHSmallKokkos<DeviceType>::operator()(TagRigidNHSetXV<TRICLINIC,NEIGHFLAG,EVFLAG>,
                                                    const int &i, EV_FLOAT &ev) const
{
  const int ibody = d_atom2body(i);
  if (ibody < 0) return;

  const BodyKokkos &bk = d_body[ibody];

  const int xbox = (d_xcmimage(i) & IMGMASK) - IMGMAX;
  const int ybox = (d_xcmimage(i) >> IMGBITS & IMGMASK) - IMGMAX;
  const int zbox = (d_xcmimage(i) >> IMG2BITS) - IMGMAX;

  const double deltax = xbox*d_prd[0] + (TRICLINIC ? ybox*d_h[5] + zbox*d_h[4] : 0.0);
  const double deltay = ybox*d_prd[1] + (TRICLINIC ? zbox*d_h[3] : 0.0);
  const double deltaz = zbox*d_prd[2];

  double x0 = 0.0, x1 = 0.0, x2 = 0.0, vx = 0.0, vy = 0.0, vz = 0.0;
  if constexpr (EVFLAG) {
    x0 = d_x(i,0) + deltax;
    x1 = d_x(i,1) + deltay;
    x2 = d_x(i,2) + deltaz;
    vx = d_v(i,0);
    vy = d_v(i,1);
    vz = d_v(i,2);
  }

  KK_FLOAT disp[3] = {(KK_FLOAT)d_displace(i,0), (KK_FLOAT)d_displace(i,1), (KK_FLOAT)d_displace(i,2)};
  KK_FLOAT xnew[3];
  MathExtraKokkos::matvec(bk.ex_space, bk.ey_space, bk.ez_space, disp, xnew);

  d_v(i,0) = bk.omega[1]*xnew[2] - bk.omega[2]*xnew[1] + bk.vcm[0];
  d_v(i,1) = bk.omega[2]*xnew[0] - bk.omega[0]*xnew[2] + bk.vcm[1];
  d_v(i,2) = bk.omega[0]*xnew[1] - bk.omega[1]*xnew[0] + bk.vcm[2];

  d_x(i,0) = xnew[0] + bk.xcm[0] - deltax;
  d_x(i,1) = xnew[1] + bk.xcm[1] - deltay;
  d_x(i,2) = xnew[2] + bk.xcm[2] - deltaz;

  if constexpr (EVFLAG) {
    double massone;
    if (d_rmass.data()) massone = d_rmass(i);
    else massone = d_mass(d_type(i));

    const double fc0 = 0.5*(massone*(d_v(i,0) - vx)/dtf - d_f(i,0));
    const double fc1 = 0.5*(massone*(d_v(i,1) - vy)/dtf - d_f(i,1));
    const double fc2 = 0.5*(massone*(d_v(i,2) - vz)/dtf - d_f(i,2));

    const KK_ACC_FLOAT vd00 = static_cast<KK_ACC_FLOAT>(x0*fc0);
    const KK_ACC_FLOAT vd11 = static_cast<KK_ACC_FLOAT>(x1*fc1);
    const KK_ACC_FLOAT vd22 = static_cast<KK_ACC_FLOAT>(x2*fc2);
    const KK_ACC_FLOAT vd01 = static_cast<KK_ACC_FLOAT>(x0*fc1);
    const KK_ACC_FLOAT vd02 = static_cast<KK_ACC_FLOAT>(x0*fc2);
    const KK_ACC_FLOAT vd12 = static_cast<KK_ACC_FLOAT>(x1*fc2);

    if (vflag_global) {
      ev.v[0] += vd00;
      ev.v[1] += vd11;
      ev.v[2] += vd22;
      ev.v[3] += vd01;
      ev.v[4] += vd02;
      ev.v[5] += vd12;
    }
    if (vflag_atom) {
      auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,
        decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom, ndup_vatom);
      auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
      a_vatom(i,0) += vd00;
      a_vatom(i,1) += vd11;
      a_vatom(i,2) += vd22;
      a_vatom(i,3) += vd01;
      a_vatom(i,4) += vd02;
      a_vatom(i,5) += vd12;
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int TRICLINIC, int EVFLAG>
void FixRigidNHSmallKokkos<DeviceType>::set_v_kokkos()
{
  int nlocal = atomKK->nlocal;

  copymode = 1;
  if constexpr (!EVFLAG) {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,
      TagRigidNHSetV<TRICLINIC,HALF,0>>(0, nlocal), *this);
  } else {
    int neighflag = lmp->kokkos->neighflag;
    if (neighflag == FULL) {
      neighflag =
        (lmp->kokkos->nthreads > 1 || lmp->kokkos->ngpus > 0) ? HALFTHREAD : HALF;
    }

    const bool need_dup =
      (neighflag != HALF) &&
      std::is_same_v<NeedDup_v<HALFTHREAD, DeviceType>,
                     Kokkos::Experimental::ScatterDuplicated>;

    if (vflag_atom) {
      if (need_dup)
        dup_vatom = Kokkos::Experimental::create_scatter_view<
          Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
      else
        ndup_vatom = Kokkos::Experimental::create_scatter_view<
          Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
    }

    EV_FLOAT ev{};
    if (neighflag == HALF) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,
        TagRigidNHSetV<TRICLINIC,HALF,1>>(0, nlocal), *this, ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,
        TagRigidNHSetV<TRICLINIC,HALFTHREAD,1>>(0, nlocal), *this, ev);
    }

    if (vflag_global) {
      virial[0] += static_cast<double>(ev.v[0]);
      virial[1] += static_cast<double>(ev.v[1]);
      virial[2] += static_cast<double>(ev.v[2]);
      virial[3] += static_cast<double>(ev.v[3]);
      virial[4] += static_cast<double>(ev.v[4]);
      virial[5] += static_cast<double>(ev.v[5]);
    }

    if (vflag_atom && need_dup) Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    if (vflag_atom && need_dup) dup_vatom = {};
  }
  copymode = 0;
}

template<class DeviceType>
template<int TRICLINIC, int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void FixRigidNHSmallKokkos<DeviceType>::operator()(TagRigidNHSetV<TRICLINIC,NEIGHFLAG,EVFLAG>,
                                                    const int &i) const
{
  EV_FLOAT ev;
  this->template operator()(TagRigidNHSetV<TRICLINIC,NEIGHFLAG,EVFLAG>(), i, ev);
}

template<class DeviceType>
template<int TRICLINIC, int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void FixRigidNHSmallKokkos<DeviceType>::operator()(TagRigidNHSetV<TRICLINIC,NEIGHFLAG,EVFLAG>,
                                                    const int &i, EV_FLOAT &ev) const
{
  const int ibody = d_atom2body(i);
  if (ibody < 0) return;

  const BodyKokkos &bk = d_body[ibody];

  KK_FLOAT disp[3] = {(KK_FLOAT)d_displace(i,0), (KK_FLOAT)d_displace(i,1), (KK_FLOAT)d_displace(i,2)};
  KK_FLOAT delta[3];
  MathExtraKokkos::matvec(bk.ex_space, bk.ey_space, bk.ez_space, disp, delta);

  double vx = 0.0, vy = 0.0, vz = 0.0;
  if constexpr (EVFLAG) {
    vx = d_v(i,0);
    vy = d_v(i,1);
    vz = d_v(i,2);
  }

  d_v(i,0) = bk.omega[1]*delta[2] - bk.omega[2]*delta[1] + bk.vcm[0];
  d_v(i,1) = bk.omega[2]*delta[0] - bk.omega[0]*delta[2] + bk.vcm[1];
  d_v(i,2) = bk.omega[0]*delta[1] - bk.omega[1]*delta[0] + bk.vcm[2];

  if constexpr (EVFLAG) {
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

    const KK_ACC_FLOAT vd00 = static_cast<KK_ACC_FLOAT>(x0*fc0);
    const KK_ACC_FLOAT vd11 = static_cast<KK_ACC_FLOAT>(x1*fc1);
    const KK_ACC_FLOAT vd22 = static_cast<KK_ACC_FLOAT>(x2*fc2);
    const KK_ACC_FLOAT vd01 = static_cast<KK_ACC_FLOAT>(x0*fc1);
    const KK_ACC_FLOAT vd02 = static_cast<KK_ACC_FLOAT>(x0*fc2);
    const KK_ACC_FLOAT vd12 = static_cast<KK_ACC_FLOAT>(x1*fc2);

    if (vflag_global) {
      ev.v[0] += vd00;
      ev.v[1] += vd11;
      ev.v[2] += vd22;
      ev.v[3] += vd01;
      ev.v[4] += vd02;
      ev.v[5] += vd12;
    }
    if (vflag_atom) {
      auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,
        decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom, ndup_vatom);
      auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
      a_vatom(i,0) += vd00;
      a_vatom(i,1) += vd11;
      a_vatom(i,2) += vd22;
      a_vatom(i,3) += vd01;
      a_vatom(i,4) += vd02;
      a_vatom(i,5) += vd12;
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::setup_pre_neighbor()
{
  atomKK->sync(Host, ALL_MASK);
  // Avoid sync_host() on fix-owned DualViews before setup_bodies_static():
  // Kokkos device buffers are uninitialized and would overwrite valid host data.

  FixRigidSmall::setup_pre_neighbor();

  k_body.modify_host();
  k_bodyown.modify_host();
  k_bodytag.modify_host();
  k_atom2body.modify_host();
  k_xcmimage.modify_host();
  k_displace.modify_host();
  if (extended) k_eflags.modify_host();
  atomKK->modified(Host, X_MASK | IMAGE_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::pre_neighbor()
{
  // Device-side body[] is uninitialized until integration; do not pull it onto host.

  for (int ibody = 0; ibody < nlocal_body; ibody++) {
    Body *b = &body[ibody];
    domain->remap(b->xcm, b->image);
  }
  k_body.modify_host();

  nghost_body = 0;
  commflag = FULL_BODY;
  comm->forward_comm(this);

  k_body.modify_host();
  k_body.sync_device();
  k_bodyown.sync_host();
  k_bodytag.sync_host();
  k_atom2body.sync_host();

  reset_atom2body();
  k_atom2body.modify_host();
  image_shift_kokkos();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::grow_arrays(int nmax)
{
  // Do not sync from device before grow: uninitialized device data must not
  // overwrite the host arrays this routine is about to resize.

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
void FixRigidNHSmallKokkos<DeviceType>::grow_body()
{
  nmax_body += DELTA_BODY;
  k_body.resize(nmax_body);
  body = k_body.view_host().data();
  d_body = k_body.view_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::copy_arrays(int i, int j, int delflag)
{
  k_bodyown.sync_host(); k_bodytag.sync_host(); k_atom2body.sync_host();
  k_xcmimage.sync_host(); k_displace.sync_host(); k_body.sync_host();
  if (extended) k_eflags.sync_host();

  FixRigidSmall::copy_arrays(i, j, delflag);

  k_bodyown.modify_host(); k_bodytag.modify_host(); k_atom2body.modify_host();
  k_xcmimage.modify_host(); k_displace.modify_host(); k_body.modify_host();
  if (extended) k_eflags.modify_host();
}

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::set_arrays(int i)
{
  k_bodyown.sync_host(); k_bodytag.sync_host(); k_atom2body.sync_host();
  k_xcmimage.sync_host(); k_displace.sync_host();

  FixRigidSmall::set_arrays(i);

  k_bodyown.modify_host(); k_bodytag.modify_host(); k_atom2body.modify_host();
  k_xcmimage.modify_host(); k_displace.modify_host();
}

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::set_molecule(int nlocalprev, tagint tagprev,
                                                      int imol, double *xgeom,
                                                      double *vcm, double *quat)
{
  k_bodyown.sync_host(); k_bodytag.sync_host(); k_displace.sync_host(); k_body.sync_host();
  if (extended) k_eflags.sync_host();

  FixRigidSmall::set_molecule(nlocalprev, tagprev, imol, xgeom, vcm, quat);

  k_bodyown.modify_host(); k_bodytag.modify_host(); k_displace.modify_host(); k_body.modify_host();
  if (extended) k_eflags.modify_host();
}

template<class DeviceType>
int FixRigidNHSmallKokkos<DeviceType>::pack_exchange(int i, double *buf)
{
  k_bodytag.sync_host(); k_xcmimage.sync_host(); k_displace.sync_host();
  k_bodyown.sync_host(); k_body.sync_host();
  if (extended) k_eflags.sync_host();
  return FixRigidSmall::pack_exchange(i, buf);
}

template<class DeviceType>
int FixRigidNHSmallKokkos<DeviceType>::unpack_exchange(int nlocal, double *buf)
{
  int result = FixRigidSmall::unpack_exchange(nlocal, buf);
  k_bodytag.modify_host(); k_xcmimage.modify_host(); k_displace.modify_host();
  k_bodyown.modify_host(); k_body.modify_host();
  if (extended) k_eflags.modify_host();
  return result;
}

template<class DeviceType>
int FixRigidNHSmallKokkos<DeviceType>::pack_forward_comm(int n, int *list,
                                                          double *buf, int pbc_flag, int *pbc)
{
  k_bodyown.sync_host();
  return FixRigidSmall::pack_forward_comm(n, list, buf, pbc_flag, pbc);
}

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  FixRigidSmall::unpack_forward_comm(n, first, buf);
  k_body.modify_host(); k_bodyown.modify_host();
}

template<class DeviceType>
int FixRigidNHSmallKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  k_body.sync_host();
  return FixRigidSmall::pack_reverse_comm(n, first, buf);
}

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  FixRigidSmall::unpack_reverse_comm(n, list, buf);
  k_body.modify_host();
}

template<class DeviceType>
double FixRigidNHSmallKokkos<DeviceType>::compute_scalar()
{
  k_body.sync_host();
  return FixRigidNHSmall::compute_scalar();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::zero_momentum()
{
  copymode = 1;
  auto l_body = d_body;
  Kokkos::parallel_for(nlocal_body + nghost_body, KOKKOS_LAMBDA(const int &ibody) {
    BodyKokkos &bk = l_body[ibody];
    bk.vcm[0] = bk.vcm[1] = bk.vcm[2] = 0.0;
  });
  copymode = 0;
  k_body.modify_device();
  k_body.sync_host();
  // forward communicate of omega to all ghost copies
  commflag = FINAL;
  comm->forward_comm(this,10);
  k_body.modify_host();
  k_body.sync_device();

  // set velocity of atoms in rigid bodues
  if (triclinic) set_v_kokkos<1,0>();
  else set_v_kokkos<0,0>();

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::zero_rotation()
{
  copymode = 1;
  auto l_body = d_body;
  Kokkos::parallel_for(nlocal_body + nghost_body, KOKKOS_LAMBDA(const int &ibody) {
    BodyKokkos &bk = l_body[ibody];
    bk.angmom[0] = bk.angmom[1] = bk.angmom[2] = 0.0;
    bk.omega[0] = bk.omega[1] = bk.omega[2] = 0.0;
  });
  copymode = 0;
  k_body.modify_device();
  k_body.sync_host();

  // forward communicate of omega to all ghost copies
  commflag = FINAL;
  comm->forward_comm(this,10);

  // set velocity of atoms in rigid bodues
  if (triclinic) set_v_kokkos<1,0>();
  else set_v_kokkos<0,0>();

}


/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixRigidNHSmallKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixRigidNHSmallKokkos<LMPHostType>;
#endif
}
