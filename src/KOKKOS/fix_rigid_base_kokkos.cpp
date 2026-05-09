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

#include "fix_rigid_base_kokkos.h"

#include "fix.h"
#include "fix_rigid.h"
#include "fix_rigid_small.h"
#include "fix_rigid_nh.h"
#include "fix_rigid_nh_small.h"
#include "fix_rigid_small_kokkos.h"

#include "atom_kokkos.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "compute.h"
#include "domain_kokkos.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "math_extra_kokkos.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "molecule.h"
#include "rigid_const.h"
#include "update.h"

#include <type_traits>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace RigidConst;

using MathExtraKokkos::angmom_to_omega;
using MathExtraKokkos::exyz_to_q;
using MathExtraKokkos::invquatvec;
using MathExtraKokkos::matvec;
using MathExtraKokkos::no_squish_rotate;
using MathExtraKokkos::q_to_exyz;
using MathExtraKokkos::quat_to_mat;
using MathExtraKokkos::quatvec;
using MathExtraKokkos::richardson;
using MathExtraKokkos::transpose_matvec;

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
FixRigidBaseKokkos<DeviceType,FixRigidBase>::
  FixRigidBaseKokkos(LAMMPS *lmp, int narg, char **arg) :
    FixRigidBase(lmp, narg, arg), KokkosBase()
{
  kokkosable = 1;
  atomKK = static_cast<AtomKokkos*>(atom);
  domainKK = static_cast<DomainKokkos*>(domain);
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  forward_comm_device = 1;
  reverse_comm_device = 1;
  exchange_comm_device = 1;
  sort_device = 1;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
FixRigidBaseKokkos<DeviceType,FixRigidBase>::~FixRigidBaseKokkos()
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
#ifdef LMP_KOKKOS_DEBUG_RNG
  rand_pool.destroy();
#endif
}

/* ----------------------------------------------------------------------
   FIX METHODS
------------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::post_constructor()
{
  const int nmax = atomKK->nmax;
  const int nlocal = atomKK->nlocal;
  // save bodytag and bodyown filled by the base constructor's create_bodies()
  int *old_bodyown = bodyown;
  bodyown = nullptr;
  tagint *old_bodytag = bodytag;
  bodytag = nullptr;
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
    k_displace = TransformView<KK_FLOAT**, double**, Kokkos::LayoutRight, DeviceType>("rigid/small:displace", nmax, 3);
    double *dh = const_cast<double *>(k_displace.view_host().data());
    memcpy(dh, displace_backup.data(), displace_backup.size() * sizeof(double));
    const bigint nbytes = ((bigint) sizeof(double *)) * nmax;
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
  if (extended) {
    memoryKK->destroy_kokkos(k_eflags, eflags);
    memoryKK->create_kokkos(k_eflags, eflags, nmax, "rigid/small:eflags");
  }
  k_body = TransformView<BodyKokkos*, Body*, Kokkos::LayoutRight, DeviceType>("rigid/small:body", nmax_body);
  if (nmax_body > 0 && body != nullptr) {
    memcpy(k_body.view_host().data(), body, (bigint) nmax_body * sizeof(Body));
    memory->sfree(body);
    body = k_body.view_host().data();
    k_body.modify_host();
    k_body.sync_device();
  }

#ifdef LMP_KOKKOS_DEBUG_RNG
  this->rand_pool = Kokkos::Random_XorShift64_Pool<DeviceType>(seed + comm->me, lmp);
#else
  this->rand_pool = Kokkos::Random_XorShift64_Pool<DeviceType>(seed + comm->me);
#endif

}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::init()
{
  FixRigidBase::init();
  atomKK->k_mass.modify_host();
  atomKK->k_mass.template sync<DeviceType>();
#ifdef LMP_KOKKOS_DEBUG_RNG
  rand_pool.init(random,seed + comm->me);
#endif
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::pre_neighbor()
{
  k_body.sync_device();
  auto l_body = k_body.template view<DeviceType>();
  auto l_xperiodic = domainKK->xperiodic;
  auto l_yperiodic = domainKK->yperiodic;
  auto l_zperiodic = domainKK->zperiodic;
  auto l_lo0 = static_cast<KK_FLOAT>(domainKK->boxlo[0]);
  auto l_lo1 = static_cast<KK_FLOAT>(domainKK->boxlo[1]);
  auto l_lo2 = static_cast<KK_FLOAT>(domainKK->boxlo[2]);
  auto l_hi0 = static_cast<KK_FLOAT>(domainKK->boxhi[0]);
  auto l_hi1 = static_cast<KK_FLOAT>(domainKK->boxhi[1]);
  auto l_hi2 = static_cast<KK_FLOAT>(domainKK->boxhi[2]);
  auto l_prd0 = static_cast<KK_FLOAT>(domainKK->prd[0]);
  auto l_prd1 = static_cast<KK_FLOAT>(domainKK->prd[1]);
  auto l_prd2 = static_cast<KK_FLOAT>(domainKK->prd[2]);

  copymode = 1;
  Kokkos::parallel_for("rigid/small:pre_neighbor",
    Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      imageint idim, otherdims;

      if (l_xperiodic) {
        while (bk.xcm[0] < l_lo0) {
          bk.xcm[0] += l_prd0;
          idim = bk.image & IMGMASK;
          otherdims = bk.image ^ idim;
          idim--;
          idim &= IMGMASK;
          bk.image = otherdims | idim;
        }
        while (bk.xcm[0] >= l_hi0) {
          bk.xcm[0] -= l_prd0;
          idim = bk.image & IMGMASK;
          otherdims = bk.image ^ idim;
          idim++;
          idim &= IMGMASK;
          bk.image = otherdims | idim;
        }
        bk.xcm[0] = MAX(bk.xcm[0], l_lo0);
      }
      if (l_yperiodic) {
        while (bk.xcm[1] < l_lo1) {
          bk.xcm[1] += l_prd1;
          idim = (bk.image >> IMGBITS) & IMGMASK;
          otherdims = bk.image ^ (idim << IMGBITS);
          idim--;
          idim &= IMGMASK;
          bk.image = otherdims | (idim << IMGBITS);
        }
        while (bk.xcm[1] >= l_hi1) {
          bk.xcm[1] -= l_prd1;
          idim = (bk.image >> IMGBITS) & IMGMASK;
          otherdims = bk.image ^ (idim << IMGBITS);
          idim++;
          idim &= IMGMASK;
          bk.image = otherdims | (idim << IMGBITS);
        }
        bk.xcm[1] = MAX(bk.xcm[1], l_lo1);
      }
      if (l_zperiodic) {
        while (bk.xcm[2] < l_lo2) {
          bk.xcm[2] += l_prd2;
          idim = bk.image >> IMG2BITS;
          otherdims = bk.image ^ (idim << IMG2BITS);
          idim--;
          idim &= IMGMASK;
          bk.image = otherdims | (idim << IMG2BITS);
        }
        while (bk.xcm[2] >= l_hi2) {
          bk.xcm[2] -= l_prd2;
          idim = bk.image >> IMG2BITS;
          otherdims = bk.image ^ (idim << IMG2BITS);
          idim++;
          idim &= IMGMASK;
          bk.image = otherdims | (idim << IMG2BITS);
        }
        bk.xcm[2] = MAX(bk.xcm[2], l_lo2);
      }
    }
  );
  copymode = 0;
  k_body.modify_device();

  nghost_body = 0;
  commflag = FULL_BODY;
  comm->forward_comm(this);

  reset_atom2body();
  image_shift();
}


/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::setup(int vflag)
{

  const int nlocal = atomKK->nlocal;

  // error if maxextent > comm->cutghost
  // NOTE: could just warn if an override flag set
  // NOTE: this could fail for comm multi mode if user sets a wrong cutoff
  //       for atom types in rigid bodies - need a more careful test
  // must check here, not in init, b/c neigh/comm values set after fix init

  double cutghost = MAX(neighbor->cutneighmax, comm->cutghostuser);
  if (maxextent > cutghost)
    error->all(FLERR, Error::NOLASTLINE,
               "Rigid body extent {} > ghost atom cutoff - use comm_modify cutoff", maxextent);

  if (langflag && (nlocal_body > maxlang)) {
    memoryKK->destroy_kokkos(k_langextra, langextra);
    maxlang = nbody_total();
    memoryKK->create_kokkos(k_langextra, langextra, 6, "rigid/small:langextra");
  }

  compute_forces_and_torques();
  // enforce 2d body forces and torques
  if (domainKK->dimension == 2) enforce2d();

  // virial setup before call to set_v
  v_init(vflag);

  // compute and forward communicate vcm and omega of all bodies
  k_body.template sync<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();

  copymode = 1;
  Kokkos::parallel_for("rigid/small:setup",
    Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      angmom_to_omega(bk.angmom, bk.ex_space, bk.ey_space,
                                       bk.ez_space, bk.inertia, bk.omega);
    }
  );
  copymode = 0;
  k_body.template modify<DeviceType>();

  commflag = FINAL;
  comm->forward_comm(this,10);

  // set velocity/rotation of atoms in rigid bodues
  if (evflag) set_v<true>();
  else set_v<false>();

  // guesstimate virial as 2x the set_v contribution
  if (vflag_global) {
    for (int n = 0; n < 6; n++) virial[n] *= 2.0;
  }
  if (vflag_atom) {
    for (int i = 0; i < nlocal; i++) {
      for (int n = 0; n < 6; n++) vatom[i][n] *= 2.0;
    }
  }

  if constexpr(is_nh) {

  compute_dof();

  copymode = 1;
  k_body.sync_device();
  auto l_tstat_flag = tstat_flag;
  auto l_pstat_flag = pstat_flag;
  KK_ACC_FLOAT ke[2], keall[2];
  Kokkos::parallel_reduce("rigid/small:setup_nh",
    Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
    KOKKOS_LAMBDA(const int &ibody, KK_ACC_FLOAT &l_akin_t, KK_ACC_FLOAT &l_akin_r ) {
      BodyKokkos &bk = l_body(ibody);
      KK_FLOAT mbody[3];
      transpose_matvec(bk.ex_space, bk.ey_space, bk.ez_space,
                                        bk.angmom, mbody);
      quatvec(bk.quat, mbody, bk.conjqm);
      bk.conjqm[0] *= 2.0;
      bk.conjqm[1] *= 2.0;
      bk.conjqm[2] *= 2.0;
      bk.conjqm[3] *= 2.0;
      if (l_tstat_flag || l_pstat_flag) {
        l_akin_t += bk.mass * (bk.vcm[0] * bk.vcm[0] + bk.vcm[1] * bk.vcm[1]
                               + bk.vcm[2] * bk.vcm[2]);
        l_akin_r += bk.angmom[0] * bk.omega[0] + bk.angmom[1] * bk.omega[1]
                    + bk.angmom[2] * bk.omega[2];
      }
    }, ke[0], ke[1]
  );
  copymode = 0;
  k_body.modify_device();
  if (l_tstat_flag || l_pstat_flag) {
    MPI_Allreduce(ke, keall, 2, MPI_KK_ACC_FLOAT, MPI_SUM, world);
    FixRigidBase::akin_t = keall[0];
    FixRigidBase::akin_r = keall[1];
  }
  if (l_tstat_flag) FixRigidBase::compute_temp_target();
  else if (l_pstat_flag) {
    FixRigidBase::t0 = FixRigidBase::temperature->compute_scalar();
    if (FixRigidBase::t0 == 0.0) {
      if (strcmp(update->unit_style, "lj") == 0) FixRigidBase::t0 = 1.0;
      else FixRigidBase::t0 = 300.0;
    }
    FixRigidBase::t_target = FixRigidBase::t0;
  }
  if (l_pstat_flag) {
    FixRigidBase::compute_press_target();
    if (pstyle == ISO) {
      FixRigidBase::temperature->compute_scalar();
      FixRigidBase::pressure->compute_scalar();
    } else {
      FixRigidBase::temperature->compute_vector();
      FixRigidBase::pressure->compute_vector();
    }
    FixRigidBase::couple();
    FixRigidBase::pressure->addstep(update->ntimestep+1);
  }
  double t_mass, tb_mass;
  const double kt = FixRigidBase::boltz * FixRigidBase::t_target;
  if (l_tstat_flag) {
    t_mass = kt / (FixRigidBase::t_freq * FixRigidBase::t_freq);
    FixRigidBase::q_t[0] = FixRigidBase::nf_t * t_mass;
    FixRigidBase::q_r[0] = FixRigidBase::nf_r * t_mass;
    for (int i = 1; i < t_chain; i++)
      FixRigidBase::q_t[i] = FixRigidBase::q_r[i] = t_mass;
    for (int i = 1; i < t_chain; i++) {
      FixRigidBase::f_eta_t[i] = (FixRigidBase::q_t[i-1] * FixRigidBase::eta_dot_t[i-1] * FixRigidBase::eta_dot_t[i-1] - kt)/FixRigidBase::q_t[i];
      FixRigidBase::f_eta_r[i] = (FixRigidBase::q_r[i-1] * FixRigidBase::eta_dot_r[i-1] * FixRigidBase::eta_dot_r[i-1] - kt)/FixRigidBase::q_r[i];
    }
  }
  const int dimension = domainKK->dimension;
  if (l_pstat_flag) {
    for (int i = 0; i < 3; i++)
      if (p_flag[i]) {
        const auto p_freq_i_sq = (p_freq[i]) * (p_freq[i]);
        FixRigidBase::epsilon_mass[i] = (FixRigidBase::g_f + dimension) * kt / p_freq_i_sq;
        FixRigidBase::epsilon[i] = log(FixRigidBase::vol0)/dimension;
      }
    tb_mass = kt / (FixRigidBase::p_freq_max * FixRigidBase::p_freq_max);
    FixRigidBase::q_b[0] = dimension * dimension * tb_mass;
    for (int i = 1; i < p_chain; i++) {
      FixRigidBase::q_b[i] = tb_mass;
      const auto eta_dot_sq = (FixRigidBase::eta_dot_b[i-1]) * (FixRigidBase::eta_dot_b[i-1]);
      FixRigidBase::f_eta_b[i] = (FixRigidBase::q_b[i] * eta_dot_sq - kt)/FixRigidBase::q_b[i];
    }
  }
  if (l_tstat_flag || l_pstat_flag) {
    for (int i = 0; i < t_order; i++) {
      FixRigidBase::wdti1[i] = FixRigidBase::w[i] * dtv / t_iter;
      FixRigidBase::wdti2[i] = FixRigidBase::wdti1[i] / 2.0;
      FixRigidBase::wdti4[i] = FixRigidBase::wdti1[i] / 4.0;
    }
  }
  if (l_pstat_flag) {
    FixRigidBase::compute_press_target();
    FixRigidBase::nh_epsilon_dot();
  }

  }
}

/* ---------------------------------------------------------------------- */


template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::initial_integrate(int vflag)
{
  auto lambda = [&]<bool NH, bool TSTAT, bool PSTAT>() {

    // kokkos view
    k_body.template sync<DeviceType>();
    auto l_body = k_body.template view<DeviceType>();

    KK_FLOAT l_dtf = static_cast<KK_FLOAT>(dtf);
    KK_FLOAT l_dtf2 = l_dtf * KK_TWO;
    KK_FLOAT l_dtv = static_cast<KK_FLOAT>(dtv);
    KK_FLOAT l_dtq = static_cast<KK_FLOAT>(dtq);

    KK_FLOAT l_scale_t0, l_scale_t1, l_scale_t2;
    KK_FLOAT l_scale_v0, l_scale_v1, l_scale_v2, l_scale_r;

    if constexpr (NH) {
      KK_FLOAT tmp;
      l_scale_t0 = l_scale_t1 = l_scale_t2 = KK_ONE;
      l_scale_v0 = l_scale_v1 = l_scale_v2 = KK_ONE;
      l_scale_r = KK_ONE;
      if constexpr(TSTAT) {
        tmp = exp(-l_dtq * static_cast<KK_FLOAT>(nh()->eta_dot_t[0]));
        l_scale_t0 = l_scale_t1 = l_scale_t2 = tmp;
        tmp = exp(-l_dtq * static_cast<KK_FLOAT>(nh()->eta_dot_r[0]));
        l_scale_r = tmp;
      }
      if constexpr (PSTAT) {

        auto maclaurin = [&](KK_FLOAT x) -> KK_FLOAT {
          const KK_FLOAT y = x * x;
          static constexpr KK_FLOAT c1 = KK_FLOAT(1.0 / 6.0);
          static constexpr KK_FLOAT c2 = KK_FLOAT(1.0 / 120.0);
          static constexpr KK_FLOAT c3 = KK_FLOAT(1.0 / 5040.0);
          static constexpr KK_FLOAT c4 = KK_FLOAT(1.0 / 362880.0);
          // horner form with fused multiply add
          // P(y) = 1 + y (1/6 + y(1/120 + y(1/5040 + y(1/362880))))
          return fma(y, fma(y, fma(y, fma(y, c4, c3), c2), c1), KK_ONE);
        };

        auto epsilon_dot = nh()->epsilon_dot;
        KK_FLOAT l_mtk_term2 = static_cast<KK_FLOAT>(nh()->mtk_term2);
        l_scale_t0 *= exp(-l_dtq * (static_cast<KK_FLOAT>(epsilon_dot[0]) + l_mtk_term2));
        l_scale_t1 *= exp(-l_dtq * (static_cast<KK_FLOAT>(epsilon_dot[1]) + l_mtk_term2));
        l_scale_t2 *= exp(-l_dtq * (static_cast<KK_FLOAT>(epsilon_dot[2]) + l_mtk_term2));
        l_scale_r *= exp(-l_dtq * (static_cast<KK_FLOAT>(nh()->pdim) * l_mtk_term2));
        tmp = l_dtq * static_cast<KK_FLOAT>(epsilon_dot[0]);
        l_scale_v0 = l_dtv * exp(tmp) * maclaurin(tmp);
        tmp = l_dtq * static_cast<KK_FLOAT>(epsilon_dot[1]);
        l_scale_v1 = l_dtv * exp(tmp) * maclaurin(tmp);
        tmp = l_dtq * static_cast<KK_FLOAT>(epsilon_dot[2]);
        l_scale_v2 = l_dtv * exp(tmp) * maclaurin(tmp);
      }
    }

    Kokkos::parallel_for("rigid/small:initial_integrate",
      Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
      KOKKOS_LAMBDA(const int &ibody){
        BodyKokkos &bk = l_body(ibody);
        // update vcm by 1/2 step
        const KK_FLOAT dtfm = l_dtf / bk.mass;
        bk.vcm[0] = fma(dtfm, bk.fcm[0], bk.vcm[0]);
        bk.vcm[1] = fma(dtfm, bk.fcm[1], bk.vcm[1]);
        bk.vcm[2] = fma(dtfm, bk.fcm[2], bk.vcm[2]);
        if constexpr(TSTAT || PSTAT) {
          bk.vcm[0] *= l_scale_t0;
          bk.vcm[1] *= l_scale_t1;
          bk.vcm[2] *= l_scale_t2;
        }
        // update xcm by full step
        if constexpr(!PSTAT) {
          bk.xcm[0] = fma(l_dtv, bk.vcm[0], bk.xcm[0]);
          bk.xcm[1] = fma(l_dtv, bk.vcm[1], bk.xcm[1]);
          bk.xcm[2] = fma(l_dtv, bk.vcm[2], bk.xcm[2]);
        } else {
          bk.xcm[0] = fma(l_scale_v0, bk.vcm[0], bk.xcm[0]);
          bk.xcm[1] = fma(l_scale_v1, bk.vcm[1], bk.xcm[1]);
          bk.xcm[2] = fma(l_scale_v2, bk.vcm[2], bk.xcm[2]);
        }
        if constexpr(!NH) {
          // update angular momentum by 1/2 step
          bk.angmom[0] = fma(l_dtf, bk.torque[0], bk.angmom[0]);
          bk.angmom[1] = fma(l_dtf, bk.torque[1], bk.angmom[1]);
          bk.angmom[2] = fma(l_dtf, bk.torque[2], bk.angmom[2]);
          // compute omega at 1/2 step from angmom at 1/2 step and current q
          // update quaternion a full step via Richardson iteration
          // returns new normalized quaternion, also updated omega at 1/2 step
          // update ex,ey,ez to reflect new quaternion
          angmom_to_omega(bk.angmom, bk.ex_space, bk.ey_space,
                          bk.ez_space, bk.inertia, bk.omega);
          richardson(bk.quat, bk.angmom, bk.omega, bk.inertia, l_dtq);
          q_to_exyz(bk.quat, bk.ex_space, bk.ey_space, bk.ez_space);
        } else {
          // apply torque to quaternion momentum
          KK_FLOAT mbody[3], tbody[3], fquat[4];
          transpose_matvec(bk.ex_space, bk.ey_space, bk.ez_space,
                                      bk.torque, tbody);
          quatvec(bk.quat, tbody, fquat);
          bk.conjqm[0] = fma(l_dtf2, fquat[0], bk.conjqm[0]);
          bk.conjqm[1] = fma(l_dtf2, fquat[1], bk.conjqm[1]);
          bk.conjqm[2] = fma(l_dtf2, fquat[2], bk.conjqm[2]);
          bk.conjqm[3] = fma(l_dtf2, fquat[3], bk.conjqm[3]);
          if constexpr (TSTAT || PSTAT) {
            bk.conjqm[0] *= l_scale_r;
            bk.conjqm[1] *= l_scale_r;
            bk.conjqm[2] *= l_scale_r;
            bk.conjqm[3] *= l_scale_r;
          }
          // step 1.4 to 1.13 - no_squish_rotate
          no_squish_rotate(3, bk.conjqm, bk.quat, bk.inertia, l_dtq);
          no_squish_rotate(2, bk.conjqm, bk.quat, bk.inertia, l_dtq);
          no_squish_rotate(1, bk.conjqm, bk.quat, bk.inertia, l_dtv);
          no_squish_rotate(2, bk.conjqm, bk.quat, bk.inertia, l_dtq);
          no_squish_rotate(3, bk.conjqm, bk.quat, bk.inertia, l_dtq);
          // update exyz_space, transform p back to angmom, update omega
          q_to_exyz(bk.quat, bk.ex_space, bk.ey_space, bk.ez_space);
          invquatvec(bk.quat, bk.conjqm, mbody);
          matvec(bk.ex_space, bk.ey_space, bk.ez_space, mbody, bk.angmom);
          bk.angmom[0] *= KK_HALF;
          bk.angmom[1] *= KK_HALF;
          bk.angmom[2] *= KK_HALF;
          angmom_to_omega(bk.angmom, bk.ex_space, bk.ey_space,
                          bk.ez_space, bk.inertia, bk.omega);
        }

      }
    );
    k_body.modify_device();
  };

  copymode = 1;
  if constexpr (is_nh) {
    const bool is_t = tstat_flag, is_p = pstat_flag;
    if      ( is_t &&  is_p) lambda.template operator()<true,true,true>();   // NPT
    else if ( is_t && !is_p) lambda.template operator()<true,true,false>();  // NVT
    else if (!is_t &&  is_p) lambda.template operator()<true,false,true>();  // NPH
    else                     lambda.template operator()<true,false,false>(); // NVE
  } else                     lambda.template operator()<false,false,false>();
  copymode = 0;

  commflag = INITIAL;
  comm->forward_comm(this, 29);

  if constexpr(is_nh) {
  // accumulate kinetic energies
  if (tstat_flag || pstat_flag) {

    // kokkos view
    k_body.template sync<DeviceType>();
    auto l_body = k_body.template view<DeviceType>();
    KK_ACC_FLOAT ke[2], keall[2];

    copymode = 1;
    Kokkos::parallel_reduce("rigid/small:initial_integrate_nh",
      Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
      KOKKOS_LAMBDA(const int &ibody, KK_ACC_FLOAT &l_akin_t, KK_ACC_FLOAT &l_akin_r) {
        BodyKokkos &bk = l_body(ibody);
        l_akin_t += bk.mass*(bk.vcm[0]*bk.vcm[0] + bk.vcm[1]*bk.vcm[1] +
                    bk.vcm[2]*bk.vcm[2]);
        l_akin_r += bk.angmom[0]*bk.omega[0] + bk.angmom[1]*bk.omega[1] +
                    bk.angmom[2]*bk.omega[2];
      }, ke[0], ke[1]);
    copymode = 0;
    MPI_Allreduce(ke, keall, 2, MPI_KK_ACC_FLOAT, MPI_SUM, world);
    nh()->akin_t = keall[0];
    nh()->akin_r = keall[1];
  }
  if (tstat_flag) {
    FixRigidBase::compute_temp_target();
    if (dynamic) compute_dof();
    FixRigidBase::nhc_temp_integrate();
  }
  if (pstat_flag) {
    FixRigidBase::nhc_press_integrate();
    remap();
  }
  }

  // virial setup
  v_init(vflag);
  if (vflag_atom) {
    if (atomKK->nmax > (int) d_vatom.extent(0)) {
      memoryKK->destroy_kokkos(k_vatom, vatom);
      memoryKK->create_kokkos(k_vatom, vatom, atomKK->nmax, 6, "fix:vatom");
      d_vatom = k_vatom.template view<DeviceType>();
    } else {
      k_vatom.template sync<DeviceType>();
    }
    Kokkos::deep_copy(d_vatom, 0.0);
  }

  if (evflag) set_xv<true>();
  else set_xv<false>();

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  if (extended) {
    // not implemented
  }

  if constexpr(is_nh) {
    // remap again + kspace
    if (pstat_flag) {
      remap();
      if (nh()->kspace_flag) {
        atomKK->sync(Host, X_MASK | V_MASK);
        force->kspace->setup();
      }
    }
  }

}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::final_integrate()
{

  if (!earlyflag) compute_forces_and_torques();
  if (domainKK->dimension == 2) enforce2d();

  // kokkos views

  k_body.template sync<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();

  auto lambda = [&]<bool NH, bool TSTAT, bool PSTAT>() {
    KK_FLOAT l_dtf = static_cast<KK_FLOAT>(dtf);
    KK_FLOAT l_dtf2 = l_dtf * KK_TWO;
    KK_FLOAT l_scale_t0, l_scale_t1, l_scale_t2, l_scale_r;
    if constexpr (NH) {
      KK_FLOAT l_dtq = static_cast<KK_FLOAT>(nh()->dtq);
      l_scale_t0 = l_scale_t1 = l_scale_t2 = l_scale_r = KK_ONE;
      if constexpr (TSTAT) {
        const KK_FLOAT tmp = exp(-KK_ONE * l_dtq * static_cast<KK_FLOAT>(nh()->eta_dot_t[0]));
        l_scale_t0 = l_scale_t1 = l_scale_t2 = tmp;
        l_scale_r = exp(-KK_ONE * l_dtq * static_cast<KK_FLOAT>(nh()->eta_dot_r[0]));
      }
      if constexpr (PSTAT) {
        auto epsilon_dot = nh()->epsilon_dot;
        KK_FLOAT l_mtk_term2 = static_cast<KK_FLOAT>(nh()->mtk_term2);
        l_scale_t0 *= exp(-l_dtq * (static_cast<KK_FLOAT>(epsilon_dot[0]) + l_mtk_term2));
        l_scale_t1 *= exp(-l_dtq * (static_cast<KK_FLOAT>(epsilon_dot[1]) + l_mtk_term2));
        l_scale_t2 *= exp(-l_dtq * (static_cast<KK_FLOAT>(epsilon_dot[2]) + l_mtk_term2));
        l_scale_r *= exp(-l_dtq * (static_cast<KK_FLOAT>(nh()->pdim) * l_mtk_term2));
      }
    }

    Kokkos::parallel_for("rigid/small:final_integrate",
      Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
      KOKKOS_LAMBDA(const int &ibody) {
        BodyKokkos &bk = l_body(ibody);
        // update vcm by 1/2 step
        const KK_FLOAT dtfm = l_dtf / bk.mass;
        if constexpr (TSTAT || PSTAT) {
          bk.vcm[0] *= l_scale_t0;
          bk.vcm[1] *= l_scale_t1;
          bk.vcm[2] *= l_scale_t2;
        }
        bk.vcm[0] = fma(dtfm, bk.fcm[0], bk.vcm[0]);
        bk.vcm[1] = fma(dtfm, bk.fcm[1], bk.vcm[1]);
        bk.vcm[2] = fma(dtfm, bk.fcm[2], bk.vcm[2]);

        if constexpr (!NH) {
          // update angular momentum by 1/2 step
          bk.angmom[0] = fma(l_dtf, bk.torque[0], bk.angmom[0]);
          bk.angmom[1] = fma(l_dtf, bk.torque[1], bk.angmom[1]);
          bk.angmom[2] = fma(l_dtf, bk.torque[2], bk.angmom[2]);
          angmom_to_omega(bk.angmom, bk.ex_space, bk.ey_space,
                                   bk.ez_space, bk.inertia, bk.omega);
        } else {
          KK_FLOAT mbody[3], tbody[3], fquat[4];
          transpose_matvec(bk.ex_space, bk.ey_space,
                                    bk.ez_space, bk.torque, tbody);
          quatvec(bk.quat, tbody, fquat);
          if constexpr (TSTAT || PSTAT) {
            bk.conjqm[0] = fma(l_scale_r, bk.conjqm[0], l_dtf2 * fquat[0]);
            bk.conjqm[1] = fma(l_scale_r, bk.conjqm[1], l_dtf2 * fquat[1]);
            bk.conjqm[2] = fma(l_scale_r, bk.conjqm[2], l_dtf2 * fquat[2]);
            bk.conjqm[3] = fma(l_scale_r, bk.conjqm[3], l_dtf2 * fquat[3]);
          } else {
            bk.conjqm[0] = Kokkos::fma(l_dtf2, fquat[0], bk.conjqm[0]);
            bk.conjqm[1] = Kokkos::fma(l_dtf2, fquat[1], bk.conjqm[1]);
            bk.conjqm[2] = Kokkos::fma(l_dtf2, fquat[2], bk.conjqm[2]);
            bk.conjqm[3] = Kokkos::fma(l_dtf2, fquat[3], bk.conjqm[3]);
          }
          invquatvec(bk.quat, bk.conjqm, mbody);
          matvec(bk.ex_space, bk.ey_space, bk.ez_space, mbody, bk.angmom);
          bk.angmom[0] *= KK_HALF;
          bk.angmom[1] *= KK_HALF;
          bk.angmom[2] *= KK_HALF;
          angmom_to_omega(bk.angmom, bk.ex_space, bk.ey_space,
                                   bk.ez_space, bk.inertia, bk.omega);
        }
      }
    );
    k_body.modify_device();
  };

  copymode = 1;
  if constexpr (is_nh) {
    const bool is_t = tstat_flag, is_p = pstat_flag;
    if      ( is_t &&  is_p) lambda.template operator()<true,true,true>();   // NPT
    else if ( is_t && !is_p) lambda.template operator()<true,true,false>();  // NVT
    else if (!is_t &&  is_p) lambda.template operator()<true,false,true>();  // NPH
    else                     lambda.template operator()<true,false,false>(); // NVE
  } else                     lambda.template operator()<false,false,false>();
  copymode = 0;

  commflag = FINAL;
  comm->forward_comm(this, 10);

  if (vflag_atom) k_vatom.template sync<DeviceType>();

  if (evflag) set_v<true>();
  else set_v<false>();

  if (extended) {
    // not implemented
  }
  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  if constexpr (is_nh) {
    auto *temperature = nh()->temperature;
    // compute temperature
    if (FixRigidBase::tcomputeflag) {
      atomKK->sync(temperature->execution_space, temperature->datamask_read);
      nh()->t_current = temperature->compute_scalar();
      atomKK->modified(temperature->execution_space, temperature->datamask_modify);
      atomKK->sync(execution_space, temperature->datamask_modify);
    }
    // pressure
    if (pstat_flag) {
      auto *pressure = nh()->pressure;
      // accumulate kinetic energies for pstat
      copymode = 1;
      k_body.sync_device();
      KK_ACC_FLOAT ke[2], keall[2];
      Kokkos::parallel_reduce("rigid/small:final_integrate_nh",
        Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
        KOKKOS_LAMBDA(const int ibody, KK_ACC_FLOAT &l_akin_t, KK_ACC_FLOAT &l_akin_r ) {
          BodyKokkos &bk = l_body(ibody);
          l_akin_t += bk.mass*(bk.vcm[0]*bk.vcm[0] + bk.vcm[1]*bk.vcm[1] +
                    bk.vcm[2]*bk.vcm[2]);
          l_akin_r += bk.angmom[0]*bk.omega[0] + bk.angmom[1]*bk.omega[1] +
                    bk.angmom[2]*bk.omega[2];
        }, ke[0], ke[1]
      );
    copymode = 0;
    MPI_Allreduce(ke, keall, 2, MPI_KK_ACC_FLOAT, MPI_SUM, world);
    nh()->akin_t = keall[0];
    nh()->akin_r = keall[1];

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
    FixRigidBase::couple();
    pressure->addstep(update->ntimestep+1);
    FixRigidBase::compute_press_target();
    FixRigidBase::nh_epsilon_dot();
  }
  }
}

/* ----------------------------------------------------------------------
   count # of DOF removed by rigid bodies for atoms in igroup
   return total count of DOF
------------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
bigint FixRigidBaseKokkos<DeviceType,FixRigidBase>::dof(int tgroup)
{

  // cannot count DOF correctly unless setup_bodies_static() has been called
  if (!setupflag) {
    if (comm->me == 0)
      error->warning(FLERR,"Cannot count rigid body degrees-of-freedom "
                     "before bodies are fully initialized");
    return 0;
  }

  // kokkos views

  atomKK->sync(execution_space, MASK_MASK);
  auto l_mask = atomKK->k_mask.template view<DeviceType>();

  k_body.template sync<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();

  k_atom2body.template sync<DeviceType>();
  auto l_atom2body = k_atom2body.template view<DeviceType>();

  k_eflags.template sync<DeviceType>();
  auto l_eflags = k_eflags.template view<DeviceType>();

  // counts = 3 values per rigid body I own
  // 0 = # of point particles in rigid body and in temperature group
  // 1 = # of finite-size particles in rigid body and in temperature group
  // 2 = # of particles in rigid body, disregarding temperature group

  memoryKK->create_kokkos(k_counts, counts, nbody_total(), 3, "rigid/small:counts");
  auto l_counts = k_counts.template view<DeviceType>();
  deep_copy(k_counts.view_device(), 0);

  // tally counts from my owned atoms
  // 0 = # of point particles in rigid body and in temperature group
  // 1 = # of finite-size particles in rigid body and in temperature group
  // 2 = # of particles in rigid body, disregarding temperature group

  int l_tgroupbit = group->bitmask[tgroup];
  auto l_extended = extended;

  Kokkos::parallel_for("rigid/small:dof_count",
    Kokkos::RangePolicy<DeviceType>(0, nlocal()),
    KOKKOS_LAMBDA(const int i) {
      const int ibody = l_atom2body(i);
      if (ibody < 0) return;
      Kokkos::atomic_inc(&l_counts(ibody,2));
      if (l_mask(i) & l_tgroupbit) {
        if (l_extended && (l_eflags(i) & ~(POINT | DIPOLE)))
          Kokkos::atomic_inc(&l_counts(ibody,1));
        else
          Kokkos::atomic_inc(&l_counts(ibody,0));
    }
  });
  k_counts.modify_device();

  commflag = DOF;
  comm->reverse_comm(this, 3);

  // nall = count0 = # of point particles in each rigid body
  // mall = count1 = # of finite-size particles in each rigid body
  // warn if nall+mall != nrigid for any body included in temperature group

  int flag = 0;
  Kokkos::parallel_reduce("rigid/small:dof_check",
    Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
    KOKKOS_LAMBDA(const int ibody, int &l_flag) {
      const int l_counts01 = l_counts(ibody,0) + l_counts(ibody,1);
      if (l_counts01 > 0 && l_counts01 != l_counts(ibody,2)) l_flag = 1;
      else l_flag = 0;
    }, Kokkos::Max<int>(flag)
  );
  int flag_all;
  MPI_Allreduce(&flag, &flag_all, 1, MPI_INT, MPI_MAX, world);
  if (flag_all && comm->me == 0)
    error->warning(FLERR,"Computing temperature of portions of rigid bodies");

  // remove appropriate DOFs for each rigid body wholly in temperature group
  // N = # of point particles in body
  // M = # of finite-size particles in body
  // 3d body has 3N + 6M dof to start with
  // 2d body has 2N + 3M dof to start with
  // 3d point-particle body with all non-zero I should have 6 dof, remove 3N-6
  // 3d point-particle body (linear) with a 0 I should have 5 dof, remove 3N-5
  // 2d point-particle body should have 3 dof, remove 2N-3
  // 3d body with any finite-size M should have 6 dof, remove (3N+6M) - 6
  // 2d body with any finite-size M should have 3 dof, remove (2N+3M) - 3

  bigint n = 0;
  nlinear = 0;
  if (domainKK->dimension == 3) {
    Kokkos::parallel_reduce("rigid/small:dof_remove_3d",
      Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
      KOKKOS_LAMBDA(const int ibody, bigint &l_n, int &l_nlinear) {
        if (l_counts(ibody,0) + l_counts(ibody,1) == l_counts(ibody,2)) {
          l_n += 3*l_counts(ibody,0) + 6*l_counts(ibody,1) - 6;
          auto inertia = l_body(ibody).inertia;
          if (inertia[0] == KK_ZERO || inertia[1] == KK_ZERO || inertia[2] == KK_ZERO) {
            l_n++;
            l_nlinear++;
          }
        }
      }, n, nlinear
    );
  } else if (domainKK->dimension == 2) {
    Kokkos::parallel_reduce("rigid/small:dof_remove_2d",
      Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
      KOKKOS_LAMBDA(const int ibody, bigint &l_n) {
        if (l_counts(ibody,0) + l_counts(ibody,1) == l_counts(ibody,2)) {
          l_n += 2*l_counts(ibody,0) + 3*l_counts(ibody,1) - 3;
        }
      }, n
    );
  }

  bigint n_all;
  MPI_Allreduce(&n, &n_all, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  return n_all;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
double FixRigidBaseKokkos<DeviceType,FixRigidBase>::compute_scalar()
{

  // kokkos view
  k_body.template sync<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();

  auto lambda = [&]<bool NH>(const int &ibody, KK_ACC_FLOAT &l_ke) {
    const BodyKokkos &bk = l_body(ibody);
    if constexpr (!NH) {
      l_ke += bk.mass * fma(bk.vcm[0] , bk.vcm[0],
                        fma(bk.vcm[1] , bk.vcm[1],
                            bk.vcm[2] * bk.vcm[2]));
      // for Iw^2 rotational term, need wbody = angular velocity in body frame
      // not omega = angular velocity in space frame
      KK_FLOAT wbody[3], rot[3][3];
      quat_to_mat(bk.quat, rot);
      transpose_matvec(rot, bk.angmom, wbody);
      if (bk.inertia[0] == KK_ZERO) wbody[0] = KK_ZERO;
      else wbody[0] /= bk.inertia[0];
      if (bk.inertia[1] == KK_ZERO) wbody[1] = KK_ZERO;
      else wbody[1] /= bk.inertia[1];
      if (bk.inertia[2] == KK_ZERO) wbody[2] = KK_ZERO;
      else wbody[2] /= bk.inertia[2];
      l_ke += fma(bk.inertia[0] , wbody[0] * wbody[0],
              fma(bk.inertia[1] , wbody[1] * wbody[1],
                  bk.inertia[2] * wbody[2] * wbody[2]));
    } else {
      l_ke += bk.mass * fma(bk.vcm[0] , bk.vcm[0],
                        fma(bk.vcm[1] , bk.vcm[1],
                            bk.vcm[2] * bk.vcm[2]));
      KK_FLOAT Pkq[4];
      for (int k = 1; k < 4; k++) {
        if (k == 1) {
          Pkq[0] = -bk.quat[1];
          Pkq[1] =  bk.quat[0];
          Pkq[2] =  bk.quat[3];
          Pkq[3] = -bk.quat[2];
        } else if (k == 2) {
          Pkq[0] = -bk.quat[2];
          Pkq[1] = -bk.quat[3];
          Pkq[2] =  bk.quat[0];
          Pkq[3] =  bk.quat[1];
        } else if (k == 3) {
          Pkq[0] = -bk.quat[3];
          Pkq[1] =  bk.quat[2];
          Pkq[2] = -bk.quat[1];
          Pkq[3] =  bk.quat[0];
        }
        KK_ACC_FLOAT tmp = static_cast<KK_ACC_FLOAT>(
          fma(bk.conjqm[0], Pkq[0],
          fma(bk.conjqm[1], Pkq[1],
          fma(bk.conjqm[2], Pkq[2],
              bk.conjqm[3] * Pkq[3])))
        );
        tmp *= tmp;
        if (Kokkos::fabs(bk.inertia[k-1]) < 1e-6) tmp = KK_ACC_FLOAT(0.0);
        else tmp /= (KK_ACC_FLOAT(8.0) * bk.inertia[k-1]);
        l_ke += tmp;
      }
    }
  };

  KK_ACC_FLOAT ke, ke_all;
  copymode = 1;
  Kokkos::parallel_reduce("rigid/small:compute_scalar",
    Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
    KOKKOS_LAMBDA(const int &ibody, KK_ACC_FLOAT &l_ke) {
      lambda.template operator()<is_nh>(ibody, l_ke);
    }, ke
  );
  copymode = 0;
  MPI_Allreduce(&ke, &ke_all, 1, MPI_KK_ACC_FLOAT, MPI_SUM, world);

  if constexpr (!is_nh) {
    KK_ACC_FLOAT tfactor = force->mvv2e / ((6.0*nbody - nlinear) * force->boltz);
    return static_cast<double>(ke_all * tfactor);
  } else {
    // compute the kinetic parts of H_NVE in Kameraj et al (JCP 2005, pp 224114)
    // translational and rotational kinetic energy
    KK_ACC_FLOAT energy = ke_all * FixRigidBase::mvv2e;
    KK_ACC_FLOAT kt = force->boltz * nh()->t_target;
    if (tstat_flag) {
      // thermostat chain energy: from equation 12 in Kameraj et al (JCP 2005)
      energy += kt * (nh()->nf_t * nh()->eta_t[0] + nh()->nf_r * nh()->eta_r[0]);
      for (int i = 1; i < t_chain; i++)
        energy += kt * (nh()->eta_t[i] + nh()->eta_r[i]);
      for (int i = 0; i < t_chain; i++) {
        energy += 0.5 * nh()->q_t[i] * (nh()->eta_dot_t[i] * nh()->eta_dot_t[i]);
        energy += 0.5 * nh()->q_r[i] * (nh()->eta_dot_r[i] * nh()->eta_dot_r[i]);
      }
    }
    if (pstat_flag) {
      // using equation 22 in Kameraj et al for H_NPT
      KK_ACC_FLOAT e = KK_ACC_FLOAT(0.0);
      for (int i = 0; i < 3; i++) {
        if (p_flag[i])
          e += nh()->epsilon_mass[i] * nh()->epsilon_dot[i] * nh()->epsilon_dot[i];
      }
      energy += e*(0.5/nh()->pdim);
      double vol = domainKK->xprd * domainKK->yprd;
      if (domainKK->dimension == 3) vol *= domainKK->zprd;
      double p0 = (nh()->p_target[0] + nh()->p_target[1] + nh()->p_target[2]) / 3.0;
      energy += p0 * vol / nh()->nktv2p;
      for (int i = 0;  i < p_chain; i++) {
        energy += kt * nh()->eta_b[i];
        energy += 0.5 * nh()->q_b[i] * (nh()->eta_dot_b[i] * nh()->eta_dot_b[i]);
      }
    }
    return energy;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::zero_momentum()
{

  // kokkos views
  k_body.template sync<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();

  copymode = 1;
  Kokkos::parallel_for("rigid/small:zero_momentum",
    Kokkos::RangePolicy<DeviceType>(0, nbody_total()),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      bk.vcm[0] = bk.vcm[1] = bk.vcm[2] = KK_ZERO;
    }
  );
  copymode = 0;
  k_body.modify_device();

  // forward communicate of omega to all ghost copies
  commflag = FINAL;
  comm->forward_comm(this, 10);
  // set velocity of atoms in rigid bodues
  set_v<false>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::zero_rotation()
{

  // kokkos views
  k_body.template sync<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();

  copymode = 1;
  Kokkos::parallel_for("rigid/small:zero_rotation",
    Kokkos::RangePolicy<DeviceType>(0, nbody_total()),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      bk.angmom[0] = bk.angmom[1] = bk.angmom[2] = KK_ZERO;
      bk.omega[0] = bk.omega[1] = bk.omega[2] = KK_ZERO;
    }
  );
  copymode = 0;
  k_body.modify_device();

  // forward communicate of omega to all ghost copies
  commflag = FINAL;
  comm->forward_comm(this, 10);

  // set velocity of atoms in rigid bodues
  set_v<false>();

}

/* ----------------------------------------------------------------------
   FixRigidSmall PROTECTED METHODS
------------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::setup_bodies_static()
{

  // FIXME: -------- setup_bodies_static() --------

  // extended = 1 if any particle in a rigid body is finite size
  //              or has a dipole moment

  extended = orientflag = dorientflag = 0;

  AtomVecEllipsoid::Bonus *ebonus;
  if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
  AtomVecLine::Bonus *lbonus;
  if (avec_line) lbonus = avec_line->bonus;
  AtomVecTri::Bonus *tbonus;
  if (avec_tri) tbonus = avec_tri->bonus;
  double **mu = atomKK->mu;
  double *radius = atomKK->radius;
  double *rmass = atomKK->rmass;
  double *mass = atomKK->mass;
  int *ellipsoid = atomKK->ellipsoid;
  int *line = atomKK->line;
  int *tri = atomKK->tri;
  int *type = atomKK->type;

  if (atomKK->radius_flag || atomKK->ellipsoid_flag || atomKK->line_flag ||
      atomKK->tri_flag || atomKK->mu_flag) {
    int flag = 0;
    for (int i = 0; i < nlocal(); i++) {
      if (bodytag[i] == 0) continue;
      if (radius && radius[i] > 0.0) flag = 1;
      if (ellipsoid && ellipsoid[i] >= 0) flag = 1;
      if (line && line[i] >= 0) flag = 1;
      if (tri && tri[i] >= 0) flag = 1;
      if (mu && mu[i][3] > 0.0) flag = 1;
    }
    MPI_Allreduce(&flag, &extended, 1, MPI_INT, MPI_MAX, world);
  }

  // extended = 1 if using molecule template with finite-size particles
  // require all molecules in template to have consistent radiusflag

  if (onemols) {
    int radiusflag = onemols[0]->radiusflag;
    for (int i = 1; i < nmol; i++) {
      if (onemols[i]->radiusflag != radiusflag)
        error->all(FLERR, Error::NOLASTLINE, "Inconsistent use of finite-size particles "
                   "by molecule template molecules");
    }
    if (radiusflag) extended = 1;
  }

  if (extended) {
    error->all(FLERR, Error::NOLASTLINE,
               "Fix {}: extended particles not implemented in KOKKOS", style);
  }

  // grow extended arrays and set extended flags for each particle
  // orientflag = 4 if any particle stores ellipsoid or tri orientation or quat
  // orientflag = 1 if any particle stores line orientation
  // dorientflag = 1 if any particle stores dipole orientation

  if (extended) {
    // not implemented
  }

  // kokkos views

  atomKK->sync(execution_space, X_MASK | IMAGE_MASK | RMASS_MASK | TYPE_MASK );
  auto l_x = atomKK->k_x.template view<DeviceType>();
  auto l_image = atomKK->k_image.template view<DeviceType>();
  auto l_rmass = atomKK->k_rmass.template view<DeviceType>();
  auto l_type = atomKK->k_type.template view<DeviceType>();
  auto l_mass = atomKK->k_mass.template view<DeviceType>();

  k_atom2body.template sync<DeviceType>();
  auto l_atom2body = k_atom2body.template view<DeviceType>();

  k_bodytag.template sync<DeviceType>();
  auto l_bodytag = k_bodytag.template view<DeviceType>();

  k_xcmimage.template sync<DeviceType>();
  auto l_xcmimage = k_xcmimage.template view<DeviceType>();

  k_displace.template sync<DeviceType>();
  auto l_displace = k_displace.template view<DeviceType>();

  // domainKK

  auto l_triclinic = triclinic;
  auto l_periodic = Few<bool,3>(domainKK->periodicity);
  auto l_prd = Few<KK_FLOAT,3>(domainKK->prd);
  auto l_h = Few<KK_FLOAT,6>(domainKK->h);
  auto l_invprd = Few<KK_FLOAT,3>(
    KK_ONE / domainKK->prd[0],
    KK_ONE / domainKK->prd[1],
    KK_ONE / domainKK->prd[2]
  );

  // set body xcmimage flags = true image flags
  copymode = 1;
  Kokkos::parallel_for("rigid/small:setup_bodies_static_xcmimage",
    Kokkos::RangePolicy<DeviceType>(0, nlocal()),
    KOKKOS_LAMBDA(const int &i) {
      if (l_bodytag(i) >= 0) l_xcmimage(i) = l_image(i);
      else l_xcmimage(i) = 0;
    }
  );
  copymode = 0;
  k_xcmimage.template modify<DeviceType>();

  // acquire ghost bodies via forward comm
  // set atom2body for ghost atoms via forward comm
  // set atom2body for other owned atoms via reset_atom2body()
  nghost_body = 0;
  commflag = FULL_BODY;
  comm->forward_comm(this);
  reset_atom2body();

  k_body.template sync<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();

  // compute mass & center-of-mass of each rigid body
  k_body.sync_device();
  copymode = 1;
  Kokkos::parallel_for("rigid/small:setup_bodies_static_zero",
    Kokkos::RangePolicy<DeviceType>(0, nbody_total()),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      bk.xcm[0] = bk.xcm[1] = bk.xcm[2] = KK_ZERO;
      bk.xgc[0] = bk.xgc[1] = bk.xgc[2] = KK_ZERO;
      bk.mass = KK_ZERO;
      bk.natoms = 0;
    }
  );
  auto lambda_xcm_mass = [&]<bool RMASS>(const int &i) {
    const int ibody = l_atom2body(i);
    if (ibody < 0) return;
    BodyKokkos &bk = l_body(ibody);
    KK_FLOAT massone;
    if constexpr (RMASS) massone = l_rmass(i);
    else massone = l_mass(l_type(i));
    auto unwrap = DomainKokkos::unmap(l_prd, l_h, l_triclinic, &l_x(i,0), l_xcmimage(i));
    Kokkos::atomic_add(&bk.xcm[0], unwrap[0] * massone);
    Kokkos::atomic_add(&bk.xcm[1], unwrap[1] * massone);
    Kokkos::atomic_add(&bk.xcm[2], unwrap[2] * massone);
    Kokkos::atomic_add(&bk.xgc[0], unwrap[0]);
    Kokkos::atomic_add(&bk.xgc[1], unwrap[1]);
    Kokkos::atomic_add(&bk.xgc[2], unwrap[2]);
    Kokkos::atomic_add(&bk.mass, massone);
    Kokkos::atomic_add(&bk.natoms, 1);
  };
  if (l_rmass.data()) {
    Kokkos::parallel_for("rigid/small:setup_bodies_static_xcm_mass_rmass",
    Kokkos::RangePolicy<DeviceType>(0, nlocal()),
    KOKKOS_LAMBDA(const int i) {
      lambda_xcm_mass.template operator()<true>(i);
    });
  } else {
    Kokkos::parallel_for("rigid/small:setup_bodies_static_xcm_mass",
    Kokkos::RangePolicy<DeviceType>(0, nlocal()),
    KOKKOS_LAMBDA(const int i) {
      lambda_xcm_mass.template operator()<false>(i);
    });
  }
  copymode = 0;
  k_body.modify_device();

  // reverse communicate xcm, mass of all bodies
  commflag = XCM_MASS;
  comm->reverse_comm(this, 8);

  k_body.sync_device();
  copymode = 1;
  Kokkos::parallel_for("rigid/small:setup_bodies_static_xcm_mass_normalize",
  Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
  KOKKOS_LAMBDA(const int &ibody) {
    BodyKokkos &bk = l_body(ibody);
    bk.xcm[0] /= bk.mass;
    bk.xcm[1] /= bk.mass;
    bk.xcm[2] /= bk.mass;
    bk.xgc[0] /= bk.natoms;
    bk.xgc[1] /= bk.natoms;
    bk.xgc[2] /= bk.natoms;
    // set vcm, angmom = 0.0 in case inpfile is used
    // and doesn't overwrite all body's values
    // since setup_bodies_dynamic() will not be called
    // set rigid body image flags to default values
    bk.vcm[0] = bk.vcm[1] = bk.vcm[2] = KK_ZERO;
    bk.angmom[0] = bk.angmom[1] = bk.angmom[2] = KK_ZERO;
    // set rigid body image flags to default values
    bk.image = ((imageint) IMGMAX << IMG2BITS) |
               ((imageint) IMGMAX << IMGBITS) | IMGMAX;

  });
  copymode = 0;
  k_body.modify_device();

  // overwrite masstotal, center-of-mass, image flags with file values
  // inbody[i] = 0/1 if Ith rigid body is initialized by file
  DAT::tdual_int_1d k_inbody("rigid/small:inbody", nlocal_body);
  deep_copy(k_inbody.view_host(), 0);

  auto l_inbody = k_inbody.template view<DeviceType>();
  if (inpfile) {
    // must call it here so it doesn't override read in data but
    // initialize bodies whose dynamic settings not set in inpfile
    setup_bodies_dynamic();
    readfile(0, nullptr, k_inbody.view_host().data());
    k_inbody.modify_host();
  }

  // remap the xcm of each body back into simulation box
  //   and reset body and atom xcmimage flags via pre_neighbor()
  pre_neighbor();

  // compute 6 moments of inertia of each body in Cartesian reference frame
  // dx,dy,dz = coords relative to center-of-mass
  // symmetric 3x3 inertia tensor stored in Voigt notation as 6-vector
  //base()->memory->create(base()->itensor, nbody_total(), 6, "rigid/small:itensor");
  memoryKK->create_kokkos(k_itensor, itensor, nbody_total(), 6, "rigid/small:itensor");
  deep_copy(k_itensor.view_device(), KK_ZERO);

  auto l_itensor = k_itensor.template view<DeviceType>();
  auto lambda_itensor = [&]<bool RMASS, bool UNWRAP>(const int &i) {
    const int ibody = l_atom2body(i);
    if (ibody < 0) return;
    BodyKokkos &bk = l_body(ibody);
    KK_FLOAT massone;
    if constexpr (RMASS) massone = l_rmass(i);
    else massone = l_mass(l_type(i));
    KK_FLOAT dx, dy, dz;
    if constexpr (UNWRAP) {
      auto unwrap = DomainKokkos::unmap(l_prd, l_h, l_triclinic, &l_x(i,0), l_xcmimage(i));
      dx = unwrap[0] - bk.xcm[0];
      dy = unwrap[1] - bk.xcm[1];
      dz = unwrap[2] - bk.xcm[2];
    } else {
      dx = l_displace(i,0);
      dy = l_displace(i,1);
      dz = l_displace(i,2);
    }
    Kokkos::atomic_add(&l_itensor(ibody,0), massone * (dy*dy + dz*dz));
    Kokkos::atomic_add(&l_itensor(ibody,1), massone * (dx*dx + dz*dz));
    Kokkos::atomic_add(&l_itensor(ibody,2), massone * (dx*dx + dy*dy));
    Kokkos::atomic_sub(&l_itensor(ibody,3), massone * dy*dz);
    Kokkos::atomic_sub(&l_itensor(ibody,4), massone * dx*dz);
    Kokkos::atomic_sub(&l_itensor(ibody,5), massone * dx*dy);

  };
  if (l_rmass.data()) {
    Kokkos::parallel_for("rigid/small:setup_bodies_static_itensor_rmass",
    Kokkos::RangePolicy<DeviceType>(0, nlocal()),
    KOKKOS_LAMBDA(const int i) {
      lambda_itensor.template operator()<true,true>(i);
    });
  } else {
    Kokkos::parallel_for("rigid/small:setup_bodies_static_itensor",
    Kokkos::RangePolicy<DeviceType>(0, nlocal()),
    KOKKOS_LAMBDA(const int i) {
      lambda_itensor.template operator()<false,true>(i);
    });
  }
  copymode = 0;
  k_itensor.template modify<DeviceType>();

  // extended particles may contribute extra terms to moments of inertia
  if (extended) {
    // not implemented
  }

  // reverse communicate inertia tensor of all bodies
  commflag = ITENSOR;
  comm->reverse_comm(this, 6);

  // overwrite Cartesian inertia tensor with file values
  if (inpfile) {
    k_itensor.sync_host();
    k_inbody.sync_host();
    readfile(1, itensor, k_inbody.view_host().data());
    k_itensor.modify_host();
    k_inbody.modify_host();
  }

  // diagonalize inertia tensor for each body via Jacobi rotations
  // inertia = 3 eigenvalues = principal moments of inertia
  //   request that jacobi3() returns them in ascending order,
  //   so that in 2d last evector is z-axis
  // evectors and exzy_space = 3 evectors = principal axes of rigid body
  auto lambda_initial = [&]<bool DIMENSION2, bool TRICLINIC>(const int &ibody, KK_FLOAT &l_dflag) {
    BodyKokkos &bk = l_body(ibody);
    KK_FLOAT evectors[3][3], tensor[3][3] = {
      { l_itensor(ibody,0), l_itensor(ibody,5), l_itensor(ibody,4) },
      { l_itensor(ibody,5), l_itensor(ibody,1), l_itensor(ibody,3) },
      { l_itensor(ibody,4), l_itensor(ibody,3), l_itensor(ibody,2) }
    };
    // secondary benefit of switching from Jacobi to Cardano:
    // analytic methods do not have convergence failures
    const int ignored = MathExtraKokkos::sym3x3_eigen(tensor, bk.inertia, evectors, 1);
    bk.ex_space[0] = evectors[0][0];
    bk.ex_space[1] = evectors[1][0];
    bk.ex_space[2] = evectors[2][0];
    bk.ey_space[0] = evectors[0][1];
    bk.ey_space[1] = evectors[1][1];
    bk.ey_space[2] = evectors[2][1];
    bk.ez_space[0] = evectors[0][2];
    bk.ez_space[1] = evectors[1][2];
    bk.ez_space[2] = evectors[2][2];
    // for 2d, ensure that evector along z axis is last
    // necessary so that quaternion is a simple rotation around +z axis
    //   or a 180 degree rotation for a -z axis
    // otherwise richardson() method for a body with a tiny evalue (near-linear)
    //  may not preserve the correct z-aligned quat and associated evectors
    //  over time due to round-off accumulation
    if constexpr (DIMENSION2) {
      if (fabs(bk.ez_space[0]) > EPSILON || fabs(bk.ez_space[1]) > EPSILON) {
        Kokkos::kokkos_swap(bk.inertia[1], bk.inertia[2]);
        Kokkos::kokkos_swap(bk.ey_space[0], bk.ez_space[0]);
        Kokkos::kokkos_swap(bk.ey_space[1], bk.ez_space[1]);
        Kokkos::kokkos_swap(bk.ey_space[2], bk.ez_space[2]);
      }
    }
    // if any principal moment < scaled EPSILON, set to 0.0
    const KK_FLOAT max = Kokkos::max({bk.inertia[0], bk.inertia[1], bk.inertia[2]});
    if (bk.inertia[0] < EPSILON*max) bk.inertia[0] = KK_ZERO;
    if (bk.inertia[1] < EPSILON*max) bk.inertia[1] = KK_ZERO;
    if (bk.inertia[2] < EPSILON*max) bk.inertia[2] = KK_ZERO;
    // enforce 3 evectors as a right-handed coordinate system
    // flip 3rd vector if needed
    KK_FLOAT cross[3];
    MathExtraKokkos::cross3(bk.ex_space, bk.ey_space, cross);
    if (MathExtraKokkos::dot3(cross, bk.ez_space) < 0.0) MathExtraKokkos::negate3(bk.ez_space);
    // create initial quaternion
    exyz_to_q(bk.ex_space, bk.ey_space, bk.ez_space, bk.quat);
    // convert geometric center position to principal axis coordinates
    // xcm is wrapped, but xgc is not initially
    Few<KK_FLOAT,3> delta(
      bk.xgc[0] - bk.xcm[0],
      bk.xgc[1] - bk.xcm[1],
      bk.xgc[2] - bk.xcm[2]
    );
    DomainKokkos::minimum_image_big<TRICLINIC>(l_periodic, l_prd, l_invprd, l_h, delta, l_dflag);
    MathExtraKokkos::transpose_matvec(bk.ex_space, bk.ey_space, bk.ez_space, delta.data(), bk.xgc_body);
    bk.xgc[0] = bk.xcm[0] + delta[0];
    bk.xgc[1] = bk.xcm[1] + delta[1];
    bk.xgc[2] = bk.xcm[2] + delta[2];
  };

  k_itensor.template sync<DeviceType>();
  k_body.template sync<DeviceType>();
  copymode = 1;
  KK_FLOAT dflag = KK_ZERO;
  if (domainKK->dimension == 2) {
    if (l_triclinic) {
      Kokkos::parallel_reduce("rigid/small:setup_bodies_static_diagonalize_2d_triclinic",
      Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
        KOKKOS_LAMBDA(const int ibody, KK_FLOAT &l_dflag) {
        lambda_initial.template operator()<true,true>(ibody, l_dflag);
      }, Kokkos::Max<KK_FLOAT>(dflag));
    } else {
      Kokkos::parallel_reduce("rigid/small:setup_bodies_static_diagonalize_2d",
      Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
        KOKKOS_LAMBDA(const int ibody, KK_FLOAT &l_dflag) {
        lambda_initial.template operator()<true,false>(ibody, l_dflag);
      }, Kokkos::Max<KK_FLOAT>(dflag));
    }
  } else {
    if (l_triclinic) {
      Kokkos::parallel_reduce("rigid/small:setup_bodies_static_diagonalize_3d_triclinic",
      Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
        KOKKOS_LAMBDA(const int ibody, KK_FLOAT &l_dflag) {
        lambda_initial.template operator()<false,true>(ibody, l_dflag);
      }, Kokkos::Max<KK_FLOAT>(dflag));
    } else {
      Kokkos::parallel_reduce("rigid/small:setup_bodies_static_diagonalize_3d",
      Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
        KOKKOS_LAMBDA(const int ibody, KK_FLOAT &l_dflag) {
        lambda_initial.template operator()<false,false>(ibody, l_dflag);
      }, Kokkos::Max<KK_FLOAT>(dflag));
    }
  }
  copymode = 0;
  k_body.template modify<DeviceType>();

  if( dflag > KK_ZERO)
    error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dflag);

  // forward communicate updated info of all bodies
  commflag = INITIAL;
  comm->forward_comm(this, 29);

  // displace = initial atom coords in basis of principal axes
  // set displace = 0.0 for atoms not in any rigid body
  // for extended particles, set their orientation wrt to rigid body

  auto lambda_displace = [&]<bool EXTENDED>(const int i) {
    const int ibody = l_atom2body(i);
    if (ibody < 0) {
      l_displace(i,0) = l_displace(i,1) = l_displace(i,2) = KK_ZERO;
      return;
    }
    BodyKokkos &bk = l_body(ibody);
    auto unwrap = DomainKokkos::unmap(l_prd, l_h, l_triclinic, &l_x(i,0), l_xcmimage(i));
    auto delta3 = unwrap - Few<KK_FLOAT,3>(bk.xcm);
    transpose_matvec(bk.ex_space,bk.ey_space,bk.ez_space, delta3.data(), &l_displace(i,0));
    if constexpr (EXTENDED) {
      // not implemented
    }
  };

  k_displace.template sync<DeviceType>();
  copymode = 1;
  if (extended) {
    Kokkos::parallel_for("rigid/small:setup_bodies_static_displace_extended",
      Kokkos::RangePolicy<DeviceType>(0, nlocal()),
    KOKKOS_LAMBDA(const int i) {
      lambda_displace.template operator()<true>(i);
    });
  } else {
    Kokkos::parallel_for("rigid/small:setup_bodies_static_displace",
    Kokkos::RangePolicy<DeviceType>(0, nlocal()),
    KOKKOS_LAMBDA(const int i) {
      lambda_displace.template operator()<false>(i);
    });
  }
  copymode = 0;
  k_displace.template modify<DeviceType>();

  // test for valid principal moments & axes
  // recompute moments of inertia around new axes
  // 3 diagonal moments should equal principal moments
  // 3 off-diagonal moments should be 0.0
  // extended particles may contribute extra terms to moments of inertia
  deep_copy(k_itensor.view_device(), KK_ZERO);
  if (l_rmass.data()) {
    Kokkos::parallel_for("rigid/small:setup_bodies_static_itensor_rmass",
    Kokkos::RangePolicy<DeviceType>(0, nlocal()),
    KOKKOS_LAMBDA(const int i) {
      lambda_itensor.template operator()<true,false>(i);
    });
  } else {
    Kokkos::parallel_for("rigid/small:setup_bodies_static_itensor",
    Kokkos::RangePolicy<DeviceType>(0, nlocal()),
    KOKKOS_LAMBDA(const int i) {
      lambda_itensor.template operator()<false,false>(i);
    });
  }
  if (extended) {
    // not implemented
  }
  copymode = 0;
  k_itensor.modify_device();

  // reverse communicate inertia tensor of all bodies
  commflag = ITENSOR;
  comm->reverse_comm(this, 6);

  // error check that re-computed moments of inertia match diagonalized ones
  // do not do test for bodies with params read from inpfile

  int flag = 0;
  bool l_inpfile = (inpfile != nullptr);
  k_inbody.template sync<DeviceType>();
  copymode = 1;
  Kokkos::parallel_reduce("rigid/small:setup_bodies_static_check",
  Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
  KOKKOS_LAMBDA(const int ibody, int& l_flag) {
    if (l_inpfile && l_inbody(ibody)) return;
    const BodyKokkos &bk = l_body(ibody);
    KK_FLOAT norm = KK_ZERO;
    l_flag = 0;
    for( int i=0; i<3 ; i++ ) {
      const auto a = bk.inertia[i];
      const auto b = l_itensor(ibody,i);
      if (a == 0.0 && Kokkos::fabs(b) > TOLERANCE) l_flag = 1;
      if (a != 0.0 && Kokkos::fabs((b-a)/a) > TOLERANCE) l_flag = 1;
      norm += a;
    }
    norm *= KK_THIRD;
    if (Kokkos::fabs(l_itensor(ibody,3) / norm) > TOLERANCE ||
        Kokkos::fabs(l_itensor(ibody,4) / norm) > TOLERANCE ||
        Kokkos::fabs(l_itensor(ibody,5) / norm) > TOLERANCE) {
          l_flag = 1;
      }
  }, Kokkos::Max<int>(flag));
  copymode = 0;
  if (flag) error->all(FLERR, Error::NOLASTLINE, "Fix {}: Bad principal moments", style);

  // cleanup
  memoryKK->destroy_kokkos(k_itensor, itensor);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::setup_bodies_dynamic()
{

  // kokkos views

  atomKK->sync(execution_space, X_MASK | V_MASK | RMASS_MASK | TYPE_MASK );
  auto l_x = atomKK->k_x.template view<DeviceType>();
  auto l_v = atomKK->k_v.template view<DeviceType>();
  auto l_rmass = atomKK->k_rmass.template view<DeviceType>();
  auto l_type = atomKK->k_type.template view<DeviceType>();
  // FIXME: what *_MASK to sync for k_mass ???
  auto l_mass = atomKK->k_mass.template view<DeviceType>();

  k_atom2body.template sync<DeviceType>();
  auto l_atom2body = k_atom2body.template view<DeviceType>();

  k_body.template sync<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();

  k_xcmimage.template sync<DeviceType>();
  auto l_xcmimage = k_xcmimage.template view<DeviceType>();

  // domainKK

  auto l_prd = Few<KK_FLOAT,3>(domainKK->prd);
  auto l_h = Few<KK_FLOAT,6>(domainKK->h);
  auto l_triclinic = triclinic;

  // sum vcm, angmom across all rigid bodies
  // vcm = velocity of COM
  // angmom = angular momentum around COM

  copymode = 1;
  Kokkos::parallel_for("rigid/small:setup_bodies_dynamic_zero",
    Kokkos::RangePolicy<DeviceType>(0, nbody_total()),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      bk.vcm[0] = bk.vcm[1] = bk.vcm[2] = KK_ZERO;
      bk.angmom[0] = bk.angmom[1] = bk.angmom[2] = KK_ZERO;
    }
  );

  auto lambda = [&]<bool RMASS>(const int &i) {
    const int ibody = l_atom2body(i);
    if (ibody < 0) return;
    BodyKokkos &bk = l_body(ibody);
    KK_FLOAT massone;
    if constexpr (RMASS) massone = l_rmass(i);
    else massone = l_mass(l_type(i));
    Kokkos::atomic_add(&bk.vcm[0], l_v(i,0) * massone);
    Kokkos::atomic_add(&bk.vcm[1], l_v(i,1) * massone);
    Kokkos::atomic_add(&bk.vcm[2], l_v(i,2) * massone);
    auto unwrap = DomainKokkos::unmap(l_prd, l_h, l_triclinic, &l_x(i,0), l_xcmimage(i));
    const KK_FLOAT dx = unwrap[0] - bk.xcm[0];
    const KK_FLOAT dy = unwrap[1] - bk.xcm[1];
    const KK_FLOAT dz = unwrap[2] - bk.xcm[2];
    Kokkos::atomic_add(&bk.angmom[0], massone * (dy * l_v(i,2) - dz * l_v(i,1)));
    Kokkos::atomic_add(&bk.angmom[1], massone * (dz * l_v(i,0) - dx * l_v(i,2)));
    Kokkos::atomic_add(&bk.angmom[2], massone * (dx * l_v(i,1) - dy * l_v(i,0)));
  };
  auto policy = Kokkos::RangePolicy<DeviceType>(0, atomKK->nlocal);
  if (l_rmass.data()) {
    Kokkos::parallel_for("rigid/small:setup_bodies_dynamic_vcm_angmom_rmass",
      policy, KOKKOS_LAMBDA(const int i) { lambda.template operator()<true>(i); }
    );
  } else {
    Kokkos::parallel_for("rigid/small:setup_bodies_dynamic_vcm_angmom",
      policy, KOKKOS_LAMBDA(const int i) { lambda.template operator()<false>(i); }
    );
  }
  copymode = 0;
  k_body.modify_device();

  // extended particles add their rotation to angmom of body
  if (extended) {
    // not implemented
  }

  // reverse communicate vcm, angmom of all bodies
  commflag = VCM_ANGMOM;
  comm->reverse_comm(this, 6);

  // normalize velocity of COM
  k_body.template sync<DeviceType>();
  copymode = 1;
  Kokkos::parallel_for("rigid/small:setup_bodies_dynamic_normalize",
    Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      bk.vcm[0] /= bk.mass;
      bk.vcm[1] /= bk.mass;
      bk.vcm[2] /= bk.mass;
    }
  );
  copymode = 0;
  k_body.template modify<DeviceType>();

}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::apply_langevin_thermostat()
{
  // grow langextra if needed
  if (nlocal_body > maxlang) {
    maxlang = nlocal_body + nghost_body;
    memoryKK->destroy_kokkos(k_langextra, langextra);
    memoryKK->create_kokkos(k_langextra, langextra, 6, "rigid/small:langextra");
  }

  KK_FLOAT delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  const KK_FLOAT l_t_target = t_start + delta * (t_stop-t_start);
  const KK_FLOAT l_tsqrt = sqrt(l_t_target);
  const KK_FLOAT l_t_period = static_cast<KK_FLOAT>(t_period);

  const KK_FLOAT l_ftm2v = force->ftm2v;
  const KK_FLOAT l_langfactor = sqrt(24.0 * force->boltz / t_period
    / update->dt / force->mvv2e) / l_ftm2v;

  copymode = 1;
  k_body.sync_device();
  auto l_body = k_body.template view<DeviceType>();
  auto l_rand_pool = rand_pool;
  auto l_langextra = k_langextra.template view<DeviceType>();
  Kokkos::parallel_for("rigid/small:langevin",
    Kokkos::RangePolicy<DeviceType>(0, nbody_total()),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      KK_FLOAT gamma1 = -bk.mass / l_t_period / l_ftm2v;
      KK_FLOAT gamma2 = sqrt(bk.mass) * l_tsqrt * l_langfactor;

      rand_type rand_gen = l_rand_pool.get_state();
#if defined (LMP_KOKKOS_SINGLE_SINGLE) || defined (LMP_KOKKOS_SINGLE_DOUBLE)
      const float rnd1 = rand_gen.frand() - 0.5f;
      const float rnd2 = rand_gen.frand() - 0.5f;
      const float rnd3 = rand_gen.frand() - 0.5f;
      const float rnd4 = rand_gen.frand() - 0.5f;
      const float rnd5 = rand_gen.frand() - 0.5f;
      const float rnd6 = rand_gen.frand() - 0.5f;
#elif defined (LMP_KOKKOS_DOUBLE_DOUBLE)
      const double rnd1 = rand_gen.drand() - 0.5;
      const double rnd2 = rand_gen.drand() - 0.5;
      const double rnd3 = rand_gen.drand() - 0.5;
      const double rnd4 = rand_gen.drand() - 0.5;
      const double rnd5 = rand_gen.drand() - 0.5;
      const double rnd6 = rand_gen.drand() - 0.5;
#endif
      l_rand_pool.free_state(rand_gen);

      l_langextra(ibody, 0) = gamma1 * bk.vcm[0] + gamma2 * rnd1;
      l_langextra(ibody, 1) = gamma1 * bk.vcm[1] + gamma2 * rnd2;
      l_langextra(ibody, 2) = gamma1 * bk.vcm[2] + gamma2 * rnd3;
      gamma1 = -KK_ONE / l_t_period / l_ftm2v;
      gamma2 = l_tsqrt * l_langfactor;
      KK_FLOAT wbody[3], tbody[3];
      // convert omega from space frame to body frame
      transpose_matvec(bk.ex_space,bk.ey_space,bk.ez_space,bk.omega,wbody);
      // compute langevin torques in the body frame
      tbody[0] = bk.inertia[0] * gamma1 * wbody[0]
                 + sqrt(bk.inertia[0]) * gamma2 * rnd4;
      tbody[1] = bk.inertia[1] * gamma1 * wbody[1]
                 + sqrt(bk.inertia[1]) * gamma2 * rnd5;
      tbody[2] = bk.inertia[2] * gamma1 * wbody[2] +
                 + sqrt(bk.inertia[2]) * gamma2 * rnd6;
      // convert langevin torques from body frame back to space frame
      matvec(bk.ex_space,bk.ey_space,bk.ez_space,tbody,
                              &l_langextra(ibody, 3));
    }
  );
  copymode = 0;
  k_langextra.modify_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::remap()
{
  if constexpr (is_nh) {
  int nlocal = atomKK->nlocal;
  // epsilon is not used, except for book-keeping
  for (int i = 0; i < 3; i++) nh()->epsilon[i] += dtq * nh()->epsilon_dot[i];
  // convert pertinent atoms and rigid bodies to lamda coords
  if (allremap) domainKK->x2lamda(nlocal);
  else domainKK->x2lamda(nlocal, dilate_group_bit);
  for (auto &ifix : nh()->rfix) ifix->deform(0);
  // reset global and local box to new size/shape
  for (int i = 0; i < 3; i++) {
    if (p_flag[i]) {
      const double oldlo = domainKK->boxlo[i];
      const double oldhi = domainKK->boxhi[i];
      const double ctr = 0.5 * (oldlo + oldhi);
      const double expfac = exp(dtq * nh()->epsilon_dot[i]);
      domainKK->boxlo[i] = (oldlo-ctr)*expfac + ctr;
      domainKK->boxhi[i] = (oldhi-ctr)*expfac + ctr;
    }
  }
  domainKK->set_global_box();
  domainKK->set_local_box();
  // convert pertinent atoms and rigid bodies back to box coords
  if (allremap) domainKK->lamda2x(nlocal);
  else domainKK->lamda2x(nlocal, dilate_group_bit);
  for (auto &ifix : nh()->rfix) ifix->deform(1);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::deform(int flag)
{
  k_body.sync_device();
  auto l_body = k_body.template view<DeviceType>();
  copymode = 1;
  auto l_lo0 = domainKK->boxlo[0];
  auto l_lo1 = domainKK->boxlo[1];
  auto l_lo2 = domainKK->boxlo[2];
  if (flag == 0) {
    auto l_h_inv0 = domainKK->h_inv[0];
    auto l_h_inv1 = domainKK->h_inv[1];
    auto l_h_inv2 = domainKK->h_inv[2];
    auto l_h_inv3 = domainKK->h_inv[3];
    auto l_h_inv4 = domainKK->h_inv[4];
    auto l_h_inv5 = domainKK->h_inv[5];
    Kokkos::parallel_for("rigid/small:x2lamba",
      Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
      KOKKOS_LAMBDA(const int ibody) {
        BodyKokkos &bk = l_body(ibody);
        const KK_FLOAT delta0 = bk.xcm[0] - l_lo0;
        const KK_FLOAT delta1 = bk.xcm[1] - l_lo1;
        const KK_FLOAT delta2 = bk.xcm[2] - l_lo2;
        bk.xcm[0] = l_h_inv0*delta0 + l_h_inv5*delta1 + l_h_inv4*delta2;
        bk.xcm[1] = l_h_inv1*delta1 + l_h_inv3*delta2;
        bk.xcm[2] = l_h_inv2*delta2;
      }
    );
  }
  else {
    auto l_h0 = domainKK->h[0];
    auto l_h1 = domainKK->h[1];
    auto l_h2 = domainKK->h[2];
    auto l_h3 = domainKK->h[3];
    auto l_h4 = domainKK->h[4];
    auto l_h5 = domainKK->h[5];
    Kokkos::parallel_for("rigid/small:deform_lamda2x",
      Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
      KOKKOS_LAMBDA(const int ibody) {
        BodyKokkos &bk = l_body(ibody);
        const KK_FLOAT xi1 = bk.xcm[1];
        const KK_FLOAT xi2 = bk.xcm[2];
        bk.xcm[0] = l_h0*bk.xcm[0] + l_h5*xi1 + l_h4*xi2 + l_lo0;
        bk.xcm[1] = l_h1*xi1 + l_h3*xi2 + l_lo1;
        bk.xcm[2] = l_h2*xi2 + l_lo2;
      }
    );
  }
  copymode = 0;
  k_body.modify_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
template<bool EVFLAG>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::set_xv()
{

  // kokkos views

  atomKK->sync(execution_space, X_MASK | V_MASK | F_MASK | TYPE_MASK | RMASS_MASK);
  auto l_x = atomKK->k_x.template view<DeviceType>();
  auto l_v = atomKK->k_v.template view<DeviceType>();
  auto l_f = atomKK->k_f.template view<DeviceType>();
  auto l_rmass = atomKK->k_rmass.template view<DeviceType>();
  auto l_mass = atomKK->k_mass.template view<DeviceType>();
  auto l_type = atomKK->k_type.template view<DeviceType>();

  k_atom2body.template sync<DeviceType>();
  auto l_atom2body = k_atom2body.template view<DeviceType>();

  k_body.template sync<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();

  k_xcmimage.template sync<DeviceType>();
  auto l_xcmimage = k_xcmimage.template view<DeviceType>();

  k_displace.template sync<DeviceType>();
  auto l_displace = k_displace.template view<DeviceType>();

  // domainKK

  auto l_h0 = domainKK->h[0];
  auto l_h1 = domainKK->h[1];
  auto l_h2 = domainKK->h[2];
  auto l_h3 = domainKK->h[3];
  auto l_h4 = domainKK->h[4];
  auto l_h5 = domainKK->h[5];
  auto l_prd0 = static_cast<KK_FLOAT>(domainKK->prd[0]);
  auto l_prd1 = static_cast<KK_FLOAT>(domainKK->prd[1]);
  auto l_prd2 = static_cast<KK_FLOAT>(domainKK->prd[2]);

  KK_ACC_FLOAT l_half_dt = static_cast<KK_ACC_FLOAT>(0.5 / dtf);
  auto l_vflag_global = vflag_global;
  auto l_vflag_atom = vflag_atom;

  auto lambda = [&]<bool TRICLINIC, int NEIGHFLAG>(const int &i, EV_FLOAT &ev) {
    const int ibody = l_atom2body(i);
    if (ibody < 0) return;
    const BodyKokkos &bk = l_body(ibody);

    const KK_FLOAT xbox = static_cast<KK_FLOAT>((l_xcmimage(i) & IMGMASK) - IMGMAX);
    const KK_FLOAT ybox = static_cast<KK_FLOAT>((l_xcmimage(i) >> IMGBITS & IMGMASK) - IMGMAX);
    const KK_FLOAT zbox = static_cast<KK_FLOAT>((l_xcmimage(i) >> IMG2BITS) - IMGMAX);

    KK_FLOAT deltax, deltay;
    if constexpr(TRICLINIC) {
      deltax = fma(xbox, l_prd0, fma(ybox, l_h5, zbox * l_h4));
      deltay = fma(ybox, l_prd1, zbox * l_h3);
    } else {
      deltax = xbox * l_prd0;
      deltay = ybox * l_prd1;
    }
    const KK_FLOAT deltaz = zbox * l_prd2;

    KK_FLOAT x0 = KK_ZERO, x1 = KK_ZERO, x2 = KK_ZERO;
    KK_FLOAT vx = KK_ZERO, vy = KK_ZERO, vz = KK_ZERO;
    if constexpr (EVFLAG) {
      x0 = l_x(i,0) + deltax;
      x1 = l_x(i,1) + deltay;
      x2 = l_x(i,2) + deltaz;
      vx = l_v(i,0);
      vy = l_v(i,1);
      vz = l_v(i,2);
    }

    KK_FLOAT xnew[3];
    matvec(bk.ex_space, bk.ey_space, bk.ez_space, &l_displace(i,0), xnew);

    if constexpr (EVFLAG) {
      // Compute v_new in KK_ACC_FLOAT before truncating to KK_FLOAT for storage,
      // so the pre-truncation value can be used for the constraint-force virial.
      const KK_ACC_FLOAT vx = fma(KK_ACC_FLOAT(bk.omega[1]), KK_ACC_FLOAT(xnew[2]),
                               fma(KK_ACC_FLOAT(-bk.omega[2]), KK_ACC_FLOAT(xnew[1]),
                               KK_ACC_FLOAT(bk.vcm[0])));
          const KK_ACC_FLOAT vy = fma(KK_ACC_FLOAT(bk.omega[2]), KK_ACC_FLOAT(xnew[0]),
                               fma(KK_ACC_FLOAT(-bk.omega[0]), KK_ACC_FLOAT(xnew[2]),
                               KK_ACC_FLOAT(bk.vcm[1])));
          const KK_ACC_FLOAT vz = fma(KK_ACC_FLOAT(bk.omega[0]), KK_ACC_FLOAT(xnew[1]),
                               fma(KK_ACC_FLOAT(-bk.omega[1]), KK_ACC_FLOAT(xnew[0]),
                               KK_ACC_FLOAT(bk.vcm[2])));
          const KK_ACC_FLOAT dvx = vx - KK_ACC_FLOAT(l_v(i,0));
          const KK_ACC_FLOAT dvy = vy - KK_ACC_FLOAT(l_v(i,1));
          const KK_ACC_FLOAT dvz = vz - KK_ACC_FLOAT(l_v(i,2));
          l_v(i,0) = KK_FLOAT(vx);
          l_v(i,1) = KK_FLOAT(vy);
          l_v(i,2) = KK_FLOAT(vz);
          l_x(i,0) = xnew[0] + bk.xcm[0] - deltax;
          l_x(i,1) = xnew[1] + bk.xcm[1] - deltay;
          l_x(i,2) = xnew[2] + bk.xcm[2] - deltaz;

          KK_ACC_FLOAT massone;
          if (l_rmass.data()) massone = l_rmass(i);
          else massone = l_mass(l_type(i));
          const KK_ACC_FLOAT half_m_dt = l_half_dt * massone;
          const KK_ACC_FLOAT fc0 = fma(half_m_dt, dvx, KK_ACC_FLOAT(-0.5)*KK_ACC_FLOAT(l_f(i,0)));
          const KK_ACC_FLOAT fc1 = fma(half_m_dt, dvy, KK_ACC_FLOAT(-0.5)*KK_ACC_FLOAT(l_f(i,1)));
          const KK_ACC_FLOAT fc2 = fma(half_m_dt, dvz, KK_ACC_FLOAT(-0.5)*KK_ACC_FLOAT(l_f(i,2)));

          const KK_ACC_FLOAT vd00 = x0*fc0;
          const KK_ACC_FLOAT vd11 = x1*fc1;
          const KK_ACC_FLOAT vd22 = x2*fc2;
          const KK_ACC_FLOAT vd01 = x0*fc1;
          const KK_ACC_FLOAT vd02 = x0*fc2;
          const KK_ACC_FLOAT vd12 = x1*fc2;

          if (l_vflag_global) {
            ev.v[0] += vd00;
            ev.v[1] += vd11;
            ev.v[2] += vd22;
            ev.v[3] += vd01;
            ev.v[4] += vd02;
            ev.v[5] += vd12;
          }
          if (l_vflag_atom) {
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
        } else {
          l_v(i,0) = fma(bk.omega[1], xnew[2], fma(-bk.omega[2], xnew[1], bk.vcm[0]));
          l_v(i,1) = fma(bk.omega[2], xnew[0], fma(-bk.omega[0], xnew[2], bk.vcm[1]));
          l_v(i,2) = fma(bk.omega[0], xnew[1], fma(-bk.omega[1], xnew[0], bk.vcm[2]));
          l_x(i,0) = xnew[0] + bk.xcm[0] - deltax;
          l_x(i,1) = xnew[1] + bk.xcm[1] - deltay;
          l_x(i,2) = xnew[2] + bk.xcm[2] - deltaz;
        }
      };

  const int nlocal = atomKK->nlocal;
  copymode = 1;
  if constexpr (!EVFLAG) {
    // TagRigidSetXV<TRICLINIC,HALF,0>>(0, nlocal),
    EV_FLOAT ev;
    if (triclinic) {
      Kokkos::parallel_for("rigid/small:set_xv_triclinic",
        Kokkos::RangePolicy<DeviceType>(0, nlocal),
        KOKKOS_LAMBDA(const int &i) {
          EV_FLOAT ev;
          lambda.template operator()<true,HALF>(i, ev);
        }
      );
    } else {
      Kokkos::parallel_for("rigid/small:set_xv",
        Kokkos::RangePolicy<DeviceType>(0, nlocal),
        KOKKOS_LAMBDA(const int &i) {
          EV_FLOAT ev;
          lambda.template operator()<false,HALF>(i, ev);
        }
      );
    }
  } else {
    auto kokkos = lmp->kokkos;
    int neighflag = kokkos->neighflag;
    if (neighflag == FULL) {
      neighflag =
        (kokkos->nthreads > 1 || kokkos->ngpus > 0) ? HALFTHREAD : HALF;
    }
    const bool need_dup =
      (neighflag != HALF) &&
      std::is_same_v<NeedDup_v<HALFTHREAD, DeviceType>,
                     Kokkos::Experimental::ScatterDuplicated>;
    if (l_vflag_atom) {
      if (need_dup)
        dup_vatom = Kokkos::Experimental::create_scatter_view<
          Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
      else
        ndup_vatom = Kokkos::Experimental::create_scatter_view<
          Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
    }
    EV_FLOAT ev{};
    if (neighflag == HALF) {
      // TagRigidSetXV<TRICLINIC,HALF>>(0, nlocal),
      if (triclinic) {
        Kokkos::parallel_reduce("rigid/small:set_xv_evflag_half_triclinic",
          Kokkos::RangePolicy<DeviceType>(0, nlocal),
          KOKKOS_LAMBDA(const int &i, EV_FLOAT &ev_) {
            lambda.template operator()<true,HALF>(i, ev_);
          }, ev
        );
      } else {
        Kokkos::parallel_reduce("rigid/small:set_xv_evflag_half",
          Kokkos::RangePolicy<DeviceType>(0, nlocal),
          KOKKOS_LAMBDA(const int &i, EV_FLOAT &ev_) {
            lambda.template operator()<false,HALF>(i, ev_);
          }, ev
        );
      }
    } else {
      // TagRigidSetXV<TRICLINIC,HALFTHREAD>>(0, nlocal),
      if (triclinic) {
        Kokkos::parallel_reduce("rigid/small:set_xv_evflag_halfthread_triclinic",
          Kokkos::RangePolicy<DeviceType>(0, nlocal),
          KOKKOS_LAMBDA(const int &i, EV_FLOAT &ev_) {
            lambda.template operator()<true,HALFTHREAD>(i, ev_);
          }, ev
        );
      } else {
        Kokkos::parallel_reduce("rigid/small:set_xv_evflag_halfthread",
          Kokkos::RangePolicy<DeviceType>(0, nlocal),
          KOKKOS_LAMBDA(const int &i, EV_FLOAT &ev_) {
            lambda.template operator()<false,HALFTHREAD>(i, ev_);
          }, ev
        );
      }
    }
    if (l_vflag_global) {
      virial[0] += static_cast<double>(ev.v[0]);
      virial[1] += static_cast<double>(ev.v[1]);
      virial[2] += static_cast<double>(ev.v[2]);
      virial[3] += static_cast<double>(ev.v[3]);
      virial[4] += static_cast<double>(ev.v[4]);
      virial[5] += static_cast<double>(ev.v[5]);
    }
    if (l_vflag_atom && need_dup) {
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
      dup_vatom = {};
    }
  }
  // update geometric center of bodies
  Kokkos::parallel_for("rigid/small:set_xv_xgc",
    Kokkos::RangePolicy<DeviceType>(0, nbody_total()),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      KK_FLOAT xgc_tmp[3];
      matvec(bk.ex_space, bk.ey_space, bk.ez_space, bk.xgc_body, xgc_tmp);
      bk.xgc[0] = xgc_tmp[0] + bk.xcm[0];
      bk.xgc[1] = xgc_tmp[1] + bk.xcm[1];
      bk.xgc[2] = xgc_tmp[2] + bk.xcm[2];
    }
  );
  copymode = 0;
  atomKK->modified(execution_space, X_MASK | V_MASK);
  k_body.modify_device();

  if (extended) {
    // not implemented
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
template<bool EVFLAG>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::set_v()
{

  atomKK->sync(execution_space, X_MASK | V_MASK | F_MASK | TYPE_MASK | RMASS_MASK);
  auto l_x = atomKK->k_x.template view<DeviceType>();
  auto l_v = atomKK->k_v.template view<DeviceType>();
  auto l_f = atomKK->k_f.template view<DeviceType>();
  auto l_rmass = atomKK->k_rmass.template view<DeviceType>();
  auto l_mass = atomKK->k_mass.template view<DeviceType>();
  auto l_type = atomKK->k_type.template view<DeviceType>();

  k_atom2body.sync_device();
  k_body.sync_device();
  k_xcmimage.sync_device();
  k_displace.sync_device();
  auto l_atom2body = k_atom2body.template view<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();
  auto l_xcmimage = k_xcmimage.template view<DeviceType>();
  auto l_displace = k_displace.template view<DeviceType>();

  auto l_h0 = domainKK->h[0];
  auto l_h1 = domainKK->h[1];
  auto l_h2 = domainKK->h[2];
  auto l_h3 = domainKK->h[3];
  auto l_h4 = domainKK->h[4];
  auto l_h5 = domainKK->h[5];
  auto l_prd0 = static_cast<KK_FLOAT>(domainKK->prd[0]);
  auto l_prd1 = static_cast<KK_FLOAT>(domainKK->prd[1]);
  auto l_prd2 = static_cast<KK_FLOAT>(domainKK->prd[2]);

  KK_ACC_FLOAT l_half_dt = static_cast<KK_ACC_FLOAT>(0.5 / dtf);
  auto l_vflag_global = vflag_global;
  auto l_vflag_atom = vflag_atom;

  auto lambda = [&]<bool TRICLINIC, int NEIGHFLAG>(const int &i, EV_FLOAT &ev) {
    const int ibody = l_atom2body(i);
    if (ibody < 0) return;
    const BodyKokkos &bk = l_body(ibody);

    KK_FLOAT delta[3];
    matvec(bk.ex_space, bk.ey_space, bk.ez_space, &l_displace(i, 0), delta);

    if constexpr (EVFLAG) {
      // Compute v_new in KK_ACC_FLOAT before truncating to KK_FLOAT for storage,
      // so the pre-truncation value can be used for the constraint-force virial.
      const KK_ACC_FLOAT vx = fma(KK_ACC_FLOAT(bk.omega[1]), KK_ACC_FLOAT(delta[2]),
                               fma(KK_ACC_FLOAT(-bk.omega[2]), KK_ACC_FLOAT(delta[1]),
                               KK_ACC_FLOAT(bk.vcm[0])));
      const KK_ACC_FLOAT vy = fma(KK_ACC_FLOAT(bk.omega[2]), KK_ACC_FLOAT(delta[0]),
                               fma(KK_ACC_FLOAT(-bk.omega[0]), KK_ACC_FLOAT(delta[2]),
                               KK_ACC_FLOAT(bk.vcm[1])));
      const KK_ACC_FLOAT vz = fma(KK_ACC_FLOAT(bk.omega[0]), KK_ACC_FLOAT(delta[1]),
                               fma(KK_ACC_FLOAT(-bk.omega[1]), KK_ACC_FLOAT(delta[0]),
                               KK_ACC_FLOAT(bk.vcm[2])));
      const KK_ACC_FLOAT dvx = vx - KK_ACC_FLOAT(l_v(i,0));
      const KK_ACC_FLOAT dvy = vy - KK_ACC_FLOAT(l_v(i,1));
      const KK_ACC_FLOAT dvz = vz - KK_ACC_FLOAT(l_v(i,2));
      l_v(i,0) = KK_FLOAT(vx);
      l_v(i,1) = KK_FLOAT(vy);
      l_v(i,2) = KK_FLOAT(vz);
      KK_ACC_FLOAT massone;
      if (l_rmass.data()) massone = l_rmass(i);
      else massone = l_mass(l_type(i));
      const KK_ACC_FLOAT half_m_dt = l_half_dt * massone;
      const KK_ACC_FLOAT fc0 = fma(half_m_dt, dvx, KK_ACC_FLOAT(-0.5)*KK_ACC_FLOAT(l_f(i,0)));
      const KK_ACC_FLOAT fc1 = fma(half_m_dt, dvy, KK_ACC_FLOAT(-0.5)*KK_ACC_FLOAT(l_f(i,1)));
      const KK_ACC_FLOAT fc2 = fma(half_m_dt, dvz, KK_ACC_FLOAT(-0.5)*KK_ACC_FLOAT(l_f(i,2)));

      const KK_FLOAT xbox = static_cast<KK_FLOAT>((l_xcmimage(i) & IMGMASK) - IMGMAX);
      const KK_FLOAT ybox = static_cast<KK_FLOAT>((l_xcmimage(i) >> IMGBITS & IMGMASK) - IMGMAX);
      const KK_FLOAT zbox = static_cast<KK_FLOAT>((l_xcmimage(i) >> IMG2BITS) - IMGMAX);

      KK_ACC_FLOAT x0, x1;
      if constexpr(TRICLINIC) {
        x0 = fma(xbox, l_prd0, fma(ybox, l_h5, fma(zbox, l_h4, l_x(i,0))));
        x1 = fma(ybox, l_prd1, fma(zbox, l_h3, l_x(i,1)));
      } else {
        x0 = fma(xbox, l_prd0, l_x(i,0));
        x1 = fma(ybox, l_prd1, l_x(i,1));
      }
      const KK_ACC_FLOAT x2 = fma(zbox, l_prd2, l_x(i,2));

      const KK_ACC_FLOAT vd00 = x0*fc0;
      const KK_ACC_FLOAT vd11 = x1*fc1;
      const KK_ACC_FLOAT vd22 = x2*fc2;
      const KK_ACC_FLOAT vd01 = x0*fc1;
      const KK_ACC_FLOAT vd02 = x0*fc2;
      const KK_ACC_FLOAT vd12 = x1*fc2;

      if (l_vflag_global) {
        ev.v[0] += vd00;
        ev.v[1] += vd11;
        ev.v[2] += vd22;
        ev.v[3] += vd01;
        ev.v[4] += vd02;
        ev.v[5] += vd12;
      }
      if (l_vflag_atom) {
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
    } else {
      l_v(i,0) = fma(bk.omega[1], delta[2], fma(-bk.omega[2], delta[1], bk.vcm[0]));
      l_v(i,1) = fma(bk.omega[2], delta[0], fma(-bk.omega[0], delta[2], bk.vcm[1]));
      l_v(i,2) = fma(bk.omega[0], delta[1], fma(-bk.omega[1], delta[0], bk.vcm[2]));
    }
  };

  const int nlocal = atomKK->nlocal;
  copymode = 1;
  if constexpr (!EVFLAG) {
    // TagRigidSetV<TRICLINIC,HALF,0>>(0, nlocal),
    if (triclinic) {
      Kokkos::parallel_for("rigid/small:set_v_half_triclinic",
        Kokkos::RangePolicy<DeviceType>(0, nlocal),
        KOKKOS_LAMBDA(const int &i) {
          EV_FLOAT ev;
          lambda.template operator()<true,HALF>(i, ev);
        }
      );
    } else {
      Kokkos::parallel_for("rigid/small:set_v_half",
        Kokkos::RangePolicy<DeviceType>(0, nlocal),
        KOKKOS_LAMBDA(const int &i) {
          EV_FLOAT ev;
          lambda.template operator()<false,HALF>(i, ev);
        }
      );
    }
  } else {
    auto kokkos = lmp->kokkos;
    int neighflag = kokkos->neighflag;
    if (neighflag == FULL) {
      neighflag =
        (kokkos->nthreads > 1 || kokkos->ngpus > 0) ? HALFTHREAD : HALF;
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
      // TagRigidSetV<TRICLINIC,HALF>>(0, nlocal),
      if (triclinic) {
        Kokkos::parallel_reduce("rigid/small:set_v_evflag_half_triclinic",
          Kokkos::RangePolicy<DeviceType>(0, nlocal),
          KOKKOS_LAMBDA(const int &i, EV_FLOAT &ev_) {
            lambda.template operator()<true,HALF>(i, ev_);
          }, ev
        );
      } else {
        Kokkos::parallel_reduce("rigid/small:set_v_evflag_half",
          Kokkos::RangePolicy<DeviceType>(0, nlocal),
          KOKKOS_LAMBDA(const int &i, EV_FLOAT &ev_) {
            lambda.template operator()<false,HALF>(i, ev_);
          }, ev
        );
      }
    } else {
      // TagRigidSetV<TRICLINIC,HALFTHREAD>>(0, nlocal),
      if (triclinic) {
        Kokkos::parallel_reduce("rigid/small:set_v_evflag_halfthread_triclinic",
          Kokkos::RangePolicy<DeviceType>(0, nlocal),
          KOKKOS_LAMBDA(const int &i, EV_FLOAT &ev_) {
            lambda.template operator()<true,HALFTHREAD>(i, ev_);
          }, ev
        );
      } else {
        Kokkos::parallel_reduce("rigid/small:set_v_evflag_halfthread",
          Kokkos::RangePolicy<DeviceType>(0, nlocal),
          KOKKOS_LAMBDA(const int &i, EV_FLOAT &ev_) {
            lambda.template operator()<false,HALFTHREAD>(i, ev_);
          }, ev
        );
      }
    }
    if (l_vflag_global) {
      virial[0] += static_cast<double>(ev.v[0]);
      virial[1] += static_cast<double>(ev.v[1]);
      virial[2] += static_cast<double>(ev.v[2]);
      virial[3] += static_cast<double>(ev.v[3]);
      virial[4] += static_cast<double>(ev.v[4]);
      virial[5] += static_cast<double>(ev.v[5]);
    }
    if (l_vflag_atom && need_dup) {
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
      dup_vatom = {};
    }
  }
  atomKK->modified(execution_space, V_MASK);
  k_body.modify_device();
  copymode = 0;

  if (extended) {
    // not implemented
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::compute_forces_and_torques()
{

  const int nlocal = atomKK->nlocal;

  // kokkos views

  atomKK->sync(execution_space, X_MASK | F_MASK );
  auto l_x = atomKK->k_x.template view<DeviceType>();
  auto l_f = atomKK->k_f.template view<DeviceType>();

  k_atom2body.template sync<DeviceType>();
  auto l_atom2body = k_atom2body.template view<DeviceType>();

  k_body.template sync<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();

  k_xcmimage.template sync<DeviceType>();
  auto l_xcmimage = k_xcmimage.template view<DeviceType>();

  // domainKK

  auto l_prd = Few<KK_FLOAT,3>(domainKK->prd);
  auto l_h = Few<KK_FLOAT,6>(domainKK->h);
  auto l_triclinic = triclinic;

  copymode = 1;
  Kokkos::parallel_for("rigid/small:compute_forces_and_torques_zero",
    Kokkos::RangePolicy<DeviceType>(0, nbody_total()),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      bk.fcm[0] = bk.fcm[1] = bk.fcm[2] = KK_ZERO;
      bk.torque[0] = bk.torque[1] = bk.torque[2] = KK_ZERO;
    }
  );
  Kokkos::parallel_for("rigid/small:compute_forces_and_torques",
    Kokkos::RangePolicy<DeviceType>(0, nlocal),
    KOKKOS_LAMBDA(const int &i) {
      const int ibody = l_atom2body(i);
      if (ibody < 0) return;
      BodyKokkos &bk = l_body(ibody);
      Kokkos::atomic_add(&bk.fcm[0], l_f(i,0));
      Kokkos::atomic_add(&bk.fcm[1], l_f(i,1));
      Kokkos::atomic_add(&bk.fcm[2], l_f(i,2));
      auto unwrap = DomainKokkos::unmap(l_prd, l_h, l_triclinic, &l_x(i,0), l_xcmimage(i));
      const KK_FLOAT dx = unwrap[0] - bk.xcm[0];
      const KK_FLOAT dy = unwrap[1] - bk.xcm[1];
      const KK_FLOAT dz = unwrap[2] - bk.xcm[2];
      Kokkos::atomic_add(&bk.torque[0], fma(dy, l_f(i,2), -dz*l_f(i,1)));
      Kokkos::atomic_add(&bk.torque[1], fma(dz, l_f(i,0), -dx*l_f(i,2)));
      Kokkos::atomic_add(&bk.torque[2], fma(dx, l_f(i,1), -dy*l_f(i,0)));
    }
  );
  copymode = 0;
  k_body.template modify<DeviceType>();

  if (extended) {
    // not implemented
  }

  commflag = FORCE_TORQUE;
  comm->reverse_comm(this, 6);

  if (langflag) {
    k_body.sync_device();
    auto l_body = k_body.template view<DeviceType>();
    auto l_langextra = k_langextra.template view<DeviceType>();
    copymode = 1;
    Kokkos::parallel_for("rigid/small:compute_forces_and_torques_langflag",
      Kokkos::RangePolicy<DeviceType>(0, nbody_total()),
      KOKKOS_LAMBDA(const int &ibody) {
        BodyKokkos &bk = l_body(ibody);
        bk.fcm[0] += l_langextra(ibody, 0);
        bk.fcm[1] += l_langextra(ibody, 1);
        bk.fcm[2] += l_langextra(ibody, 2);
        bk.torque[0] += l_langextra(ibody, 3);
        bk.torque[1] += l_langextra(ibody, 4);
        bk.torque[2] += l_langextra(ibody, 5);
      }
    );
    copymode = 0;
  }
  if (id_gravity) {
    k_body.sync_device();
    auto l_body = k_body.template view<DeviceType>();
    auto l_gvec0 = gvec[0];
    auto l_gvec1 = gvec[1];
    auto l_gvec2 = gvec[2];
    copymode = 1;
    Kokkos::parallel_for("rigid/small:compute_forces_and_torques_gravity",
      Kokkos::RangePolicy<DeviceType>(0, nbody_total()),
      KOKKOS_LAMBDA(const int &ibody) {
        BodyKokkos &bk = l_body(ibody);
        bk.fcm[0] += l_gvec0 * bk.mass;
        bk.fcm[1] += l_gvec1 * bk.mass;
        bk.fcm[2] += l_gvec2 * bk.mass;
      }
    );
    copymode = 0;
  }
  if (langflag || id_gravity) k_body.modify_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::enforce2d()
{
  k_body.sync_device();
  auto l_body = k_body.template view<DeviceType>();
  copymode = 1;
  Kokkos::parallel_for("rigid/small:enforce2d",
    Kokkos::RangePolicy<DeviceType>(0, nlocal_body),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      bk.xcm[2] = bk.vcm[2] = bk.fcm[2] = bk.xgc[2] = KK_ZERO;
      bk.torque[0] = bk.torque[1] = KK_ZERO;
      bk.angmom[0] = bk.angmom[1] = KK_ZERO;
      bk.omega[0] = bk.omega[1] = KK_ZERO;
    }
  );
  copymode = 0;
  k_body.modify_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::image_shift()
{
  copymode = 1;
  atomKK->sync(execution_space, IMAGE_MASK);
  k_atom2body.template sync<DeviceType>();
  k_xcmimage.template sync<DeviceType>();
  k_body.sync_device();
  auto l_image = atomKK->k_image.template view<DeviceType>();
  auto l_atom2body = k_atom2body.template view<DeviceType>();
  auto l_xcmimage = k_xcmimage.template view<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();
  Kokkos::parallel_for("rigid/small:image_shift",
    Kokkos::RangePolicy<DeviceType>(0, nlocal()),
    KOKKOS_LAMBDA(const int &i) {
      const int ibody = l_atom2body(i);
      if (ibody < 0) return;
      const BodyKokkos &bk = l_body(ibody);
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
    }
  );
  k_xcmimage.template modify<DeviceType>();
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::reset_atom2body()
{
  auto l_map_style = atomKK->map_style;
  auto l_map_array = atomKK->k_map_array;
  auto l_map_hash = atomKK->k_map_hash;
  if (l_map_style == Atom::MAP_ARRAY) l_map_array.template sync<DeviceType>();
  else if (l_map_style == Atom::MAP_HASH) l_map_hash.template sync<DeviceType>();

  atomKK->sync(execution_space, TAG_MASK);
  k_atom2body.template sync<DeviceType>();
  k_bodytag.template sync<DeviceType>();
  k_bodyown.template sync<DeviceType>();
  auto l_tag = atomKK->k_tag.template view<DeviceType>();
  auto l_atom2body = k_atom2body.template view<DeviceType>();
  auto l_bodytag = k_bodytag.template view<DeviceType>();
  auto l_bodyown = k_bodyown.template view<DeviceType>();
  auto l_comm_me = comm->me;
  auto l_ntimestep = update->ntimestep;

  copymode = 1;
  Kokkos::parallel_for("rigid/small:reset_atom2body",
    Kokkos::RangePolicy<DeviceType>(0, atomKK->nlocal),
    KOKKOS_LAMBDA(const int &i) {
      l_atom2body(i) = -1;
      if (l_bodytag(i)) {
        const int iowner = AtomKokkos::map_kokkos<DeviceType>(
          l_bodytag(i), l_map_style, l_map_array, l_map_hash);
        if (iowner == -1) {
          Kokkos::printf("Rigid body atoms %lld %lld missing on proc %d at step %lld\n",
                         (long long) l_tag(i), (long long) l_bodytag(i),
                         l_comm_me, (long long) l_ntimestep);
          Kokkos::abort("Rigid body atom missing");
        }
        l_atom2body(i) = l_bodyown(iowner);
      }
    }
  );
  copymode = 0;
  k_atom2body.modify_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::compute_dof()
{

  // kokkos view
  k_body.template sync<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();

  // total translational and rotational degrees of freedom
  const int l_dimension = domainKK->dimension;
  const int l_nlocal_body = nlocal_body;
  double nf[2] = {static_cast<double>(l_dimension * l_nlocal_body)};
  if (l_dimension == 3) {
    int nf_r = 0;
    Kokkos::parallel_reduce("rigid/small:compute_dof",
      Kokkos::RangePolicy<DeviceType>(0, l_nlocal_body),
      KOKKOS_LAMBDA(const int ibody, int &l_nf_r) {
        auto inertia = l_body(ibody).inertia;
        l_nf_r += 3;
        if (Kokkos::fabs(inertia[0]) < EPSILON) l_nf_r--;
        if (Kokkos::fabs(inertia[1]) < EPSILON) l_nf_r--;
        if (Kokkos::fabs(inertia[2]) < EPSILON) l_nf_r--;
      }, nf_r
    );
    nf[1] = static_cast<double>(nf_r);
  } else if (l_dimension == 2) nf[1] = l_nlocal_body;

  double nf_all[2];
  MPI_Allreduce(nf, nf_all, 2, MPI_DOUBLE, MPI_SUM, world);
  nh()->nf_t = static_cast<int>(nf_all[0]);
  nh()->nf_r = static_cast<int>(nf_all[1]);
  nh()->g_f = nh()->nf_t + nh()->nf_r;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::grow_arrays(int nmax)
{
  // Do not sync from device before grow: uninitialized device data must not
  // overwrite the host arrays this routine is about to resize.

  memoryKK->grow_kokkos(k_bodyown, bodyown, nmax, "rigid/small:bodyown");
  memoryKK->grow_kokkos(k_bodytag, bodytag, nmax, "rigid/small:bodytag");
  memoryKK->grow_kokkos(k_atom2body, atom2body, nmax, "rigid/small:atom2body");
  memoryKK->grow_kokkos(k_xcmimage, xcmimage, nmax, "rigid/small:xcmimage");

  // k_displace is a TransformView: grow_kokkos cannot maintain the displace[i]
  // pointer array into the view's host buffer, so resize explicitly and
  // re-point each displace[i] into the new host allocation.
  k_displace.sync_host();
  k_displace.resize(nmax, 3);
  bigint nbytes = ((bigint) sizeof(double *)) * nmax;
  displace = (double **) memory->srealloc(displace, nbytes, "rigid/small:displace");
  for (int i = 0; i < nmax; i++)
    displace[i] = (k_displace.extent_int(1) > 0) ? &k_displace.view_host()(i, 0) : nullptr;
  k_displace.modify_host();

  if (extended) {
    k_eflags.template sync<DeviceType>();
    memoryKK->grow_kokkos(k_eflags, eflags, nmax, "rigid/small:eflags");
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

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::grow_body()
{
  grow_body(nmax_body + DELTA_BODY);
}

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::grow_body(int nmax_body)
{
  nmax_body = nmax_body;
  k_body.sync_host();
  k_body.resize(nmax_body);
  body = k_body.view_host().data();
  k_body.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::copy_arrays(int i, int j, int delflag)
{
  sync_host();
  FixRigidBase::copy_arrays(i, j, delflag);
  modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::set_arrays(int i)
{
  sync_host();
  FixRigidBase::set_arrays(i);
  modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::set_molecule(int nlocalprev, tagint tagprev,
                                                    int imol, double *xgeom,
                                                    double *vcm, double *quat)
{
  sync_host();
  FixRigidBase::set_molecule(nlocalprev, tagprev, imol, xgeom, vcm, quat);
  modify_host();
}

/* ----------------------------------------------------------------------
  HOST COMM
------------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
int FixRigidBaseKokkos<DeviceType,FixRigidBase>::pack_exchange(int i, double *buf)
{

  k_body.template sync<LMPHostType>();
  k_bodyown.template sync<LMPHostType>();
  k_bodytag.template sync<LMPHostType>();
  k_atom2body.template sync<LMPHostType>();
  k_xcmimage.template sync<LMPHostType>();
  k_displace.template sync<LMPHostType>();
  if (extended) k_eflags.template sync<LMPHostType>();

  return FixRigidBase::pack_exchange(i, buf);

}

template<class DeviceType, class FixRigidBase>
int FixRigidBaseKokkos<DeviceType,FixRigidBase>::unpack_exchange(int nlocal, double *buf)
{

  int result = FixRigidBase::unpack_exchange(nlocal, buf);

  k_body.template modify<LMPHostType>();
  k_bodyown.template modify<LMPHostType>();
  k_bodytag.template modify<LMPHostType>();
  k_atom2body.template modify<LMPHostType>();
  k_xcmimage.template modify<LMPHostType>();
  k_displace.template modify<LMPHostType>();
  if (extended) k_eflags.template modify<LMPHostType>();

  return result;

}

template<class DeviceType, class FixRigidBase>
int FixRigidBaseKokkos<DeviceType,FixRigidBase>::pack_forward_comm(int n, int *list,
                                                       double *buf, int pbc_flag, int *pbc)
{

  k_bodyown.template sync<LMPHostType>();
  k_body.template sync<LMPHostType>();

  return FixRigidBase::pack_forward_comm(n, list, buf, pbc_flag, pbc);

}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::unpack_forward_comm( int n, int first, double *buf)
{

  FixRigidBase::unpack_forward_comm(n, first, buf);

  if (commflag == FULL_BODY) k_bodyown.template sync<LMPHostType>();
  k_body.template sync<LMPHostType>();

}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
int FixRigidBaseKokkos<DeviceType,FixRigidBase>::pack_reverse_comm(int n, int first, double *buf)
{

  k_bodyown.template sync<LMPHostType>();
  if (commflag == FORCE_TORQUE || commflag == VCM_ANGMOM || commflag == XCM_MASS )
    k_body.template sync<LMPHostType>();
  else if (commflag == ITENSOR)
    k_itensor.template sync<LMPHostType>();
  else if (commflag == DOF)
    k_counts.template sync<LMPHostType>();

  return FixRigidBase::pack_reverse_comm(n, first, buf);

}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::unpack_reverse_comm(int n, int *list, double *buf)
{

  FixRigidBase::unpack_reverse_comm(n, list, buf);

  k_bodyown.template modify<LMPHostType>();
  if (commflag == FORCE_TORQUE || commflag == VCM_ANGMOM || commflag == XCM_MASS )
    k_body.template modify<LMPHostType>();
  else if (commflag == ITENSOR)
    k_itensor.template modify<LMPHostType>();
  else if (commflag == DOF)
    k_counts.template modify<LMPHostType>();

}


/* ----------------------------------------------------------------------
  KOKKOS BASE
------------------------------------------------------------------------- */

inline constexpr bool is_kk_fp32 = std::is_same_v<KK_FLOAT, float>;

template<class DeviceType, class FixRigidBase>
int FixRigidBaseKokkos<DeviceType,FixRigidBase>::pack_forward_comm_kokkos(
    int n, DAT::tdual_int_1d k_sendlist, DAT::tdual_double_1d &k_buf,
    int pbc_flag, int *pbc)
{

  // kokkos views

  k_sendlist.template sync<DeviceType>();
  auto l_sendlist = k_sendlist.template view<DeviceType>();

  k_buf.template sync<DeviceType>();
  auto l_buf = k_buf.template view<DeviceType>();

  k_bodyown.template sync<DeviceType>();
  auto l_bodyown = k_bodyown.template view<DeviceType>();

  k_body.template sync<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();

  auto l_comm_me = comm->me;
  int result = 0;
  copymode = 1;
  if (commflag == FULL_BODY) {
    Kokkos::parallel_scan("rigid/small:pack_forward_full_body",
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int &i, int &offset, const bool &final) {
        const int ibody = l_bodyown(l_sendlist(i));
        if (!final) {
          if (ibody >= 0) offset += 47;
          return;
        }
        if (ibody < 0) {
          l_buf[i] = d_ubuf(0).d; // flag
          return;
        }
        l_buf[i] = d_ubuf(1).d; // flag
        int m = n + offset;
        offset += 47;
        const BodyKokkos &bk = l_body(ibody);
        l_buf[m++] = d_ubuf(bk.natoms).d;
        l_buf[m++] = static_cast<double>(bk.mass);
        l_buf[m++] = static_cast<double>(bk.xcm[0]);
        l_buf[m++] = static_cast<double>(bk.xcm[1]);
        l_buf[m++] = static_cast<double>(bk.xcm[2]);
        l_buf[m++] = static_cast<double>(bk.xgc[0]);
        l_buf[m++] = static_cast<double>(bk.xgc[1]);
        l_buf[m++] = static_cast<double>(bk.xgc[2]);
        l_buf[m++] = static_cast<double>(bk.vcm[0]);
        l_buf[m++] = static_cast<double>(bk.vcm[1]);
        l_buf[m++] = static_cast<double>(bk.vcm[2]);
        l_buf[m++] = static_cast<double>(bk.fcm[0]);
        l_buf[m++] = static_cast<double>(bk.fcm[1]);
        l_buf[m++] = static_cast<double>(bk.fcm[2]);
        l_buf[m++] = static_cast<double>(bk.torque[0]);
        l_buf[m++] = static_cast<double>(bk.torque[1]);
        l_buf[m++] = static_cast<double>(bk.torque[2]);
        l_buf[m++] = static_cast<double>(bk.quat[0]);
        l_buf[m++] = static_cast<double>(bk.quat[1]);
        l_buf[m++] = static_cast<double>(bk.quat[2]);
        l_buf[m++] = static_cast<double>(bk.quat[3]);
        l_buf[m++] = static_cast<double>(bk.inertia[0]);
        l_buf[m++] = static_cast<double>(bk.inertia[1]);
        l_buf[m++] = static_cast<double>(bk.inertia[2]);
        l_buf[m++] = static_cast<double>(bk.ex_space[0]);
        l_buf[m++] = static_cast<double>(bk.ex_space[1]);
        l_buf[m++] = static_cast<double>(bk.ex_space[2]);
        l_buf[m++] = static_cast<double>(bk.ey_space[0]);
        l_buf[m++] = static_cast<double>(bk.ey_space[1]);
        l_buf[m++] = static_cast<double>(bk.ey_space[2]);
        l_buf[m++] = static_cast<double>(bk.ez_space[0]);
        l_buf[m++] = static_cast<double>(bk.ez_space[1]);
        l_buf[m++] = static_cast<double>(bk.ez_space[2]);
        l_buf[m++] = static_cast<double>(bk.xgc_body[0]);
        l_buf[m++] = static_cast<double>(bk.xgc_body[1]);
        l_buf[m++] = static_cast<double>(bk.xgc_body[2]);
        l_buf[m++] = static_cast<double>(bk.angmom[0]);
        l_buf[m++] = static_cast<double>(bk.angmom[1]);
        l_buf[m++] = static_cast<double>(bk.angmom[2]);
        l_buf[m++] = static_cast<double>(bk.omega[0]);
        l_buf[m++] = static_cast<double>(bk.omega[1]);
        l_buf[m++] = static_cast<double>(bk.omega[2]);
        l_buf[m++] = static_cast<double>(bk.conjqm[0]);
        l_buf[m++] = static_cast<double>(bk.conjqm[1]);
        l_buf[m++] = static_cast<double>(bk.conjqm[2]);
        l_buf[m++] = static_cast<double>(bk.conjqm[3]);
        l_buf[m++] = d_ubuf(bk.image).d;
      }, result
    );
    return n + result;
  } else if (commflag == INITIAL) {
    Kokkos::parallel_scan("rigid/small:pack_forward_initial",
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(l_sendlist(i));
        if (ibody < 0) return;
        if (final) {
          const BodyKokkos &bk = l_body(ibody);
          l_buf[m++] = static_cast<double>(bk.xcm[0]);
          l_buf[m++] = static_cast<double>(bk.xcm[1]);
          l_buf[m++] = static_cast<double>(bk.xcm[2]);
          l_buf[m++] = static_cast<double>(bk.xgc[0]);
          l_buf[m++] = static_cast<double>(bk.xgc[1]);
          l_buf[m++] = static_cast<double>(bk.xgc[2]);
          l_buf[m++] = static_cast<double>(bk.vcm[0]);
          l_buf[m++] = static_cast<double>(bk.vcm[1]);
          l_buf[m++] = static_cast<double>(bk.vcm[2]);
          l_buf[m++] = static_cast<double>(bk.quat[0]);
          l_buf[m++] = static_cast<double>(bk.quat[1]);
          l_buf[m++] = static_cast<double>(bk.quat[2]);
          l_buf[m++] = static_cast<double>(bk.quat[3]);
          l_buf[m++] = static_cast<double>(bk.omega[0]);
          l_buf[m++] = static_cast<double>(bk.omega[1]);
          l_buf[m++] = static_cast<double>(bk.omega[2]);
          l_buf[m++] = static_cast<double>(bk.ex_space[0]);
          l_buf[m++] = static_cast<double>(bk.ex_space[1]);
          l_buf[m++] = static_cast<double>(bk.ex_space[2]);
          l_buf[m++] = static_cast<double>(bk.ey_space[0]);
          l_buf[m++] = static_cast<double>(bk.ey_space[1]);
          l_buf[m++] = static_cast<double>(bk.ey_space[2]);
          l_buf[m++] = static_cast<double>(bk.ez_space[0]);
          l_buf[m++] = static_cast<double>(bk.ez_space[1]);
          l_buf[m++] = static_cast<double>(bk.ez_space[2]);
          l_buf[m++] = static_cast<double>(bk.conjqm[0]);
          l_buf[m++] = static_cast<double>(bk.conjqm[1]);
          l_buf[m++] = static_cast<double>(bk.conjqm[2]);
          l_buf[m++] = static_cast<double>(bk.conjqm[3]);
        } else m += 29;
      }, result
    );
  } else if (commflag == FINAL) {
    Kokkos::parallel_scan("rigid/small:pack_forward_final",
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(l_sendlist(i));
        if (ibody < 0) return;
        if (final) {
          const BodyKokkos &bk = l_body(ibody);
          l_buf[m++] = static_cast<double>(bk.vcm[0]);
          l_buf[m++] = static_cast<double>(bk.vcm[1]);
          l_buf[m++] = static_cast<double>(bk.vcm[2]);
          l_buf[m++] = static_cast<double>(bk.omega[0]);
          l_buf[m++] = static_cast<double>(bk.omega[1]);
          l_buf[m++] = static_cast<double>(bk.omega[2]);
          l_buf[m++] = static_cast<double>(bk.conjqm[0]);
          l_buf[m++] = static_cast<double>(bk.conjqm[1]);
          l_buf[m++] = static_cast<double>(bk.conjqm[2]);
          l_buf[m++] = static_cast<double>(bk.conjqm[3]);
        } else m += 10;
      }, result
    );
  }
  copymode = 0;
  return result;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::unpack_forward_comm_kokkos(
    int n, int first, DAT::tdual_double_1d &k_buf)
{

  // grow body views if needed
  // receiving n atoms so at most n bodies

  if (commflag == FULL_BODY && nlocal_body + n > nmax_body)
    grow_body(nlocal_body + n);

  // kokkos views

  k_buf.template sync<DeviceType>();
  auto l_buf = k_buf.template view<DeviceType>();

  k_bodyown.template sync<DeviceType>();
  auto l_bodyown = k_bodyown.template view<DeviceType>();

  k_body.template sync<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();

  copymode = 1;
  if (commflag == FULL_BODY) {

    auto l_comm_me = comm->me;
    auto l_nlocal_body = nlocal_body;
    int nbody_recv;

    Kokkos::parallel_scan("rigid/small:unpack_forward_full_body",
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int &i, int &l_nbody_recv, const bool &final) {
        const int flag = d_ubuf(l_buf[i]).i;
        if (!final) {
          if( flag == 1 ) l_nbody_recv++;
          return;
        }
        if( flag == 0 ) {
          l_bodyown(first + i) = -1;
          return;
        }
        const int j = l_nlocal_body + l_nbody_recv;
        int m = n + l_nbody_recv * 47;
        l_body(j).ilocal = first + i;
        l_bodyown(first + i) = j;
        BodyKokkos &bk = l_body(j);
        bk.natoms = d_ubuf(l_buf[m++]).i;
        bk.mass = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.xcm[0] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.xcm[1] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.xcm[2] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.xgc[0] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.xgc[1] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.xgc[2] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.vcm[0] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.vcm[1] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.vcm[2] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.fcm[0] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.fcm[1] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.fcm[2] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.torque[0] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.torque[1] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.torque[2] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.quat[0] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.quat[1] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.quat[2] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.quat[3] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.inertia[0] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.inertia[1] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.inertia[2] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.ex_space[0] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.ex_space[1] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.ex_space[2] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.ey_space[0] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.ey_space[1] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.ey_space[2] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.ez_space[0] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.ez_space[1] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.ez_space[2] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.xgc_body[0] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.xgc_body[1] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.xgc_body[2] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.angmom[0] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.angmom[1] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.angmom[2] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.omega[0] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.omega[1] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.omega[2] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.conjqm[0] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.conjqm[1] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.conjqm[2] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.conjqm[3] = static_cast<KK_FLOAT>(l_buf[m++]);
        bk.image = d_ubuf(l_buf[m++]).i;
        l_nbody_recv++;
      }, nbody_recv
    );

    nghost_body += nbody_recv;
    k_bodyown.modify_device();

  } else if (commflag == INITIAL) {
    Kokkos::parallel_scan("rigid/small:unpack_forward_initial",
      Kokkos::RangePolicy<DeviceType>(first, first+n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(i);
        if (ibody < 0) return;
        if (final) {
          BodyKokkos &bk = l_body(ibody);
          bk.xcm[0] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.xcm[1] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.xcm[2] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.xgc[0] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.xgc[1] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.xgc[2] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.vcm[0] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.vcm[1] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.vcm[2] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.quat[0] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.quat[1] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.quat[2] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.quat[3] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.omega[0] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.omega[1] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.omega[2] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.ex_space[0] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.ex_space[1] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.ex_space[2] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.ey_space[0] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.ey_space[1] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.ey_space[2] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.ez_space[0] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.ez_space[1] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.ez_space[2] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.conjqm[0] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.conjqm[1] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.conjqm[2] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.conjqm[3] = static_cast<KK_FLOAT>(l_buf[m++]);
        } else m += 29;
      }
    );
  } else if (commflag == FINAL) {
    Kokkos::parallel_scan("rigid/small:unpack_forward_final",
      Kokkos::RangePolicy<DeviceType>(first, first+n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(i);
        if (ibody < 0) return;
        if (final) {
          BodyKokkos &bk = l_body(ibody);
          bk.vcm[0] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.vcm[1] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.vcm[2] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.omega[0] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.omega[1] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.omega[2] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.conjqm[0] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.conjqm[1] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.conjqm[2] = static_cast<KK_FLOAT>(l_buf[m++]);
          bk.conjqm[3] = static_cast<KK_FLOAT>(l_buf[m++]);
        } else m += 10;
      }
    );
  }
  copymode = 0;
  k_body.modify_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
int FixRigidBaseKokkos<DeviceType,FixRigidBase>::pack_reverse_comm_kokkos(
    int n, int first, DAT::tdual_double_1d &k_buf)
{

  // kokkos views

  k_buf.template sync<DeviceType>();
  auto l_buf = k_buf.template view<DeviceType>();

  k_bodyown.template sync<DeviceType>();
  auto l_bodyown = k_bodyown.template view<DeviceType>();

  int result = 0;
  copymode = 1;
  if (commflag == FORCE_TORQUE) {

    // kokkos views
    k_body.template sync<DeviceType>();
    auto l_body = k_body.template view<DeviceType>();

    Kokkos::parallel_scan("rigid/small:pack_reverse_force_torque",
      Kokkos::RangePolicy<DeviceType>(first, first+n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(i);
        if (ibody < 0) return;
        if (final) {
          const BodyKokkos &bk = l_body(ibody);
          l_buf[m++] = static_cast<double>(bk.fcm[0]);
          l_buf[m++] = static_cast<double>(bk.fcm[1]);
          l_buf[m++] = static_cast<double>(bk.fcm[2]);
          l_buf[m++] = static_cast<double>(bk.torque[0]);
          l_buf[m++] = static_cast<double>(bk.torque[1]);
          l_buf[m++] = static_cast<double>(bk.torque[2]);
        } else m += 6;
      }, result
    );

  } else if (commflag == VCM_ANGMOM) {

    // kokkos views
    k_body.template sync<DeviceType>();
    auto l_body = k_body.template view<DeviceType>();

    Kokkos::parallel_scan("rigid/small:vcm_angmom",
      Kokkos::RangePolicy<DeviceType>(first, first+n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(i);
        if (ibody < 0) return;
        if (final) {
          const BodyKokkos &bk = l_body(ibody);
          l_buf[m++] = static_cast<double>(bk.vcm[0]);
          l_buf[m++] = static_cast<double>(bk.vcm[1]);
          l_buf[m++] = static_cast<double>(bk.vcm[2]);
          l_buf[m++] = static_cast<double>(bk.angmom[0]);
          l_buf[m++] = static_cast<double>(bk.angmom[1]);
          l_buf[m++] = static_cast<double>(bk.angmom[2]);
        } else m += 6;
      }, result
    );

  } else if (commflag == XCM_MASS) {

    // kokkos views
    k_body.template sync<DeviceType>();
    auto l_body = k_body.template view<DeviceType>();

    Kokkos::parallel_scan("rigid/small:pack_reverse_xcm_mass",
      Kokkos::RangePolicy<DeviceType>(first, first+n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(i);
        if (ibody < 0) return;
        if (final) {
          const BodyKokkos &bk = l_body(ibody);
          l_buf[m++] = static_cast<double>(bk.xcm[0]);
          l_buf[m++] = static_cast<double>(bk.xcm[1]);
          l_buf[m++] = static_cast<double>(bk.xcm[2]);
          l_buf[m++] = static_cast<double>(bk.xgc[0]);
          l_buf[m++] = static_cast<double>(bk.xgc[1]);
          l_buf[m++] = static_cast<double>(bk.xgc[2]);
          l_buf[m++] = static_cast<double>(bk.mass);
          l_buf[m++] = static_cast<double>(bk.natoms);
        } else m += 8;
      }, result
    );

  } else if (commflag == ITENSOR) {

    // kokkos views
    k_itensor.template sync<DeviceType>();
    auto l_itensor = k_itensor.template view<DeviceType>();

    Kokkos::parallel_scan("rigid/small:pack_reverse_itensor",
      Kokkos::RangePolicy<DeviceType>(first, first+n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(i);
        if (ibody < 0) return;
        if (final) {
          l_buf[m++] = static_cast<double>(l_itensor(ibody,0));
          l_buf[m++] = static_cast<double>(l_itensor(ibody,1));
          l_buf[m++] = static_cast<double>(l_itensor(ibody,2));
          l_buf[m++] = static_cast<double>(l_itensor(ibody,3));
          l_buf[m++] = static_cast<double>(l_itensor(ibody,4));
          l_buf[m++] = static_cast<double>(l_itensor(ibody,5));
        } else m += 6;
      }, result
    );

  } else if (commflag == DOF) {

    // kokkos views
    k_counts.template sync<DeviceType>();
    auto l_counts = k_counts.template view<DeviceType>();

    Kokkos::parallel_scan("rigid/small:pack_reverse_dof",
      Kokkos::RangePolicy<DeviceType>(first, first+n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(i);
        if (ibody < 0) return;
        if (final) {
          l_buf[m++] = static_cast<double>(l_counts(ibody,0));
          l_buf[m++] = static_cast<double>(l_counts(ibody,1));
          l_buf[m++] = static_cast<double>(l_counts(ibody,2));
        } else m += 3;
      }, result
    );

  }

  copymode = 0;
  k_buf.modify_device();
  return result;

}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::unpack_reverse_comm_kokkos(
    int n, DAT::tdual_int_1d k_sendlist, DAT::tdual_double_1d &k_buf)
{

  // kokkos views

  k_sendlist.template sync<DeviceType>();
  auto l_sendlist = k_sendlist.template view<DeviceType>();

  k_buf.template sync<DeviceType>();
  auto l_buf = k_buf.template view<DeviceType>();

  k_bodyown.template sync<DeviceType>();
  auto l_bodyown = k_bodyown.template view<DeviceType>();

  copymode = 1;
  if (commflag == FORCE_TORQUE) {

    // kokkos views
    k_body.template sync<DeviceType>();
    auto l_body = k_body.template view<DeviceType>();

    Kokkos::parallel_scan("rigid/small:unpack_reverse_force_torque",
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(l_sendlist(i));
        if (ibody < 0) return;
        if (final) {
          BodyKokkos &bk = l_body(ibody);
          Kokkos::atomic_add(&bk.fcm[0], static_cast<KK_FLOAT>(l_buf[m++]));
          Kokkos::atomic_add(&bk.fcm[1], static_cast<KK_FLOAT>(l_buf[m++]));
          Kokkos::atomic_add(&bk.fcm[2], static_cast<KK_FLOAT>(l_buf[m++]));
          Kokkos::atomic_add(&bk.torque[0], static_cast<KK_FLOAT>(l_buf[m++]));
          Kokkos::atomic_add(&bk.torque[1], static_cast<KK_FLOAT>(l_buf[m++]));
          Kokkos::atomic_add(&bk.torque[2], static_cast<KK_FLOAT>(l_buf[m++]));
        } else m += 6;
      }
    );
    k_body.template modify<DeviceType>();

  } else if (commflag == VCM_ANGMOM) {

    // kokkos views
    k_body.template sync<DeviceType>();
    auto l_body = k_body.template view<DeviceType>();

    Kokkos::parallel_scan("rigid/small:unpack_reverse_vcm_angmom",
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(l_sendlist(i));
        if (ibody < 0) return;
        if (!final) {
          m += 6;
          return;
        }
        BodyKokkos &bk = l_body(ibody);
        Kokkos::atomic_add(&bk.vcm[0], static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&bk.vcm[1], static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&bk.vcm[2], static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&bk.angmom[0], static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&bk.angmom[1], static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&bk.angmom[2], static_cast<KK_FLOAT>(l_buf[m++]));
      }
    );
    k_body.template modify<DeviceType>();

  } else if (commflag == XCM_MASS) {

    // kokkos views
    k_body.template sync<DeviceType>();
    auto l_body = k_body.template view<DeviceType>();

    Kokkos::parallel_scan("rigid/small:unpack_reverse_xcm_mass",
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(l_sendlist(i));
        if (ibody < 0) return;
        if (final) {
          BodyKokkos &bk = l_body(ibody);
          Kokkos::atomic_add(&bk.xcm[0], static_cast<KK_FLOAT>(l_buf[m++]));
          Kokkos::atomic_add(&bk.xcm[1], static_cast<KK_FLOAT>(l_buf[m++]));
          Kokkos::atomic_add(&bk.xcm[2], static_cast<KK_FLOAT>(l_buf[m++]));
          Kokkos::atomic_add(&bk.xgc[0], static_cast<KK_FLOAT>(l_buf[m++]));
          Kokkos::atomic_add(&bk.xgc[1], static_cast<KK_FLOAT>(l_buf[m++]));
          Kokkos::atomic_add(&bk.xgc[2], static_cast<KK_FLOAT>(l_buf[m++]));
          Kokkos::atomic_add(&bk.mass, static_cast<KK_FLOAT>(l_buf[m++]));
          Kokkos::atomic_add(&bk.natoms, static_cast<int>(l_buf[m++]));
        } else m += 8;
      }
    );
    k_body.template modify<DeviceType>();

  } else if (commflag == ITENSOR) {

    // kokkos views
    k_itensor.template sync<DeviceType>();
    auto l_itensor = k_itensor.template view<DeviceType>();

    Kokkos::parallel_scan("rigid/small:unpack_reverse_itensor",
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(l_sendlist(i));
        if (ibody < 0) return;
        if (final) {
          Kokkos::atomic_add(&l_itensor(ibody,0), static_cast<KK_FLOAT>(l_buf[m++]));
          Kokkos::atomic_add(&l_itensor(ibody,1), static_cast<KK_FLOAT>(l_buf[m++]));
          Kokkos::atomic_add(&l_itensor(ibody,2), static_cast<KK_FLOAT>(l_buf[m++]));
          Kokkos::atomic_add(&l_itensor(ibody,3), static_cast<KK_FLOAT>(l_buf[m++]));
          Kokkos::atomic_add(&l_itensor(ibody,4), static_cast<KK_FLOAT>(l_buf[m++]));
          Kokkos::atomic_add(&l_itensor(ibody,5), static_cast<KK_FLOAT>(l_buf[m++]));
        } else m += 6;
      }
    );
    k_itensor.template modify<DeviceType>();

  } else if (commflag == DOF) {

    // kokkos views
    k_counts.template sync<DeviceType>();
    auto l_counts = k_counts.template view<DeviceType>();

    Kokkos::parallel_scan("rigid/small:unpack_reverse_dof",
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(l_sendlist(i));
        if (ibody < 0) return;
        if (final) {
          Kokkos::atomic_add(&l_counts(ibody,0), static_cast<int>(l_buf[m++]));
          Kokkos::atomic_add(&l_counts(ibody,1), static_cast<int>(l_buf[m++]));
          Kokkos::atomic_add(&l_counts(ibody,2), static_cast<int>(l_buf[m++]));
        } else m += 3;
      }
    );
    k_counts.template modify<DeviceType>();

  }

  copymode = 0;

}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
int FixRigidBaseKokkos<DeviceType,FixRigidBase>::pack_exchange_kokkos(
    const int &nsend, DAT::tdual_double_2d_lr &k_buf,
    DAT::tdual_int_1d k_sendlist, DAT::tdual_int_1d k_copylist,
    ExecutionSpace space)
{
  int result = 0;
  auto l_buf = typename AT::t_double_1d_um(
      k_buf.template view<DeviceType>().data(),
      k_buf.extent(0)*k_buf.extent(1)
  );
  auto l_sendlist = k_sendlist.view<DeviceType>();
  auto l_copylist = k_copylist.view<DeviceType>();
  auto l_bodytag = k_bodytag.template view<DeviceType>();
  auto l_xcmimage = k_xcmimage.template view<DeviceType>();
  auto l_displace = k_displace.template view<DeviceType>();
  auto l_bodyown = k_bodyown.template view<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();
  auto l_vatom = k_vatom.view<DeviceType>();
  auto l_vflag_atom = vflag_atom;
  int l_bodysize_kk = sizeof(BodyKokkos)/sizeof(double);
  if (l_bodysize_kk * sizeof(double) != sizeof(BodyKokkos)) l_bodysize_kk++;

  copymode = 1;
  Kokkos::parallel_scan("rigid/small:pack_exchange",
    Kokkos::RangePolicy<DeviceType>(0, nsend),
    KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
      const int local_idx = l_sendlist(i);
        if (!final) {
          m += 5;
          if (l_bodytag(local_idx)) {
            if (l_vflag_atom) m += 6;
            m += 1;
            if (l_bodyown(local_idx) >= 0) m += l_bodysize_kk;
        }
        return;
      }
      l_buf(m++) = d_ubuf(l_bodytag(local_idx)).d;
      l_buf(m++) = d_ubuf(l_xcmimage(local_idx)).d;
      l_buf(m++) = l_displace(local_idx, 0);
      l_buf(m++) = l_displace(local_idx, 1);
      l_buf(m++) = l_displace(local_idx, 2);

      /*
      // extended attribute info
      if (l_extended) {
        l_buf(m++) = l_eflags(local_idx);
        for (int j = 0; j < l_orientflag; j++)
        l_buf[m++] = l_orient(local_idx,j);
        if (l_dorientflag) {
          l_buf(m++) = l_dorient(local_idx, 0);
          l_buf(m++) = l_dorient(local_idx, 1);
          l_buf(m++) = l_dorient(local_idx, 2);
        }
      }
      */

      // atom not in a rigid body
      if (!l_bodytag(local_idx)) return;

      // must also pack vatom if per-atom virial calculated on this timestep
      // since vatom is calculated before and after atom migration
      if (l_vflag_atom) {
        for (int k = 0; k < 6; k++) l_buf(m++) = l_vatom(local_idx, k);
      }

      // atom does not own its rigid body
      if (l_bodyown(local_idx) < 0) {
        l_buf(m++) = 0;
        return;
      }

      // body info for atom that owns a rigid body
      l_buf(m++) = 1;
      memcpy(&l_buf(m),&l_body(l_bodyown(local_idx)),sizeof(BodyKokkos));
      m += l_bodysize_kk;
    }, result
  );
  copymode = 0;
  k_buf.modify_device();
  return result;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::unpack_exchange_kokkos(
    DAT::tdual_double_2d_lr &k_buf, DAT::tdual_int_1d &k_indices,
    int nrecv, int nrecv1, int nextrarecv1, ExecutionSpace space)
{
  auto l_buf = typename AT::t_double_1d_um(
      k_buf.template view<DeviceType>().data(),
      k_buf.extent(0)*k_buf.extent(1)
  );
  auto l_bodytag = k_bodytag.template view<DeviceType>();
  auto l_indices = k_indices.template view<DeviceType>();
  auto l_xcmimage = k_xcmimage.template view<DeviceType>();
  auto l_displace = k_displace.template view<DeviceType>();
  auto l_bodyown = k_bodyown.template view<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();
  auto l_vatom = k_vatom.view<DeviceType>();
  auto l_vflag_atom = vflag_atom;
  int l_bodysize_kk = sizeof(BodyKokkos)/sizeof(double);
  if (l_bodysize_kk * sizeof(double) != sizeof(BodyKokkos)) l_bodysize_kk++;
  auto l_nlocal_body = nlocal_body;
  int nbody_recv = 0;

  copymode = 1;
  Kokkos::parallel_scan("rigid/small:unpack_exchange",
    Kokkos::RangePolicy<DeviceType>(0, nrecv),
    KOKKOS_LAMBDA(const int &i, int &l_nbody_recv, const bool &final) {

      const int local_idx = l_indices(i);
      if (local_idx < 0) return;
      int m = static_cast<int>(ubuf(l_buf(i)).i);
      if (i >= nrecv1)
        m = nextrarecv1 + static_cast<int>(ubuf(l_buf(nextrarecv1 + i - nrecv1)).i);

      tagint tag = static_cast<tagint>(ubuf(l_buf(m++)).i);
      if (!final) {
        if (tag) {
          // Look ahead in the buffer to check if this atom owns the body
          // Skip xcmimage (1), displace (3), and vatom (6 if active)
          const int bodyown_offset = m + 4 + (l_vflag_atom ? 6 : 0);
          const int bodyown_val = static_cast<int>(l_buf(bodyown_offset));
          if (bodyown_val != 0) l_nbody_recv++;
        }
        return;
      }

      l_bodytag(local_idx) = tag;
      l_xcmimage(local_idx) = static_cast<imageint>(ubuf(l_buf(m++)).i);
      l_displace(local_idx, 0) = l_buf(m++);
      l_displace(local_idx, 1) = l_buf(m++);
      l_displace(local_idx, 2) = l_buf(m++);

      /*
      // extended attribute info
      if (l_extended) {
        l_eflags(nlocal) = static_cast<int>(l_buf(m++));
        for (int j = 0; j < l_orientflag; j++)
          l_orient[nlocal][j] = l_buf(m++);
        if (l_dorientflag) {
          l_dorient[nlocal][0] = l_buf(m++);
          l_dorient[nlocal][1] = l_buf(m++);
          l_dorient[nlocal][2] = l_buf(m++);
        }
      }
      */

      // atom not in a rigid body
      if (!l_bodytag(local_idx)) {
        l_bodyown(local_idx) = -1;
        return;
      }

      // must also unpack vatom if per-atom virial calculated on this timestep
      // since vatom is calculated before and after atom migration

      if (l_vflag_atom) {
        for (int k = 0; k < 6; k++) l_vatom(local_idx, k) = l_buf(m++);
      }

      // atom does not own its rigid body
      l_bodyown(local_idx) = static_cast<int>(l_buf(m++));
      if (l_bodyown(local_idx) == 0) {
        l_bodyown(local_idx) = -1;
        return;
      }

      // body info for atom that owns a rigid body
      const int l_nbody_idx = l_nlocal_body + l_nbody_recv;
      memcpy(&l_body(l_nbody_idx), &l_buf(m), sizeof(BodyKokkos));
      m += l_bodysize_kk;
      l_body(l_nbody_idx).ilocal = local_idx;
      l_bodyown(local_idx) = l_nbody_idx;
      l_nbody_recv++;
    }, nbody_recv
  );
  nlocal_body += nbody_recv;
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  // always sort on the device

  k_body.sync_device();
  k_bodytag.sync_device();
  k_xcmimage.sync_device();
  k_displace.sync_device();
  k_bodyown.sync_device();
  if (extended) k_eflags.sync_device();

  Sorter.sort(LMPDeviceType(), k_body.view_device());
  Sorter.sort(LMPDeviceType(), k_bodytag.view_device());
  Sorter.sort(LMPDeviceType(), k_xcmimage.view_device());
  Sorter.sort(LMPDeviceType(), k_displace.view_device());
  Sorter.sort(LMPDeviceType(), k_bodyown.view_device());
  if (extended) Sorter.sort(LMPDeviceType(), k_eflags.view_device());

  k_body.modify_device();
  k_bodytag.modify_device();
  k_xcmimage.modify_device();
  k_displace.modify_device();
  k_bodyown.modify_device();
  if (extended) k_eflags.modify_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::modify_host()
{
  k_body.modify_host();
  k_bodyown.modify_host();
  k_bodytag.modify_host();
  k_atom2body.modify_host();
  k_xcmimage.modify_host();
  k_displace.modify_host();
  if (extended) k_eflags.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::sync_host()
{
  k_body.sync_host();
  k_bodyown.sync_host();
  k_bodytag.sync_host();
  k_atom2body.sync_host();
  k_xcmimage.sync_host();
  k_displace.sync_host();
  if (extended) k_eflags.sync_host();
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {

  // fix rigid/small/kk
  template class FixRigidBaseKokkos<LMPDeviceType, FixRigidSmall>;

  // fix rigid/(nve|nvt|nph|npt)/small/kk
  template class FixRigidBaseKokkos<LMPDeviceType, FixRigidNHSmall>;

#ifdef LMP_KOKKOS_GPU

    // fix rigid/small/kk/host
    template class FixRigidBaseKokkos<LMPHostType, FixRigidSmall>;

    // fix rigid/(nve|nvt|nph|npt)/small/kk/host
    template class FixRigidBaseKokkos<LMPHostType, FixRigidNHSmall>;

  #endif // LMP_KOKKOS_GPU

}

