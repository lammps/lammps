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
#include "domain_kokkos.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "math_extra_kokkos.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "molecule.h"
#include "rigid_const.h"
#include "update.h"

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
FixRigidBaseKokkos<DeviceType,FixRigidBase>::FixRigidBaseKokkos(Atom* atom, Domain* domain) :
    KokkosBase(),
    atomKK(static_cast<AtomKokkos*>(atom)), domainKK(static_cast<DomainKokkos*>(domain))
{
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
FixRigidBaseKokkos<DeviceType,FixRigidBase>::~FixRigidBaseKokkos()
{
  if (base_ptr->copymode) return;
  base_ptr->memoryKK->destroy_kokkos(k_bodyown, base_ptr->bodyown);
  base_ptr->memoryKK->destroy_kokkos(k_bodytag, base_ptr->bodytag);
  base_ptr->memoryKK->destroy_kokkos(k_atom2body, base_ptr->atom2body);
  base_ptr->memoryKK->destroy_kokkos(k_xcmimage, base_ptr->xcmimage);
  if (base_ptr->displace) {
    base_ptr->memory->sfree(base_ptr->displace);
    base_ptr->displace = nullptr;
  }
  base_ptr->memoryKK->destroy_kokkos(k_displace);
  if (base_ptr->extended) base_ptr->memoryKK->destroy_kokkos(k_eflags, base_ptr->eflags);
  base_ptr->body = nullptr;
  base_ptr->bodyown = nullptr;
  base_ptr->bodytag = nullptr;
  base_ptr->atom2body = nullptr;
  base_ptr->xcmimage = nullptr;
  base_ptr->eflags = nullptr;
#ifdef LMP_KOKKOS_DEBUG_RNG
  rand_pool.destroy();
#endif
}

/* ----------------------------------------------------------------------
   FIX METHODS
------------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::post_constructor_base()
{
  base_ptr = base(); // keep the pointer for destructor later
  const int nmax = atomKK->nmax;
  const int nlocal = atomKK->nlocal;
  // save bodytag and bodyown filled by the base constructor's create_bodies()
  int *old_bodyown = base()->bodyown;
  base()->bodyown = nullptr;
  tagint *old_bodytag = base()->bodytag;
  base()->bodytag = nullptr;
  base()->memoryKK->destroy_kokkos(k_bodyown, base()->bodyown);
  base()->memoryKK->destroy_kokkos(k_bodytag, base()->bodytag);
  base()->memoryKK->destroy_kokkos(k_atom2body, base()->atom2body);
  base()->memoryKK->destroy_kokkos(k_xcmimage, base()->xcmimage);
  {
    double **old_displace = base()->displace;
    std::vector<double> displace_backup((bigint) nmax * 3);
    for (int i = 0; i < nmax; i++)
      for (int j = 0; j < 3; j++) displace_backup[(bigint) i * 3 + j] = old_displace[i][j];
    base()->memory->destroy(base()->displace);
    k_displace = TransformView<KK_FLOAT**, double**, Kokkos::LayoutRight, DeviceType>("rigid/small:displace", nmax, 3);
    double *dh = const_cast<double *>(k_displace.view_host().data());
    memcpy(dh, displace_backup.data(), displace_backup.size() * sizeof(double));
    const bigint nbytes = ((bigint) sizeof(double *)) * nmax;
    base()->displace = (double **) base()->memory->smalloc(nbytes, "rigid/small:displace");
    for (int i = 0; i < nmax; i++) base()->displace[i] = &k_displace.view_host()(i, 0);
    k_displace.modify_host();
    k_displace.sync_device();
  }
  base()->memoryKK->create_kokkos(k_bodyown, base()->bodyown, nmax, "rigid/small:bodyown");
  base()->memoryKK->create_kokkos(k_bodytag, base()->bodytag, nmax, "rigid/small:bodytag");
  base()->memoryKK->create_kokkos(k_atom2body, base()->atom2body, nmax, "rigid/small:atom2body");
  base()->memoryKK->create_kokkos(k_xcmimage, base()->xcmimage, nmax, "rigid/small:xcmimage");
  if (nlocal > 0) {
    memcpy(base()->bodyown, old_bodyown, nlocal * sizeof(int));
    memcpy(base()->bodytag, old_bodytag, nlocal * sizeof(tagint));
    k_bodyown.modify_host();
    k_bodytag.modify_host();
  }
  base()->memory->sfree(old_bodyown);
  base()->memory->sfree(old_bodytag);
  if (base()->extended) {
    base()->memoryKK->destroy_kokkos(k_eflags, base()->eflags);
    base()->memoryKK->create_kokkos(k_eflags, base()->eflags, nmax, "rigid/small:eflags");
  }
  k_body = TransformView<BodyKokkos*, Body*, Kokkos::LayoutRight, DeviceType>("rigid/small:body", base()->nmax_body);
  if (base()->nmax_body > 0 && base()->body != nullptr) {
    memcpy(k_body.view_host().data(), base()->body, (bigint) base()->nmax_body * sizeof(Body));
    base()->memory->sfree(base()->body);
    base()->body = k_body.view_host().data();
    k_body.modify_host();
    k_body.sync_device();
  }

#ifdef LMP_KOKKOS_DEBUG_RNG
  this->rand_pool = Kokkos::Random_XorShift64_Pool<DeviceType>(base()->seed + base()->comm->me, lmp);
#else
  this->rand_pool = Kokkos::Random_XorShift64_Pool<DeviceType>(base()->seed + base()->comm->me);
#endif

}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::init_base()
{
  base()->FixRigidBase::init();
  atomKK->k_mass.modify_host();
  atomKK->k_mass.template sync<DeviceType>();
#ifdef LMP_KOKKOS_DEBUG_RNG
  rand_pool.init(random,seed + comm->me);
#endif
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::setup_pre_neighbor_base()
{
  if (base()->reinitflag || !base()->setupflag)
    setup_bodies_static_base();
  else
    pre_neighbor_base();
  if ((base()->reinitflag || !base()->setupflag) && !base()->inpfile)
    setup_bodies_dynamic_base();
  base()->setupflag = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::pre_neighbor_base()
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

  utils::logmesg(base()->lmp, "*** pre_neighbor() {} {} {}\n{} {} {}\n{} {} {}\n{} {} {}\n", l_xperiodic, l_yperiodic, l_zperiodic, l_lo0, l_lo1, l_lo2, l_hi0, l_hi1, l_hi2,l_prd0, l_prd1, l_prd2);

  base()->copymode = 1;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<DeviceType>(0, base()->nlocal_body),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      imageint idim, otherdims;
      Kokkos::printf("*** ibody %i bk.xcm %f %f %f l_prd0 %f bk.image %i\n", ibody, bk.xcm[0], bk.xcm[1], bk.xcm[2], l_prd0, bk.image);

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
  base()->copymode = 0;
  k_body.modify_device();

  base()->nghost_body = 0;
  base()->commflag = FULL_BODY;
  base()->comm->forward_comm(fix_base());

  reset_atom2body_base();
  image_shift_base();
}


/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::setup_base(int vflag)
{

  const int nlocal = atomKK->nlocal;

  // error if maxextent > comm->cutghost
  // NOTE: could just warn if an override flag set
  // NOTE: this could fail for comm multi mode if user sets a wrong cutoff
  //       for atom types in rigid bodies - need a more careful test
  // must check here, not in init, b/c neigh/comm values set after fix init

  double cutghost = MAX(base()->neighbor->cutneighmax, base()->comm->cutghostuser);
  if (base()->maxextent > cutghost)
    base()->error->all(FLERR, Error::NOLASTLINE,
               "Rigid body extent {} > ghost atom cutoff - use comm_modify cutoff", base()->maxextent);

  if (base()->langflag && (base()->nlocal_body > base()->maxlang)) {
    base()->memoryKK->destroy_kokkos(k_langextra, base()->langextra);
    base()->maxlang = nbody_total();
    base()->memoryKK->create_kokkos(k_langextra, base()->langextra, 6, "rigid/small:langextra");
  }

  compute_forces_and_torques_base();
  // enforce 2d body forces and torques
  if (domainKK->dimension == 2) enforce2d_base();

  // virial setup before call to set_v
  base()->v_init(vflag);

  // compute and forward communicate vcm and omega of all bodies
  base()->copymode = 1;
  k_body.sync_device();
  auto l_body = k_body.template view<DeviceType>();
  Kokkos::parallel_for(
    Kokkos::RangePolicy<DeviceType>(0, base()->nlocal_body),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      angmom_to_omega(bk.angmom, bk.ex_space, bk.ey_space,
                                       bk.ez_space, bk.inertia, bk.omega);
    }
  );
  k_body.modify_device();
  base()->copymode = 0;

  base()->commflag = FINAL;
  base()->comm->forward_comm(fix_base(),10);

  // set velocity/rotation of atoms in rigid bodues
  if (base()->evflag) set_v_base<true>();
  else set_v_base<false>();

  // guesstimate virial as 2x the set_v contribution
  if (base()->vflag_global) {
    for (int n = 0; n < 6; n++) base()->virial[n] *= 2.0;
  }
  if (base()->vflag_atom) {
    for (int i = 0; i < nlocal; i++) {
      for (int n = 0; n < 6; n++) base()->vatom[i][n] *= 2.0;
    }
  }

  if constexpr(is_nh) {

  // FIXME
  base()->compute_dof();

  base()->copymode = 1;
  k_body.sync_device();
  auto l_tstat_flag = nh_base()->tstat_flag;
  auto l_pstat_flag = nh_base()->pstat_flag;
  KK_ACC_FLOAT ke[2], keall[2];
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<DeviceType>(0, nh_base()->nlocal_body),
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
  base()->copymode = 0;
  k_body.modify_device();
  if (l_tstat_flag || l_pstat_flag) {
    MPI_Allreduce(ke, keall, 2, MPI_KK_ACC_FLOAT, MPI_SUM, base()->world);
    nh_base()->akin_t = keall[0];
    nh_base()->akin_r = keall[1];
  }
  auto temperature = nh_base()->temperature;
  if (nh_base()->tstat_flag) nh_base()->compute_temp_target();
  else if (nh_base()->pstat_flag) {
    auto pressure = nh_base()->pressure;
    nh_base()->t0 = temperature->compute_scalar();
    if (nh_base()->t0 == 0.0) {
      if (strcmp(base()->update->unit_style, "lj") == 0) nh_base()->t0 = 1.0;
      else nh_base()->t0 = 300.0;
    }
    nh_base()->t_target = nh_base()->t0;
    nh_base()->compute_press_target();
    if (nh_base()->pstyle == ISO) {
      temperature->compute_scalar();
      pressure->compute_scalar();
    } else {
      temperature->compute_vector();
      pressure->compute_vector();
    }
    nh_base()->couple();
    pressure->addstep(base()->update->ntimestep+1);
  }
  double t_mass, tb_mass;
  const double kt = nh_base()->boltz * nh_base()->t_target;
  if (nh_base()->tstat_flag) {
    t_mass = kt / (nh_base()->t_freq * nh_base()->t_freq);
    nh_base()->q_t[0] = nh_base()->nf_t * t_mass;
    nh_base()->q_r[0] = nh_base()->nf_r * t_mass;
    for (int i = 1; i < nh_base()->t_chain; i++)
      nh_base()->q_t[i] = nh_base()->q_r[i] = t_mass;
    for (int i = 1; i < nh_base()->t_chain; i++) {
      nh_base()->f_eta_t[i] = (nh_base()->q_t[i-1] * nh_base()->eta_dot_t[i-1] * nh_base()->eta_dot_t[i-1] - kt)/nh_base()->q_t[i];
      nh_base()->f_eta_r[i] = (nh_base()->q_r[i-1] * nh_base()->eta_dot_r[i-1] * nh_base()->eta_dot_r[i-1] - kt)/nh_base()->q_r[i];
    }
  }
  const int dimension = domainKK->dimension;
  if (nh_base()->pstat_flag) {
    for (int i = 0; i < 3; i++)
      if (nh_base()->p_flag[i]) {
        const auto p_freq_i_sq = (nh_base()->p_freq[i]) * (nh_base()->p_freq[i]);
        nh_base()->epsilon_mass[i] = (nh_base()->g_f + dimension) * kt / p_freq_i_sq;
        nh_base()->epsilon[i] = log(nh_base()->vol0)/dimension;
      }
    tb_mass = kt / (nh_base()->p_freq_max * nh_base()->p_freq_max);
    nh_base()->q_b[0] = dimension * dimension * tb_mass;
    for (int i = 1; i < nh_base()->p_chain; i++) {
      nh_base()->q_b[i] = tb_mass;
      const auto eta_dot_sq = (nh_base()->eta_dot_b[i-1]) * (nh_base()->eta_dot_b[i-1]);
      nh_base()->f_eta_b[i] = (nh_base()->q_b[i] * eta_dot_sq - kt)/nh_base()->q_b[i];
    }
  }
  if (nh_base()->tstat_flag || nh_base()->pstat_flag) {
    for (int i = 0; i < nh_base()->t_order; i++) {
      nh_base()->wdti1[i] = nh_base()->w[i] * nh_base()->dtv / nh_base()->t_iter;
      nh_base()->wdti2[i] = nh_base()->wdti1[i] / 2.0;
      nh_base()->wdti4[i] = nh_base()->wdti1[i] / 4.0;
    }
  }
  if (nh_base()->pstat_flag) {
    nh_base()->compute_press_target();
    nh_base()->nh_epsilon_dot();
  }

  }
}

/* ---------------------------------------------------------------------- */


template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::initial_integrate_base(int vflag)
{

  auto lambda = [&]<bool NH, bool TSTAT, bool PSTAT>() {

    k_body.template sync<DeviceType>();
    auto l_body = k_body.template view<DeviceType>();

    KK_FLOAT l_dtf = static_cast<KK_FLOAT>(base()->dtf);
    KK_FLOAT l_dtf2 = l_dtf * KK_FLOAT(2.0);
    KK_FLOAT l_dtv = static_cast<KK_FLOAT>(base()->dtv);
    KK_FLOAT l_dtq = static_cast<KK_FLOAT>(base()->dtq);

    KK_FLOAT l_scale_t0, l_scale_t1, l_scale_t2;
    KK_FLOAT l_scale_v0, l_scale_v1, l_scale_v2, l_scale_r;

    if constexpr (NH) {
      KK_FLOAT tmp;
      l_scale_t0 = l_scale_t1 = l_scale_t2 = KK_FLOAT(1.0);
      l_scale_v0 = l_scale_v1 = l_scale_v2 = KK_FLOAT(1.0);
      l_scale_r = KK_FLOAT(1.0);
      if constexpr(TSTAT) {
        tmp = exp(-l_dtq * static_cast<KK_FLOAT>(nh_base()->eta_dot_t[0]));
        l_scale_t0 = l_scale_t1 = l_scale_t2 = tmp;
        tmp = exp(-l_dtq * static_cast<KK_FLOAT>(nh_base()->eta_dot_r[0]));
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
    return fma(y, fma(y, fma(y, fma(y, c4, c3), c2), c1), 1.0);
  };

        KK_FLOAT mtk_term2 = static_cast<KK_FLOAT>(nh_base()->mtk_term2);
        l_scale_t0 *= exp(-l_dtq * (static_cast<KK_FLOAT>(nh_base()->epsilon_dot[0]) + mtk_term2));
        l_scale_t1 *= exp(-l_dtq * (static_cast<KK_FLOAT>(nh_base()->epsilon_dot[1]) + mtk_term2));
        l_scale_t2 *= exp(-l_dtq * (static_cast<KK_FLOAT>(nh_base()->epsilon_dot[2]) + mtk_term2));
        l_scale_r *= exp(-l_dtq * (static_cast<KK_FLOAT>(nh_base()->pdim) * mtk_term2));
        tmp = l_dtq * static_cast<KK_FLOAT>(nh_base()->epsilon_dot[0]);
        l_scale_v0 = l_dtv * exp(tmp) * maclaurin(tmp);
        tmp = l_dtq * static_cast<KK_FLOAT>(nh_base()->epsilon_dot[1]);
        l_scale_v1 = l_dtv * exp(tmp) * maclaurin(tmp);
        tmp = l_dtq * static_cast<KK_FLOAT>(nh_base()->epsilon_dot[2]);
        l_scale_v2 = l_dtv * exp(tmp) * maclaurin(tmp);
      }
    }

    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType>(0, base()->nlocal_body),
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
          richardson(bk.quat,bk.angmom,bk.omega,bk.inertia,l_dtq);
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
          bk.angmom[0] *= 0.5;
          bk.angmom[1] *= 0.5;
          bk.angmom[2] *= 0.5;
          angmom_to_omega(bk.angmom, bk.ex_space, bk.ey_space,
                                     bk.ez_space, bk.inertia, bk.omega);
        }
      }
    );
    k_body.modify_device();
  };

  base()->copymode = 1;
  if constexpr (is_nh) {
    const bool is_t = base()->tstat_flag, is_p = base()->pstat_flag;
    if      ( is_t &&  is_p) lambda.template operator()<true,true,true>();   // NPT
    else if ( is_t && !is_p) lambda.template operator()<true,true,false>();  // NVT
    else if (!is_t &&  is_p) lambda.template operator()<true,false,true>();  // NPH
    else                     lambda.template operator()<true,false,false>(); // NVE
  } else                     lambda.template operator()<false,false,false>();
  base()->copymode = 0;

  base()->commflag = INITIAL;
  base()->comm->forward_comm(fix_base(), 29);

  // virial setup
  base()->v_init(vflag);
  if (base()->vflag_atom) {
    if (atomKK->nmax > (int) d_vatom.extent(0)) {
      base()->memoryKK->destroy_kokkos(k_vatom, base()->vatom);
      base()->memoryKK->create_kokkos(k_vatom, base()->vatom, atomKK->nmax, 6, "fix:vatom");
      d_vatom = k_vatom.template view<DeviceType>();
    } else {
      k_vatom.template sync<DeviceType>();
    }
    Kokkos::deep_copy(d_vatom, 0.0);
  }

  if (base()->evflag) set_xv_base<true>();
  else set_xv_base<false>();

  if (base()->vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  if (base()->extended) {
    // not implemented
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::post_force_base()
{
  if (base()->langflag) apply_langevin_thermostat_base();
  if (base()->earlyflag) compute_forces_and_torques_base();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::final_integrate_base()
{

  if (!base()->earlyflag) compute_forces_and_torques_base();
  if (domainKK->dimension == 2) enforce2d_base();

  auto l_body = k_body.template view<DeviceType>();

  auto lambda = [&]<bool NH, bool TSTAT, bool PSTAT>() {
    KK_FLOAT l_dtf = static_cast<KK_FLOAT>(base()->dtf);
    KK_FLOAT l_dtf2 = l_dtf * KK_FLOAT(2.0);
    KK_FLOAT l_scale_t0, l_scale_t1, l_scale_t2, l_scale_r;
    if constexpr (NH) {
      KK_FLOAT l_dtq = static_cast<KK_FLOAT>(base()->dtq);
      l_scale_t0 = l_scale_t1 = l_scale_t2 = l_scale_r = KK_FLOAT(1.0);
      if constexpr (TSTAT) {
        const KK_FLOAT tmp = exp(-1.0 * l_dtq * static_cast<KK_FLOAT>(nh_base()->eta_dot_t[0]));
        l_scale_t0 = l_scale_t1 = l_scale_t2 = tmp;
        l_scale_r = exp(-1.0 * l_dtq * static_cast<KK_FLOAT>(nh_base()->eta_dot_r[0]));
      }
      if constexpr (PSTAT) {
        KK_FLOAT mtk_term2 = static_cast<KK_FLOAT>(nh_base()->mtk_term2);
        l_scale_t0 *= exp(-l_dtq * (static_cast<KK_FLOAT>(nh_base()->epsilon_dot[0]) + mtk_term2));
        l_scale_t1 *= exp(-l_dtq * (static_cast<KK_FLOAT>(nh_base()->epsilon_dot[1]) + mtk_term2));
        l_scale_t2 *= exp(-l_dtq * (static_cast<KK_FLOAT>(nh_base()->epsilon_dot[2]) + mtk_term2));
        l_scale_r *= exp(-l_dtq * (static_cast<KK_FLOAT>(nh_base()->pdim) * mtk_term2));
      }
    }

    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType>(0, base()->nlocal_body),
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
          bk.angmom[0] *= 0.5;
          bk.angmom[1] *= 0.5;
          bk.angmom[2] *= 0.5;
          angmom_to_omega(bk.angmom, bk.ex_space, bk.ey_space,
                                   bk.ez_space, bk.inertia, bk.omega);
        }
      }
    );
  };

  k_body.sync_device();
  base()->copymode = 1;
  if constexpr (is_nh) {
    const bool is_t = base()->tstat_flag, is_p = base()->pstat_flag;
    if      ( is_t &&  is_p) lambda.template operator()<true,true,true>();   // NPT
    else if ( is_t && !is_p) lambda.template operator()<true,true,false>();  // NVT
    else if (!is_t &&  is_p) lambda.template operator()<true,false,true>();  // NPH
    else                     lambda.template operator()<true,false,false>(); // NVE
  } else                     lambda.template operator()<false,false,false>();
  base()->copymode = 0;
  k_body.modify_device();

  base()->commflag = FINAL;
  base()->comm->forward_comm(fix_base(), 10);

  if (base()->vflag_atom) k_vatom.template sync<DeviceType>();

  if (base()->evflag) set_v_base<true>();
  else set_v_base<false>();

  atomKK->modified(execution_space, V_MASK);
  if (base()->extended) {
    // not implemented
  }
  if (base()->vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  if constexpr (is_nh) {
    auto temperature = nh_base()->temperature;
    // compute temperature
    if (nh_base()->tcomputeflag) {
      atomKK->sync(temperature->execution_space, temperature->datamask_read);
      base()->t_current = temperature->compute_scalar();
      atomKK->modified(temperature->execution_space, temperature->datamask_modify);
      atomKK->sync(execution_space, temperature->datamask_modify);
    }
    // pressure
    if (base()->pstat_flag) {
      auto pressure = base()->pressure;
      // accumulate kinetic energies for pstat
      base()->copymode = 1;
      k_body.sync_device();
      KK_ACC_FLOAT ke[2], keall[2];
      Kokkos::parallel_reduce(
        Kokkos::RangePolicy<DeviceType>(0, base()->nlocal_body),
        KOKKOS_LAMBDA(const int ibody, KK_ACC_FLOAT &l_akin_t, KK_ACC_FLOAT &l_akin_r ) {
          BodyKokkos &bk = l_body(ibody);
          l_akin_t += bk.mass*(bk.vcm[0]*bk.vcm[0] + bk.vcm[1]*bk.vcm[1] +
                    bk.vcm[2]*bk.vcm[2]);
          l_akin_r += bk.angmom[0]*bk.omega[0] + bk.angmom[1]*bk.omega[1] +
                    bk.angmom[2]*bk.omega[2];
        }, ke[0], ke[1]
      );
    base()->copymode = 0;
    MPI_Allreduce(ke, keall, 2, MPI_KK_ACC_FLOAT, MPI_SUM, base()->world);
    nh_base()->akin_t = keall[0];
    nh_base()->akin_r = keall[1];

    if (nh_base()->pstyle == ISO) {
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
    nh_base()->couple();
    pressure->addstep(base()->update->ntimestep+1);
    nh_base()->compute_press_target();
    nh_base()->nh_epsilon_dot();
  }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
double FixRigidBaseKokkos<DeviceType,FixRigidBase>::compute_scalar_base()
{
  KK_ACC_FLOAT t, t_all;
  k_body.sync_device();
  auto l_body = k_body.template view<DeviceType>();
  base()->copymode = 1;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<DeviceType>(0, base()->nlocal_body),
    KOKKOS_LAMBDA(const int &ibody, KK_ACC_FLOAT &l_t) {
      BodyKokkos &bk = l_body(ibody);
      l_t += bk.mass * (bk.vcm[0]*bk.vcm[0] + bk.vcm[1]*bk.vcm[1] + bk.vcm[2]*bk.vcm[2]);
      // for Iw^2 rotational term, need wbody = angular velocity in body frame
      // not omega = angular velocity in space frame
      KK_FLOAT wbody[3], rot[3][3];
      quat_to_mat(bk.quat, rot);
      transpose_matvec(rot, bk.angmom, wbody);
      if (bk.inertia[0] == 0.0) wbody[0] = 0.0;
      else wbody[0] /= bk.inertia[0];
      if (bk.inertia[1] == 0.0) wbody[1] = 0.0;
      else wbody[1] /= bk.inertia[1];
      if (bk.inertia[2] == 0.0) wbody[2] = 0.0;
      else wbody[2] /= bk.inertia[2];
      l_t += bk.inertia[0] * wbody[0] * wbody[0]
             + bk.inertia[1] * wbody[1] * wbody[1]
             + bk.inertia[2] * wbody[2] * wbody[2];
    }, t
  );
  base()->copymode = 0;
  MPI_Allreduce(&t, &t_all, 1, MPI_KK_ACC_FLOAT, MPI_SUM, base()->world);
  KK_ACC_FLOAT tfactor = base()->force->mvv2e / ((6.0*base()->nbody - base()->nlinear) * base()->force->boltz);
  t_all *= tfactor;
  return static_cast<double>(t_all);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::zero_momentum_base()
{
  base()->copymode = 1;
  k_body.sync_device();
  auto l_body = k_body.template view<DeviceType>();
  Kokkos::parallel_for(
    Kokkos::RangePolicy<DeviceType>(0, nbody_total()),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      bk.vcm[0] = bk.vcm[1] = bk.vcm[2] = 0.0;
    }
  );
  k_body.modify_device();
  base()->copymode = 0;
  // forward communicate of omega to all ghost copies
  base()->commflag = FINAL;
  base()->comm->forward_comm(fix_base(),10);
  // set velocity of atoms in rigid bodues
  set_v_base<false>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::zero_rotation_base()
{
  k_body.sync_device();
  auto l_body = k_body.template view<DeviceType>();
  base()->copymode = 1;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<DeviceType>(0, nbody_total()),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      bk.angmom[0] = bk.angmom[1] = bk.angmom[2] = 0.0;
      bk.omega[0] = bk.omega[1] = bk.omega[2] = 0.0;
    }
  );
  base()->copymode = 0;
  k_body.modify_device();
  // forward communicate of omega to all ghost copies
  base()->commflag = FINAL;
  base()->comm->forward_comm(fix_base(),10);
  // set velocity of atoms in rigid bodues
  set_v_base<false>();
}

/* ----------------------------------------------------------------------
   PROTECTED METHODS
------------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::setup_bodies_static_base()
{

  const int nlocal = atomKK->nlocal;
  const int nlocal_body = base()->nlocal_body;
  auto policy_nlocal = Kokkos::RangePolicy<DeviceType>(0, nlocal);
  auto policy_nlocal_body = Kokkos::RangePolicy<DeviceType>(0, nlocal_body);
  auto policy_nbody_total = Kokkos::RangePolicy<DeviceType>(0, nbody_total());

  // extended = 1 if any particle in a rigid body is finite size
  //              or has a dipole moment

  base()->extended = base()->orientflag = base()->dorientflag = 0;

  AtomVecEllipsoid::Bonus *ebonus;
  if (base()->avec_ellipsoid) ebonus = base()->avec_ellipsoid->bonus;
  AtomVecLine::Bonus *lbonus;
  if (base()->avec_line) lbonus = base()->avec_line->bonus;
  AtomVecTri::Bonus *tbonus;
  if (base()->avec_tri) tbonus = base()->avec_tri->bonus;
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
    for (int i = 0; i < nlocal; i++) {
      if (base()->bodytag[i] == 0) continue;
      if (radius && radius[i] > 0.0) flag = 1;
      if (ellipsoid && ellipsoid[i] >= 0) flag = 1;
      if (line && line[i] >= 0) flag = 1;
      if (tri && tri[i] >= 0) flag = 1;
      if (mu && mu[i][3] > 0.0) flag = 1;
    }
    MPI_Allreduce(&flag, &base()->extended, 1, MPI_INT, MPI_MAX, base()->world);
  }

  // extended = 1 if using molecule template with finite-size particles
  // require all molecules in template to have consistent radiusflag

  if (base()->onemols) {
    int radiusflag = base()->onemols[0]->radiusflag;
    for (int i = 1; i < base()->nmol; i++) {
      if (base()->onemols[i]->radiusflag != radiusflag)
        base()->error->all(FLERR, Error::NOLASTLINE, "Inconsistent use of finite-size particles "
                   "by molecule template molecules");
    }
    if (radiusflag) base()->extended = 1;
  }

  if (base()->extended) {
    base()->error->all(FLERR, Error::NOLASTLINE,
               "Fix {}: extended particles not implemented in KOKKOS", base()->style);
  }

  // grow extended arrays and set extended flags for each particle
  // orientflag = 4 if any particle stores ellipsoid or tri orientation or quat
  // orientflag = 1 if any particle stores line orientation
  // dorientflag = 1 if any particle stores dipole orientation

  if (base()->extended) {
    // not implemented
  }

  // kokkos views

  //atomKK->sync(execution_space, X_MASK | F_MASK | IMAGE_MASK | RMASS_MASK | TYPE_MASK );
  atomKK->sync(execution_space, ALL_MASK );
  auto l_x = atomKK->k_x.template view<DeviceType>();
  auto l_f = atomKK->k_f.template view<DeviceType>();
  auto l_image = atomKK->k_image.template view<DeviceType>();
  auto l_rmass = atomKK->k_rmass.template view<DeviceType>();
  auto l_type = atomKK->k_type.template view<DeviceType>();
  auto l_mass = atomKK->k_mass.template view<DeviceType>();

  auto l_atom2body = k_atom2body.template view<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();
  auto l_bodytag = k_bodytag.template view<DeviceType>();
  auto l_xcmimage = k_xcmimage.template view<DeviceType>();
  auto l_displace = k_displace.template view<DeviceType>();

  k_atom2body.template sync<DeviceType>();
  k_body.template sync<DeviceType>();
  k_bodytag.template sync<DeviceType>();
  k_xcmimage.template sync<DeviceType>();
  k_displace.template sync<DeviceType>();

  // domainKK

  auto l_triclinic = base()->triclinic;
  auto l_periodic = Few<bool,3>(domainKK->periodicity);
  auto l_prd = Few<KK_FLOAT,3>(domainKK->prd);
  auto l_h = Few<KK_FLOAT,6>(domainKK->h);
  auto l_invprd = Few<KK_FLOAT,3>(
    KK_FLOAT(1.0) / domainKK->prd[0],
    KK_FLOAT(1.0) / domainKK->prd[1],
    KK_FLOAT(1.0) / domainKK->prd[2]
  );

  // set body xcmimage flags = true image flags
  k_xcmimage.template sync<DeviceType>();
  base()->copymode = 1;
  Kokkos::parallel_for(policy_nlocal, KOKKOS_LAMBDA(const int &i) {
    if (l_bodytag(i) >= 0) l_xcmimage(i) = l_image(i);
    else l_xcmimage(i) = 0;
    //Kokkos::printf("*** i %i l_bodytag %i l_xcmimage %i l_image %i\n", i, l_bodytag(i), l_xcmimage(i), l_image(i));

  });
  base()->copymode = 0;
  k_xcmimage.template modify<DeviceType>();

  // acquire ghost bodies via forward comm
  // set atom2body for ghost atoms via forward comm
  // set atom2body for other owned atoms via reset_atom2body()
  base()->nghost_body = 0;
  base()->commflag = FULL_BODY;
  base()->comm->forward_comm(fix_base());
  reset_atom2body_base();

  // compute mass & center-of-mass of each rigid body
  k_body.sync_device();
  base()->copymode = 1;
  Kokkos::parallel_for(policy_nbody_total, KOKKOS_LAMBDA(const int &ibody) {
    BodyKokkos &bk = l_body(ibody);
    bk.xcm[0] = bk.xcm[1] = bk.xcm[2] = KK_FLOAT(0.0);
    bk.xgc[0] = bk.xgc[1] = bk.xgc[2] = KK_FLOAT(0.0);
    bk.mass = KK_FLOAT(0.0);
    bk.natoms = 0;
  });
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
    Kokkos::printf("*** i %i bk.mass %f massone %f bk.natoms %i unwrap %f %f %f bk.xcm %f %f %f\n", i, bk.mass, massone, bk.natoms, unwrap[0], unwrap[1], unwrap[2], bk.xcm[0], bk.xcm[1], bk.xcm[2]);
  };
  if (l_rmass.data()) {
    Kokkos::parallel_for(policy_nlocal, KOKKOS_LAMBDA(const int i) {
      lambda_xcm_mass.template operator()<true>(i);
    });
  } else {
    Kokkos::parallel_for(policy_nlocal, KOKKOS_LAMBDA(const int i) {
      lambda_xcm_mass.template operator()<false>(i);
    });
  }
  base()->copymode = 0;
  k_body.modify_device();

  // reverse communicate xcm, mass of all bodies
  base()->commflag = XCM_MASS;
  base()->comm->reverse_comm(fix_base(),8);

  k_body.sync_device();
  base()->copymode = 1;
  Kokkos::parallel_for(policy_nlocal_body, KOKKOS_LAMBDA(const int &ibody) {
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
    bk.vcm[0] = bk.vcm[1] = bk.vcm[2] = 0.0;
    bk.angmom[0] = bk.angmom[1] = bk.angmom[2] = 0.0;
    // set rigid body image flags to default values
    bk.image = ((imageint) IMGMAX << IMG2BITS) |
               ((imageint) IMGMAX << IMGBITS) | IMGMAX;

    Kokkos::printf("*** ibody %i bk.mass %f bk.natoms %i bk.xcm %f %f %f bk.image %i\n", ibody, bk.mass, bk.natoms, bk.xcm[0], bk.xcm[1], bk.xcm[2], bk.image);
  });
  base()->copymode = 0;
  k_body.modify_device();

  // overwrite masstotal, center-of-mass, image flags with file values
  // inbody[i] = 0/1 if Ith rigid body is initialized by file
  DAT::tdual_int_1d k_inbody("rigid/small:inbody", nlocal_body);
  deep_copy(k_inbody.view_host(), 0);

  auto l_inbody = k_inbody.view_device();
  if (base()->inpfile) {
    // must call it here so it doesn't override read in data but
    // initialize bodies whose dynamic settings not set in inpfile
    setup_bodies_dynamic_base();
    base()->readfile(0, nullptr, k_inbody.view_host().data());
  }
  k_inbody.modify_host();

  // remap the xcm of each body back into simulation box
  //   and reset body and atom xcmimage flags via pre_neighbor()
  pre_neighbor_base();

  // compute 6 moments of inertia of each body in Cartesian reference frame
  // dx,dy,dz = coords relative to center-of-mass
  // symmetric 3x3 inertia tensor stored in Voigt notation as 6-vector
  //base()->memory->create(base()->itensor, nbody_total(), 6, "rigid/small:itensor");
  base()->memoryKK->create_kokkos(k_itensor, base()->itensor, nbody_total(), 6, "rigid/small:itensor");

  deep_copy(k_itensor.view_device(), KK_FLOAT(0.0));

  auto l_itensor = k_itensor.view_device();
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
    Kokkos::parallel_for(policy_nlocal, KOKKOS_LAMBDA(const int i) {
      lambda_itensor.template operator()<true,true>(i);
    });
  } else {
    Kokkos::parallel_for(policy_nlocal, KOKKOS_LAMBDA(const int i) {
      lambda_itensor.template operator()<false,true>(i);
    });
  }
  base()->copymode = 0;

  // extended particles may contribute extra terms to moments of inertia
  if (base()->extended) {
    // not implemented
  }

  // reverse communicate inertia tensor of all bodies
  base()->commflag = ITENSOR;
  base()->comm->reverse_comm(fix_base(),6);

  // overwrite Cartesian inertia tensor with file values
  if (base()->inpfile) {
    k_itensor.sync_host();
    k_inbody.sync_host();
    base()->readfile(1, base()->itensor, k_inbody.view_host().data());
    k_itensor.modify_host();
    k_inbody.modify_host();
    k_itensor.sync_device();
    k_inbody.sync_device();
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
    if (bk.inertia[0] < EPSILON*max) bk.inertia[0] = 0.0;
    if (bk.inertia[1] < EPSILON*max) bk.inertia[1] = 0.0;
    if (bk.inertia[2] < EPSILON*max) bk.inertia[2] = 0.0;
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
  base()->copymode = 1;
  KK_FLOAT dflag = KK_FLOAT(0.0);
  if (domainKK->dimension == 2) {
    if (l_triclinic) {
      Kokkos::parallel_reduce(policy_nlocal_body,
        KOKKOS_LAMBDA(const int ibody, KK_FLOAT &l_dflag) {
        lambda_initial.template operator()<true,true>(ibody, l_dflag);
      }, Kokkos::Max<KK_FLOAT>(dflag));
    } else {
      Kokkos::parallel_reduce(policy_nlocal_body,
        KOKKOS_LAMBDA(const int ibody, KK_FLOAT &l_dflag) {
        lambda_initial.template operator()<true,false>(ibody, l_dflag);
      }, Kokkos::Max<KK_FLOAT>(dflag));
    }
  } else {
    if (l_triclinic) {
      Kokkos::parallel_reduce(policy_nlocal_body,
        KOKKOS_LAMBDA(const int ibody, KK_FLOAT &l_dflag) {
        lambda_initial.template operator()<false,true>(ibody, l_dflag);
      }, Kokkos::Max<KK_FLOAT>(dflag));
    } else {
      Kokkos::parallel_reduce(policy_nlocal_body,
        KOKKOS_LAMBDA(const int ibody, KK_FLOAT &l_dflag) {
        lambda_initial.template operator()<false,false>(ibody, l_dflag);
      }, Kokkos::Max<KK_FLOAT>(dflag));
    }
  }
  base()->copymode = 0;
  if( dflag > KK_FLOAT(0.0))
    base()->error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dflag);

  // forward communicate updated info of all bodies
  base()->commflag = INITIAL;
  base()->comm->forward_comm(base(),29);

  // displace = initial atom coords in basis of principal axes
  // set displace = 0.0 for atoms not in any rigid body
  // for extended particles, set their orientation wrt to rigid body

  auto lambda_displace = [&]<bool EXTENDED>(const int i) {
    const int ibody = l_atom2body(i);
    if (ibody < 0) {
      l_displace(i,0) = l_displace(i,1) = l_displace(i,2) = 0.0;
      return;
    }
    BodyKokkos &bk = l_body(ibody);
    //Kokkos::printf("*** ex %f %f %f ey %f %f %f ez %f %f %f\n", bk.ex_space[0], bk.ex_space[1], bk.ex_space[2], bk.ey_space[0], bk.ey_space[1], bk.ey_space[2], bk.ez_space[0], bk.ez_space[1], bk.ez_space[2]);
    auto unwrap = DomainKokkos::unmap(l_prd, l_h, l_triclinic, &l_x(i,0), l_xcmimage(i));
    auto delta3 = unwrap - Few<KK_FLOAT,3>(bk.xcm);
    transpose_matvec(bk.ex_space,bk.ey_space,bk.ez_space, delta3.data(), &l_displace(i,0));
    //Kokkos::printf("*** unwrap %f %f %f delta3 %f %f %f l_displace(%i) %f %f %f\n", unwrap[0], unwrap[1], unwrap[2], delta3[0], delta3[1], delta3[2], i, l_displace(i,0), l_displace(i,1), l_displace(i,2));
    if constexpr (EXTENDED) {
      // not implemented
    }
  };
  base()->copymode = 1;
  if (base()->extended) {
    Kokkos::parallel_for(policy_nlocal, KOKKOS_LAMBDA(const int i) {
      lambda_displace.template operator()<true>(i);
    });
  } else {
    Kokkos::parallel_for(policy_nlocal, KOKKOS_LAMBDA(const int i) {
      lambda_displace.template operator()<false>(i);
    });
  }
  k_displace.modify_device();

  // test for valid principal moments & axes
  // recompute moments of inertia around new axes
  // 3 diagonal moments should equal principal moments
  // 3 off-diagonal moments should be 0.0
  // extended particles may contribute extra terms to moments of inertia
  deep_copy(k_itensor.view_device(), KK_FLOAT(0.0));
  if (l_rmass.data()) {
    Kokkos::parallel_for(policy_nlocal, KOKKOS_LAMBDA(const int i) {
      lambda_itensor.template operator()<true,false>(i);
    });
  } else {
    Kokkos::parallel_for(policy_nlocal, KOKKOS_LAMBDA(const int i) {
      lambda_itensor.template operator()<false,false>(i);
    });
  }
  k_itensor.modify_device();
  if (base()->extended) {
    // not implemented
  }
  base()->copymode = 0;

  // reverse communicate inertia tensor of all bodies
  base()->commflag = ITENSOR;
  base()->comm->reverse_comm(fix_base(),6);

  // error check that re-computed moments of inertia match diagonalized ones
  // do not do test for bodies with params read from inpfile

  int flag = 0;
  bool l_inpfile = (base()->inpfile != nullptr);
  k_inbody.sync_device();
  base()->copymode = 1;
  Kokkos::parallel_reduce(policy_nlocal_body, KOKKOS_LAMBDA(const int ibody, int& l_flag) {
    if (l_inpfile && l_inbody(ibody)) return;
    BodyKokkos &bk = l_body(ibody);
    KK_FLOAT norm = KK_FLOAT(0.0);
    l_flag = 0;
    for( int i=0; i<3 ; i++ ) {
      const auto a = bk.inertia[i];
      const auto b = l_itensor(ibody,i);
      //Kokkos::printf("*** ibody %i a %.16g b %.16g TOLERANCE %g l_flag %i\n", ibody, a, b, TOLERANCE, l_flag);
      if (a == 0.0 && Kokkos::fabs(b) > TOLERANCE) l_flag = 1;
      if (a != 0.0 && Kokkos::fabs((b-a)/a) > TOLERANCE) l_flag = 1;
      //if (l_flag) Kokkos::printf("*** ibody %i a %f b %f\n", ibody, a, b);
      norm += a;
    }
    norm *= KK_FLOAT(1.0/3.0);
    if (Kokkos::fabs(l_itensor(ibody,3) / norm) > TOLERANCE ||
        Kokkos::fabs(l_itensor(ibody,4) / norm) > TOLERANCE ||
        Kokkos::fabs(l_itensor(ibody,5) / norm) > TOLERANCE) {
          l_flag = 1;
          //Kokkos::printf("*** ibody %i l_itensor(ibody,3,4,5) %f %f %f norm %f\n", ibody, l_itensor(ibody,3), l_itensor(ibody,4), l_itensor(ibody,5), norm);
      }
  }, Kokkos::Max<int>(flag));
  base()->copymode = 0;
  if (flag) base()->error->all(FLERR, Error::NOLASTLINE, "Fix {}: Bad principal moments", base()->style);

  // cleanup
  base()->memoryKK->destroy_kokkos(k_itensor, base()->itensor);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::setup_bodies_dynamic_base()
{

  // sum vcm, angmom across all rigid bodies
  // vcm = velocity of COM
  // angmom = angular momentum around COM

  atomKK->sync(execution_space, X_MASK | V_MASK | F_MASK | RMASS_MASK | TYPE_MASK );
  k_atom2body.template sync<DeviceType>();
  k_body.template sync<DeviceType>();
  k_xcmimage.template sync<DeviceType>();
  auto l_x = atomKK->k_x.template view<DeviceType>();
  auto l_v = atomKK->k_v.template view<DeviceType>();
  auto l_f = atomKK->k_f.template view<DeviceType>();
  auto l_rmass = atomKK->k_rmass.template view<DeviceType>();
  auto l_type = atomKK->k_type.template view<DeviceType>();

  // FIXME: what *_MASK to sync for k_mass ???
  auto l_mass = atomKK->k_mass.template view<DeviceType>();

  auto l_atom2body = k_atom2body.template view<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();
  auto l_xcmimage = k_xcmimage.template view<DeviceType>();
  auto l_prd = Few<KK_FLOAT,3>(domainKK->prd);
  auto l_h = Few<KK_FLOAT,6>(domainKK->h);
  auto l_triclinic = base()->triclinic;

  base()->copymode = 1;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<DeviceType>(0, nbody_total()),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      bk.vcm[0] = bk.vcm[1] = bk.vcm[2] = 0.0;
      bk.angmom[0] = bk.angmom[1] = bk.angmom[2] = 0.0;
    }
  );

  auto lambda = [&]<bool RMASS>(const int &i) {
    const int ibody = l_atom2body(i);
    if (ibody < 0) return;
    BodyKokkos &bk = l_body(ibody);
    KK_FLOAT massone;
    if constexpr (RMASS) massone = l_rmass(i);
    else massone = l_mass(l_type(i));
    bk.vcm[0] += l_v(i,0) * massone;
    bk.vcm[1] += l_v(i,1) * massone;
    bk.vcm[2] += l_v(i,2) * massone;
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
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int i) {
      lambda.template operator()<true>(i);
    });
  } else {
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int i) {
      lambda.template operator()<false>(i);
    });
  }
  base()->copymode = 0;
  k_body.modify_device();

  // extended particles add their rotation to angmom of body
  if (base()->extended) {
    // not implemented
  }

  // reverse communicate vcm, angmom of all bodies
  base()->commflag = VCM_ANGMOM;
  base()->comm->reverse_comm(base(),6);

  // normalize velocity of COM
  k_body.template sync<DeviceType>();
  base()->copymode = 1;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<DeviceType>(0, base()->nlocal_body),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      bk.vcm[0] /= bk.mass;
      bk.vcm[1] /= bk.mass;
      bk.vcm[2] /= bk.mass;
    }
  );
  base()->copymode = 0;
  k_body.modify_device();

}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::apply_langevin_thermostat_base()
{
  // grow langextra if needed
  if (base()->nlocal_body > base()->maxlang) {
    base()->maxlang = base()->nlocal_body + base()->nghost_body;
    base()->memoryKK->destroy_kokkos(k_langextra, base()->langextra);
    base()->memoryKK->create_kokkos(k_langextra, base()->langextra, 6, "rigid/small:langextra");
  }

  auto update = base()->update;
  auto force = base()->force;
  KK_FLOAT delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  const KK_FLOAT l_t_target = base()->t_start + delta * (base()->t_stop-base()->t_start);
  const KK_FLOAT l_tsqrt = sqrt(l_t_target);
  const KK_FLOAT l_t_period = static_cast<KK_FLOAT>(base()->t_period);

  const KK_FLOAT l_ftm2v = force->ftm2v;
  const KK_FLOAT l_langfactor = sqrt(24.0 * force->boltz / base()->t_period
    / update->dt / force->mvv2e) / l_ftm2v;

  base()->copymode = 1;
  k_body.sync_device();
  auto l_body = k_body.template view<DeviceType>();
  auto l_rand_pool = rand_pool;
  auto l_langextra = k_langextra.template view<DeviceType>();
  Kokkos::parallel_for(
    Kokkos::RangePolicy<DeviceType>(0, nbody_total()),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      KK_FLOAT gamma1 = -bk.mass / l_t_period / l_ftm2v;
      KK_FLOAT gamma2 = sqrt(bk.mass) * l_tsqrt * l_langfactor;

      rand_type rand_gen = l_rand_pool.get_state();
#if defined (LMP_KOKKOS_SINGLE_SINGLE) || defined (LMP_KOKKOS_SINGLE_DOUBLE)
      const float rnd1 = rand_gen.frand() - 0.5;
      const float rnd2 = rand_gen.frand() - 0.5;
      const float rnd3 = rand_gen.frand() - 0.5;
      const float rnd4 = rand_gen.frand() - 0.5;
      const float rnd5 = rand_gen.frand() - 0.5;
      const float rnd6 = rand_gen.frand() - 0.5;
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
      gamma1 = -1.0 / l_t_period / l_ftm2v;
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
  base()->copymode = 0;
  k_langextra.modify_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::remap_base()
{
  if constexpr (is_nh) {
  int nlocal = atomKK->nlocal;
  // epsilon is not used, except for book-keeping
  for (int i = 0; i < 3; i++) nh_base()->epsilon[i] += base()->dtq * nh_base()->epsilon_dot[i];
  // convert pertinent atoms and rigid bodies to lamda coords
  if (base()->allremap) domainKK->x2lamda(nlocal);
  else domainKK->x2lamda(nlocal, base()->dilate_group_bit);
  for (auto &ifix : nh_base()->rfix) ifix->deform_base(0);
  // reset global and local box to new size/shape
  for (int i = 0; i < 3; i++) {
    if (base()->p_flag[i]) {
      const double oldlo = domainKK->boxlo[i];
      const double oldhi = domainKK->boxhi[i];
      const double ctr = 0.5 * (oldlo + oldhi);
      const double expfac = exp(base()->dtq * nh_base()->epsilon_dot[i]);
      domainKK->boxlo[i] = (oldlo-ctr)*expfac + ctr;
      domainKK->boxhi[i] = (oldhi-ctr)*expfac + ctr;
    }
  }
  domainKK->set_global_box();
  domainKK->set_local_box();
  // convert pertinent atoms and rigid bodies back to box coords
  if (base()->allremap) domainKK->lamda2x(nlocal);
  else domainKK->lamda2x(nlocal, base()->dilate_group_bit);
  for (auto &ifix : nh_base()->rfix) ifix->deform_base(1);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::deform_base(int flag)
{
  k_body.sync_device();
  auto l_body = k_body.template view<DeviceType>();
  base()->copymode = 1;
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
    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType>(0, base()->nlocal_body),
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
    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType>(0, base()->nlocal_body),
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
  base()->copymode = 0;
  k_body.modify_device();
}




/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
template<bool EVFLAG>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::set_xv_base()
{
  auto l_x = atomKK->k_x.template view<DeviceType>();
  auto l_v = atomKK->k_v.template view<DeviceType>();
  auto l_f = atomKK->k_f.template view<DeviceType>();
  auto l_rmass = atomKK->k_rmass.template view<DeviceType>();
  auto l_mass = atomKK->k_mass.template view<DeviceType>();
  auto l_type = atomKK->k_type.template view<DeviceType>();
  atomKK->sync(execution_space, X_MASK | V_MASK | F_MASK | TYPE_MASK | RMASS_MASK);

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

  KK_FLOAT l_dtf = static_cast<KK_FLOAT>(base()->dtf);
  auto l_vflag_global = base()->vflag_global;
  auto l_vflag_atom = base()->vflag_atom;

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

    KK_FLOAT x0 = 0.0, x1 = 0.0, x2 = 0.0, vx = 0.0, vy = 0.0, vz = 0.0;
    if constexpr (EVFLAG) {
      x0 = l_x(i,0) + deltax;
      x1 = l_x(i,1) + deltay;
      x2 = l_x(i,2) + deltaz;
      vx = l_v(i,0);
      vy = l_v(i,1);
      vz = l_v(i,2);
    }

    KK_FLOAT xnew[3];
    matvec(bk.ex_space, bk.ey_space, bk.ez_space, &l_displace(i, 0), xnew);

    if constexpr (EVFLAG) {
      // Compute v_new in KK_ACC_FLOAT before truncating to KK_FLOAT for storage,
      // so the pre-truncation value can be used for the constraint-force virial.
      const KK_ACC_FLOAT vnew0 = fma(KK_ACC_FLOAT(bk.omega[1]), KK_ACC_FLOAT(xnew[2]),
                               fma(KK_ACC_FLOAT(-bk.omega[2]), KK_ACC_FLOAT(xnew[1]),
                               KK_ACC_FLOAT(bk.vcm[0])));
          const KK_ACC_FLOAT vnew1 = fma(KK_ACC_FLOAT(bk.omega[2]), KK_ACC_FLOAT(xnew[0]),
                               fma(KK_ACC_FLOAT(-bk.omega[0]), KK_ACC_FLOAT(xnew[2]),
                               KK_ACC_FLOAT(bk.vcm[1])));
          const KK_ACC_FLOAT vnew2 = fma(KK_ACC_FLOAT(bk.omega[0]), KK_ACC_FLOAT(xnew[1]),
                               fma(KK_ACC_FLOAT(-bk.omega[1]), KK_ACC_FLOAT(xnew[0]),
                               KK_ACC_FLOAT(bk.vcm[2])));
          l_v(i,0) = KK_FLOAT(vnew0);
          l_v(i,1) = KK_FLOAT(vnew1);
          l_v(i,2) = KK_FLOAT(vnew2);
          l_x(i,0) = xnew[0] + bk.xcm[0] - deltax;
          l_x(i,1) = xnew[1] + bk.xcm[1] - deltay;
          l_x(i,2) = xnew[2] + bk.xcm[2] - deltaz;

          double massone;
          if (l_rmass.data()) massone = l_rmass(i);
          else massone = l_mass(l_type(i));
          const double half_m_dt = 0.5 * massone / l_dtf;
          const KK_ACC_FLOAT fc0 = fma(half_m_dt, vnew0 - vx, -0.5*l_f(i,0));
          const KK_ACC_FLOAT fc1 = fma(half_m_dt, vnew1 - vy, -0.5*l_f(i,1));
          const KK_ACC_FLOAT fc2 = fma(half_m_dt, vnew2 - vz, -0.5*l_f(i,2));

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
  base()->copymode = 1;
  if constexpr (!EVFLAG) {
    // TagRigidSetXV<TRICLINIC,HALF,0>>(0, nlocal),
    EV_FLOAT ev;
    if (base()->triclinic) {
      Kokkos::parallel_for(
        Kokkos::RangePolicy<DeviceType>(0, nlocal),
        KOKKOS_LAMBDA(const int &i) {
          EV_FLOAT ev;
          lambda.template operator()<true,HALF>(i, ev);
        }
      );
    } else {
      Kokkos::parallel_for(
        Kokkos::RangePolicy<DeviceType>(0, nlocal),
        KOKKOS_LAMBDA(const int &i) {
          EV_FLOAT ev;
          lambda.template operator()<false,HALF>(i, ev);
        }
      );
    }
  } else {
    auto kokkos = base()->lmp->kokkos;
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
      if (base()->triclinic) {
        Kokkos::parallel_reduce(
          Kokkos::RangePolicy<DeviceType>(0, nlocal),
          KOKKOS_LAMBDA(const int &i, EV_FLOAT &ev_) {
            lambda.template operator()<true,HALF>(i, ev_);
          }, ev
        );
      } else {
        Kokkos::parallel_reduce(
          Kokkos::RangePolicy<DeviceType>(0, nlocal),
          KOKKOS_LAMBDA(const int &i, EV_FLOAT &ev_) {
            lambda.template operator()<false,HALF>(i, ev_);
          }, ev
        );
      }
    } else {
      // TagRigidSetXV<TRICLINIC,HALFTHREAD>>(0, nlocal),
      if (base()->triclinic) {
        Kokkos::parallel_reduce(
          Kokkos::RangePolicy<DeviceType>(0, nlocal),
          KOKKOS_LAMBDA(const int &i, EV_FLOAT &ev_) {
            lambda.template operator()<true,HALFTHREAD>(i, ev_);
          }, ev
        );
      } else {
        Kokkos::parallel_reduce(
          Kokkos::RangePolicy<DeviceType>(0, nlocal),
          KOKKOS_LAMBDA(const int &i, EV_FLOAT &ev_) {
            lambda.template operator()<false,HALFTHREAD>(i, ev_);
          }, ev
        );
      }
    }
    if (l_vflag_global) {
      base()->virial[0] += static_cast<double>(ev.v[0]);
      base()->virial[1] += static_cast<double>(ev.v[1]);
      base()->virial[2] += static_cast<double>(ev.v[2]);
      base()->virial[3] += static_cast<double>(ev.v[3]);
      base()->virial[4] += static_cast<double>(ev.v[4]);
      base()->virial[5] += static_cast<double>(ev.v[5]);
    }
    if (l_vflag_atom && need_dup) {
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
      dup_vatom = {};
    }
  }
  // update geometric center of bodies
  Kokkos::parallel_for(
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
  base()->copymode = 0;
  atomKK->modified(execution_space, X_MASK | V_MASK);
  k_body.modify_device();

  if (base()->extended) {
    // not implemented
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
template<bool EVFLAG>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::set_v_base()
{

  auto l_x = atomKK->k_x.template view<DeviceType>();
  auto l_v = atomKK->k_v.template view<DeviceType>();
  auto l_f = atomKK->k_f.template view<DeviceType>();
  auto l_rmass = atomKK->k_rmass.template view<DeviceType>();
  auto l_mass = atomKK->k_mass.template view<DeviceType>();
  auto l_type = atomKK->k_type.template view<DeviceType>();
  atomKK->sync(execution_space, X_MASK | V_MASK | F_MASK | TYPE_MASK | RMASS_MASK);

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

  KK_FLOAT l_dtf = static_cast<KK_FLOAT>(base()->dtf);
  auto l_vflag_global = base()->vflag_global;
  auto l_vflag_atom = base()->vflag_atom;

  auto lambda = [&]<bool TRICLINIC, int NEIGHFLAG>(const int &i, EV_FLOAT &ev) {
    const int ibody = l_atom2body(i);
    if (ibody < 0) return;
    const BodyKokkos &bk = l_body(ibody);

    KK_FLOAT delta[3];
    matvec(bk.ex_space, bk.ey_space, bk.ez_space, &l_displace(i, 0), delta);

    if constexpr (EVFLAG) {
      const double vx = l_v(i,0);
      const double vy = l_v(i,1);
      const double vz = l_v(i,2);
      // Compute v_new in KK_ACC_FLOAT before truncating to KK_FLOAT for storage,
      // so the pre-truncation value can be used for the constraint-force virial.
      const KK_ACC_FLOAT vnew0 = fma(KK_ACC_FLOAT(bk.omega[1]), KK_ACC_FLOAT(delta[2]),
                               fma(KK_ACC_FLOAT(-bk.omega[2]), KK_ACC_FLOAT(delta[1]),
                               KK_ACC_FLOAT(bk.vcm[0])));
      const KK_ACC_FLOAT vnew1 = fma(KK_ACC_FLOAT(bk.omega[2]), KK_ACC_FLOAT(delta[0]),
                               fma(KK_ACC_FLOAT(-bk.omega[0]), KK_ACC_FLOAT(delta[2]),
                               KK_ACC_FLOAT(bk.vcm[1])));
      const KK_ACC_FLOAT vnew2 = fma(KK_ACC_FLOAT(bk.omega[0]), KK_ACC_FLOAT(delta[1]),
                               fma(KK_ACC_FLOAT(-bk.omega[1]), KK_ACC_FLOAT(delta[0]),
                               KK_ACC_FLOAT(bk.vcm[2])));
      l_v(i,0) = KK_FLOAT(vnew0);
      l_v(i,1) = KK_FLOAT(vnew1);
      l_v(i,2) = KK_FLOAT(vnew2);
      double massone;
      if (l_rmass.data()) massone = l_rmass(i);
      else massone = l_mass(l_type(i));
      const double half_m_dt = 0.5 * massone / l_dtf;
      const KK_ACC_FLOAT fc0 = fma(half_m_dt, vnew0 - vx, -0.5*l_f(i,0));
      const KK_ACC_FLOAT fc1 = fma(half_m_dt, vnew1 - vy, -0.5*l_f(i,1));
      const KK_ACC_FLOAT fc2 = fma(half_m_dt, vnew2 - vz, -0.5*l_f(i,2));

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
  base()->copymode = 1;
  if constexpr (!EVFLAG) {
    // TagRigidSetV<TRICLINIC,HALF,0>>(0, nlocal),
    if (base()->triclinic) {
      Kokkos::parallel_for(
        Kokkos::RangePolicy<DeviceType>(0, nlocal),
        KOKKOS_LAMBDA(const int &i) {
          EV_FLOAT ev;
          lambda.template operator()<true,HALF>(i, ev);
        }
      );
    } else {
      Kokkos::parallel_for(
        Kokkos::RangePolicy<DeviceType>(0, nlocal),
        KOKKOS_LAMBDA(const int &i) {
          EV_FLOAT ev;
          lambda.template operator()<false,HALF>(i, ev);
        }
      );
    }
  } else {
    auto kokkos = base()->lmp->kokkos;
    int neighflag = kokkos->neighflag;
    if (neighflag == FULL) {
      neighflag =
        (kokkos->nthreads > 1 || kokkos->ngpus > 0) ? HALFTHREAD : HALF;
    }
    const bool need_dup =
      (neighflag != HALF) &&
      std::is_same_v<NeedDup_v<HALFTHREAD, DeviceType>,
                     Kokkos::Experimental::ScatterDuplicated>;
    if (base()->vflag_atom) {
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
      if (base()->triclinic) {
        Kokkos::parallel_reduce(
          Kokkos::RangePolicy<DeviceType>(0, nlocal),
          KOKKOS_LAMBDA(const int &i, EV_FLOAT &ev_) {
            lambda.template operator()<true,HALF>(i, ev_);
          }, ev
        );
      } else {
        Kokkos::parallel_reduce(
          Kokkos::RangePolicy<DeviceType>(0, nlocal),
          KOKKOS_LAMBDA(const int &i, EV_FLOAT &ev_) {
            lambda.template operator()<false,HALF>(i, ev_);
          }, ev
        );
      }
    } else {
      // TagRigidSetV<TRICLINIC,HALFTHREAD>>(0, nlocal),
      if (base()->triclinic) {
        Kokkos::parallel_reduce(
          Kokkos::RangePolicy<DeviceType>(0, nlocal),
          KOKKOS_LAMBDA(const int &i, EV_FLOAT &ev_) {
            lambda.template operator()<true,HALFTHREAD>(i, ev_);
          }, ev
        );
      } else {
        Kokkos::parallel_reduce(
          Kokkos::RangePolicy<DeviceType>(0, nlocal),
          KOKKOS_LAMBDA(const int &i, EV_FLOAT &ev_) {
            lambda.template operator()<false,HALFTHREAD>(i, ev_);
          }, ev
        );
      }
    }
    if (l_vflag_global) {
      base()->virial[0] += static_cast<double>(ev.v[0]);
      base()->virial[1] += static_cast<double>(ev.v[1]);
      base()->virial[2] += static_cast<double>(ev.v[2]);
      base()->virial[3] += static_cast<double>(ev.v[3]);
      base()->virial[4] += static_cast<double>(ev.v[4]);
      base()->virial[5] += static_cast<double>(ev.v[5]);
    }
    if (l_vflag_atom && need_dup) {
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
      dup_vatom = {};
    }
  }
  atomKK->modified(execution_space, V_MASK);
  k_body.modify_device();
  base()->copymode = 0;

  if (base()->extended) {
    // not implemented
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::compute_forces_and_torques_base()
{

  const int nlocal = atomKK->nlocal;
  base()->copymode = 1;
  k_body.template sync<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();
  Kokkos::parallel_for(
    Kokkos::RangePolicy<DeviceType>(0, nbody_total()),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      bk.fcm[0] = bk.fcm[1] = bk.fcm[2] = KK_FLOAT(0.0);
      bk.torque[0] = bk.torque[1] = bk.torque[2] = KK_FLOAT(0.0);
    }
  );

  atomKK->sync(execution_space, X_MASK | F_MASK );
  k_atom2body.template sync<DeviceType>();
  k_xcmimage.template sync<DeviceType>();
  auto l_atom2body = k_atom2body.template view<DeviceType>();
  auto l_xcmimage = k_xcmimage.template view<DeviceType>();
  auto l_x = atomKK->k_x.template view<DeviceType>();
  auto l_f = atomKK->k_f.template view<DeviceType>();
  auto l_prd = Few<KK_FLOAT,3>(domainKK->prd);
  auto l_h = Few<KK_FLOAT,6>(domainKK->h);
  auto l_triclinic = base()->triclinic;
  Kokkos::parallel_for(
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
  k_body.modify_device();
  base()->copymode = 0;

  if (base()->extended) {
    // not implemented
  } else {
    k_body.template modify<DeviceType>();
    k_body.sync_host();
  }

  base()->commflag = FORCE_TORQUE;
  base()->comm->reverse_comm(fix_base(), 6);

  if (base()->langflag) {
    k_body.sync_device();
    auto l_body = k_body.template view<DeviceType>();
    auto l_langextra = k_langextra.template view<DeviceType>();
    base()->copymode = 1;
    Kokkos::parallel_for(
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
    base()->copymode = 0;
  }
  if (base()->id_gravity) {
    k_body.sync_device();
    auto l_body = k_body.template view<DeviceType>();
    auto l_gvec0 = base()->gvec[0];
    auto l_gvec1 = base()->gvec[1];
    auto l_gvec2 = base()->gvec[2];
    base()->copymode = 1;
    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType>(0, nbody_total()),
      KOKKOS_LAMBDA(const int &ibody) {
        BodyKokkos &bk = l_body(ibody);
        bk.fcm[0] += l_gvec0 * bk.mass;
        bk.fcm[1] += l_gvec1 * bk.mass;
        bk.fcm[2] += l_gvec2 * bk.mass;
      }
    );
    base()->copymode = 0;
  }
  if (base()->langflag || base()->id_gravity) k_body.modify_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::enforce2d_base()
{
  k_body.sync_device();
  auto l_body = k_body.template view<DeviceType>();
  base()->copymode = 1;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<DeviceType>(0, base()->nlocal_body),
    KOKKOS_LAMBDA(const int &ibody) {
      BodyKokkos &bk = l_body(ibody);
      bk.xcm[2] = bk.vcm[2] = bk.fcm[2] = bk.xgc[2] = 0.0;
      bk.torque[0] = bk.torque[1] = 0.0;
      bk.angmom[0] = bk.angmom[1] = 0.0;
      bk.omega[0] = bk.omega[1] = 0.0;
    }
  );
  base()->copymode = 0;
  k_body.modify_device();
}

/* ----------------------------------------------------------------------
   count # of DOF removed by rigid bodies for atoms in igroup
   return total count of DOF
------------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
bigint FixRigidBaseKokkos<DeviceType,FixRigidBase>::dof_base(int tgroup)
{

  // cannot count DOF correctly unless setup_bodies_static() has been called
  if (!base()->setupflag) {
    if (base()->comm->me == 0)
      base()->error->warning(FLERR,"Cannot count rigid body degrees-of-freedom "
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

  k_counts = DAT::tdual_int_2d("rigid/small:counts", nbody_total(), 3);
  auto l_counts = k_counts.template view<DeviceType>();
  deep_copy(k_counts.view_device(), 0);

  // tally counts from my owned atoms
  // 0 = # of point particles in rigid body and in temperature group
  // 1 = # of finite-size particles in rigid body and in temperature group
  // 2 = # of particles in rigid body, disregarding temperature group

  int l_tgroupbit = base()->group->bitmask[tgroup];
  auto l_extended = base()->extended;

  Kokkos::parallel_for(
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

  base()->commflag = DOF;
  base()->comm->reverse_comm(fix_base(), 3);

  // nall = count0 = # of point particles in each rigid body
  // mall = count1 = # of finite-size particles in each rigid body
  // warn if nall+mall != nrigid for any body included in temperature group

  int flag = 0;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<DeviceType>(0, nlocal_body()),
    KOKKOS_LAMBDA(const int ibody, int &l_flag) {
      const int l_counts01 = l_counts(ibody,0) + l_counts(ibody,1);
      if (l_counts01 > 0 && l_counts01 != l_counts(ibody,2)) l_flag = 1;
    }, Kokkos::Max<int>(flag)
  );
  int flag_all;
  MPI_Allreduce(&flag, &flag_all, 1, MPI_INT, MPI_MAX, world());
  if (flag_all && base()->comm->me == 0)
    base()->error->warning(FLERR,"Computing temperature of portions of rigid bodies");

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
  base()->nlinear = 0;
  if (domainKK->dimension == 3) {
    Kokkos::parallel_reduce(
      Kokkos::RangePolicy<DeviceType>(0, nlocal_body()),
      KOKKOS_LAMBDA(const int ibody, bigint &l_n, int &l_nlinear) {
        if (l_counts(ibody,0) + l_counts(ibody,1) == l_counts(ibody,2)) {
          l_n += 3*l_counts(ibody,0) + 6*l_counts(ibody,1) - 6;
          auto inertia = l_body(ibody).inertia;
          if (inertia[0] == 0.0 || inertia[1] == 0.0 || inertia[2] == 0.0) {
            l_n++;
            l_nlinear++;
          }
        }
      }, n, base()->nlinear
    );
  } else if (domainKK->dimension == 2) {
    Kokkos::parallel_reduce(
      Kokkos::RangePolicy<DeviceType>(0, nlocal_body()),
      KOKKOS_LAMBDA(const int ibody, bigint &l_n) {
        if (l_counts(ibody,0) + l_counts(ibody,1) == l_counts(ibody,2)) {
          l_n += 2*l_counts(ibody,0) + 3*l_counts(ibody,1) - 3;
        }
      }, n
    );
  }

  bigint n_all;
  MPI_Allreduce(&n, &n_all, 1, MPI_LMP_BIGINT, MPI_SUM, base()->world);
  return n_all;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::image_shift_base()
{
  base()->copymode = 1;
  atomKK->sync(execution_space, IMAGE_MASK);
  k_atom2body.template sync<DeviceType>();
  k_xcmimage.template sync<DeviceType>();
  k_body.sync_device();
  auto l_image = atomKK->k_image.template view<DeviceType>();
  auto l_atom2body = k_atom2body.template view<DeviceType>();
  auto l_xcmimage = k_xcmimage.template view<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();
  Kokkos::parallel_for(
    Kokkos::RangePolicy<DeviceType>(0, atomKK->nlocal),
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
      Kokkos::printf("*** ibody %i l_image(%i) %i bk.image %i l_xcmimage(%i) %i\n", ibody, i, l_image(i), bk.image, i, l_xcmimage(i));
    }
  );
  k_xcmimage.template modify<DeviceType>();
  base()->copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::reset_atom2body_base()
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
  auto l_comm_me = base()->comm->me;
  auto l_ntimestep = base()->update->ntimestep;

  base()->copymode = 1;
  Kokkos::parallel_for(
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
  base()->copymode = 0;
  k_atom2body.modify_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::grow_arrays_base(int nmax)
{
  // Do not sync from device before grow: uninitialized device data must not
  // overwrite the host arrays this routine is about to resize.

  base()->memoryKK->grow_kokkos(k_bodyown, base()->bodyown, nmax, "rigid/small:bodyown");
  base()->memoryKK->grow_kokkos(k_bodytag, base()->bodytag, nmax, "rigid/small:bodytag");
  base()->memoryKK->grow_kokkos(k_atom2body, base()->atom2body, nmax, "rigid/small:atom2body");
  base()->memoryKK->grow_kokkos(k_xcmimage, base()->xcmimage, nmax, "rigid/small:xcmimage");

  // k_displace is a TransformView: grow_kokkos cannot maintain the displace[i]
  // pointer array into the view's host buffer, so resize explicitly and
  // re-point each displace[i] into the new host allocation.
  k_displace.sync_host();
  k_displace.resize(nmax, 3);
  bigint nbytes = ((bigint) sizeof(double *)) * nmax;
  base()->displace = (double **) base()->memory->srealloc(base()->displace, nbytes, "rigid/small:displace");
  for (int i = 0; i < nmax; i++)
    base()->displace[i] = (k_displace.extent_int(1) > 0) ? &k_displace.view_host()(i, 0) : nullptr;
  k_displace.modify_host();

  if (base()->extended) {
    k_eflags.template sync<DeviceType>();
    base()->memoryKK->grow_kokkos(k_eflags, base()->eflags, nmax, "rigid/small:eflags");
    if (base()->orientflag) base()->memory->grow(base()->orient, nmax, base()->orientflag, "rigid/small:orient");
    if (base()->dorientflag) base()->memory->grow(base()->dorient, nmax, 3, "rigid/small:dorient");
  }

  if (nmax > base()->maxvatom) {
    base()->maxvatom = nmax;
    base()->memoryKK->destroy_kokkos(k_vatom, base()->vatom);
    base()->memoryKK->create_kokkos(k_vatom, base()->vatom, base()->maxvatom, 6, "fix:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::grow_body_base()
{
  grow_body_base(base()->nmax_body + DELTA_BODY);
}

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::grow_body_base(int nmax_body)
{
  base()->nmax_body = nmax_body;
  k_body.resize(nmax_body);
  memcpy(k_body.view_host().data(), base()->body, (bigint) nmax_body * sizeof(Body));
  base()->memory->sfree(base()->body);
  base()->body = k_body.view_host().data();
  k_body.modify_host();
  k_body.sync_device();
}

/* ----------------------------------------------------------------------
  KOKKOS BASE
------------------------------------------------------------------------- */

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

  int result = 0;
  base()->copymode = 1;
  if (base()->commflag == FULL_BODY) {
    Kokkos::parallel_scan(
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(l_sendlist(i));
        if (!final) {
          if (ibody < 0) m++;
          else m += 48;
          return;
        }
        if (ibody < 0) {
          l_buf[m++] = d_ubuf(0).d; // flag
          return;
        }
        const BodyKokkos &bk = l_body(ibody);
        l_buf[m++] = d_ubuf(1).d; // flag
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
  } else if (base()->commflag == INITIAL) {
    Kokkos::parallel_scan(
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(l_sendlist(i));
        if (ibody < 0) return;
        if (!final) {
          m += 29;
          return;
        }
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
      }, result
    );
  } else if (base()->commflag == FINAL) {
    Kokkos::parallel_scan(
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(l_sendlist(i));
        if (ibody < 0) return;
        if (!final) {
          m += 10;
          return;
        }
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
      }, result
    );
  }
  base()->copymode = 0;
  k_buf.modify_device();
  return result;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::unpack_forward_comm_kokkos(
    int n, int first, DAT::tdual_double_1d &k_buf)
{

  // grow body views if needed
  // receiving n atoms so at most n bodies

  if (base()->commflag == FULL_BODY && nlocal_body() + n > nmax_body())
    grow_body_base(nlocal_body() + n);

  // kokkos views

  k_buf.template sync<DeviceType>();
  auto l_buf = k_buf.template view<DeviceType>();
  
  k_bodyown.template sync<DeviceType>();
  auto l_bodyown = k_bodyown.template view<DeviceType>();

  k_body.template sync<DeviceType>();
  auto l_body = k_body.template view<DeviceType>();

  base()->copymode = 1;
  if (base()->commflag == FULL_BODY) {

    auto l_nlocal_body = nlocal_body();
    Few<int,2> recv(0,0);

    Kokkos::parallel_scan(
      Kokkos::RangePolicy<DeviceType>(first, first+n),
      KOKKOS_LAMBDA(const int &i, Few<int,2> &l_recv, const bool &final) {
        const int flag = d_ubuf(l_buf[l_recv[0]++]).i;
        if (!final) {
          if (flag == 1) { l_recv[0] += 47; l_recv[1]++; }
          return;
        }
        if (flag == 0) {
          l_bodyown(i) = -1;
          return;
        }
        const int j = l_nlocal_body + (l_recv[1]++);
        int m = l_recv[0];
        l_body(j).ilocal = i;
        l_bodyown(i) = j;
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
        l_recv[0] += 47;

        Kokkos::printf("*** unpack_forward_comm_kokkos i %i j %i mass %f natoms %i xcm %f %f %f\n", i, j, bk.mass, bk.natoms, bk.xcm[0], bk.xcm[1], bk.xcm[2]);

      }, recv
    );

    base()->nghost_body += recv[1];
    k_bodyown.modify_device();

  } else if (base()->commflag == INITIAL) {
    Kokkos::parallel_scan(
      Kokkos::RangePolicy<DeviceType>(first, first+n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(i);
        if (ibody < 0) return;
        if (!final) {
          m += 29;
          return;
        }
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
      }
    );
  } else if (base()->commflag == FINAL) {
    Kokkos::parallel_scan(
      Kokkos::RangePolicy<DeviceType>(first, first+n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(i);
        if (ibody < 0) return;
        if (!final) {
          m += 10;
          return;
        }
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
      }
    );
  }
  base()->copymode = 0;
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
  base()->copymode = 1;
  if (base()->commflag == FORCE_TORQUE) {

    // kokkos views
    k_body.template sync<DeviceType>();
    auto l_body = k_body.template view<DeviceType>();

    Kokkos::parallel_scan(
      Kokkos::RangePolicy<DeviceType>(first, first+n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(i);
        if (ibody < 0) return;
        if (!final) {
          m += 6;
          return;
        }
        const BodyKokkos &bk = l_body(ibody);
        l_buf[m++] = static_cast<double>(bk.fcm[0]);
        l_buf[m++] = static_cast<double>(bk.fcm[1]);
        l_buf[m++] = static_cast<double>(bk.fcm[2]);
        l_buf[m++] = static_cast<double>(bk.torque[0]);
        l_buf[m++] = static_cast<double>(bk.torque[1]);
        l_buf[m++] = static_cast<double>(bk.torque[2]);
      }, result
    );

  } else if (base()->commflag == VCM_ANGMOM) {

    // kokkos views
    k_body.template sync<DeviceType>();
    auto l_body = k_body.template view<DeviceType>();

    Kokkos::parallel_scan(
      Kokkos::RangePolicy<DeviceType>(first, first+n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(i);
        if (ibody < 0) return;
        if (!final) {
          m += 6;
          return;
        }
        const BodyKokkos &bk = l_body(ibody);
        l_buf[m++] = static_cast<double>(bk.vcm[0]);
        l_buf[m++] = static_cast<double>(bk.vcm[1]);
        l_buf[m++] = static_cast<double>(bk.vcm[2]);
        l_buf[m++] = static_cast<double>(bk.angmom[0]);
        l_buf[m++] = static_cast<double>(bk.angmom[1]);
        l_buf[m++] = static_cast<double>(bk.angmom[2]);
      }, result
    );

  } else if (base()->commflag == XCM_MASS) {

    // kokkos views
    k_body.template sync<DeviceType>();
    auto l_body = k_body.template view<DeviceType>();

    Kokkos::parallel_scan(
      Kokkos::RangePolicy<DeviceType>(first, first+n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(i);
        if (ibody < 0) return;
        if (!final) {
          m += 8;
          return;
        }
        const BodyKokkos &bk = l_body(ibody);
        l_buf[m++] = static_cast<double>(bk.xcm[0]);
        l_buf[m++] = static_cast<double>(bk.xcm[1]);
        l_buf[m++] = static_cast<double>(bk.xcm[2]);
        l_buf[m++] = static_cast<double>(bk.xgc[0]);
        l_buf[m++] = static_cast<double>(bk.xgc[1]);
        l_buf[m++] = static_cast<double>(bk.xgc[2]);
        l_buf[m++] = static_cast<double>(bk.mass);
        l_buf[m++] = static_cast<double>(bk.natoms);
        Kokkos::printf("*** pack_reverse_comm_kokkos i %i ibody %i bk.xcm %f %f %f natoms %i\n", i, ibody, bk.xcm[0], bk.xcm[1], bk.xcm[2], bk.natoms);
      }, result
    );

  } else if (base()->commflag == ITENSOR) {

    // kokkos views
    k_itensor.template sync<DeviceType>();
    auto l_itensor = k_itensor.template view<DeviceType>();

    Kokkos::parallel_scan(
      Kokkos::RangePolicy<DeviceType>(first, first+n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(i);
        if (ibody < 0) return;
        if (!final) {
          m += 6;
          return;
        }
        l_buf[m++] = static_cast<double>(l_itensor(ibody,0));
        l_buf[m++] = static_cast<double>(l_itensor(ibody,1));
        l_buf[m++] = static_cast<double>(l_itensor(ibody,2));
        l_buf[m++] = static_cast<double>(l_itensor(ibody,3));
        l_buf[m++] = static_cast<double>(l_itensor(ibody,4));
        l_buf[m++] = static_cast<double>(l_itensor(ibody,5));
      }, result
    );

  }

  base()->copymode = 0;
  k_buf.modify_device();
  return result;

}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::unpack_reverse_comm_kokkos(
    int n, DAT::tdual_int_1d k_sendlist, DAT::tdual_double_1d &k_buf)
{

  //utils::logmesg(base()->lmp, "*** unpack_reverse_comm_kokkos() commflag {}\n", base()->commflag);

  // kokkos views

  k_sendlist.template sync<DeviceType>();
  auto l_sendlist = k_sendlist.template view<DeviceType>();

  k_buf.template sync<DeviceType>();
  auto l_buf = k_buf.template view<DeviceType>();

  k_bodyown.template sync<DeviceType>();
  auto l_bodyown = k_bodyown.template view<DeviceType>();

  base()->copymode = 1;
  if (base()->commflag == FORCE_TORQUE) {

    // kokkos views
    k_body.template sync<DeviceType>();
    auto l_body = k_body.template view<DeviceType>();

    Kokkos::parallel_scan(
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(l_sendlist(i));
        if (ibody < 0) return;
        if (!final) {
          m += 6;
          return;
        }
        BodyKokkos &bk = l_body(ibody);
        Kokkos::atomic_add(&bk.fcm[0], static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&bk.fcm[1], static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&bk.fcm[2], static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&bk.torque[0], static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&bk.torque[1], static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&bk.torque[2], static_cast<KK_FLOAT>(l_buf[m++]));
      }
    );
    k_body.modify_device();

  } else if (base()->commflag == VCM_ANGMOM) {

    // kokkos views
    k_body.template sync<DeviceType>();
    auto l_body = k_body.template view<DeviceType>();

    Kokkos::parallel_scan(
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
    k_body.modify_device();
  } else if (base()->commflag == XCM_MASS) {

    // kokkos views
    k_body.template sync<DeviceType>();
    auto l_body = k_body.template view<DeviceType>();

    Kokkos::parallel_scan(
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(l_sendlist(i));
        if (ibody < 0) return;
        if (!final) {
          m += 8;
          return;
        }
        BodyKokkos &bk = l_body(ibody);
        Kokkos::atomic_add(&bk.xcm[0], static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&bk.xcm[1], static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&bk.xcm[2], static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&bk.xgc[0], static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&bk.xgc[1], static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&bk.xgc[2], static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&bk.mass, static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&bk.natoms, static_cast<int>(l_buf[m++]));
        Kokkos::printf("*** unpack_reverse_comm_kokkos i %i ibody %i mass %f natoms %i xcm %f %f %f\n", i, ibody, bk.mass, bk.natoms, bk.xcm[0], bk.xcm[1], bk.xcm[2]);
      }
    );
    k_body.modify_device();

  } else if (base()->commflag == ITENSOR) {

    // kokkos views
    k_itensor.template sync<DeviceType>();
    auto l_itensor = k_itensor.template view<DeviceType>();

    Kokkos::parallel_scan(
      Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int &i, int &m, const bool &final) {
        const int ibody = l_bodyown(l_sendlist(i));
        if (ibody < 0) return;
        if (!final) {
          m += 6;
          return;
        }
        Kokkos::atomic_add(&l_itensor(ibody,0), static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&l_itensor(ibody,1), static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&l_itensor(ibody,2), static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&l_itensor(ibody,3), static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&l_itensor(ibody,4), static_cast<KK_FLOAT>(l_buf[m++]));
        Kokkos::atomic_add(&l_itensor(ibody,5), static_cast<KK_FLOAT>(l_buf[m++]));
      }
    );
    k_itensor.modify_device();

  }
  base()->copymode = 0;
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
  auto l_vflag_atom = base()->vflag_atom;
  int l_bodysize_kk = sizeof(BodyKokkos)/sizeof(double);
  if (l_bodysize_kk * sizeof(double) != sizeof(BodyKokkos)) l_bodysize_kk++;

  base()->copymode = 1;
  Kokkos::parallel_scan(
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
  base()->copymode = 0;
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
  auto l_vflag_atom = base()->vflag_atom;
  int l_bodysize_kk = sizeof(BodyKokkos)/sizeof(double);
  if (l_bodysize_kk * sizeof(double) != sizeof(BodyKokkos)) l_bodysize_kk++;
  auto l_nlocal_body = base()->nlocal_body;
  int nbody_recv = 0;

  base()->copymode = 1;
  Kokkos::parallel_scan(
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
  base()->nlocal_body += nbody_recv;
  base()->copymode = 0;
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
  if (base()->extended) k_eflags.sync_device();

  Sorter.sort(LMPDeviceType(), k_body.view_device());
  Sorter.sort(LMPDeviceType(), k_bodytag.view_device());
  Sorter.sort(LMPDeviceType(), k_xcmimage.view_device());
  Sorter.sort(LMPDeviceType(), k_displace.view_device());
  Sorter.sort(LMPDeviceType(), k_bodyown.view_device());
  if (base()->extended) Sorter.sort(LMPDeviceType(), k_eflags.view_device());

  k_body.modify_device();
  k_bodytag.modify_device();
  k_xcmimage.modify_device();
  k_displace.modify_device();
  k_bodyown.modify_device();
  if (base()->extended) k_eflags.modify_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::modify_host_base()
{
  k_body.modify_host();
  k_bodyown.modify_host();
  k_bodytag.modify_host();
  k_atom2body.modify_host();
  k_xcmimage.modify_host();
  k_displace.modify_host();
  if (base()->extended) k_eflags.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::modify_device_base()
{
  k_body.modify_device();
  k_bodyown.modify_device();
  k_bodytag.modify_device();
  k_atom2body.modify_device();
  k_xcmimage.modify_device();
  k_displace.modify_device();
  if (base()->extended) k_eflags.modify_device();
}


/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::sync_host_base()
{
  k_body.sync_host();
  k_bodyown.sync_host();
  k_bodytag.sync_host();
  k_atom2body.sync_host();
  k_xcmimage.sync_host();
  k_displace.sync_host();
  if (base()->extended) k_eflags.sync_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class FixRigidBase>
void FixRigidBaseKokkos<DeviceType,FixRigidBase>::sync_device_base()
{
  k_body.sync_device();
  k_bodyown.sync_device();
  k_bodytag.sync_device();
  k_atom2body.sync_device();
  k_xcmimage.sync_device();
  k_displace.sync_device();
  if (base()->extended) k_eflags.sync_device();
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixRigidBaseKokkos<LMPDeviceType, FixRigidSmall>;
template class FixRigidBaseKokkos<LMPDeviceType, FixRigidNHSmall>;
#ifdef LMP_KOKKOS_GPU
//template class FixRigidBaseKokkos<LMPHostType, FixRigidSmall>;
//template class FixRigidBaseKokkos<LMPHostType, FixRigidNHSmall>;
#endif
}

