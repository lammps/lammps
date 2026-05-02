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
   Contributing authors: Trung Nguyen (U Chicago)
                         Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#include "fix_efield_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain_kokkos.h"
#include "error.h"
#include "input.h"
#include "kokkos_base.h"
#include "memory_kokkos.h"
#include "modify_kokkos.h"
#include "region.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

template<class DeviceType, bool TIP4P>
FixEfieldKokkos<DeviceType,TIP4P>::FixEfieldKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixEfield(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TORQUE_MASK | Q_MASK | MU_MASK | IMAGE_MASK | MASK_MASK;
  datamask_modify = F_MASK | TORQUE_MASK;

  memory->destroy(efield);
  memoryKK->create_kokkos(k_efield,efield,maxatom,4,"efield:efield");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, bool TIP4P>
FixEfieldKokkos<DeviceType,TIP4P>::~FixEfieldKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_efield,efield);
  memoryKK->destroy_kokkos(k_vatom,vatom);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, bool TIP4P>
void FixEfieldKokkos<DeviceType,TIP4P>::init()
{
  FixEfield::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with Kokkos");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, bool TIP4P>
void FixEfieldKokkos<DeviceType,TIP4P>::post_force(int vflag)
{

  force_flag = 0;
  const int nlocal = atomKK->nlocal;

  // virial setup

  v_init(vflag);

  // reallocate per-atom arrays if necessary

  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"efield:vatom");
  }

  // update region if necessary

  typename AT::t_int_1d l_match;
  if (region) {
    if (!(utils::strmatch(region->style, "^block") || utils::strmatch(region->style, "^sphere")))
      error->all(FLERR,"Cannot (yet) use {}-style region with fix efield/kk",region->style);
    region->prematch();
    DAT::tdual_int_1d k_match = DAT::tdual_int_1d("efield:k_match",nlocal);
    KokkosBase* regionKKBase = dynamic_cast<KokkosBase*>(region);
    regionKKBase->match_all_kokkos(groupbit,k_match);
    k_match.template sync<DeviceType>();
    l_match = k_match.template view<DeviceType>();
  }

  // reallocate sforce array if necessary

  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memoryKK->destroy_kokkos(k_efield,efield);
    memoryKK->create_kokkos(k_efield,efield,maxatom,4,"efield:efield");
  }

  // kokkos views

  atomKK->sync(execution_space, datamask_read);
  auto l_x = atomKK->k_x.template view<DeviceType>();
  auto l_f = atomKK->k_f.template view<DeviceType>();
  auto l_q = atomKK->k_q.template view<DeviceType>();
  auto l_mu = atomKK->k_mu.template view<DeviceType>();
  auto l_torque = atomKK->k_torque.template view<DeviceType>();
  auto l_image = atomKK->k_image.template view<DeviceType>();
  auto l_mask = atomKK->k_mask.template view<DeviceType>();
  auto l_efield = k_efield.template view<DeviceType>();
  auto l_vatom = k_vatom.template view<DeviceType>();

  // domainKK

  auto l_prd = Few<KK_FLOAT,3>(domain->prd);
  auto l_h = Few<KK_FLOAT,6>(domain->h);
  auto l_triclinic = domain->triclinic;

  auto lambda = [&]<bool QFLAG, bool MUFLAG, bool CONSTANT_FLAG>(auto& result_) {

    auto l_groupbit = groupbit;
    auto l_region = region;
    auto l_ex = static_cast<KK_FLOAT>(ex);
    auto l_ey = static_cast<KK_FLOAT>(ey);
    auto l_ez = static_cast<KK_FLOAT>(ez);
    auto l_qe2f = static_cast<KK_FLOAT>(qe2f);
    auto l_xstyle = xstyle;
    auto l_ystyle = ystyle;
    auto l_zstyle = zstyle;
    auto l_pstyle = pstyle;
    auto l_estyle = estyle;
    auto l_evflag = evflag;
    auto l_vflag_global = vflag_global;
    auto l_vflag_atom = vflag_atom;

    copymode = 1;
    Kokkos::parallel_reduce(
      Kokkos::RangePolicy<DeviceType>(0, nlocal),
      KOKKOS_LAMBDA(const int &i, auto& l_result) {

        if( !(l_mask(i) & l_groupbit) || (l_region && !l_match(i))) return;

        if constexpr (QFLAG) {
          Few<KK_FLOAT,3> xi(l_x(i,0),l_x(i,1),l_x(i,2));
          auto u = DomainKokkos::unmap(l_prd, l_h, l_triclinic, xi, l_image(i));
          KK_FLOAT fx, fy, fz;

          if constexpr (CONSTANT_FLAG) {
            fx = l_q(i) * l_ex;
            fy = l_q(i) * l_ey;
            fz = l_q(i) * l_ez;
          } else {
            if (l_xstyle == ATOM) fx = l_qe2f * l_q(i) * l_efield(i,0);
            else fx = l_q(i) * l_ex;
            if (l_ystyle == ATOM) fy = l_qe2f * l_q(i) * l_efield(i,1);
            else fy = l_q(i) * l_ey;
            if (l_zstyle == ATOM) fz = l_qe2f * l_q(i) * l_efield(i,2);
            else fz = l_q(i) * l_ez;
          }

          l_f(i,0) += fx;
          l_f(i,1) += fy;
          l_f(i,2) += fz;

          if constexpr (CONSTANT_FLAG) {
            l_result[0] -= fma(fx, u[0], fma(fy, u[1], fz * u[2]));
          } else {
            if (l_pstyle == ATOM) l_result[0] += l_qe2f * l_q(i) * l_efield(i,3);
            else if (l_estyle == ATOM) l_result[0] += l_efield(i,3);
          }

          l_result[1] += fx;
          l_result[2] += fy;
          l_result[3] += fz;

          if constexpr (CONSTANT_FLAG) {
          if (l_evflag) {

            // tally virial into global and per-atom accumulators
            // only when field is CONSTANT
            // i = local index of atom, v = total virial for the interaction
            // this method can be used when fix computes forces in post_force()
            // and the force depends on a distance to some external object
            // eg. fix wall/lj93: compute virial only on owned atoms

            KK_ACC_FLOAT v[6] = {fx*u[0], fy*u[1], fz*u[2], fx*u[1], fx*u[2], fy*u[2]};

            if ( l_vflag_global ) {
              // increment global virial by v
              l_result[4] += v[0];
              l_result[5] += v[1];
              l_result[6] += v[2];
              l_result[7] += v[3];
              l_result[8] += v[4];
              l_result[9] += v[5];
            }

            if (l_vflag_atom) {
              // increment per-atom virial by v
              Kokkos::atomic_add(&(l_vatom(i,0)), v[0]);
              Kokkos::atomic_add(&(l_vatom(i,1)), v[1]);
              Kokkos::atomic_add(&(l_vatom(i,2)), v[2]);
              Kokkos::atomic_add(&(l_vatom(i,3)), v[3]);
              Kokkos::atomic_add(&(l_vatom(i,4)), v[4]);
              Kokkos::atomic_add(&(l_vatom(i,5)), v[5]);
            }
          }
          }
        }

        if constexpr (MUFLAG) {
          l_torque(i,0) += l_ez * l_mu(i,1) - l_ey * l_mu(i,2);
          l_torque(i,1) += l_ex * l_mu(i,2) - l_ez * l_mu(i,0);
          l_torque(i,2) += l_ey * l_mu(i,0) - l_ex * l_mu(i,1);
          if constexpr (CONSTANT_FLAG)
            l_result[0] -= l_mu(i,0) * l_ex + l_mu(i,1) * l_ey + l_mu(i,2) * l_ez;
        }

      }, result_
    );
    copymode = 0;
  };

  Few<KK_ACC_FLOAT, 10> result = {0.0};

  if (varflag == CONSTANT) {

    if (qflag && muflag)      lambda.template operator()<true,true,true>(result);
    else if(qflag && !muflag) lambda.template operator()<true,false,true>(result);
    else if(!qflag && muflag) lambda.template operator()<false,true,true>(result);
    else                      lambda.template operator()<false,false,true>(result);

  } else {

    update_efield_variables();

    if(qflag && muflag)       lambda.template operator()<true,true,false>(result);
    else if(qflag && !muflag) lambda.template operator()<true,false,false>(result);
    else if(!qflag && muflag) lambda.template operator()<false,true,false>(result);
    else                      lambda.template operator()<false,false,false>(result);

  }

  atomKK->modified(execution_space, datamask_modify);

  fsum[0] = result[0];
  fsum[1] = result[1];
  fsum[2] = result[2];
  fsum[3] = result[3];

  if (vflag_global) {
    virial[0] += result[4];
    virial[1] += result[5];
    virial[2] += result[6];
    virial[3] += result[7];
    virial[4] += result[8];
    virial[5] += result[9];
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }
}


namespace LAMMPS_NS {
template class FixEfieldKokkos<LMPDeviceType,false>;
template class FixEfieldKokkos<LMPDeviceType,true>;
#ifdef LMP_KOKKOS_GPU
//template class FixEfieldKokkos<LMPHostType>;
#endif
}
