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

#include "pair_gran_hooke_history_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "fix_neigh_history_kokkos.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairGranHookeHistoryKokkos<DeviceType>::PairGranHookeHistoryKokkos(LAMMPS *lmp) : PairGranHookeHistory(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | OMEGA_MASK | F_MASK | TORQUE_MASK | TYPE_MASK | MASK_MASK | ENERGY_MASK | VIRIAL_MASK | RMASS_MASK | RADIUS_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairGranHookeHistoryKokkos<DeviceType>::~PairGranHookeHistoryKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    eatom = nullptr;
    vatom = nullptr;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairGranHookeHistoryKokkos<DeviceType>::init_style()
{
  // if history is stored and first init, create Fix to store history
  // it replaces FixDummy, created in the constructor
  // this is so its order in the fix list is preserved

  if (history && fix_history == nullptr) {
    auto cmd = std::string("NEIGH_HISTORY_HH") + std::to_string(instance_me) + " all ";
    if (execution_space == Device)
      cmd += "NEIGH_HISTORY/KK/DEVICE 3";
    else
      cmd += "NEIGH_HISTORY/KK/HOST 3";
    fix_history = (FixNeighHistory *)
      modify->replace_fix("NEIGH_HISTORY_HH_DUMMY"+std::to_string(instance_me),cmd,1);
    fix_history->pair = this;
    fix_historyKK = (FixNeighHistoryKokkos<DeviceType> *)fix_history;
  }

  PairGranHookeHistory::init_style();

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL)
    error->all(FLERR,"Must use half neighbor list with gran/hooke/history/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairGranHookeHistoryKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  copymode = 1;

  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  int shearupdate = 1;
  if (update->setupflag) shearupdate = 0;

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK | TORQUE_MASK);

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  omega = atomKK->k_omega.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];

  int inum = list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  kt_kk = static_cast<KK_FLOAT>(kt);
  kn_kk = static_cast<KK_FLOAT>(kn);
  xmu_kk = static_cast<KK_FLOAT>(xmu);
  gammat_kk = static_cast<KK_FLOAT>(gammat);
  gamman_kk = static_cast<KK_FLOAT>(gamman);
  dt_kk = static_cast<KK_FLOAT>(dt);

  if (d_numneigh.extent(0) != d_numneigh_touch.extent(0))
    d_numneigh_touch = typename AT::t_int_1d("pair:numneigh_touch",d_numneigh.extent(0));
  if (d_neighbors.extent(0) != d_neighbors_touch.extent(0) ||
      d_neighbors.extent(1) != d_neighbors_touch.extent(1))
    d_neighbors_touch = typename AT::t_neighbors_2d("pair:neighbors_touch",d_neighbors.extent(0),d_neighbors.extent(1));

  fix_historyKK->k_firstflag.template sync<DeviceType>();
  fix_historyKK->k_firstvalue.template sync<DeviceType>();

  d_firsttouch = fix_historyKK->k_firstflag.template view<DeviceType>();
  d_firstshear = fix_historyKK->k_firstvalue.template view<DeviceType>();

  Kokkos::deep_copy(d_firsttouch,0);

  EV_FLOAT ev;

  if (neighflag == HALF) {
    if (force->newton_pair) {
      if (vflag_either) {
        if (shearupdate) {
          Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,1,1,1>>(0,inum),*this, ev);
        } else {
          Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,1,1,0>>(0,inum),*this, ev);
        }
      } else {
        if (shearupdate) {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,1,0,1>>(0,inum),*this);
        } else {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,1,0,0>>(0,inum),*this);
        }
      }
    } else {
      if (vflag_either) {
        if (shearupdate) {
          Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,0,1,1>>(0,inum),*this, ev);
        } else {
          Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,0,1,0>>(0,inum),*this, ev);
        }
      } else {
        if (shearupdate) {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,0,0,1>>(0,inum),*this);
        } else {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,0,0,0>>(0,inum),*this);
        }
      }
    }
  } else { // HALFTHREAD
    if (force->newton_pair) {
      if (vflag_either) {
        if (shearupdate) {
          Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,1,1,1>>(0,inum),*this, ev);
        } else {
          Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,1,1,0>>(0,inum),*this, ev);
        }
      } else {
        if (shearupdate) {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,1,0,1>>(0,inum),*this);
        } else {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,1,0,0>>(0,inum),*this);
        }
      }
    } else {
      if (vflag_either) {
        if (shearupdate) {
          Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,0,1,1>>(0,inum),*this, ev);
        } else {
          Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,0,1,0>>(0,inum),*this, ev);
        }
      } else {
        if (shearupdate) {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,0,0,1>>(0,inum),*this);
        } else {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,0,0,0>>(0,inum),*this);
        }
      }
    }
  }

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }

  if (vflag_global) {
    virial[0] += static_cast<double>(ev.v[0]);
    virial[1] += static_cast<double>(ev.v[1]);
    virial[2] += static_cast<double>(ev.v[2]);
    virial[3] += static_cast<double>(ev.v[3]);
    virial[4] += static_cast<double>(ev.v[4]);
    virial[5] += static_cast<double>(ev.v[5]);
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int SHEARUPDATE>
KOKKOS_INLINE_FUNCTION
void PairGranHookeHistoryKokkos<DeviceType>::operator()(TagPairGranHookeHistoryCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,SHEARUPDATE>, const int ii, EV_FLOAT &ev) const {

  // The f and torque arrays are atomic for Half/Thread neighbor style
  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;
  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_torque = torque;

  const int i = d_ilist[ii];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);
  const KK_FLOAT imass = rmass[i];
  const KK_FLOAT irad = radius[i];
  const int jnum = d_numneigh[i];
  const int mask_i = mask[i];

  const KK_FLOAT vx_i = v(i,0);
  const KK_FLOAT vy_i = v(i,1);
  const KK_FLOAT vz_i = v(i,2);

  const KK_FLOAT omegax_i = omega(i,0);
  const KK_FLOAT omegay_i = omega(i,1);
  const KK_FLOAT omegaz_i = omega(i,2);

  KK_ACC_FLOAT fx_i = 0.0;
  KK_ACC_FLOAT fy_i = 0.0;
  KK_ACC_FLOAT fz_i = 0.0;

  KK_ACC_FLOAT torquex_i = 0.0;
  KK_ACC_FLOAT torquey_i = 0.0;
  KK_ACC_FLOAT torquez_i = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    KK_FLOAT factor_lj = static_cast<KK_FLOAT>(special_lj[sbmask(j)]);
    j &= NEIGHMASK;

    if (factor_lj == 0) continue;

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    const KK_FLOAT jmass = rmass[j];
    const KK_FLOAT jrad = radius[j];
    const KK_FLOAT radsum = irad + jrad;

    // check for touching neighbors

    if (rsq >= radsum * radsum) {
      d_firstshear(i,3*jj) = 0;
      d_firstshear(i,3*jj+1) = 0;
      d_firstshear(i,3*jj+2) = 0;
      continue;
    }

    d_firsttouch(i,jj) = 1;

    const KK_FLOAT r = sqrt(rsq);
    const KK_FLOAT rinv = static_cast<KK_FLOAT>(1.0)/r;
    const KK_FLOAT rsqinv = static_cast<KK_FLOAT>(1.0)/rsq;

    // relative translational velocity

    KK_FLOAT vr1 = vx_i - v(j,0);
    KK_FLOAT vr2 = vy_i - v(j,1);
    KK_FLOAT vr3 = vz_i - v(j,2);

    // normal component

    KK_FLOAT vnnr = vr1*delx + vr2*dely + vr3*delz;
    KK_FLOAT vn1 = delx*vnnr * rsqinv;
    KK_FLOAT vn2 = dely*vnnr * rsqinv;
    KK_FLOAT vn3 = delz*vnnr * rsqinv;

    // tangential component

    KK_FLOAT vt1 = vr1 - vn1;
    KK_FLOAT vt2 = vr2 - vn2;
    KK_FLOAT vt3 = vr3 - vn3;

    // relative rotational velocity

    KK_FLOAT wr1 = (irad*omegax_i + jrad*omega(j,0)) * rinv;
    KK_FLOAT wr2 = (irad*omegay_i + jrad*omega(j,1)) * rinv;
    KK_FLOAT wr3 = (irad*omegaz_i + jrad*omega(j,2)) * rinv;

    KK_FLOAT meff = imass*jmass / (imass+jmass);
    if (mask_i & freeze_group_bit) meff = jmass;
    if (mask[j] & freeze_group_bit) meff = imass;

    KK_FLOAT damp = meff*gamman_kk*vnnr*rsqinv;
    KK_FLOAT ccel = kn_kk*(radsum-r)*rinv - damp;
    if (limit_damping && (ccel < static_cast<KK_FLOAT>(0.0))) ccel = static_cast<KK_FLOAT>(0.0);

    // relative velocities

    KK_FLOAT vtr1 = vt1 - (delz*wr2-dely*wr3);
    KK_FLOAT vtr2 = vt2 - (delx*wr3-delz*wr1);
    KK_FLOAT vtr3 = vt3 - (dely*wr1-delx*wr2);

    // shear history effects

    KK_FLOAT shear1 = d_firstshear(i,3*jj);
    KK_FLOAT shear2 = d_firstshear(i,3*jj+1);
    KK_FLOAT shear3 = d_firstshear(i,3*jj+2);

    if (SHEARUPDATE) {
      shear1 += vtr1*dt_kk;
      shear2 += vtr2*dt_kk;
      shear3 += vtr3*dt_kk;
    }
    KK_FLOAT shrmag = sqrt(shear1*shear1 + shear2*shear2 +
                          shear3*shear3);

    if (SHEARUPDATE) {
      // rotate shear displacements

      KK_FLOAT rsht = shear1*delx + shear2*dely + shear3*delz;
      rsht *= rsqinv;

      shear1 -= rsht*delx;
      shear2 -= rsht*dely;
      shear3 -= rsht*delz;
    }

    // tangential forces = shear + tangential velocity damping

    KK_FLOAT fs1 = - (kt_kk*shear1 + meff*gammat_kk*vtr1);
    KK_FLOAT fs2 = - (kt_kk*shear2 + meff*gammat_kk*vtr2);
    KK_FLOAT fs3 = - (kt_kk*shear3 + meff*gammat_kk*vtr3);

    // rescale frictional displacements and forces if needed

    KK_FLOAT fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
    KK_FLOAT fn = xmu_kk * fabs(ccel*r);

    if (fs > fn) {
      if (shrmag != static_cast<KK_FLOAT>(0.0)) {
        shear1 = (fn/fs) * (shear1 + meff*gammat_kk*vtr1/kt_kk) -
          meff*gammat_kk*vtr1/kt_kk;
        shear2 = (fn/fs) * (shear2 + meff*gammat_kk*vtr2/kt_kk) -
          meff*gammat_kk*vtr2/kt_kk;
        shear3 = (fn/fs) * (shear3 + meff*gammat_kk*vtr3/kt_kk) -
          meff*gammat_kk*vtr3/kt_kk;
        fs1 *= fn/fs;
        fs2 *= fn/fs;
        fs3 *= fn/fs;
      } else fs1 = fs2 = fs3 = 0;
    }

    if (SHEARUPDATE) {
      d_firstshear(i,3*jj) = shear1;
      d_firstshear(i,3*jj+1) = shear2;
      d_firstshear(i,3*jj+2) = shear3;
    }

    // forces & torques

    KK_FLOAT fx = delx*ccel + fs1;
    KK_FLOAT fy = dely*ccel + fs2;
    KK_FLOAT fz = delz*ccel + fs3;
    fx *= factor_lj;
    fy *= factor_lj;
    fz *= factor_lj;
    fx_i += static_cast<KK_ACC_FLOAT>(fx);
    fy_i += static_cast<KK_ACC_FLOAT>(fy);
    fz_i += static_cast<KK_ACC_FLOAT>(fz);

    KK_FLOAT tor1 = rinv * (dely*fs3 - delz*fs2);
    KK_FLOAT tor2 = rinv * (delz*fs1 - delx*fs3);
    KK_FLOAT tor3 = rinv * (delx*fs2 - dely*fs1);
    tor1 *= factor_lj;
    tor2 *= factor_lj;
    tor3 *= factor_lj;
    torquex_i -= static_cast<KK_ACC_FLOAT>(irad*tor1);
    torquey_i -= static_cast<KK_ACC_FLOAT>(irad*tor2);
    torquez_i -= static_cast<KK_ACC_FLOAT>(irad*tor3);

    if (NEWTON_PAIR || j < nlocal) {
      a_f(j,0) -= static_cast<KK_ACC_FLOAT>(fx);
      a_f(j,1) -= static_cast<KK_ACC_FLOAT>(fy);
      a_f(j,2) -= static_cast<KK_ACC_FLOAT>(fz);
      a_torque(j,0) -= static_cast<KK_ACC_FLOAT>(jrad*tor1);
      a_torque(j,1) -= static_cast<KK_ACC_FLOAT>(jrad*tor2);
      a_torque(j,2) -= static_cast<KK_ACC_FLOAT>(jrad*tor3);
    }

    if (VFLAG)
      ev_tally_xyz<NEIGHFLAG, NEWTON_PAIR>(ev, i, j, fx, fy, fz, delx, dely, delz);
  }

  a_f(i,0) += fx_i;
  a_f(i,1) += fy_i;
  a_f(i,2) += fz_i;
  a_torque(i,0) += torquex_i;
  a_torque(i,1) += torquey_i;
  a_torque(i,2) += torquez_i;
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int SHEARUPDATE>
KOKKOS_INLINE_FUNCTION
void PairGranHookeHistoryKokkos<DeviceType>::operator()(TagPairGranHookeHistoryCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,SHEARUPDATE>, const int ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,VFLAG,SHEARUPDATE>(TagPairGranHookeHistoryCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,SHEARUPDATE>(), ii, ev);
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairGranHookeHistoryKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT &ev, int i, int j,
                                                          KK_FLOAT fx, KK_FLOAT fy, KK_FLOAT fz,
                                                          KK_FLOAT delx, KK_FLOAT dely, KK_FLOAT delz) const
{
  Kokkos::View<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom = d_vatom;

  if (vflag_global || vflag_atom) {
    const KK_ACC_FLOAT v_acc[6] =
      { static_cast<KK_ACC_FLOAT>(delx*fx),
        static_cast<KK_ACC_FLOAT>(dely*fy),
        static_cast<KK_ACC_FLOAT>(delz*fz),
        static_cast<KK_ACC_FLOAT>(delx*fy),
        static_cast<KK_ACC_FLOAT>(delx*fz),
        static_cast<KK_ACC_FLOAT>(dely*fz) };

    if (vflag_global) {
      if (NEWTON_PAIR) { // neigh half, newton on
        for (int n = 0; n < 6; n++)
          ev.v[n] += v_acc[n];
      } else { // neigh half, newton off
        if (i < nlocal) {
          for (int n = 0; n < 6; n++)
            ev.v[n] += static_cast<KK_ACC_FLOAT>(0.5) * v_acc[n];
        }
        if (j < nlocal) {
          for (int n = 0; n < 6; n++)
            ev.v[n] += static_cast<KK_ACC_FLOAT>(0.5) * v_acc[n];
        }
      }
    }

    if (vflag_atom) {

      if (NEWTON_PAIR || i < nlocal) {
        for (int n = 0; n < 6; n++)
          v_vatom(i,n) += static_cast<KK_ACC_FLOAT>(0.5) * v_acc[n];
      }
      if (NEWTON_PAIR || j < nlocal) {
        for (int n = 0; n < 6; n++)
          v_vatom(j,n) += static_cast<KK_ACC_FLOAT>(0.5) * v_acc[n];
      }
    }
  }
}

namespace LAMMPS_NS {
template class PairGranHookeHistoryKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairGranHookeHistoryKokkos<LMPHostType>;
#endif
}
