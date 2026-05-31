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

#include "pair_lj_cut_tip4p_long_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kspace.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace EwaldConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCutTIP4PLongKokkos<DeviceType>::PairLJCutTIP4PLongKokkos(LAMMPS *lmp) : PairLJCutTIP4PLong(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | Q_MASK | TAG_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
  no_virial_fdotr_compute = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCutTIP4PLongKokkos<DeviceType>::~PairLJCutTIP4PLongKokkos()
{
  if (copymode) return;
  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom, eatom);
    memoryKK->destroy_kokkos(k_vatom, vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutTIP4PLongKokkos<DeviceType>::init_style()
{
  PairLJCutTIP4PLong::init_style();

  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutTIP4PLongKokkos<DeviceType>::init_tables(double cut_coul, double *cut_respa)
{
  Pair::init_tables(cut_coul,cut_respa);

  typedef typename AT::t_kkfloat_1d table_type;
  typedef HAT::t_kkfloat_1d host_table_type;

  int ntable = 1;
  for (int i = 0; i < ncoultablebits; i++) ntable *= 2;

  tabinnersq_kk = static_cast<KK_FLOAT>(tabinnersq);

  auto copy_table = [&](double *src, table_type &dst) {
    host_table_type h_table("HostTable",ntable);
    table_type d_table("DeviceTable",ntable);
    for (int i = 0; i < ntable; i++) h_table(i) = static_cast<KK_FLOAT>(src[i]);
    Kokkos::deep_copy(d_table,h_table);
    dst = d_table;
  };

  copy_table(rtable,d_rtable);   copy_table(drtable,d_drtable);
  copy_table(ftable,d_ftable);   copy_table(dftable,d_dftable);
  copy_table(ctable,d_ctable);   copy_table(dctable,d_dctable);
  copy_table(etable,d_etable);   copy_table(detable,d_detable);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutTIP4PLongKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag);

  atomKK->sync(execution_space,datamask_read);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  tag = atomKK->k_tag.view<DeviceType>();
  atomKK->k_sametag.template sync<DeviceType>();
  d_sametag = atomKK->k_sametag.view<DeviceType>();

  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  qqrd2e = force->qqrd2e;
  for (int i = 0; i < 4; i++) { special_coul[i] = force->special_coul[i]; special_lj[i] = force->special_lj[i]; }

  m_alpha = alpha;
  m_qdist = qdist;
  m_typeO = typeO;
  m_typeH = typeH;
  m_cut_coulsq = cut_coulsq;
  m_cut_coulsqplus = (cut_coul + 2.0*qdist) * (cut_coul + 2.0*qdist);
  g_ewald_kk = static_cast<KK_FLOAT>(g_ewald);
  m_ncoultablebits = ncoultablebits;
  m_ncoulmask = ncoulmask;
  m_ncoulshiftbits = ncoulshiftbits;

  const int ntp1 = atom->ntypes + 1;
  if ((int)d_lj1.extent(0) != ntp1) {
    d_lj1 = typename AT::t_kkfloat_2d("lj/cut/tip4p/long/kk:lj1",ntp1,ntp1);
    d_lj2 = typename AT::t_kkfloat_2d("lj/cut/tip4p/long/kk:lj2",ntp1,ntp1);
    d_lj3 = typename AT::t_kkfloat_2d("lj/cut/tip4p/long/kk:lj3",ntp1,ntp1);
    d_lj4 = typename AT::t_kkfloat_2d("lj/cut/tip4p/long/kk:lj4",ntp1,ntp1);
    d_offset = typename AT::t_kkfloat_2d("lj/cut/tip4p/long/kk:offset",ntp1,ntp1);
    d_cut_ljsq = typename AT::t_kkfloat_2d("lj/cut/tip4p/long/kk:cut_ljsq",ntp1,ntp1);
  }
  {
    auto h_lj1 = Kokkos::create_mirror_view(d_lj1);
    auto h_lj2 = Kokkos::create_mirror_view(d_lj2);
    auto h_lj3 = Kokkos::create_mirror_view(d_lj3);
    auto h_lj4 = Kokkos::create_mirror_view(d_lj4);
    auto h_offset = Kokkos::create_mirror_view(d_offset);
    auto h_cut_ljsq = Kokkos::create_mirror_view(d_cut_ljsq);
    for (int i = 1; i < ntp1; i++)
      for (int j = 1; j < ntp1; j++) {
        h_lj1(i,j) = lj1[i][j]; h_lj2(i,j) = lj2[i][j];
        h_lj3(i,j) = lj3[i][j]; h_lj4(i,j) = lj4[i][j];
        h_offset(i,j) = offset[i][j]; h_cut_ljsq(i,j) = cut_ljsq[i][j];
      }
    Kokkos::deep_copy(d_lj1,h_lj1); Kokkos::deep_copy(d_lj2,h_lj2);
    Kokkos::deep_copy(d_lj3,h_lj3); Kokkos::deep_copy(d_lj4,h_lj4);
    Kokkos::deep_copy(d_offset,h_offset); Kokkos::deep_copy(d_cut_ljsq,h_cut_ljsq);
  }

  map_style = atom->map_style;
  if (map_style == Atom::MAP_ARRAY) {
    k_map_array = atomKK->k_map_array;
    k_map_array.template sync<DeviceType>();
  } else if (map_style == Atom::MAP_HASH) {
    k_map_hash = atomKK->k_map_hash;
    k_map_hash.template sync<DeviceType>();
  }

  if ((int)d_newsite.extent(0) < atom->nmax) {
    d_newsite = typename AT::t_kkfloat_1d_3("lj/cut/tip4p/long/kk:newsite", atom->nmax);
    d_hneigh  = typename AT::t_int_1d_3("lj/cut/tip4p/long/kk:hneigh", atom->nmax);
  }

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom, eatom);
    memoryKK->create_kokkos(k_eatom, eatom, maxeatom, "pair:eatom");
    d_eatom = k_eatom.template view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom, vatom);
    memoryKK->create_kokkos(k_vatom, vatom, maxvatom, "pair:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  const int inum = list->inum;

  copymode = 1;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairLJCutTIP4PLongNewsite>(0,nall), *this);

  EV_FLOAT ev;
  if (eflag || vflag)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairLJCutTIP4PLongCompute<1>>(0,inum), *this, ev);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairLJCutTIP4PLongCompute<0>>(0,inum), *this, ev);

  copymode = 0;

  if (eflag_global) { eng_vdwl += ev.evdwl; eng_coul += ev.ecoul; }
  if (vflag_global) {
    virial[0] += ev.v[0]; virial[1] += ev.v[1]; virial[2] += ev.v[2];
    virial[3] += ev.v[3]; virial[4] += ev.v[4]; virial[5] += ev.v[5];
  }

  if (eflag_atom) { k_eatom.template modify<DeviceType>(); k_eatom.sync_host(); }
  if (vflag_atom) { k_vatom.template modify<DeviceType>(); k_vatom.sync_host(); }

  atomKK->modified(execution_space,datamask_modify);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
int PairLJCutTIP4PLongKokkos<DeviceType>::closest_image(const int i, int j) const
{
  if (j < 0) return j;
  const KK_FLOAT xi0 = x(i,0), xi1 = x(i,1), xi2 = x(i,2);
  int closest = j;
  KK_FLOAT delx = xi0 - x(j,0), dely = xi1 - x(j,1), delz = xi2 - x(j,2);
  KK_FLOAT rsqmin = delx*delx + dely*dely + delz*delz;
  while (d_sametag[j] >= 0) {
    j = d_sametag[j];
    delx = xi0 - x(j,0); dely = xi1 - x(j,1); delz = xi2 - x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    if (rsq < rsqmin) { rsqmin = rsq; closest = j; }
  }
  return closest;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairLJCutTIP4PLongKokkos<DeviceType>::operator()(TagPairLJCutTIP4PLongNewsite, const int &i) const
{
  if (type(i) != m_typeO) return;
  int iH1 = AtomKokkos::map_kokkos<DeviceType>(tag(i)+1,map_style,k_map_array,k_map_hash);
  int iH2 = AtomKokkos::map_kokkos<DeviceType>(tag(i)+2,map_style,k_map_array,k_map_hash);
  if (iH1 < 0 || iH2 < 0) { d_hneigh(i,0) = -1; return; }
  iH1 = closest_image(i,iH1);
  iH2 = closest_image(i,iH2);
  d_hneigh(i,0) = iH1;
  d_hneigh(i,1) = iH2;
  d_newsite(i,0) = x(i,0) + m_alpha*(KK_FLOAT)0.5*(x(iH1,0) + x(iH2,0) - (KK_FLOAT)2.0*x(i,0));
  d_newsite(i,1) = x(i,1) + m_alpha*(KK_FLOAT)0.5*(x(iH1,1) + x(iH2,1) - (KK_FLOAT)2.0*x(i,1));
  d_newsite(i,2) = x(i,2) + m_alpha*(KK_FLOAT)0.5*(x(iH1,2) + x(iH2,2) - (KK_FLOAT)2.0*x(i,2));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairLJCutTIP4PLongKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
    const KK_FLOAT &evdwl, const KK_FLOAT &fpair,
    const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const
{
  if (eflag_global) ev.evdwl += evdwl;
  if (eflag_atom) {
    Kokkos::atomic_add(&d_eatom[i], (KK_ACC_FLOAT)(0.5*evdwl));
    Kokkos::atomic_add(&d_eatom[j], (KK_ACC_FLOAT)(0.5*evdwl));
  }
  if (vflag_global || vflag_atom) {
    const KK_FLOAT v0 = delx*delx*fpair, v1 = dely*dely*fpair, v2 = delz*delz*fpair;
    const KK_FLOAT v3 = delx*dely*fpair, v4 = delx*delz*fpair, v5 = dely*delz*fpair;
    if (vflag_global) {
      ev.v[0] += v0; ev.v[1] += v1; ev.v[2] += v2;
      ev.v[3] += v3; ev.v[4] += v4; ev.v[5] += v5;
    }
    if (vflag_atom) {
      Kokkos::atomic_add(&d_vatom(i,0), (KK_ACC_FLOAT)(0.5*v0));
      Kokkos::atomic_add(&d_vatom(i,1), (KK_ACC_FLOAT)(0.5*v1));
      Kokkos::atomic_add(&d_vatom(i,2), (KK_ACC_FLOAT)(0.5*v2));
      Kokkos::atomic_add(&d_vatom(i,3), (KK_ACC_FLOAT)(0.5*v3));
      Kokkos::atomic_add(&d_vatom(i,4), (KK_ACC_FLOAT)(0.5*v4));
      Kokkos::atomic_add(&d_vatom(i,5), (KK_ACC_FLOAT)(0.5*v5));
      Kokkos::atomic_add(&d_vatom(j,0), (KK_ACC_FLOAT)(0.5*v0));
      Kokkos::atomic_add(&d_vatom(j,1), (KK_ACC_FLOAT)(0.5*v1));
      Kokkos::atomic_add(&d_vatom(j,2), (KK_ACC_FLOAT)(0.5*v2));
      Kokkos::atomic_add(&d_vatom(j,3), (KK_ACC_FLOAT)(0.5*v3));
      Kokkos::atomic_add(&d_vatom(j,4), (KK_ACC_FLOAT)(0.5*v4));
      Kokkos::atomic_add(&d_vatom(j,5), (KK_ACC_FLOAT)(0.5*v5));
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairLJCutTIP4PLongKokkos<DeviceType>::ev_tally_tip4p(EV_FLOAT &ev, const int &key,
    const int (&vlist)[6], const KK_FLOAT (&v)[6], const KK_FLOAT &ecoul) const
{
  if (eflag_global) ev.ecoul += ecoul;
  if (vflag_global)
    for (int k = 0; k < 6; k++) ev.v[k] += v[k];

  if (eflag_atom) {
    const KK_FLOAT a = m_alpha;
    if (key == 0) {
      Kokkos::atomic_add(&d_eatom[vlist[0]], (KK_ACC_FLOAT)(0.5*ecoul));
      Kokkos::atomic_add(&d_eatom[vlist[1]], (KK_ACC_FLOAT)(0.5*ecoul));
    } else if (key == 1) {
      Kokkos::atomic_add(&d_eatom[vlist[0]], (KK_ACC_FLOAT)(0.5*ecoul*(1.0-a)));
      Kokkos::atomic_add(&d_eatom[vlist[1]], (KK_ACC_FLOAT)(0.25*ecoul*a));
      Kokkos::atomic_add(&d_eatom[vlist[2]], (KK_ACC_FLOAT)(0.25*ecoul*a));
      Kokkos::atomic_add(&d_eatom[vlist[3]], (KK_ACC_FLOAT)(0.5*ecoul));
    } else if (key == 2) {
      Kokkos::atomic_add(&d_eatom[vlist[0]], (KK_ACC_FLOAT)(0.5*ecoul));
      Kokkos::atomic_add(&d_eatom[vlist[1]], (KK_ACC_FLOAT)(0.5*ecoul*(1.0-a)));
      Kokkos::atomic_add(&d_eatom[vlist[2]], (KK_ACC_FLOAT)(0.25*ecoul*a));
      Kokkos::atomic_add(&d_eatom[vlist[3]], (KK_ACC_FLOAT)(0.25*ecoul*a));
    } else {
      Kokkos::atomic_add(&d_eatom[vlist[0]], (KK_ACC_FLOAT)(0.5*ecoul*(1.0-a)));
      Kokkos::atomic_add(&d_eatom[vlist[1]], (KK_ACC_FLOAT)(0.25*ecoul*a));
      Kokkos::atomic_add(&d_eatom[vlist[2]], (KK_ACC_FLOAT)(0.25*ecoul*a));
      Kokkos::atomic_add(&d_eatom[vlist[3]], (KK_ACC_FLOAT)(0.5*ecoul*(1.0-a)));
      Kokkos::atomic_add(&d_eatom[vlist[4]], (KK_ACC_FLOAT)(0.25*ecoul*a));
      Kokkos::atomic_add(&d_eatom[vlist[5]], (KK_ACC_FLOAT)(0.25*ecoul*a));
    }
  }

  if (vflag_atom) {
    const KK_FLOAT a = m_alpha;
    for (int k = 0; k < 6; k++) {
      const KK_FLOAT vk = v[k];
      if (key == 0) {
        Kokkos::atomic_add(&d_vatom(vlist[0],k), (KK_ACC_FLOAT)(0.5*vk));
        Kokkos::atomic_add(&d_vatom(vlist[1],k), (KK_ACC_FLOAT)(0.5*vk));
      } else if (key == 1) {
        Kokkos::atomic_add(&d_vatom(vlist[0],k), (KK_ACC_FLOAT)(0.5*vk*(1.0-a)));
        Kokkos::atomic_add(&d_vatom(vlist[1],k), (KK_ACC_FLOAT)(0.25*vk*a));
        Kokkos::atomic_add(&d_vatom(vlist[2],k), (KK_ACC_FLOAT)(0.25*vk*a));
        Kokkos::atomic_add(&d_vatom(vlist[3],k), (KK_ACC_FLOAT)(0.5*vk));
      } else if (key == 2) {
        Kokkos::atomic_add(&d_vatom(vlist[0],k), (KK_ACC_FLOAT)(0.5*vk));
        Kokkos::atomic_add(&d_vatom(vlist[1],k), (KK_ACC_FLOAT)(0.5*vk*(1.0-a)));
        Kokkos::atomic_add(&d_vatom(vlist[2],k), (KK_ACC_FLOAT)(0.25*vk*a));
        Kokkos::atomic_add(&d_vatom(vlist[3],k), (KK_ACC_FLOAT)(0.25*vk*a));
      } else {
        Kokkos::atomic_add(&d_vatom(vlist[0],k), (KK_ACC_FLOAT)(0.5*vk*(1.0-a)));
        Kokkos::atomic_add(&d_vatom(vlist[1],k), (KK_ACC_FLOAT)(0.25*vk*a));
        Kokkos::atomic_add(&d_vatom(vlist[2],k), (KK_ACC_FLOAT)(0.25*vk*a));
        Kokkos::atomic_add(&d_vatom(vlist[3],k), (KK_ACC_FLOAT)(0.5*vk*(1.0-a)));
        Kokkos::atomic_add(&d_vatom(vlist[4],k), (KK_ACC_FLOAT)(0.25*vk*a));
        Kokkos::atomic_add(&d_vatom(vlist[5],k), (KK_ACC_FLOAT)(0.25*vk*a));
      }
    }
  }
}

/* ----------------------------------------------------------------------
   main kernel: LJ on atom positions + long-range Coulomb on the M sites
------------------------------------------------------------------------- */

template<class DeviceType>
template<int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairLJCutTIP4PLongKokkos<DeviceType>::operator()(TagPairLJCutTIP4PLongCompute<EVFLAG>,
                                                      const int &ii, EV_FLOAT &ev) const
{
  const int i = d_ilist(ii);
  const int itype = type(i);
  const KK_FLOAT qtmp = q(i);
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);

  KK_FLOAT x1[3];
  int iH1 = -1, iH2 = -1;
  if (itype == m_typeO) {
    iH1 = d_hneigh(i,0);
    iH2 = d_hneigh(i,1);
    x1[0] = d_newsite(i,0); x1[1] = d_newsite(i,1); x1[2] = d_newsite(i,2);
  } else {
    x1[0] = xtmp; x1[1] = ytmp; x1[2] = ztmp;
  }

  const int jnum = d_numneigh(i);

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    const KK_FLOAT factor_lj = special_lj[sbmask(j)];
    const KK_FLOAT factor_coul = special_coul[sbmask(j)];
    j &= NEIGHMASK;
    const int jtype = type(j);

    const KK_FLOAT adelx = xtmp - x(j,0);
    const KK_FLOAT adely = ytmp - x(j,1);
    const KK_FLOAT adelz = ztmp - x(j,2);
    const KK_FLOAT arsq = adelx*adelx + adely*adely + adelz*adelz;

    // ---- Lennard-Jones (real atom positions) ----
    if (arsq < d_cut_ljsq(itype,jtype)) {
      const KK_FLOAT r2inv = (KK_FLOAT)1.0 / arsq;
      const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
      const KK_FLOAT forcelj = factor_lj * r6inv *
          (d_lj1(itype,jtype)*r6inv - d_lj2(itype,jtype)) * r2inv;
      Kokkos::atomic_add(&f(i,0), (KK_ACC_FLOAT)(adelx*forcelj));
      Kokkos::atomic_add(&f(i,1), (KK_ACC_FLOAT)(adely*forcelj));
      Kokkos::atomic_add(&f(i,2), (KK_ACC_FLOAT)(adelz*forcelj));
      Kokkos::atomic_add(&f(j,0), (KK_ACC_FLOAT)(-adelx*forcelj));
      Kokkos::atomic_add(&f(j,1), (KK_ACC_FLOAT)(-adely*forcelj));
      Kokkos::atomic_add(&f(j,2), (KK_ACC_FLOAT)(-adelz*forcelj));
      if (EVFLAG) {
        KK_FLOAT evdwl = 0.0;
        if (eflag)
          evdwl = factor_lj * (r6inv*(d_lj3(itype,jtype)*r6inv - d_lj4(itype,jtype))
                               - d_offset(itype,jtype));
        ev_tally(ev,i,j,evdwl,forcelj,adelx,adely,adelz);
      }
    }

    // ---- long-range Coulomb (M sites) ----
    if (arsq >= m_cut_coulsqplus) continue;

    KK_FLOAT delx = adelx, dely = adely, delz = adelz, rsq = arsq;
    int jH1 = -1, jH2 = -1;
    if (itype == m_typeO || jtype == m_typeO) {
      KK_FLOAT x2[3];
      if (jtype == m_typeO) {
        jH1 = d_hneigh(j,0);
        jH2 = d_hneigh(j,1);
        x2[0] = d_newsite(j,0); x2[1] = d_newsite(j,1); x2[2] = d_newsite(j,2);
      } else {
        x2[0] = x(j,0); x2[1] = x(j,1); x2[2] = x(j,2);
      }
      delx = x1[0] - x2[0];
      dely = x1[1] - x2[1];
      delz = x1[2] - x2[2];
      rsq = delx*delx + dely*dely + delz*delz;
    }

    if (rsq >= m_cut_coulsq) continue;

    const KK_FLOAT r2inv = (KK_FLOAT)1.0 / rsq;
    KK_FLOAT cforce, ecoul = 0.0;
    if (!m_ncoultablebits || rsq <= tabinnersq_kk) {
      const KK_FLOAT r = Kokkos::sqrt(rsq);
      const KK_FLOAT grij = g_ewald_kk * r;
      const KK_FLOAT expm2 = Kokkos::exp(-grij*grij);
      const KK_FLOAT t = (KK_FLOAT)1.0 / ((KK_FLOAT)1.0 + (KK_FLOAT)EWALD_P*grij);
      const KK_FLOAT erfc = t*((KK_FLOAT)A1+t*((KK_FLOAT)A2+t*((KK_FLOAT)A3+
                            t*((KK_FLOAT)A4+t*(KK_FLOAT)A5))))*expm2;
      const KK_FLOAT prefactor = qqrd2e * qtmp * q(j) / r;
      KK_FLOAT forcecoul = prefactor * (erfc + (KK_FLOAT)EWALD_F*grij*expm2);
      if (factor_coul < (KK_FLOAT)1.0) forcecoul -= ((KK_FLOAT)1.0-factor_coul)*prefactor;
      cforce = forcecoul * r2inv;
      if (EVFLAG && eflag) {
        ecoul = prefactor * erfc;
        if (factor_coul < (KK_FLOAT)1.0) ecoul -= ((KK_FLOAT)1.0-factor_coul)*prefactor;
      }
    } else {
      union_int_float_t rsq_lookup;
      rsq_lookup.f = rsq;
      const int itable = (rsq_lookup.i & m_ncoulmask) >> m_ncoulshiftbits;
      const KK_FLOAT fraction = ((KK_FLOAT)rsq_lookup.f - d_rtable[itable]) * d_drtable[itable];
      const KK_FLOAT tbl = d_ftable[itable] + fraction*d_dftable[itable];
      KK_FLOAT forcecoul = qtmp * q(j) * tbl;
      KK_FLOAT prefactor = 0.0;
      if (factor_coul < (KK_FLOAT)1.0) {
        const KK_FLOAT ctbl = d_ctable[itable] + fraction*d_dctable[itable];
        prefactor = qtmp * q(j) * ctbl;
        forcecoul -= ((KK_FLOAT)1.0-factor_coul)*prefactor;
      }
      cforce = forcecoul * r2inv;
      if (EVFLAG && eflag) {
        const KK_FLOAT etbl = d_etable[itable] + fraction*d_detable[itable];
        ecoul = qtmp * q(j) * etbl;
        if (factor_coul < (KK_FLOAT)1.0) ecoul -= ((KK_FLOAT)1.0-factor_coul)*prefactor;
      }
    }

    int n = 0, key = 0, vlist[6];
    KK_FLOAT v[6] = {0,0,0,0,0,0};

    if (itype != m_typeO) {
      Kokkos::atomic_add(&f(i,0), (KK_ACC_FLOAT)(delx*cforce));
      Kokkos::atomic_add(&f(i,1), (KK_ACC_FLOAT)(dely*cforce));
      Kokkos::atomic_add(&f(i,2), (KK_ACC_FLOAT)(delz*cforce));
      if (EVFLAG && vflag) {
        v[0] = x(i,0)*delx*cforce; v[1] = x(i,1)*dely*cforce; v[2] = x(i,2)*delz*cforce;
        v[3] = x(i,0)*dely*cforce; v[4] = x(i,0)*delz*cforce; v[5] = x(i,1)*delz*cforce;
      }
      vlist[n++] = i;
    } else {
      key++;
      const KK_FLOAT fdx = delx*cforce, fdy = dely*cforce, fdz = delz*cforce;
      const KK_FLOAT fOx = fdx*(1.0-m_alpha), fOy = fdy*(1.0-m_alpha), fOz = fdz*(1.0-m_alpha);
      const KK_FLOAT fHx = 0.5*m_alpha*fdx, fHy = 0.5*m_alpha*fdy, fHz = 0.5*m_alpha*fdz;
      Kokkos::atomic_add(&f(i,0), (KK_ACC_FLOAT)fOx);
      Kokkos::atomic_add(&f(i,1), (KK_ACC_FLOAT)fOy);
      Kokkos::atomic_add(&f(i,2), (KK_ACC_FLOAT)fOz);
      Kokkos::atomic_add(&f(iH1,0), (KK_ACC_FLOAT)fHx);
      Kokkos::atomic_add(&f(iH1,1), (KK_ACC_FLOAT)fHy);
      Kokkos::atomic_add(&f(iH1,2), (KK_ACC_FLOAT)fHz);
      Kokkos::atomic_add(&f(iH2,0), (KK_ACC_FLOAT)fHx);
      Kokkos::atomic_add(&f(iH2,1), (KK_ACC_FLOAT)fHy);
      Kokkos::atomic_add(&f(iH2,2), (KK_ACC_FLOAT)fHz);
      if (EVFLAG && vflag) {
        v[0] = x(i,0)*fOx + x(iH1,0)*fHx + x(iH2,0)*fHx;
        v[1] = x(i,1)*fOy + x(iH1,1)*fHy + x(iH2,1)*fHy;
        v[2] = x(i,2)*fOz + x(iH1,2)*fHz + x(iH2,2)*fHz;
        v[3] = x(i,0)*fOy + x(iH1,0)*fHy + x(iH2,0)*fHy;
        v[4] = x(i,0)*fOz + x(iH1,0)*fHz + x(iH2,0)*fHz;
        v[5] = x(i,1)*fOz + x(iH1,1)*fHz + x(iH2,1)*fHz;
      }
      vlist[n++] = i; vlist[n++] = iH1; vlist[n++] = iH2;
    }

    if (jtype != m_typeO) {
      Kokkos::atomic_add(&f(j,0), (KK_ACC_FLOAT)(-delx*cforce));
      Kokkos::atomic_add(&f(j,1), (KK_ACC_FLOAT)(-dely*cforce));
      Kokkos::atomic_add(&f(j,2), (KK_ACC_FLOAT)(-delz*cforce));
      if (EVFLAG && vflag) {
        v[0] -= x(j,0)*delx*cforce; v[1] -= x(j,1)*dely*cforce; v[2] -= x(j,2)*delz*cforce;
        v[3] -= x(j,0)*dely*cforce; v[4] -= x(j,0)*delz*cforce; v[5] -= x(j,1)*delz*cforce;
      }
      vlist[n++] = j;
    } else {
      key += 2;
      const KK_FLOAT fdx = -delx*cforce, fdy = -dely*cforce, fdz = -delz*cforce;
      const KK_FLOAT fOx = fdx*(1.0-m_alpha), fOy = fdy*(1.0-m_alpha), fOz = fdz*(1.0-m_alpha);
      const KK_FLOAT fHx = 0.5*m_alpha*fdx, fHy = 0.5*m_alpha*fdy, fHz = 0.5*m_alpha*fdz;
      Kokkos::atomic_add(&f(j,0), (KK_ACC_FLOAT)fOx);
      Kokkos::atomic_add(&f(j,1), (KK_ACC_FLOAT)fOy);
      Kokkos::atomic_add(&f(j,2), (KK_ACC_FLOAT)fOz);
      Kokkos::atomic_add(&f(jH1,0), (KK_ACC_FLOAT)fHx);
      Kokkos::atomic_add(&f(jH1,1), (KK_ACC_FLOAT)fHy);
      Kokkos::atomic_add(&f(jH1,2), (KK_ACC_FLOAT)fHz);
      Kokkos::atomic_add(&f(jH2,0), (KK_ACC_FLOAT)fHx);
      Kokkos::atomic_add(&f(jH2,1), (KK_ACC_FLOAT)fHy);
      Kokkos::atomic_add(&f(jH2,2), (KK_ACC_FLOAT)fHz);
      if (EVFLAG && vflag) {
        v[0] += x(j,0)*fOx + x(jH1,0)*fHx + x(jH2,0)*fHx;
        v[1] += x(j,1)*fOy + x(jH1,1)*fHy + x(jH2,1)*fHy;
        v[2] += x(j,2)*fOz + x(jH1,2)*fHz + x(jH2,2)*fHz;
        v[3] += x(j,0)*fOy + x(jH1,0)*fHy + x(jH2,0)*fHy;
        v[4] += x(j,0)*fOz + x(jH1,0)*fHz + x(jH2,0)*fHz;
        v[5] += x(j,1)*fOz + x(jH1,1)*fHz + x(jH2,1)*fHz;
      }
      vlist[n++] = j; vlist[n++] = jH1; vlist[n++] = jH2;
    }

    if (EVFLAG) ev_tally_tip4p(ev,key,vlist,v,ecoul);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairLJCutTIP4PLongKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJCutTIP4PLongKokkos<LMPHostType>;
#endif
}
