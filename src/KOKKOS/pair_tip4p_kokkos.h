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

/* ----------------------------------------------------------------------
   Shared device machinery for the KOKKOS TIP4P pair styles
   (tip4p/cut, tip4p/long, lj/cut/tip4p/cut, lj/cut/tip4p/long).

   PairTIP4PKokkos<DeviceType,PairCPUBase> derives from the corresponding CPU
   base class and provides the parts every TIP4P pair style needs:
     - the per-O M (virtual charge) site pre-kernel (find the two H atoms by
       tag on device, closest_image them, compute the M-site position),
     - a device closest_image(),
     - ev_tally_tip4p() (the key-based global/per-atom energy/virial split for
       the off-atom charge site), and
     - prepare()/finalize() helpers for the common compute() setup/teardown.

   Each concrete style keeps only its own main interaction kernel (LJ on/off,
   cut vs long-range Coulomb).  The shared pair_kokkos.h is NOT modified.
------------------------------------------------------------------------- */

#ifndef LMP_PAIR_TIP4P_KOKKOS_H
#define LMP_PAIR_TIP4P_KOKKOS_H

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "force.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair_kokkos.h"
#include "update.h"
#include <Kokkos_UnorderedMap.hpp>

namespace LAMMPS_NS {

// pre-kernel: compute the M (virtual charge) site for every O atom
struct TagPairTIP4PNewsite{};

template<class DeviceType, class PairCPUBase>
class PairTIP4PKokkos : public PairCPUBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairTIP4PKokkos(class LAMMPS *lmp) : PairCPUBase(lmp)
  {
    this->kokkosable = 1;
    this->atomKK = (AtomKokkos *) this->atom;
    this->execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
    this->datamask_read = X_MASK | F_MASK | TYPE_MASK | Q_MASK | TAG_MASK | ENERGY_MASK | VIRIAL_MASK;
    this->datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
    // TIP4P tallies the virial explicitly (charge site is off-atom)
    this->no_virial_fdotr_compute = 1;
  }

  ~PairTIP4PKokkos() override
  {
    if (this->copymode) return;
    if (this->allocated) {
      this->memoryKK->destroy_kokkos(k_eatom, this->eatom);
      this->memoryKK->destroy_kokkos(k_vatom, this->vatom);
    }
  }

  void init_style() override
  {
    PairCPUBase::init_style();
    auto request = this->neighbor->find_request(this);
    request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                             !std::is_same_v<DeviceType,LMPDeviceType>);
    request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  }

  // off-atom M (charge) site for every O atom (local + ghost)
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTIP4PNewsite, const int &i) const
  {
    if (type(i) != m_typeO) return;
    int iH1 = AtomKokkos::map_kokkos<DeviceType>(tag(i)+1,map_style,k_map_array,k_map_hash);
    int iH2 = AtomKokkos::map_kokkos<DeviceType>(tag(i)+2,map_style,k_map_array,k_map_hash);
    if (iH1 < 0 || iH2 < 0) { d_hneigh(i,0) = -1; return; }
    iH1 = closest_image(i,iH1);
    iH2 = closest_image(i,iH2);
    d_hneigh(i,0) = iH1;
    d_hneigh(i,1) = iH2;
    // xM = xO + alpha * 0.5 * ((xH1 - xO) + (xH2 - xO))
    d_newsite(i,0) = x(i,0) + m_alpha*(KK_FLOAT)0.5*(x(iH1,0) + x(iH2,0) - (KK_FLOAT)2.0*x(i,0));
    d_newsite(i,1) = x(i,1) + m_alpha*(KK_FLOAT)0.5*(x(iH1,1) + x(iH2,1) - (KK_FLOAT)2.0*x(i,1));
    d_newsite(i,2) = x(i,2) + m_alpha*(KK_FLOAT)0.5*(x(iH1,2) + x(iH2,2) - (KK_FLOAT)2.0*x(i,2));
  }

 protected:
  // device-callable special-bond mask (Pair::sbmask is host-only)
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int sbmask(const int j) const { return j >> SBBITS & 3; }

  // find the periodic image of j closest to i (walks the sametag chain)
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int closest_image(const int i, int j) const
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

  // global/per-atom energy and virial tally for one M-site interaction.
  // key encodes whether i and/or j are O atoms (see Pair::ev_tally_tip4p)
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally_tip4p(EV_FLOAT &ev, const int &key, const int (&vlist)[6],
                      const KK_FLOAT (&v)[6], const KK_FLOAT &ecoul) const
  {
    if (this->eflag_global) ev.ecoul += ecoul;
    if (this->vflag_global)
      for (int k = 0; k < 6; k++) ev.v[k] += v[k];

    if (this->eflag_atom) {
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

    if (this->vflag_atom) {
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

  // common compute() setup: sync atom data, bind device views, cache TIP4P
  // parameters and the atom map, (re)allocate the M-site scratch and per-atom
  // energy/virial arrays.  Returns the number of owned atoms (inum).
  int prepare(int eflag_in, int vflag_in)
  {
    this->eflag = eflag_in;
    this->vflag = vflag_in;
    this->ev_init(this->eflag,this->vflag);

    this->atomKK->sync(this->execution_space,this->datamask_read);

    x = this->atomKK->k_x.template view<DeviceType>();
    f = this->atomKK->k_f.template view<DeviceType>();
    q = this->atomKK->k_q.template view<DeviceType>();
    type = this->atomKK->k_type.template view<DeviceType>();
    tag = this->atomKK->k_tag.template view<DeviceType>();
    this->atomKK->k_sametag.template sync<DeviceType>();
    d_sametag = this->atomKK->k_sametag.template view<DeviceType>();

    nlocal = this->atom->nlocal;
    nall = this->atom->nlocal + this->atom->nghost;
    newton_pair = this->force->newton_pair;
    qqrd2e = this->force->qqrd2e;
    for (int i = 0; i < 4; i++) {
      special_coul[i] = this->force->special_coul[i];
      special_lj[i] = this->force->special_lj[i];
    }

    m_alpha = this->alpha;
    m_typeO = this->typeO;
    m_typeH = this->typeH;
    m_cut_coulsq = this->cut_coulsq;
    m_cut_coulsqplus = (this->cut_coul + 2.0*this->qdist) * (this->cut_coul + 2.0*this->qdist);

    map_style = this->atom->map_style;
    if (map_style == Atom::MAP_ARRAY) {
      k_map_array = this->atomKK->k_map_array;
      k_map_array.template sync<DeviceType>();
    } else if (map_style == Atom::MAP_HASH) {
      k_map_hash = this->atomKK->k_map_hash;
      k_map_hash.template sync<DeviceType>();
    }

    if ((int)d_newsite.extent(0) < this->atom->nmax) {
      d_newsite = typename AT::t_kkfloat_1d_3("tip4p/kk:newsite", this->atom->nmax);
      d_hneigh  = typename AT::t_int_1d_3("tip4p/kk:hneigh", this->atom->nmax);
    }

    if (this->eflag_atom) {
      this->memoryKK->destroy_kokkos(k_eatom, this->eatom);
      this->memoryKK->create_kokkos(k_eatom, this->eatom, this->maxeatom, "pair:eatom");
      d_eatom = k_eatom.template view<DeviceType>();
    }
    if (this->vflag_atom) {
      this->memoryKK->destroy_kokkos(k_vatom, this->vatom);
      this->memoryKK->create_kokkos(k_vatom, this->vatom, this->maxvatom, "pair:vatom");
      d_vatom = k_vatom.template view<DeviceType>();
    }

    NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(this->list);
    d_numneigh = k_list->d_numneigh;
    d_neighbors = k_list->d_neighbors;
    d_ilist = k_list->d_ilist;
    return this->list->inum;
  }

  // common compute() teardown: accumulate global coulomb energy and virial,
  // sync per-atom arrays back to the host, mark modified atom data
  void finalize(const EV_FLOAT &ev)
  {
    if (this->eflag_global) this->eng_coul += ev.ecoul;
    if (this->vflag_global) {
      this->virial[0] += ev.v[0]; this->virial[1] += ev.v[1]; this->virial[2] += ev.v[2];
      this->virial[3] += ev.v[3]; this->virial[4] += ev.v[4]; this->virial[5] += ev.v[5];
    }
    if (this->eflag_atom) { k_eatom.template modify<DeviceType>(); k_eatom.sync_host(); }
    if (this->vflag_atom) { k_vatom.template modify<DeviceType>(); k_vatom.sync_host(); }
    this->atomKK->modified(this->execution_space,this->datamask_modify);
  }

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkfloat_1d_randomread q;
  typename AT::t_int_1d_randomread type;
  typename AT::t_tagint_1d_randomread tag;
  typename AT::t_int_1d d_sametag;

  typename AT::t_kkfloat_1d_3 d_newsite;
  typename AT::t_int_1d_3 d_hneigh;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int map_style;
  DAT::tdual_int_1d k_map_array;
  dual_hash_type k_map_hash;

  int newton_pair, neighflag;
  int nlocal, nall, eflag, vflag;

  KK_FLOAT special_coul[4];
  KK_FLOAT special_lj[4];
  KK_FLOAT qqrd2e;
  KK_FLOAT m_alpha;
  KK_FLOAT m_cut_coulsq, m_cut_coulsqplus;
  int m_typeO, m_typeH;
};

}    // namespace LAMMPS_NS

#endif
