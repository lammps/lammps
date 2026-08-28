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
       tag on device, closest_image them, compute the M-site position).
       A missing H is recorded as d_hneigh(i,0) = -1; the compute kernels
       check that sentinel where the M site is actually used (same semantics
       as the CPU styles' lazy find_M()) and set d_h_missing, which
       finalize() turns into the "TIP4P hydrogen is missing" error,
     - a device closest_image(),
     - apply_site_force() (per-side M-site force redistribution, virial and
       tally-list bookkeeping used by every interaction kernel),
     - ev_tally_tip4p() (the key-based global/per-atom energy/virial split for
       the off-atom charge site) and the standard pairwise ev_tally() for the
       LJ part of the lj/cut variants,
     - the long-range Coulomb machinery for the *long variants (init_tables(),
       prepare_coul_long(), device coul_long()), and
     - prepare()/finalize() helpers for the common compute() setup/teardown.

   Each concrete style keeps only its own main interaction kernel (LJ on/off,
   cut vs long-range Coulomb).  The shared pair_kokkos.h is NOT modified.
------------------------------------------------------------------------- */

#ifndef LMP_PAIR_TIP4P_KOKKOS_H
#define LMP_PAIR_TIP4P_KOKKOS_H

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "ewald_const.h"
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
    k_h_missing = DAT::tdual_int_scalar("pair:tip4p_h_missing");
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
  int closest_image(const int i, const int j) const
  {
    return AtomKokkos::closest_image_kokkos(i,j,x,d_sametag);
  }

  // apply the Coulomb force with prefactor cforce (sign folded in by the
  // caller: +cforce for the i side, -cforce for the j side) acting on the
  // charge site of atom idx along (delx,dely,delz).  For an O atom the force
  // acts on the off-atom M site and is redistributed onto O and its two H
  // atoms.  Accumulates the absolute-position virial in v and records the
  // participating atoms in vlist/n and the O-topology code in key for
  // ev_tally_tip4p().
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void apply_site_force(const int idx, const int iH1, const int iH2, const bool isO,
                        const KK_FLOAT delx, const KK_FLOAT dely, const KK_FLOAT delz,
                        const KK_FLOAT cforce, const bool do_virial, const int keyinc,
                        int &n, int &key, int (&vlist)[6], KK_ACC_FLOAT (&v)[6]) const
  {
    if (!isO) {
      Kokkos::atomic_add(&f(idx,0), (KK_ACC_FLOAT)(delx*cforce));
      Kokkos::atomic_add(&f(idx,1), (KK_ACC_FLOAT)(dely*cforce));
      Kokkos::atomic_add(&f(idx,2), (KK_ACC_FLOAT)(delz*cforce));
      if (do_virial) {
        v[0] += static_cast<KK_ACC_FLOAT>(x(idx,0)*delx*cforce);
        v[1] += static_cast<KK_ACC_FLOAT>(x(idx,1)*dely*cforce);
        v[2] += static_cast<KK_ACC_FLOAT>(x(idx,2)*delz*cforce);
        v[3] += static_cast<KK_ACC_FLOAT>(x(idx,0)*dely*cforce);
        v[4] += static_cast<KK_ACC_FLOAT>(x(idx,0)*delz*cforce);
        v[5] += static_cast<KK_ACC_FLOAT>(x(idx,1)*delz*cforce);
      }
      vlist[n++] = idx;
    } else {
      key += keyinc;
      const KK_FLOAT fdx = delx*cforce, fdy = dely*cforce, fdz = delz*cforce;
      const KK_FLOAT fOx = fdx*m_alphaO, fOy = fdy*m_alphaO, fOz = fdz*m_alphaO;
      const KK_FLOAT fHx = fdx*m_alphaH, fHy = fdy*m_alphaH, fHz = fdz*m_alphaH;
      Kokkos::atomic_add(&f(idx,0), (KK_ACC_FLOAT)fOx);
      Kokkos::atomic_add(&f(idx,1), (KK_ACC_FLOAT)fOy);
      Kokkos::atomic_add(&f(idx,2), (KK_ACC_FLOAT)fOz);
      Kokkos::atomic_add(&f(iH1,0), (KK_ACC_FLOAT)fHx);
      Kokkos::atomic_add(&f(iH1,1), (KK_ACC_FLOAT)fHy);
      Kokkos::atomic_add(&f(iH1,2), (KK_ACC_FLOAT)fHz);
      Kokkos::atomic_add(&f(iH2,0), (KK_ACC_FLOAT)fHx);
      Kokkos::atomic_add(&f(iH2,1), (KK_ACC_FLOAT)fHy);
      Kokkos::atomic_add(&f(iH2,2), (KK_ACC_FLOAT)fHz);
      if (do_virial) {
        v[0] += static_cast<KK_ACC_FLOAT>(x(idx,0)*fOx + x(iH1,0)*fHx + x(iH2,0)*fHx);
        v[1] += static_cast<KK_ACC_FLOAT>(x(idx,1)*fOy + x(iH1,1)*fHy + x(iH2,1)*fHy);
        v[2] += static_cast<KK_ACC_FLOAT>(x(idx,2)*fOz + x(iH1,2)*fHz + x(iH2,2)*fHz);
        v[3] += static_cast<KK_ACC_FLOAT>(x(idx,0)*fOy + x(iH1,0)*fHy + x(iH2,0)*fHy);
        v[4] += static_cast<KK_ACC_FLOAT>(x(idx,0)*fOz + x(iH1,0)*fHz + x(iH2,0)*fHz);
        v[5] += static_cast<KK_ACC_FLOAT>(x(idx,1)*fOz + x(iH1,1)*fHz + x(iH2,1)*fHz);
      }
      vlist[n++] = idx; vlist[n++] = iH1; vlist[n++] = iH2;
    }
  }

  // global/per-atom energy and virial tally for one M-site interaction.
  // key encodes whether i and/or j are O atoms (see Pair::ev_tally_tip4p)
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally_tip4p(EV_FLOAT &ev, const int &key, const int (&vlist)[6],
                      const KK_ACC_FLOAT (&v)[6], const KK_FLOAT &ecoul) const
  {
    if (this->eflag_global) ev.ecoul += static_cast<KK_ACC_FLOAT>(ecoul);
    if (this->vflag_global)
      for (int k = 0; k < 6; k++) ev.v[k] += v[k];

    if (this->eflag_atom) {
      const KK_FLOAT eO = (KK_FLOAT)0.5*ecoul*m_alphaO;
      const KK_FLOAT eH = (KK_FLOAT)0.5*ecoul*m_alphaH;
      const KK_FLOAT eA = (KK_FLOAT)0.5*ecoul;
      if (key == 0) {
        Kokkos::atomic_add(&d_eatom[vlist[0]], (KK_ACC_FLOAT)eA);
        Kokkos::atomic_add(&d_eatom[vlist[1]], (KK_ACC_FLOAT)eA);
      } else if (key == 1) {
        Kokkos::atomic_add(&d_eatom[vlist[0]], (KK_ACC_FLOAT)eO);
        Kokkos::atomic_add(&d_eatom[vlist[1]], (KK_ACC_FLOAT)eH);
        Kokkos::atomic_add(&d_eatom[vlist[2]], (KK_ACC_FLOAT)eH);
        Kokkos::atomic_add(&d_eatom[vlist[3]], (KK_ACC_FLOAT)eA);
      } else if (key == 2) {
        Kokkos::atomic_add(&d_eatom[vlist[0]], (KK_ACC_FLOAT)eA);
        Kokkos::atomic_add(&d_eatom[vlist[1]], (KK_ACC_FLOAT)eO);
        Kokkos::atomic_add(&d_eatom[vlist[2]], (KK_ACC_FLOAT)eH);
        Kokkos::atomic_add(&d_eatom[vlist[3]], (KK_ACC_FLOAT)eH);
      } else {
        Kokkos::atomic_add(&d_eatom[vlist[0]], (KK_ACC_FLOAT)eO);
        Kokkos::atomic_add(&d_eatom[vlist[1]], (KK_ACC_FLOAT)eH);
        Kokkos::atomic_add(&d_eatom[vlist[2]], (KK_ACC_FLOAT)eH);
        Kokkos::atomic_add(&d_eatom[vlist[3]], (KK_ACC_FLOAT)eO);
        Kokkos::atomic_add(&d_eatom[vlist[4]], (KK_ACC_FLOAT)eH);
        Kokkos::atomic_add(&d_eatom[vlist[5]], (KK_ACC_FLOAT)eH);
      }
    }

    if (this->vflag_atom) {
      for (int k = 0; k < 6; k++) {
        const KK_ACC_FLOAT vO = (KK_ACC_FLOAT)0.5*v[k]*(KK_ACC_FLOAT)m_alphaO;
        const KK_ACC_FLOAT vH = (KK_ACC_FLOAT)0.5*v[k]*(KK_ACC_FLOAT)m_alphaH;
        const KK_ACC_FLOAT vA = (KK_ACC_FLOAT)0.5*v[k];
        if (key == 0) {
          Kokkos::atomic_add(&d_vatom(vlist[0],k), vA);
          Kokkos::atomic_add(&d_vatom(vlist[1],k), vA);
        } else if (key == 1) {
          Kokkos::atomic_add(&d_vatom(vlist[0],k), vO);
          Kokkos::atomic_add(&d_vatom(vlist[1],k), vH);
          Kokkos::atomic_add(&d_vatom(vlist[2],k), vH);
          Kokkos::atomic_add(&d_vatom(vlist[3],k), vA);
        } else if (key == 2) {
          Kokkos::atomic_add(&d_vatom(vlist[0],k), vA);
          Kokkos::atomic_add(&d_vatom(vlist[1],k), vO);
          Kokkos::atomic_add(&d_vatom(vlist[2],k), vH);
          Kokkos::atomic_add(&d_vatom(vlist[3],k), vH);
        } else {
          Kokkos::atomic_add(&d_vatom(vlist[0],k), vO);
          Kokkos::atomic_add(&d_vatom(vlist[1],k), vH);
          Kokkos::atomic_add(&d_vatom(vlist[2],k), vH);
          Kokkos::atomic_add(&d_vatom(vlist[3],k), vO);
          Kokkos::atomic_add(&d_vatom(vlist[4],k), vH);
          Kokkos::atomic_add(&d_vatom(vlist[5],k), vH);
        }
      }
    }
  }

  // standard pairwise (LJ) energy/virial tally for a half neighbor list,
  // used by the styles that add LJ interactions on the real atom positions
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j, const KK_FLOAT &evdwl,
                const KK_FLOAT &fpair, const KK_FLOAT &delx, const KK_FLOAT &dely,
                const KK_FLOAT &delz) const
  {
    if (this->eflag_global) ev.evdwl += static_cast<KK_ACC_FLOAT>(evdwl);
    if (this->eflag_atom) {
      Kokkos::atomic_add(&d_eatom[i], (KK_ACC_FLOAT)((KK_FLOAT)0.5*evdwl));
      Kokkos::atomic_add(&d_eatom[j], (KK_ACC_FLOAT)((KK_FLOAT)0.5*evdwl));
    }
    if (this->vflag_global || this->vflag_atom) {
      const KK_FLOAT v0 = delx*delx*fpair;
      const KK_FLOAT v1 = dely*dely*fpair;
      const KK_FLOAT v2 = delz*delz*fpair;
      const KK_FLOAT v3 = delx*dely*fpair;
      const KK_FLOAT v4 = delx*delz*fpair;
      const KK_FLOAT v5 = dely*delz*fpair;
      if (this->vflag_global) {
        ev.v[0] += static_cast<KK_ACC_FLOAT>(v0);
        ev.v[1] += static_cast<KK_ACC_FLOAT>(v1);
        ev.v[2] += static_cast<KK_ACC_FLOAT>(v2);
        ev.v[3] += static_cast<KK_ACC_FLOAT>(v3);
        ev.v[4] += static_cast<KK_ACC_FLOAT>(v4);
        ev.v[5] += static_cast<KK_ACC_FLOAT>(v5);
      }
      if (this->vflag_atom) {
        Kokkos::atomic_add(&d_vatom(i,0), (KK_ACC_FLOAT)((KK_FLOAT)0.5*v0));
        Kokkos::atomic_add(&d_vatom(i,1), (KK_ACC_FLOAT)((KK_FLOAT)0.5*v1));
        Kokkos::atomic_add(&d_vatom(i,2), (KK_ACC_FLOAT)((KK_FLOAT)0.5*v2));
        Kokkos::atomic_add(&d_vatom(i,3), (KK_ACC_FLOAT)((KK_FLOAT)0.5*v3));
        Kokkos::atomic_add(&d_vatom(i,4), (KK_ACC_FLOAT)((KK_FLOAT)0.5*v4));
        Kokkos::atomic_add(&d_vatom(i,5), (KK_ACC_FLOAT)((KK_FLOAT)0.5*v5));
        Kokkos::atomic_add(&d_vatom(j,0), (KK_ACC_FLOAT)((KK_FLOAT)0.5*v0));
        Kokkos::atomic_add(&d_vatom(j,1), (KK_ACC_FLOAT)((KK_FLOAT)0.5*v1));
        Kokkos::atomic_add(&d_vatom(j,2), (KK_ACC_FLOAT)((KK_FLOAT)0.5*v2));
        Kokkos::atomic_add(&d_vatom(j,3), (KK_ACC_FLOAT)((KK_FLOAT)0.5*v3));
        Kokkos::atomic_add(&d_vatom(j,4), (KK_ACC_FLOAT)((KK_FLOAT)0.5*v4));
        Kokkos::atomic_add(&d_vatom(j,5), (KK_ACC_FLOAT)((KK_FLOAT)0.5*v5));
      }
    }
  }

  // ----- long-range (Ewald) Coulomb machinery, used by the *long styles only

  // copy the coulomb interpolation tables to the device
  void init_tables(double cut_coul, double *cut_respa) override
  {
    Pair::init_tables(cut_coul,cut_respa);

    typedef typename AT::t_kkfloat_1d table_type;
    typedef HAT::t_kkfloat_1d host_table_type;

    int ntable = 1;
    for (int i = 0; i < this->ncoultablebits; i++) ntable *= 2;

    tabinnersq_kk = static_cast<KK_FLOAT>(this->tabinnersq);

    auto copy_table = [&](double *src, table_type &dst) {
      host_table_type h_table("HostTable",ntable);
      table_type d_table("DeviceTable",ntable);
      for (int i = 0; i < ntable; i++) h_table(i) = static_cast<KK_FLOAT>(src[i]);
      Kokkos::deep_copy(d_table,h_table);
      dst = d_table;
    };

    copy_table(this->rtable,d_rtable);   copy_table(this->drtable,d_drtable);
    copy_table(this->ftable,d_ftable);   copy_table(this->dftable,d_dftable);
    copy_table(this->ctable,d_ctable);   copy_table(this->dctable,d_dctable);
    copy_table(this->etable,d_etable);   copy_table(this->detable,d_detable);
  }

  // compute() prologue: bind the Ewald parameters not handled by prepare().
  // only instantiated by the *long styles (references their g_ewald member)
  void prepare_coul_long()
  {
    g_ewald_kk = static_cast<KK_FLOAT>(this->g_ewald);
    m_ncoultablebits = this->ncoultablebits;
    m_ncoulmask = this->ncoulmask;
    m_ncoulshiftbits = this->ncoulshiftbits;
  }

  // long-range Coulomb force prefactor at squared distance rsq, and the
  // pairwise energy if want_ecoul is set: analytic erfc() for close pairs
  // (or when tables are disabled), coulomb interpolation tables otherwise
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT coul_long(const KK_FLOAT rsq, const KK_FLOAT qtmp, const KK_FLOAT qj,
                     const KK_FLOAT factor_coul, const bool want_ecoul, KK_FLOAT &ecoul) const
  {
    const KK_FLOAT r2inv = (KK_FLOAT)1.0 / rsq;
    KK_FLOAT cforce;
    if (!m_ncoultablebits || rsq <= tabinnersq_kk) {
      const KK_FLOAT r = Kokkos::sqrt(rsq);
      const KK_FLOAT grij = g_ewald_kk * r;
      const KK_FLOAT expm2 = Kokkos::exp(-grij*grij);
      const KK_FLOAT t = (KK_FLOAT)1.0 / ((KK_FLOAT)1.0 + (KK_FLOAT)EwaldConst::EWALD_P*grij);
      const KK_FLOAT erfc = t*((KK_FLOAT)EwaldConst::A1+t*((KK_FLOAT)EwaldConst::A2+
                            t*((KK_FLOAT)EwaldConst::A3+t*((KK_FLOAT)EwaldConst::A4+
                            t*(KK_FLOAT)EwaldConst::A5))))*expm2;
      const KK_FLOAT prefactor = qqrd2e * qtmp * qj / r;
      KK_FLOAT forcecoul = prefactor * (erfc + (KK_FLOAT)EwaldConst::EWALD_F*grij*expm2);
      if (factor_coul < (KK_FLOAT)1.0) forcecoul -= ((KK_FLOAT)1.0-factor_coul)*prefactor;
      cforce = forcecoul * r2inv;
      if (want_ecoul) {
        ecoul = prefactor * erfc;
        if (factor_coul < (KK_FLOAT)1.0) ecoul -= ((KK_FLOAT)1.0-factor_coul)*prefactor;
      }
    } else {
      typename Pair::union_int_float_t rsq_lookup;
      rsq_lookup.f = rsq;
      const int itable = (rsq_lookup.i & m_ncoulmask) >> m_ncoulshiftbits;
      const KK_FLOAT fraction = ((KK_FLOAT)rsq_lookup.f - d_rtable[itable]) * d_drtable[itable];
      const KK_FLOAT tbl = d_ftable[itable] + fraction*d_dftable[itable];
      KK_FLOAT forcecoul = qtmp * qj * tbl;
      KK_FLOAT prefactor = 0.0;
      if (factor_coul < (KK_FLOAT)1.0) {
        const KK_FLOAT ctbl = d_ctable[itable] + fraction*d_dctable[itable];
        prefactor = qtmp * qj * ctbl;
        forcecoul -= ((KK_FLOAT)1.0-factor_coul)*prefactor;
      }
      cforce = forcecoul * r2inv;
      if (want_ecoul) {
        const KK_FLOAT etbl = d_etable[itable] + fraction*d_detable[itable];
        ecoul = qtmp * qj * etbl;
        if (factor_coul < (KK_FLOAT)1.0) ecoul -= ((KK_FLOAT)1.0-factor_coul)*prefactor;
      }
    }
    return cforce;
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
    qqrd2e = static_cast<KK_FLOAT>(this->force->qqrd2e);
    for (int i = 0; i < 4; i++) {
      special_coul[i] = static_cast<KK_FLOAT>(this->force->special_coul[i]);
      special_lj[i] = static_cast<KK_FLOAT>(this->force->special_lj[i]);
    }

    m_alpha = static_cast<KK_FLOAT>(this->alpha);
    // shares of the M-site force redistributed onto O and each H
    m_alphaO = static_cast<KK_FLOAT>(1.0 - this->alpha);
    m_alphaH = static_cast<KK_FLOAT>(0.5 * this->alpha);
    m_typeO = this->typeO;
    m_typeH = this->typeH;
    m_cut_coulsq = static_cast<KK_FLOAT>(this->cut_coulsq);
    m_cut_coulsqplus = static_cast<KK_FLOAT>((this->cut_coul + 2.0*this->qdist) *
                                             (this->cut_coul + 2.0*this->qdist));

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

    // reset the missing-hydrogen flag, checked in finalize()
    k_h_missing.view_host()() = 0;
    k_h_missing.modify_host();
    k_h_missing.template sync<DeviceType>();
    d_h_missing = k_h_missing.template view<DeviceType>();

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

  // common compute() teardown: stop on missing H atoms, accumulate global
  // coulomb energy and virial, sync per-atom arrays back to the host,
  // mark modified atom data
  void finalize(const EV_FLOAT &ev)
  {
    k_h_missing.template modify<DeviceType>();
    k_h_missing.sync_host();
    if (k_h_missing.view_host()())
      this->error->one(FLERR,"TIP4P hydrogen is missing");

    if (this->eflag_global) this->eng_coul += static_cast<double>(ev.ecoul);
    if (this->vflag_global)
      for (int k = 0; k < 6; k++) this->virial[k] += static_cast<double>(ev.v[k]);
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

  // set on device when a compute kernel needs an M site with a missing H
  DAT::tdual_int_scalar k_h_missing;
  typename AT::t_int_scalar d_h_missing;

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

  int nlocal, nall, eflag, vflag;

  KK_FLOAT special_coul[4];
  KK_FLOAT special_lj[4];
  KK_FLOAT qqrd2e;
  KK_FLOAT m_alpha, m_alphaO, m_alphaH;
  KK_FLOAT m_cut_coulsq, m_cut_coulsqplus;
  int m_typeO, m_typeH;

  // Ewald real-space + optional coulomb interpolation tables (device),
  // filled by init_tables()/prepare_coul_long() for the *long styles only
  typename AT::t_kkfloat_1d d_rtable, d_drtable, d_ftable, d_dftable,
                            d_ctable, d_dctable, d_etable, d_detable;
  KK_FLOAT g_ewald_kk, tabinnersq_kk;
  int m_ncoultablebits, m_ncoulmask, m_ncoulshiftbits;
};

}    // namespace LAMMPS_NS

#endif
