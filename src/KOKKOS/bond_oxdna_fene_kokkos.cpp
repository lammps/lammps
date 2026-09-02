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
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include "bond_oxdna_fene_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neighbor_kokkos.h"
#include "update.h"

#include "fix_oxdna_lrf_kokkos.h"
#include "fix_oxdna_prime_neighs_kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
BondOxdnaFENEKokkos<DeviceType>::BondOxdnaFENEKokkos(LAMMPS *lmp) : BondOxdnaFene(lmp)
{
  kokkosable = 1;

  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  // Internal FixOxdnaLRFKokkos already syncs all read masks that do not
  // change between pair/bond styles.
  datamask_read = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;

  oxdnaflag = EnabledOXDNAFlag::OXDNA;
  fix_oxdna_prime_neighsKK = nullptr;

  d_flag = typename AT::t_int_scalar("bond:flag");
  h_flag = HAT::t_int_scalar("bond:flag_mirror");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
BondOxdnaFENEKokkos<DeviceType>::~BondOxdnaFENEKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void BondOxdnaFENEKokkos<DeviceType>::init_style()
{
  if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 || force->special_lj[3] != 1.0)
    error->all(FLERR, "Must use 'special_bonds lj 0 1 1' with bond style oxdna/fene/kk, oxdna2/fene/kk, " \
                      "oxdna3/fene/kk or oxrna2/fene/kk");

  fix_oxdna_lrfKK = nullptr;
  auto fixes = modify->get_fix_by_style("^OXDNA/LRF/kk");
  if (fixes.size() == 0) error->all(FLERR, "Fix OXDNA/LRF/kk not found. Ensure pair ox*na*/excv/kk is present");
  else fix_oxdna_lrfKK = dynamic_cast<FixOxdnaLRFKokkos<DeviceType> *>(fixes[0]);

  auto prime_fixes = modify->get_fix_by_style("^OXDNA/PRIME_NEIGHS/kk");
  if (prime_fixes.size() == 0)
    fix_oxdna_prime_neighsKK =
      dynamic_cast<FixOxdnaPrimeNeighsKokkos<DeviceType> *>(modify->add_fix("prime_neighs_kk all OXDNA/PRIME_NEIGHS/kk"));
  else
    fix_oxdna_prime_neighsKK = dynamic_cast<FixOxdnaPrimeNeighsKokkos<DeviceType> *>(prime_fixes[0]);

  if (!fix_oxdna_prime_neighsKK)
    error->all(FLERR, "Fix OXDNA/PRIME_NEIGHS/kk not found");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void BondOxdnaFENEKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"bond:eatom");
    d_eatom = k_eatom.template view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"bond:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);

  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK | TORQUE_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.template view<DeviceType>();
  atomtype = atomKK->k_type.template view<DeviceType>();

  neighborKK->k_bondlist.template sync<DeviceType>();
  bondlist = neighborKK->k_bondlist.view<DeviceType>();
  nbondlist = neighborKK->nbondlist;
  d_prime_neighs_bond = fix_oxdna_prime_neighsKK->d_prime_neighs_bond;
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;

  // d_n(x/y/z)_xtrct = extracted local unit vectors in lab frame from fix_oxdna_lrf_kokkos.
  d_nx_xtrct = fix_oxdna_lrfKK->k_nx.template view<DeviceType>();
  d_ny_xtrct = fix_oxdna_lrfKK->k_ny.template view<DeviceType>();
  d_nz_xtrct = fix_oxdna_lrfKK->k_nz.template view<DeviceType>();

  copymode = 1;

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      if (oxdnaflag == OXDNA) {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagBondOxdnaFENECompute<OXDNA,1,1> >(0,nbondlist),*this,ev);
      } else if (oxdnaflag == OXDNA2) {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagBondOxdnaFENECompute<OXDNA2,1,1> >(0,nbondlist),*this,ev);
      } else if (oxdnaflag == OXRNA2) {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagBondOxdnaFENECompute<OXRNA2,1,1> >(0,nbondlist),*this,ev);
      }
    } else {
      if (oxdnaflag == OXDNA) {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagBondOxdnaFENECompute<OXDNA,0,1> >(0,nbondlist),*this,ev);
      } else if (oxdnaflag == OXDNA2) {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagBondOxdnaFENECompute<OXDNA2,0,1> >(0,nbondlist),*this,ev);
      } else if (oxdnaflag == OXRNA2) {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagBondOxdnaFENECompute<OXRNA2,0,1> >(0,nbondlist),*this,ev);
      }
    }
  } else {
    if (newton_bond) {
      if (oxdnaflag == OXDNA) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagBondOxdnaFENECompute<OXDNA,1,0> >(0,nbondlist),*this);
      } else if (oxdnaflag == OXDNA2) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagBondOxdnaFENECompute<OXDNA2,1,0> >(0,nbondlist),*this);
      } else if (oxdnaflag == OXRNA2) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagBondOxdnaFENECompute<OXRNA2,1,0> >(0,nbondlist),*this);
      }
    } else {
      if (oxdnaflag == OXDNA) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagBondOxdnaFENECompute<OXDNA,0,0> >(0,nbondlist),*this);
      } else if (oxdnaflag == OXDNA2) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagBondOxdnaFENECompute<OXDNA2,0,0> >(0,nbondlist),*this);
      } else if (oxdnaflag == OXRNA2) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagBondOxdnaFENECompute<OXRNA2,0,0> >(0,nbondlist),*this);
      }
    }
  }

  // The "FENE bond too long" flag is a benign diagnostic: the bond force is
  // capped every step inside the kernel regardless. Copying the device flag
  // back to the host every step forces a sync point that costs more than the
  // diagnostic is worth. Only read (and warn/reset) on output/thermo steps;
  // the flag accumulates across the intervening steps, so the warning reports
  // "at or before step N".
  if (eflag || vflag) {
    Kokkos::deep_copy(h_flag,d_flag);
    if (h_flag() == 1) {
      error->warning(FLERR,"FENE bond too long at or before step {}", update->ntimestep);
      Kokkos::deep_copy(d_flag,0);
    }
  }

  if (eflag_global) energy += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  copymode = 0;
}

template<class DeviceType>
template<int OXDNAFLAG, int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void BondOxdnaFENEKokkos<DeviceType>::operator()(TagBondOxdnaFENECompute<OXDNAFLAG,NEWTON_BOND,EVFLAG>, \
  const int &in, EV_FLOAT& ev) const
{
  // The f and torque arrays are atomic
  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,\
    typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;
  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,\
    typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_torque = torque;

  // Use precomputed bond and prime neighbors.
  // NOTE: already in correct order from precompute, so directionality test: a -> b is 3' -> 5' is already satisfied
  int a = d_prime_neighs_bond(in,0);
  int b = d_prime_neighs_bond(in,1);
  const int type = bondlist(in,2);
  int a3ptype, atype, btype, b5ptype;    // tetramer types

  // determine tetramer types
  // Our prime_neighs_bond ordering (a,b,id3p[a],id5p[b]) from precompute
  // is assigned such that we preserve the vanilla oxDNA convention of:
  // 3'neighbor a - a - b - 5'neighbor b
  // throughout the rest of compute.
  int id3p_local = d_prime_neighs_bond(in,2);
  a3ptype = (id3p_local != -1) ? atomtype[id3p_local] : 0;

  atype = atomtype[a];
  btype = atomtype[b];

  int id5p_local = d_prime_neighs_bond(in,3);
  b5ptype = (id5p_local != -1) ? atomtype[id5p_local] : 0;

  KK_ACC_FLOAT delf[3], delta[3], deltb[3];    // force, torque increment
  KK_FLOAT delr_bkbk[3];                   // vector backbone site b to a
  // vectors COM-backbone site in lab frame
  KK_FLOAT ra_cbk[3], rb_cbk[3];

  // vector COM-backbone site a and b - "compute_interaction_sites" vector COM-sugar-phosphate backbone in oxDNA
  if constexpr (OXDNAFLAG==OXDNA) {
    constexpr KK_FLOAT d_cs = -0.4;
    ra_cbk[0] = d_cs * d_nx_xtrct(a,0);
    ra_cbk[1] = d_cs * d_nx_xtrct(a,1);
    ra_cbk[2] = d_cs * d_nx_xtrct(a,2);
    rb_cbk[0] = d_cs * d_nx_xtrct(b,0);
    rb_cbk[1] = d_cs * d_nx_xtrct(b,1);
    rb_cbk[2] = d_cs * d_nx_xtrct(b,2);
  } else if constexpr (OXDNAFLAG==OXDNA2) {
    constexpr KK_FLOAT d_cs_x = -0.34;
    constexpr KK_FLOAT d_cs_y = +0.3408;
    ra_cbk[0] = d_cs_x * d_nx_xtrct(a,0) + d_cs_y * d_ny_xtrct(a,0);
    ra_cbk[1] = d_cs_x * d_nx_xtrct(a,1) + d_cs_y * d_ny_xtrct(a,1);
    ra_cbk[2] = d_cs_x * d_nx_xtrct(a,2) + d_cs_y * d_ny_xtrct(a,2);
    rb_cbk[0] = d_cs_x * d_nx_xtrct(b,0) + d_cs_y * d_ny_xtrct(b,0);
    rb_cbk[1] = d_cs_x * d_nx_xtrct(b,1) + d_cs_y * d_ny_xtrct(b,1);
    rb_cbk[2] = d_cs_x * d_nx_xtrct(b,2) + d_cs_y * d_ny_xtrct(b,2);
  } else {
    // OXRNA2
    constexpr KK_FLOAT d_cs_x = -0.4;
    constexpr KK_FLOAT d_cs_z = +0.2;
    ra_cbk[0] = d_cs_x * d_nx_xtrct(a,0) + d_cs_z * d_nz_xtrct(a,0);
    ra_cbk[1] = d_cs_x * d_nx_xtrct(a,1) + d_cs_z * d_nz_xtrct(a,1);
    ra_cbk[2] = d_cs_x * d_nx_xtrct(a,2) + d_cs_z * d_nz_xtrct(a,2);
    rb_cbk[0] = d_cs_x * d_nx_xtrct(b,0) + d_cs_z * d_nz_xtrct(b,0);
    rb_cbk[1] = d_cs_x * d_nx_xtrct(b,1) + d_cs_z * d_nz_xtrct(b,1);
    rb_cbk[2] = d_cs_x * d_nx_xtrct(b,2) + d_cs_z * d_nz_xtrct(b,2);
  }

  // vector backbone site b to a
  delr_bkbk[0] = x(a,0) + ra_cbk[0] - x(b,0) - rb_cbk[0];
  delr_bkbk[1] = x(a,1) + ra_cbk[1] - x(b,1) - rb_cbk[1];
  delr_bkbk[2] = x(a,2) + ra_cbk[2] - x(b,2) - rb_cbk[2];
  const KK_FLOAT rsq = delr_bkbk[0]*delr_bkbk[0] + delr_bkbk[1]*delr_bkbk[1] + delr_bkbk[2]*delr_bkbk[2];
  const KK_FLOAT r_bkbk = sqrtf(rsq);

  KK_FLOAT rr0 = r_bkbk - d_r0(type, a3ptype, atype, btype, b5ptype);
  const KK_FLOAT rr0sq = rr0 * rr0;
  const KK_FLOAT Deltasq = d_Delta(type, a3ptype, atype, btype, b5ptype)*d_Delta(type, a3ptype, atype, btype, b5ptype);
  KK_FLOAT rlogarg = 1.0 - rr0sq/Deltasq;

  // energy

  KK_FLOAT ebond = 0.0;
  if (eflag) { ebond = -0.5*d_k[type]*log(rlogarg);}

  // switching to capped force for r-r0 -> Delta at
  // r > r_max = r0 + Delta*sqrtf(1-rlogarg) OR
  // r < r_min = r0 - Delta*sqrtf(1-rlogarg)
  if (rlogarg < 0.2) { // rlogarg_min = 0.2
    // issue warning, reset rlogarg and rr0 to cap force
    d_flag() = 1;
    rlogarg = 0.2;
    // if overstretched F(r)=F(r_max)=F_max, E(r)=E(r_max)+F_max*(r-r_max)
    if (r_bkbk > d_r0(type, a3ptype, atype, btype, b5ptype)) {
      rr0 = d_Delta(type, a3ptype, atype, btype, b5ptype)*sqrtf(1.0 - rlogarg);
      // energy
      if (eflag) {
        ebond = -0.5 * d_k(type) * log(rlogarg) + d_k(type) *
                sqrtf(1.0-rlogarg) / rlogarg / d_Delta(type, a3ptype, atype, btype, b5ptype) *
                (r_bkbk - d_r0(type, a3ptype, atype, btype, b5ptype) -
                d_Delta(type, a3ptype, atype, btype, b5ptype) * sqrtf(1.0-rlogarg));
      }
    }
    // if overcompressed F(r)=F(r_min)=F_max, E(r)=E(r_min)+F_max*(r_min-r)
    else if (r_bkbk < d_r0(type, a3ptype, atype, btype, b5ptype)) {
      rr0 = -d_Delta(type, a3ptype, atype, btype, b5ptype)*sqrtf(1.0 - rlogarg);
      // energy
      if (eflag) {
        ebond = -0.5 * d_k(type) * log(rlogarg) + d_k(type) *
                sqrtf(1.0-rlogarg) / rlogarg / d_Delta(type, a3ptype, atype, btype, b5ptype) *
                (r_bkbk - d_r0(type, a3ptype, atype, btype, b5ptype) +
                d_Delta(type, a3ptype, atype, btype, b5ptype) * sqrtf(1.0-rlogarg));
      }
    }
  }

  KK_ACC_FLOAT fbond = -d_k[type] * rr0 / rlogarg / Deltasq / r_bkbk;
  delf[0] = delr_bkbk[0] * fbond;
  delf[1] = delr_bkbk[1] * fbond;
  delf[2] = delr_bkbk[2] * fbond;

  // apply force to each of 2 atoms

  if (NEWTON_BOND || a < nlocal) {
    a_f(a,0) += delf[0];
    a_f(a,1) += delf[1];
    a_f(a,2) += delf[2];
    delta[0] = ra_cbk[1]*delf[2] - ra_cbk[2]*delf[1];
    delta[1] = ra_cbk[2]*delf[0] - ra_cbk[0]*delf[2];
    delta[2] = ra_cbk[0]*delf[1] - ra_cbk[1]*delf[0];
    a_torque(a,0) += delta[0];
    a_torque(a,1) += delta[1];
    a_torque(a,2) += delta[2];
  }

  if (NEWTON_BOND || b < nlocal) {
    a_f(b,0) -= delf[0];
    a_f(b,1) -= delf[1];
    a_f(b,2) -= delf[2];
    deltb[0] = rb_cbk[1]*delf[2] - rb_cbk[2]*delf[1];
    deltb[1] = rb_cbk[2]*delf[0] - rb_cbk[0]*delf[2];
    deltb[2] = rb_cbk[0]*delf[1] - rb_cbk[1]*delf[0];
    a_torque(b,0) -= deltb[0];
    a_torque(b,1) -= deltb[1];
    a_torque(b,2) -= deltb[2];
  }

  if (EVFLAG) { ev_tally_xyz(ev, a, b, nlocal, NEWTON_BOND, ebond, delf[0], delf[1], delf[2], \
    x(a,0)-x(b,0), x(a,1)-x(b,1), x(a,2)-x(b,2)); }

}

template<class DeviceType>
template<int OXDNAFLAG, int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void BondOxdnaFENEKokkos<DeviceType>::operator()(TagBondOxdnaFENECompute<OXDNAFLAG,NEWTON_BOND,EVFLAG>, const int &in) const {
  EV_FLOAT ev;
  this->template operator()<OXDNAFLAG,NEWTON_BOND,EVFLAG>(TagBondOxdnaFENECompute<OXDNAFLAG,NEWTON_BOND,EVFLAG>(), in, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void BondOxdnaFENEKokkos<DeviceType>::allocate()
{
  BondOxdnaFene::allocate();

  int n = atom->nbondtypes;
  int m = atom->ntypes;
  k_k = DAT::tdual_kkfloat_1d("BondOxdnaFENE::k",n+1);
  k_r0 = DAT::tdual_kkfloat_5d("BondOxdnaFENE::r0",n+1,m+1,m+1,m+1,m+1);
  k_Delta = DAT::tdual_kkfloat_5d("BondOxdnaFENE::Delta",n+1,m+1,m+1,m+1,m+1);

  d_k = k_k.template view<DeviceType>();
  d_r0 = k_r0.template view<DeviceType>();
  d_Delta = k_Delta.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

template<class DeviceType>
void BondOxdnaFENEKokkos<DeviceType>::coeff(int narg, char **arg)
{
  BondOxdnaFene::coeff(narg, arg);

  // Unlike vanilla, we don't use the bounds and assert - args have already
  // been parsed.

  int m = atom->nbondtypes;
  int n = atom->ntypes;
  for (int i = 1; i <= m; i++) {
    k_k.view_host()[i] = k[i];
    for (int n1 = 0; n1 <= n; n1++) {
      for (int n2 = 0; n2 <= n; n2++) {
        for (int n3 = 0; n3 <= n; n3++) {
          for (int n4 = 0; n4 <= n; n4++) {
            k_r0.view_host()(i,n1,n2,n3,n4) = r0[i][n1][n2][n3][n4];
            k_Delta.view_host()(i,n1,n2,n3,n4) = Delta[i][n1][n2][n3][n4];
          }
        }
      }
    }
  }

  k_k.template modify<LMPHostType>();
  k_r0.template modify<LMPHostType>();
  k_Delta.template modify<LMPHostType>();

  // Sync to device
  k_k.template sync<DeviceType>();
  k_r0.template sync<DeviceType>();
  k_Delta.template sync<DeviceType>();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void BondOxdnaFENEKokkos<DeviceType>::read_restart(FILE *fp)
{
  BondOxdnaFene::read_restart(fp);

  int n = atom->nbondtypes;
  int m = atom->ntypes;
  for (int i = 1; i <= n; i++) {
    k_k.view_host()[i] = k[i];
    for (int n1 = 0; n1 <= m; n1++) {
      for (int n2 = 0; n2 <= m; n2++) {
        for (int n3 = 0; n3 <= m; n3++) {
          for (int n4 = 0; n4 <= m; n4++) {
            k_r0.view_host()(i,n1,n2,n3,n4) = r0[i][n1][n2][n3][n4];
            k_Delta.view_host()(i,n1,n2,n3,n4) = Delta[i][n1][n2][n3][n4];
          }
        }
      }
    }
  }

  k_k.template modify<LMPHostType>();
  k_r0.template modify<LMPHostType>();
  k_Delta.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void BondOxdnaFENEKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,\
      const int &nlocal, const int &newton_bond,\
      const KK_FLOAT &ebond, const KK_ACC_FLOAT &fx, const KK_ACC_FLOAT &fy, const KK_ACC_FLOAT &fz,\
      const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const
{
  KK_ACC_FLOAT ebondhalf;
  KK_ACC_FLOAT v[6];

  // The eatom and vatom arrays are atomic
  Kokkos::View<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout, \
      typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_eatom = d_eatom;
  Kokkos::View<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout, \
      typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_vatom = d_vatom;

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += static_cast<KK_ACC_FLOAT>(ebond);
      else {
        ebondhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*ebond);
        if (i < nlocal) ev.evdwl += ebondhalf;
        if (j < nlocal) ev.evdwl += ebondhalf;
      }
    }
    if (eflag_atom) {
      ebondhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*ebond);
      if (newton_bond || i < nlocal) v_eatom[i] += ebondhalf;
      if (newton_bond || j < nlocal) v_eatom[j] += ebondhalf;
    }
  }

  if (vflag_either) {
    v[0] = static_cast<KK_ACC_FLOAT>(delx * fx);
    v[1] = static_cast<KK_ACC_FLOAT>(dely * fy);
    v[2] = static_cast<KK_ACC_FLOAT>(delz * fz);
    v[3] = static_cast<KK_ACC_FLOAT>(delx * fy);
    v[4] = static_cast<KK_ACC_FLOAT>(delx * fz);
    v[5] = static_cast<KK_ACC_FLOAT>(dely * fz);

    if (vflag_global) {
      if (newton_bond) {
        ev.v[0] += v[0];
        ev.v[1] += v[1];
        ev.v[2] += v[2];
        ev.v[3] += v[3];
        ev.v[4] += v[4];
        ev.v[5] += v[5];
      } else {
        if (i < nlocal) {
          ev.v[0] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[0]);
          ev.v[1] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[1]);
          ev.v[2] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[2]);
          ev.v[3] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[3]);
          ev.v[4] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[4]);
          ev.v[5] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[5]);
        }
        if (j < nlocal) {
          ev.v[0] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[0]);
          ev.v[1] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[1]);
          ev.v[2] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[2]);
          ev.v[3] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[3]);
          ev.v[4] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[4]);
          ev.v[5] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[5]);
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        v_vatom(i,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[0]);
        v_vatom(i,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[1]);
        v_vatom(i,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[2]);
        v_vatom(i,3) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[3]);
        v_vatom(i,4) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[4]);
        v_vatom(i,5) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[5]);
      }
      if (newton_bond || j < nlocal) {
        v_vatom(j,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[0]);
        v_vatom(j,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[1]);
        v_vatom(j,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[2]);
        v_vatom(j,3) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[3]);
        v_vatom(j,4) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[4]);
        v_vatom(j,5) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v[5]);
      }
    }
  }
}

namespace LAMMPS_NS {
template class BondOxdnaFENEKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class BondOxdnaFENEKokkos<LMPHostType>;
#endif
}
