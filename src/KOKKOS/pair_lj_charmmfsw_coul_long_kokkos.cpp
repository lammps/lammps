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
 
   Contributing authors:
 
   - Ray Shan (SNL) - original PairLJCharmmCoulLongKokkos
   
   - Mitch Murphy (alphataubio) - PairLJCharmmfswCoulLongKokkos update (2023/12)

   Based on serial kspace lj-fsw sections (force-switched) provided by
   Robert Meissner and Lucio Colombi Ciacchi of Bremen University, Germany,
   with additional assistance from Robert A. Latour, Clemson University

 ------------------------------------------------------------------------- */



/* ----------------------------------------------------------------------

 *** DRAFT VERSION 1 (lots of comments to be removed just before merge) ***
 
 (1) first draft version of PairLJCharmmfswCoulLongKokkos almost exactly
 same as PairLJCharmmCoulLongKokkos but with new class name
 
 method: track changes from serial kspace pair_lj_charmm_coul_long to
 pair_lj_charmmfsw_coul_long and apply to PairLJCharmmCoulLongKokkos
 
 ISSUES:
 
 (A) charmm denom_lj_inv cache , is it to optimize code because division
 is slower that multiplication ??
 
 

 ------------------------------------------------------------------------- */


/*
 19c23
 < #include "pair_lj_charmm_coul_long.h"
 ---
 > #include "pair_lj_charmmfsw_coul_long.h"

 */

#include "pair_lj_charmmfsw_coul_long_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;


#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

/*
 47c51
 < PairLJCharmmCoulLong::PairLJCharmmCoulLong(LAMMPS *lmp) : Pair(lmp)
 ---
 > PairLJCharmmfswCoulLong::PairLJCharmmfswCoulLong(LAMMPS *lmp) : Pair(lmp)
 55a60,72
 >
 >   // short-range/long-range flag accessed by DihedralCharmmfsw
 >
 >   dihedflag = 1;
 >
 >   // switch qqr2e from LAMMPS value to CHARMM value
 >
 >   if (strcmp(update->unit_style,"real") == 0) {
 >     if ((comm->me == 0) && (force->qqr2e != force->qqr2e_charmm_real))
 >       error->message(FLERR,"Switching to CHARMM coulomb energy"
 >                      " conversion constant");
 >     force->qqr2e = force->qqr2e_charmm_real;
 >   }

 */

// added superclass constructor to inherit from PairLJCharmmfswCoulLong

template<class DeviceType>
PairLJCharmmfswCoulLongKokkos<DeviceType>::PairLJCharmmfswCoulLongKokkos(LAMMPS *lmp):PairLJCharmmfswCoulLong(lmp)
{
  
  // pair_lj_charmmfsw_coul_long_kokkos.cpp:112:28: error: qualified reference to 'PairLJCharmmfswCoulLong' is a constructor name rather than a type in this context
  // ??? PairLJCharmmfswCoulLong::PairLJCharmmfswCoulLong(lmp);

  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

/*
 
 60c77
 < PairLJCharmmCoulLong::~PairLJCharmmCoulLong()
 ---
 > PairLJCharmmfswCoulLong::~PairLJCharmmfswCoulLong()
 61a79,87
 >   // switch qqr2e back from CHARMM value to LAMMPS value
 >
 >   if (update && strcmp(update->unit_style,"real") == 0) {
 >     if ((comm->me == 0) && (force->qqr2e == force->qqr2e_charmm_real))
 >       error->message(FLERR,"Restoring original LAMMPS coulomb energy"
 >                      " conversion constant");
 >     force->qqr2e = force->qqr2e_lammps_real;
 >   }
 >

 */

// added superclass constructor to inherit from PairLJCharmmfswCoulLong

template<class DeviceType>
PairLJCharmmfswCoulLongKokkos<DeviceType>::~PairLJCharmmfswCoulLongKokkos()
{
  
  // pair_lj_charmmfsw_coul_long_kokkos.cpp:150:28: error: qualified reference to 'PairLJCharmmfswCoulLong' is a constructor name rather than a type in this context
  // ??? PairLJCharmmfswCoulLong::PairLJCharmmfswCoulLong();
  
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->destroy_kokkos(k_cutsq,cutsq);
  }
}

/* ---------------------------------------------------------------------- */

/*
 87c112
 < void PairLJCharmmCoulLong::compute(int eflag, int vflag)
 ---
 > void PairLJCharmmfswCoulLong::compute(int eflag, int vflag)
 90c115
 <   double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
 ---
 >   double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,evdwl12,evdwl6,ecoul,fpair;
 92c117
 <   double r,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
 ---
 >   double r,rinv,r2inv,r3inv,r6inv,rsq,forcecoul,forcelj,factor_coul,factor_lj;
 94c119
 <   double philj,switch1,switch2;
 ---
 >   double switch1;
 96d120
 <   double rsq;
 174,179c198,200
 <               (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) * denom_lj_inv;
 <             switch2 = 12.0*rsq * (cut_ljsq-rsq) *
 <               (rsq-cut_lj_innersq) * denom_lj_inv;
 <             philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
 <             forcelj = forcelj*switch1 + philj*switch2;
 <           }
 ---
 >               (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
 >             forcelj = forcelj*switch1;
 >          }
 205d225
 <             evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);
 207,209c227,240
 <               switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
 <                 (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) * denom_lj_inv;
 <               evdwl *= switch1;
 ---
 >               r = sqrt(rsq);
 >               rinv = 1.0/r;
 >               r3inv = rinv*rinv*rinv;
 >               evdwl12 = lj3[itype][jtype]*cut_lj6*denom_lj12 *
 >                 (r6inv - cut_lj6inv)*(r6inv - cut_lj6inv);
 >               evdwl6 = -lj4[itype][jtype]*cut_lj3*denom_lj6 *
 >                 (r3inv - cut_lj3inv)*(r3inv - cut_lj3inv);
 >               evdwl = evdwl12 + evdwl6;
 >             } else {
 >               evdwl12 = r6inv*lj3[itype][jtype]*r6inv -
 >                 lj3[itype][jtype]*cut_lj_inner6inv*cut_lj6inv;
 >               evdwl6 = -lj4[itype][jtype]*r6inv +
 >                 lj4[itype][jtype]*cut_lj_inner3inv*cut_lj3inv;
 >               evdwl = evdwl12 + evdwl6;

 */

template<class DeviceType>
void PairLJCharmmfswCoulLongKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

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
  k_cutsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];
  special_coul[0] = force->special_coul[0];
  special_coul[1] = force->special_coul[1];
  special_coul[2] = force->special_coul[2];
  special_coul[3] = force->special_coul[3];
  qqrd2e = force->qqrd2e;
  newton_pair = force->newton_pair;

  // loop over neighbors of my atoms

  copymode = 1;

  EV_FLOAT ev;
  if (ncoultablebits)
    ev = pair_compute<PairLJCharmmfswCoulLongKokkos<DeviceType>,CoulLongTable<1> >
      (this,(NeighListKokkos<DeviceType>*)list);
  else
    ev = pair_compute<PairLJCharmmfswCoulLongKokkos<DeviceType>,CoulLongTable<0> >
      (this,(NeighListKokkos<DeviceType>*)list);


  if (eflag) {
    eng_vdwl += ev.evdwl;
    eng_coul += ev.ecoul;
  }
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

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;
}

/* ----------------------------------------------------------------------
   compute LJ CHARMM pair force between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairLJCharmmfswCoulLongKokkos<DeviceType>::
compute_fpair(const F_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const {
  const F_FLOAT r2inv = 1.0/rsq;
  const F_FLOAT r6inv = r2inv*r2inv*r2inv;
  F_FLOAT forcelj, switch1, switch2, englj;

  forcelj = r6inv *
    ((STACKPARAMS?m_params[itype][jtype].lj1:params(itype,jtype).lj1)*r6inv -
     (STACKPARAMS?m_params[itype][jtype].lj2:params(itype,jtype).lj2));

  if (rsq > cut_lj_innersq) {
    switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
              (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
    switch2 = 12.0*rsq * (cut_ljsq-rsq) * (rsq-cut_lj_innersq) / denom_lj;
    englj = r6inv *
            ((STACKPARAMS?m_params[itype][jtype].lj3:params(itype,jtype).lj3)*r6inv -
             (STACKPARAMS?m_params[itype][jtype].lj4:params(itype,jtype).lj4));
    forcelj = forcelj*switch1 + englj*switch2;
  }

  return forcelj*r2inv;
}

/* ----------------------------------------------------------------------
   compute LJ CHARMM pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairLJCharmmfswCoulLongKokkos<DeviceType>::
compute_evdwl(const F_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const {
  const F_FLOAT r2inv = 1.0/rsq;
  const F_FLOAT r6inv = r2inv*r2inv*r2inv;
  F_FLOAT englj, switch1;

  englj = r6inv *
    ((STACKPARAMS?m_params[itype][jtype].lj3:params(itype,jtype).lj3)*r6inv -
     (STACKPARAMS?m_params[itype][jtype].lj4:params(itype,jtype).lj4));

  if (rsq > cut_lj_innersq) {
    switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
      (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
    englj *= switch1;
  }

  return englj;

}

/* ----------------------------------------------------------------------
   compute coulomb pair force between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS,  class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairLJCharmmfswCoulLongKokkos<DeviceType>::
compute_fcoul(const F_FLOAT& rsq, const int& /*i*/, const int&j,
              const int& /*itype*/, const int& /*jtype*/,
              const F_FLOAT& factor_coul, const F_FLOAT& qtmp) const {
  if (Specialisation::DoTable && rsq > tabinnersq) {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    const int itable = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
    const F_FLOAT fraction = (rsq_lookup.f - d_rtable[itable]) * d_drtable[itable];
    const F_FLOAT table = d_ftable[itable] + fraction*d_dftable[itable];
    F_FLOAT forcecoul = qtmp*q[j] * table;
    if (factor_coul < 1.0) {
      const F_FLOAT table = d_ctable[itable] + fraction*d_dctable[itable];
      const F_FLOAT prefactor = qtmp*q[j] * table;
      forcecoul -= (1.0-factor_coul)*prefactor;
    }
    return forcecoul/rsq;
  } else {
    const F_FLOAT r = sqrt(rsq);
    const F_FLOAT grij = g_ewald * r;
    const F_FLOAT expm2 = exp(-grij*grij);
    const F_FLOAT t = 1.0 / (1.0 + EWALD_P*grij);
    const F_FLOAT rinv = 1.0/r;
    const F_FLOAT erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
    const F_FLOAT prefactor = qqrd2e * qtmp*q[j]*rinv;
    F_FLOAT forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
    if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;

    return forcecoul*rinv*rinv;
  }
}

/* ----------------------------------------------------------------------
   compute coulomb pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairLJCharmmfswCoulLongKokkos<DeviceType>::
compute_ecoul(const F_FLOAT& rsq, const int& /*i*/, const int&j,
              const int& /*itype*/, const int& /*jtype*/, const F_FLOAT& factor_coul, const F_FLOAT& qtmp) const {
  if (Specialisation::DoTable && rsq > tabinnersq) {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    const int itable = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
    const F_FLOAT fraction = (rsq_lookup.f - d_rtable[itable]) * d_drtable[itable];
    const F_FLOAT table = d_etable[itable] + fraction*d_detable[itable];
    F_FLOAT ecoul = qtmp*q[j] * table;
    if (factor_coul < 1.0) {
      const F_FLOAT table = d_ctable[itable] + fraction*d_dctable[itable];
      const F_FLOAT prefactor = qtmp*q[j] * table;
      ecoul -= (1.0-factor_coul)*prefactor;
    }
    return ecoul;
  } else {
    const F_FLOAT r = sqrt(rsq);
    const F_FLOAT grij = g_ewald * r;
    const F_FLOAT expm2 = exp(-grij*grij);
    const F_FLOAT t = 1.0 / (1.0 + EWALD_P*grij);
    const F_FLOAT erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
    const F_FLOAT prefactor = qqrd2e * qtmp*q[j]/r;
    F_FLOAT ecoul = prefactor * erfc;
    if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
    return ecoul;
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCharmmfswCoulLongKokkos<DeviceType>::allocate()
{
  PairLJCharmmfswCoulLong::allocate();

  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  d_cut_ljsq = typename AT::t_ffloat_2d("pair:cut_ljsq",n+1,n+1);

  d_cut_coulsq = typename AT::t_ffloat_2d("pair:cut_coulsq",n+1,n+1);

  k_params = Kokkos::DualView<params_lj_coul**,Kokkos::LayoutRight,DeviceType>("PairLJCharmmCoulLong::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

template<class DeviceType>
void PairLJCharmmfswCoulLongKokkos<DeviceType>::init_tables(double cut_coul, double *cut_respa)
{
  Pair::init_tables(cut_coul,cut_respa);

  typedef typename ArrayTypes<DeviceType>::t_ffloat_1d table_type;
  typedef typename ArrayTypes<LMPHostType>::t_ffloat_1d host_table_type;

  int ntable = 1;
  for (int i = 0; i < ncoultablebits; i++) ntable *= 2;


  // Copy rtable and drtable
  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);
  for (int i = 0; i < ntable; i++) {
    h_table(i) = rtable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_rtable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);
  for (int i = 0; i < ntable; i++) {
    h_table(i) = drtable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_drtable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  // Copy ftable and dftable
  for (int i = 0; i < ntable; i++) {
    h_table(i) = ftable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_ftable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  for (int i = 0; i < ntable; i++) {
    h_table(i) = dftable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_dftable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  // Copy ctable and dctable
  for (int i = 0; i < ntable; i++) {
    h_table(i) = ctable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_ctable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  for (int i = 0; i < ntable; i++) {
    h_table(i) = dctable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_dctable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  // Copy etable and detable
  for (int i = 0; i < ntable; i++) {
    h_table(i) = etable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_etable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  for (int i = 0; i < ntable; i++) {
    h_table(i) = detable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_detable = d_table;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

/*
 682c733
 < void PairLJCharmmCoulLong::init_style()
 ---
 > void PairLJCharmmfswCoulLong::init_style()
 686c737
 <                "Pair style lj/charmm/coul/long requires atom attribute q");
 ---
 >                "Pair style lj/charmmfsw/coul/long requires atom attribute q");
 688c739
 <   // request regular or rRESPA neighbor list
 ---
 >   // request regular or rRESPA neighbor lists
 705a757,766
 >   cut_ljinv = 1.0/cut_lj;
 >   cut_lj_innerinv = 1.0/cut_lj_inner;
 >   cut_lj3 = cut_lj * cut_lj * cut_lj;
 >   cut_lj3inv = cut_ljinv * cut_ljinv * cut_ljinv;
 >   cut_lj_inner3inv = cut_lj_innerinv * cut_lj_innerinv * cut_lj_innerinv;
 >   cut_lj_inner3 = cut_lj_inner * cut_lj_inner * cut_lj_inner;
 >   cut_lj6 = cut_ljsq * cut_ljsq * cut_ljsq;
 >   cut_lj6inv = cut_lj3inv * cut_lj3inv;
 >   cut_lj_inner6inv = cut_lj_inner3inv * cut_lj_inner3inv;
 >   cut_lj_inner6 = cut_lj_innersq * cut_lj_innersq * cut_lj_innersq;
 709,711c770,773
 <   denom_lj = ( (cut_ljsq-cut_lj_innersq) * (cut_ljsq-cut_lj_innersq) *
 <                (cut_ljsq-cut_lj_innersq) );
 <   denom_lj_inv = 1.0 / denom_lj;
 ---
 >   denom_lj = (cut_ljsq-cut_lj_innersq) * (cut_ljsq-cut_lj_innersq) *
 >     (cut_ljsq-cut_lj_innersq);
 >   denom_lj12 = 1.0/(cut_lj6 - cut_lj_inner6);
 >   denom_lj6 = 1.0/(cut_lj3 - cut_lj_inner3);
 718,730d779
 <     cut_in_off = cut_respa[0];
 <     cut_in_on = cut_respa[1];
 <     cut_out_on = cut_respa[2];
 <     cut_out_off = cut_respa[3];
 <
 <     cut_in_diff = cut_in_on - cut_in_off;
 <     cut_out_diff = cut_out_off - cut_out_on;
 <     cut_in_diff_inv = 1.0 / (cut_in_diff);
 <     cut_out_diff_inv = 1.0 / (cut_out_diff);
 <     cut_in_off_sq = cut_in_off*cut_in_off;
 <     cut_in_on_sq = cut_in_on*cut_in_on;
 <     cut_out_on_sq = cut_out_on*cut_out_on;
 <     cut_out_off_sq = cut_out_off*cut_out_off;

 */

template<class DeviceType>
void PairLJCharmmfswCoulLongKokkos<DeviceType>::init_style()
{
  PairLJCharmmfswCoulLong::init_style();

  Kokkos::deep_copy(d_cut_ljsq,cut_ljsq);
  Kokkos::deep_copy(d_cut_coulsq,cut_coulsq);

  // error if rRESPA with inner levels

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa)
      error->all(FLERR,"Cannot use Kokkos pair style with rRESPA inner/middle");
  }

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairLJCharmmfswCoulLongKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairLJCharmmfswCoulLong::init_one(i,j);

  k_params.h_view(i,j).lj1 = lj1[i][j];
  k_params.h_view(i,j).lj2 = lj2[i][j];
  k_params.h_view(i,j).lj3 = lj3[i][j];
  k_params.h_view(i,j).lj4 = lj4[i][j];
  //k_params.h_view(i,j).offset = offset[i][j];
  k_params.h_view(i,j).cut_ljsq = cut_ljsq;
  k_params.h_view(i,j).cut_coulsq = cut_coulsq;

  k_params.h_view(j,i) = k_params.h_view(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.h_view(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
    m_cut_ljsq[j][i] = m_cut_ljsq[i][j] = cut_ljsq;
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = cut_coulsq;
  }

  k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutone*cutone;
  k_cutsq.template modify<LMPHostType>();
  k_params.template modify<LMPHostType>();

  return cutone;
}

namespace LAMMPS_NS {
template class PairLJCharmmfswCoulLongKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJCharmmfswCoulLongKokkos<LMPHostType>;
#endif
}




/*
 80d105
 <     memory->destroy(offset);
 598c650
 < void PairLJCharmmCoulLong::allocate()
 ---
 > void PairLJCharmmfswCoulLong::allocate()
 622d673
 <   memory->create(offset,n+1,n+1,"pair:offset");
 631c682
 < void PairLJCharmmCoulLong::settings(int narg, char **arg)
 ---
 > void PairLJCharmmfswCoulLong::settings(int narg, char **arg)
 645c696
 < void PairLJCharmmCoulLong::coeff(int narg, char **arg)
 ---
 > void PairLJCharmmfswCoulLong::coeff(int narg, char **arg)
 752c801
 < double PairLJCharmmCoulLong::init_one(int i, int j)
 ---
 > double PairLJCharmmfswCoulLong::init_one(int i, int j)
 790c839
 < void PairLJCharmmCoulLong::write_restart(FILE *fp)
 ---
 > void PairLJCharmmfswCoulLong::write_restart(FILE *fp)
 811c860
 < void PairLJCharmmCoulLong::read_restart(FILE *fp)
 ---
 > void PairLJCharmmfswCoulLong::read_restart(FILE *fp)
 842c891
 < void PairLJCharmmCoulLong::write_restart_settings(FILE *fp)
 ---
 > void PairLJCharmmfswCoulLong::write_restart_settings(FILE *fp)
 857c906
 < void PairLJCharmmCoulLong::read_restart_settings(FILE *fp)
 ---
 > void PairLJCharmmfswCoulLong::read_restart_settings(FILE *fp)
 882c931
 < void PairLJCharmmCoulLong::write_data(FILE *fp)
 ---
 > void PairLJCharmmfswCoulLong::write_data(FILE *fp)
 893c942
 < void PairLJCharmmCoulLong::write_data_all(FILE *fp)
 ---
 > void PairLJCharmmfswCoulLong::write_data_all(FILE *fp)
 903c952
 < double PairLJCharmmCoulLong::single(int i, int j, int itype, int jtype,
 ---
 > double PairLJCharmmfswCoulLong::single(int i, int j, int itype, int jtype,
 908,909c957,958
 <   double r2inv,r6inv,r,grij,expm2,t,erfc,prefactor;
 <   double switch1,switch2,fraction,table,forcecoul,forcelj,phicoul,philj;
 ---
 >   double r,rinv,r2inv,r3inv,r6inv,grij,expm2,t,erfc,prefactor;
 >   double switch1,fraction,table,forcecoul,forcelj,phicoul,philj,philj12,philj6;
 911a961,962
 >   r = sqrt(rsq);
 >   rinv = 1.0/r;
 939c990,991
 <     r6inv = r2inv*r2inv*r2inv;
 ---
 >     r3inv = rinv*rinv*rinv;
 >     r6inv = r3inv*r3inv;
 943,947c995,996
 <         (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) * denom_lj_inv;
 <       switch2 = 12.0*rsq * (cut_ljsq-rsq) *
 <         (rsq-cut_lj_innersq) * denom_lj_inv;
 <       philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
 <       forcelj = forcelj*switch1 + philj*switch2;
 ---
 >               (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
 >             forcelj = forcelj*switch1;
 965d1013
 <     philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);
 967,969c1015,1025
 <       switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
 <         (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) * denom_lj_inv;
 <       philj *= switch1;
 ---
 >       philj12 = lj3[itype][jtype]*cut_lj6*denom_lj12 *
 >         (r6inv - cut_lj6inv)*(r6inv - cut_lj6inv);
 >       philj6 = -lj4[itype][jtype]*cut_lj3*denom_lj6 *
 >         (r3inv - cut_lj3inv)*(r3inv - cut_lj3inv);
 >       philj = philj12 + philj6;
 >     } else {
 >       philj12 = r6inv*lj3[itype][jtype]*r6inv -
 >         lj3[itype][jtype]*cut_lj_inner6inv*cut_lj6inv;
 >       philj6 = -lj4[itype][jtype]*r6inv +
 >         lj4[itype][jtype]*cut_lj_inner3inv*cut_lj3inv;
 >       philj = philj12 + philj6;
 979c1035
 < void *PairLJCharmmCoulLong::extract(const char *str, int &dim)
 ---
 > void *PairLJCharmmfswCoulLong::extract(const char *str, int &dim)
 988a1045,1047
 >
 >   // info extracted by dihedral_charmmfsw
 >
 989a1049,1051
 >   if (strcmp(str,"cut_lj_inner") == 0) return (void *) &cut_lj_inner;
 >   if (strcmp(str,"cut_lj") == 0) return (void *) &cut_lj;
 >   if (strcmp(str,"dihedflag") == 0) return (void *) &dihedflag;

 
 */

// nothing to do for all these, inherited from PairLJCharmmfswCoulLong




/*
 
 226c257
 < void PairLJCharmmCoulLong::compute_inner()
 ---
 > void PairLJCharmmfswCoulLong::compute_inner()
 248a280,286
 >   double cut_out_on = cut_respa[0];
 >   double cut_out_off = cut_respa[1];
 >
 >   double cut_out_diff = cut_out_off - cut_out_on;
 >   double cut_out_on_sq = cut_out_on*cut_out_on;
 >   double cut_out_off_sq = cut_out_off*cut_out_off;
 >
 284c322
 <           rsw = (sqrt(rsq) - cut_out_on)*cut_out_diff_inv;
 ---
 >           rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
 303c341
 < void PairLJCharmmCoulLong::compute_middle()
 ---
 > void PairLJCharmmfswCoulLong::compute_middle()
 308c346
 <   double philj,switch1,switch2;
 ---
 >   double switch1;
 326a365,376
 >   double cut_in_off = cut_respa[0];
 >   double cut_in_on = cut_respa[1];
 >   double cut_out_on = cut_respa[2];
 >   double cut_out_off = cut_respa[3];
 >
 >   double cut_in_diff = cut_in_on - cut_in_off;
 >   double cut_out_diff = cut_out_off - cut_out_on;
 >   double cut_in_off_sq = cut_in_off*cut_in_off;
 >   double cut_in_on_sq = cut_in_on*cut_in_on;
 >   double cut_out_on_sq = cut_out_on*cut_out_on;
 >   double cut_out_off_sq = cut_out_off*cut_out_off;
 >
 361,365c411,412
 <             (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) * denom_lj_inv;
 <           switch2 = 12.0*rsq * (cut_ljsq-rsq) *
 <             (rsq-cut_lj_innersq) * denom_lj_inv;
 <           philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
 <           forcelj = forcelj*switch1 + philj*switch2;
 ---
 >               (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
 >           forcelj = forcelj*switch1;
 370c417
 <           rsw = (sqrt(rsq) - cut_in_off)*cut_in_diff_inv;
 ---
 >           rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff;
 374c421
 <           rsw = (sqrt(rsq) - cut_out_on)*cut_out_diff_inv;
 ---
 >           rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
 393c440
 < void PairLJCharmmCoulLong::compute_outer(int eflag, int vflag)
 ---
 > void PairLJCharmmfswCoulLong::compute_outer(int eflag, int vflag)
 396c443
 <   double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
 ---
 >   double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,evdwl6,evdwl12,ecoul,fpair;
 398c445
 <   double r,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
 ---
 >   double r,rinv,r2inv,r3inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
 400c447
 <   double philj,switch1,switch2;
 ---
 >   double switch1;
 422a470,476
 >   double cut_in_off = cut_respa[2];
 >   double cut_in_on = cut_respa[3];
 >
 >   double cut_in_diff = cut_in_on - cut_in_off;
 >   double cut_in_off_sq = cut_in_off*cut_in_off;
 >   double cut_in_on_sq = cut_in_on*cut_in_on;
 >
 448a503
 >         r6inv = r2inv*r2inv*r2inv;
 489d543
 <           r6inv = r2inv*r2inv*r2inv;
 493,497c547,548
 <               (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) * denom_lj_inv;
 <             switch2 = 12.0*rsq * (cut_ljsq-rsq) *
 <               (rsq-cut_lj_innersq) * denom_lj_inv;
 <             philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
 <             forcelj = forcelj*switch1 + philj*switch2;
 ---
 >               (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
 >             forcelj = forcelj*switch1;
 533d583
 <             r6inv = r2inv*r2inv*r2inv;
 536,538c586,598
 <               switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
 <                 (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) * denom_lj_inv;
 <               evdwl *= switch1;
 ---
 >               rinv = sqrt(r2inv);
 >               r3inv = r2inv*rinv;
 >               evdwl12 = lj3[itype][jtype]*cut_lj6*denom_lj12 *
 >                 (r6inv - cut_lj6inv)*(r6inv - cut_lj6inv);
 >               evdwl6 = -lj4[itype][jtype]*cut_lj3*denom_lj6 *
 >                 (r3inv - cut_lj3inv)*(r3inv - cut_lj3inv);
 >               evdwl = evdwl12 + evdwl6;
 >             } else {
 >               evdwl12 = r6inv*lj3[itype][jtype]*r6inv -
 >                 lj3[itype][jtype]*cut_lj_inner6inv*cut_lj6inv;
 >               evdwl6 = -lj4[itype][jtype]*r6inv +
 >                 lj4[itype][jtype]*cut_lj_inner3inv*cut_lj3inv;
 >               evdwl = evdwl12 + evdwl6;
 561d620
 <             r6inv = r2inv*r2inv*r2inv;
 565,569c624,625
 <                 (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) * denom_lj_inv;
 <               switch2 = 12.0*rsq * (cut_ljsq-rsq) *
 <                 (rsq-cut_lj_innersq) * denom_lj_inv;
 <               philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
 <               forcelj = forcelj*switch1 + philj*switch2;
 ---
 >                (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
 >               forcelj = forcelj*switch1;
 572d627
 <             r6inv = r2inv*r2inv*r2inv;
 576,580c631,632
 <                 (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) * denom_lj_inv;
 <               switch2 = 12.0*rsq * (cut_ljsq-rsq) *
 <                 (rsq-cut_lj_innersq) * denom_lj_inv;
 <               philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
 <               forcelj = forcelj*switch1 + philj*switch2;
 ---
 >                 (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
 >               forcelj = forcelj*switch1;

 */

// kokkos doesnt support respa, so ignore compute_inner / compute_middle / compute_outer
