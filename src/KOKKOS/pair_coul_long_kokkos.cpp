/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Ray Shan (SNL)
------------------------------------------------------------------------- */

#include "pair_coul_long_kokkos.h"
#include <cmath>
#include <cstring>
#include "kokkos.h"
#include "atom_kokkos.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "respa.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

#define KOKKOS_CUDA_MAX_THREADS 256
#define KOKKOS_CUDA_MIN_BLOCKS 8


#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
PairCoulLongKokkos<Space>::PairCoulLongKokkos(LAMMPS *lmp):PairCoulLong(lmp)
{
  respa_enable = 0;

  atomKK = (AtomKokkos *) atom;
  execution_space = Space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
  cutsq = NULL;
  cut_coulsq = 0.0;

}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
PairCoulLongKokkos<Space>::~PairCoulLongKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    k_cutsq = DAT::tdual_float_2d();
    k_cut_coulsq = DAT::tdual_float_2d();
    memory->sfree(cutsq);
    eatom = NULL;
    vatom = NULL;
    cutsq = NULL;
  }
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairCoulLongKokkos<Space>::cleanup_copy() {
  allocated = 0;
  cutsq = NULL;
  eatom = NULL;
  vatom = NULL;
  ftable = NULL;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairCoulLongKokkos<Space>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = DualViewHelper<Space>::view(k_eatom);
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = DualViewHelper<Space>::view(k_vatom);
  }

  atomKK->sync(execution_space,datamask_read);
  DualViewHelper<Space>::sync(k_cutsq);
  DualViewHelper<Space>::sync(k_cut_coulsq);
  DualViewHelper<Space>::sync(k_params);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = DualViewHelper<Space>::view(atomKK->k_x);
  c_x = DualViewHelper<Space>::view(atomKK->k_x);
  f = DualViewHelper<Space>::view(atomKK->k_f);
  q = DualViewHelper<Space>::view(atomKK->k_q);
  type = DualViewHelper<Space>::view(atomKK->k_type);
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  qqrd2e = force->qqrd2e;
  newton_pair = force->newton_pair;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];
  special_coul[0] = force->special_coul[0];
  special_coul[1] = force->special_coul[1];
  special_coul[2] = force->special_coul[2];
  special_coul[3] = force->special_coul[3];

  // loop over neighbors of my atoms

  copymode = 1;

  EV_FLOAT ev;
  if(ncoultablebits)
    ev = pair_compute<Space,PairCoulLongKokkos<Space>,CoulLongTable<1> >
      (this,(NeighListKokkos<Space>*)list);
  else
    ev = pair_compute<Space,PairCoulLongKokkos<Space>,CoulLongTable<0> >
      (this,(NeighListKokkos<Space>*)list);


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
    DualViewHelper<Space>::modify(k_eatom);
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    DualViewHelper<Space>::modify(k_vatom);
    k_vatom.sync_host();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute<Space>(this);

  copymode = 0;
}

/* ----------------------------------------------------------------------
   compute coulomb pair force between atoms i and j
   ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
template<bool STACKPARAMS,  class Specialisation>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairCoulLongKokkos<Space>::
compute_fcoul(const KK_FLOAT& rsq, const int& i, const int&j,
              const int& itype, const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
  if(Specialisation::DoTable && rsq > tabinnersq) {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    const int itable = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
    const KK_FLOAT fraction = (rsq_lookup.f - d_rtable[itable]) * d_drtable[itable];
    const KK_FLOAT table = d_ftable[itable] + fraction*d_dftable[itable];
    KK_FLOAT forcecoul = qtmp*q[j] * table;
    if (factor_coul < 1.0) {
      const KK_FLOAT table = d_ctable[itable] + fraction*d_dctable[itable];
      const KK_FLOAT prefactor = qtmp*q[j] * table;
      forcecoul -= (1.0-factor_coul)*prefactor;
    }
    return forcecoul/rsq;
  } else {
    const KK_FLOAT r = sqrt(rsq);
    const KK_FLOAT grij = g_ewald * r;
    const KK_FLOAT expm2 = exp(-grij*grij);
    const KK_FLOAT t = 1.0 / (1.0 + EWALD_P*grij);
    const KK_FLOAT rinv = 1.0/r;
    const KK_FLOAT erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
    const KK_FLOAT prefactor = qqrd2e * qtmp*q[j]*rinv;
    KK_FLOAT forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
    if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;

    return forcecoul*rinv*rinv;
  }
}

/* ----------------------------------------------------------------------
   compute coulomb pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairCoulLongKokkos<Space>::
compute_ecoul(const KK_FLOAT& rsq, const int& i, const int&j,
              const int& itype, const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
  if(Specialisation::DoTable && rsq > tabinnersq) {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    const int itable = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
    const KK_FLOAT fraction = (rsq_lookup.f - d_rtable[itable]) * d_drtable[itable];
    const KK_FLOAT table = d_etable[itable] + fraction*d_detable[itable];
    KK_FLOAT ecoul = qtmp*q[j] * table;
    if (factor_coul < 1.0) {
      const KK_FLOAT table = d_ctable[itable] + fraction*d_dctable[itable];
      const KK_FLOAT prefactor = qtmp*q[j] * table;
      ecoul -= (1.0-factor_coul)*prefactor;
    }
    return ecoul;
  } else {
    const KK_FLOAT r = sqrt(rsq);
    const KK_FLOAT grij = g_ewald * r;
    const KK_FLOAT expm2 = exp(-grij*grij);
    const KK_FLOAT t = 1.0 / (1.0 + EWALD_P*grij);
    const KK_FLOAT erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
    const KK_FLOAT prefactor = qqrd2e * qtmp*q[j]/r;
    KK_FLOAT ecoul = prefactor * erfc;
    if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
    return ecoul;
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairCoulLongKokkos<Space>::allocate()
{
  PairCoulLong::allocate();

  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = DualViewHelper<Space>::view(k_cutsq);

  memoryKK->create_kokkos(k_cut_coulsq,n+1,n+1,"pair:cut_coulsq");
  d_cut_coulsq = DualViewHelper<Space>::view(k_cut_coulsq);

  k_params = Kokkos::DualView<params_coul**,Kokkos::LayoutRight,DeviceType>("PairCoulLong::params",n+1,n+1);
  params = DualViewHelper<Space>::view(k_params);
}

template<ExecutionSpace Space>
void PairCoulLongKokkos<Space>::init_tables(double cut_coul, double *cut_respa)
{
  Pair::init_tables(cut_coul,cut_respa);

  typedef typename AT::t_float_1d table_type;
  typedef HAT::t_float_1d host_table_type;

  int ntable = 1;
  for (int i = 0; i < ncoultablebits; i++) ntable *= 2;


  // Copy rtable and drtable
  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);
  for(int i = 0; i < ntable; i++) {
    h_table(i) = rtable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_rtable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);
  for(int i = 0; i < ntable; i++) {
    h_table(i) = drtable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_drtable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  // Copy ftable and dftable
  for(int i = 0; i < ntable; i++) {
    h_table(i) = ftable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_ftable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  for(int i = 0; i < ntable; i++) {
    h_table(i) = dftable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_dftable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  // Copy ctable and dctable
  for(int i = 0; i < ntable; i++) {
    h_table(i) = ctable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_ctable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  for(int i = 0; i < ntable; i++) {
    h_table(i) = dctable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_dctable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  // Copy etable and detable
  for(int i = 0; i < ntable; i++) {
    h_table(i) = etable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_etable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  for(int i = 0; i < ntable; i++) {
    h_table(i) = detable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_detable = d_table;
  }
}


/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairCoulLongKokkos<Space>::settings(int narg, char **arg)
{
  if (narg > 2) error->all(FLERR,"Illegal pair_style command");

  PairCoulLong::settings(narg,arg);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairCoulLongKokkos<Space>::init_style()
{
  PairCoulLong::init_style();

  // error if rRESPA with inner levels

  if (update->whichflag == 1 && strstr(update->integrate_style,"respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa)
      error->all(FLERR,"Cannot use Kokkos pair style with rRESPA inner/middle");
  }

  // irequest = neigh request made by parent class

  neighflag = lmp->kokkos->neighflag;
  int irequest = neighbor->nrequest - 1;

  neighbor->requests[irequest]->
    kokkos_host = (Space == Host) &&
    !(Space == Device);
  neighbor->requests[irequest]->
    kokkos_device = (Space == Device);

  if (neighflag == FULL) {
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
  } else if (neighflag == HALF || neighflag == HALFTHREAD) {
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 1;
  } else {
    error->all(FLERR,"Cannot use chosen neighbor list style with buck/coul/long/kk");
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
double PairCoulLongKokkos<Space>::init_one(int i, int j)
{
  KK_FLOAT cutone = PairCoulLong::init_one(i,j);
  KK_FLOAT cut_coulsqm = cut_coulsq;

  k_params.h_view(i,j).cut_coulsq = cut_coulsqm;

  k_params.h_view(j,i) = k_params.h_view(i,j);
  if(i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.h_view(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = cut_coulsqm;
  }

  k_cutsq.h_view(i,j) = cutone*cutone;
  k_cutsq.modify_host();
  k_cut_coulsq.h_view(i,j) = cut_coulsqm;
  k_cut_coulsq.modify_host();
  k_params.modify_host();

  return cutone;

}


namespace LAMMPS_NS {
template class PairCoulLongKokkos<Device>;
template class PairCoulLongKokkos<Host>;
}

