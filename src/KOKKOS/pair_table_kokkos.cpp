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
   Contributing author: Christian Trott (SNL)
------------------------------------------------------------------------- */

#include "pair_table_kokkos.h"
#include <cstring>
#include "kokkos.h"
#include "atom.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
PairTableKokkos<Space>::PairTableKokkos(LAMMPS *lmp) : PairTable(lmp)
{
  update_table = 0;
  atomKK = (AtomKokkos *) atom;
  execution_space = Space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
  k_table = new TableDual();
  h_table = new TableHost(k_table);
  d_table = new TableDevice(k_table);
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
PairTableKokkos<Space>::~PairTableKokkos()
{
  if (copymode) return;

  delete k_table;
  k_table = nullptr;
  delete h_table;
  h_table = nullptr;
  delete d_table;
  d_table = nullptr;
  copymode = true; //prevents base class destructor from running
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairTableKokkos<Space>::compute(int eflag_in, int vflag_in)
{
  if(update_table)
    create_kokkos_tables();
  if(tabstyle == LOOKUP)
    compute_style<LOOKUP>(eflag_in,vflag_in);
  if(tabstyle == LINEAR)
    compute_style<LINEAR>(eflag_in,vflag_in);
  if(tabstyle == SPLINE)
    compute_style<SPLINE>(eflag_in,vflag_in);
  if(tabstyle == BITMAP)
    compute_style<BITMAP>(eflag_in,vflag_in);
}

template<ExecutionSpace Space>
template<int TABSTYLE>
void PairTableKokkos<Space>::compute_style(int eflag_in, int vflag_in)
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
  //DualViewHelper<Space>::sync(k_cutsq);
  //DualViewHelper<Space>::sync(k_params);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = c_x = DualViewHelper<Space>::view(atomKK->k_x);
  f = DualViewHelper<Space>::view(atomKK->k_f);
  type = DualViewHelper<Space>::view(atomKK->k_type);
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];
  newton_pair = force->newton_pair;
  d_cutsq = d_table->cutsq;
  // loop over neighbors of my atoms

  EV_FLOAT ev;
  if(atom->ntypes > MAX_TYPES_STACKPARAMS) {
    if (neighflag == FULL) {
      PairComputeFunctor<Space,PairTableKokkos<Space>,FULL,false,S_TableCompute<Space,TABSTYLE> >
        ff(this,(NeighListKokkos<Space>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,ff,ev);
      else Kokkos::parallel_for(list->inum,ff);
      ff.contribute();
    } else if (neighflag == HALFTHREAD) {
      PairComputeFunctor<Space,PairTableKokkos<Space>,HALFTHREAD,false,S_TableCompute<Space,TABSTYLE> >
        ff(this,(NeighListKokkos<Space>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,ff,ev);
      else Kokkos::parallel_for(list->inum,ff);
      ff.contribute();
    } else if (neighflag == HALF) {
      PairComputeFunctor<Space,PairTableKokkos<Space>,HALF,false,S_TableCompute<Space,TABSTYLE> >
        f(this,(NeighListKokkos<Space>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,f,ev);
      else Kokkos::parallel_for(list->inum,f);
      f.contribute();
    }
  } else {
    if (neighflag == FULL) {
      PairComputeFunctor<Space,PairTableKokkos<Space>,FULL,true,S_TableCompute<Space,TABSTYLE> >
        f(this,(NeighListKokkos<Space>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,f,ev);
      else Kokkos::parallel_for(list->inum,f);
      f.contribute();
    } else if (neighflag == HALFTHREAD) {
      PairComputeFunctor<Space,PairTableKokkos<Space>,HALFTHREAD,true,S_TableCompute<Space,TABSTYLE> >
        f(this,(NeighListKokkos<Space>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,f,ev);
      else Kokkos::parallel_for(list->inum,f);
      f.contribute();
    } else if (neighflag == HALF) {
      PairComputeFunctor<Space,PairTableKokkos<Space>,HALF,true,S_TableCompute<Space,TABSTYLE> >
        f(this,(NeighListKokkos<Space>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,f,ev);
      else Kokkos::parallel_for(list->inum,f);
      f.contribute();
    }
  }

  if (eflag) eng_vdwl += ev.evdwl;
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

}

template<ExecutionSpace Space>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairTableKokkos<Space>::
compute_fpair(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
  (void) i;
  (void) j;
  union_int_float_t rsq_lookup;
  KK_FLOAT fpair;
  const int tidx = d_table_const.tabindex(itype,jtype);
  if (Specialisation::TabStyle == LOOKUP) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    fpair = d_table_const.f(tidx,itable);
  } else if (Specialisation::TabStyle == LINEAR) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    const KK_FLOAT fraction = (rsq - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx);
    fpair = d_table_const.f(tidx,itable) + fraction*d_table_const.df(tidx,itable);
  } else if (Specialisation::TabStyle == SPLINE) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    const KK_FLOAT b = (rsq - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx);
    const KK_FLOAT a = 1.0 - b;
    fpair = a * d_table_const.f(tidx,itable) + b * d_table_const.f(tidx,itable+1) +
      ((a*a*a-a)*d_table_const.f2(tidx,itable) + (b*b*b-b)*d_table_const.f2(tidx,itable+1)) *
      d_table_const.deltasq6(tidx);
  } else {
    rsq_lookup.f = rsq;
    int itable = rsq_lookup.i & d_table_const.nmask(tidx);
    itable >>= d_table_const.nshiftbits(tidx);
    const KK_FLOAT fraction = (rsq_lookup.f - d_table_const.rsq(tidx,itable)) * d_table_const.drsq(tidx,itable);
    fpair = d_table_const.f(tidx,itable) + fraction*d_table_const.df(tidx,itable);
  }
  return fpair;
}

template<ExecutionSpace Space>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairTableKokkos<Space>::
compute_evdwl(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
  (void) i;
  (void) j;
  KK_FLOAT evdwl;
  union_int_float_t rsq_lookup;
  const int tidx = d_table_const.tabindex(itype,jtype);
  if (Specialisation::TabStyle == LOOKUP) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    evdwl = d_table_const.e(tidx,itable);
  } else if (Specialisation::TabStyle == LINEAR) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    const KK_FLOAT fraction = (rsq - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx);
    evdwl = d_table_const.e(tidx,itable) + fraction*d_table_const.de(tidx,itable);
  } else if (Specialisation::TabStyle == SPLINE) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    const KK_FLOAT b = (rsq - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx);
    const KK_FLOAT a = 1.0 - b;
    evdwl = a * d_table_const.e(tidx,itable) + b * d_table_const.e(tidx,itable+1) +
        ((a*a*a-a)*d_table_const.e2(tidx,itable) + (b*b*b-b)*d_table_const.e2(tidx,itable+1)) *
        d_table_const.deltasq6(tidx);
  } else {
    rsq_lookup.f = rsq;
    int itable = rsq_lookup.i & d_table_const.nmask(tidx);
    itable >>= d_table_const.nshiftbits(tidx);
    const KK_FLOAT fraction = (rsq_lookup.f - d_table_const.rsq(tidx,itable)) * d_table_const.drsq(tidx,itable);
    evdwl = d_table_const.e(tidx,itable) + fraction*d_table_const.de(tidx,itable);
  }
  return evdwl;
}

template<ExecutionSpace Space>
void PairTableKokkos<Space>::create_kokkos_tables()
{
  const int tlm1 = tablength-1;

  memoryKK->create_kokkos(k_table->k_nshiftbits,ntables,"Table::nshiftbits");
  memoryKK->create_kokkos(k_table->k_nmask,ntables,"Table::nmask");
  memoryKK->create_kokkos(k_table->k_innersq,ntables,"Table::innersq");
  memoryKK->create_kokkos(k_table->k_invdelta,ntables,"Table::invdelta");
  memoryKK->create_kokkos(k_table->k_deltasq6,ntables,"Table::deltasq6");

  if(tabstyle == LOOKUP) {
    memoryKK->create_kokkos(k_table->k_e,ntables,tlm1,"Table::e");
    memoryKK->create_kokkos(k_table->k_f,ntables,tlm1,"Table::f");
  }

  if(tabstyle == LINEAR) {
    memoryKK->create_kokkos(k_table->k_rsq,ntables,tablength,"Table::rsq");
    memoryKK->create_kokkos(k_table->k_e,ntables,tablength,"Table::e");
    memoryKK->create_kokkos(k_table->k_f,ntables,tablength,"Table::f");
    memoryKK->create_kokkos(k_table->k_de,ntables,tlm1,"Table::de");
    memoryKK->create_kokkos(k_table->k_df,ntables,tlm1,"Table::df");
  }

  if(tabstyle == SPLINE) {
    memoryKK->create_kokkos(k_table->k_rsq,ntables,tablength,"Table::rsq");
    memoryKK->create_kokkos(k_table->k_e,ntables,tablength,"Table::e");
    memoryKK->create_kokkos(k_table->k_f,ntables,tablength,"Table::f");
    memoryKK->create_kokkos(k_table->k_e2,ntables,tablength,"Table::e2");
    memoryKK->create_kokkos(k_table->k_f2,ntables,tablength,"Table::f2");
  }

  if(tabstyle == BITMAP) {
    int ntable = 1 << tablength;
    memoryKK->create_kokkos(k_table->k_rsq,ntables,ntable,"Table::rsq");
    memoryKK->create_kokkos(k_table->k_e,ntables,ntable,"Table::e");
    memoryKK->create_kokkos(k_table->k_f,ntables,ntable,"Table::f");
    memoryKK->create_kokkos(k_table->k_de,ntables,ntable,"Table::de");
    memoryKK->create_kokkos(k_table->k_df,ntables,ntable,"Table::df");
    memoryKK->create_kokkos(k_table->k_drsq,ntables,ntable,"Table::drsq");
  }

  for(int i=0; i < ntables; i++) {
    Table* tb = &tables[i];

    h_table->nshiftbits[i] = tb->nshiftbits;
    h_table->nmask[i] = tb->nmask;
    h_table->innersq[i] = tb->innersq;
    h_table->invdelta[i] = tb->invdelta;
    h_table->deltasq6[i] = tb->deltasq6;

    for(int j = 0; j<h_table->rsq.extent(1); j++)
      h_table->rsq(i,j) = tb->rsq[j];
    for(int j = 0; j<h_table->drsq.extent(1); j++)
      h_table->drsq(i,j) = tb->drsq[j];
    for(int j = 0; j<h_table->e.extent(1); j++)
      h_table->e(i,j) = tb->e[j];
    for(int j = 0; j<h_table->de.extent(1); j++)
      h_table->de(i,j) = tb->de[j];
    for(int j = 0; j<h_table->f.extent(1); j++)
      h_table->f(i,j) = tb->f[j];
    for(int j = 0; j<h_table->df.extent(1); j++)
      h_table->df(i,j) = tb->df[j];
    for(int j = 0; j<h_table->e2.extent(1); j++)
      h_table->e2(i,j) = tb->e2[j];
    for(int j = 0; j<h_table->f2.extent(1); j++)
      h_table->f2(i,j) = tb->f2[j];
  }

  k_table->k_nshiftbits.modify_host();
  k_table->k_nshiftbits.sync_device();
  d_table_const.nshiftbits = d_table->nshiftbits;
  k_table->k_nmask.modify_host();
  k_table->k_nmask.sync_device();
  d_table_const.nmask = d_table->nmask;
  k_table->k_innersq.modify_host();
  k_table->k_innersq.sync_device();
  d_table_const.innersq = d_table->innersq;
  k_table->k_invdelta.modify_host();
  k_table->k_invdelta.sync_device();
  d_table_const.invdelta = d_table->invdelta;
  k_table->k_deltasq6.modify_host();
  k_table->k_deltasq6.sync_device();
  d_table_const.deltasq6 = d_table->deltasq6;

  if(tabstyle == LOOKUP) {
    k_table->k_e.modify_host();
    k_table->k_e.sync_device();
    d_table_const.e = d_table->e;
    k_table->k_f.modify_host();
    k_table->k_f.sync_device();
    d_table_const.f = d_table->f;
  }

  if(tabstyle == LINEAR) {
    k_table->k_rsq.modify_host();
    k_table->k_rsq.sync_device();
    d_table_const.rsq = d_table->rsq;
    k_table->k_e.modify_host();
    k_table->k_e.sync_device();
    d_table_const.e = d_table->e;
    k_table->k_f.modify_host();
    k_table->k_f.sync_device();
    d_table_const.f = d_table->f;
    k_table->k_de.modify_host();
    k_table->k_de.sync_device();
    d_table_const.de = d_table->de;
    k_table->k_df.modify_host();
    k_table->k_df.sync_device();
    d_table_const.df = d_table->df;
  }

  if(tabstyle == SPLINE) {
    k_table->k_rsq.modify_host();
    k_table->k_rsq.sync_device();
    d_table_const.rsq = d_table->rsq;
    k_table->k_e.modify_host();
    k_table->k_e.sync_device();
    d_table_const.e = d_table->e;
    k_table->k_f.modify_host();
    k_table->k_f.sync_device();
    d_table_const.f = d_table->f;
    k_table->k_e2.modify_host();
    k_table->k_e2.sync_device();
    d_table_const.e2 = d_table->e2;
    k_table->k_f2.modify_host();
    k_table->k_f2.sync_device();
    d_table_const.f2 = d_table->f2;
  }

  if(tabstyle == BITMAP) {
    k_table->k_rsq.modify_host();
    k_table->k_rsq.sync_device();
    d_table_const.rsq = d_table->rsq;
    k_table->k_e.modify_host();
    k_table->k_e.sync_device();
    d_table_const.e = d_table->e;
    k_table->k_f.modify_host();
    k_table->k_f.sync_device();
    d_table_const.f = d_table->f;
    k_table->k_de.modify_host();
    k_table->k_de.sync_device();
    d_table_const.de = d_table->de;
    k_table->k_df.modify_host();
    k_table->k_df.sync_device();
    d_table_const.df = d_table->df;
    k_table->k_drsq.modify_host();
    k_table->k_drsq.sync_device();
    d_table_const.drsq = d_table->drsq;
  }

  k_table->k_cutsq.modify_host();
  k_table->k_cutsq.sync_device();
  d_table_const.cutsq = d_table->cutsq;
  k_table->k_tabindex.modify_host();
  k_table->k_tabindex.sync_device();
  d_table_const.tabindex = d_table->tabindex;

  update_table = 0;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairTableKokkos<Space>::allocate()
{
  allocated = 1;
  const int nt = atom->ntypes + 1;

  memory->create(setflag,nt,nt,"pair:setflag");
  memoryKK->create_kokkos(k_table->k_cutsq,cutsq,nt,nt,"pair:cutsq");
  memoryKK->create_kokkos(k_table->k_tabindex,tabindex,nt,nt,"pair:tabindex");
  d_table_const.cutsq = d_table->cutsq;
  d_table_const.tabindex = d_table->tabindex;

  memset(&setflag[0][0],0,nt*nt*sizeof(int));
  memset(&cutsq[0][0],0,nt*nt*sizeof(KK_FLOAT));
  memset(&tabindex[0][0],0,nt*nt*sizeof(int));
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairTableKokkos<Space>::settings(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal pair_style command");

  // new settings

  if (strcmp(arg[0],"lookup") == 0) tabstyle = LOOKUP;
  else if (strcmp(arg[0],"linear") == 0) tabstyle = LINEAR;
  else if (strcmp(arg[0],"spline") == 0) tabstyle = SPLINE;
  else if (strcmp(arg[0],"bitmap") == 0) tabstyle = BITMAP;
  else error->all(FLERR,"Unknown table style in pair_style command");

  tablength = force->inumeric(FLERR,arg[1]);
  if (tablength < 2) error->all(FLERR,"Illegal number of pair table entries");

  // optional keywords
  // assert the tabulation is compatible with a specific long-range solver

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"ewald") == 0) ewaldflag = 1;
    else if (strcmp(arg[iarg],"pppm") == 0) pppmflag = 1;
    else if (strcmp(arg[iarg],"msm") == 0) msmflag = 1;
    else if (strcmp(arg[iarg],"dispersion") == 0) dispersionflag = 1;
    else if (strcmp(arg[iarg],"tip4p") == 0) tip4pflag = 1;
    else error->all(FLERR,"Illegal pair_style command");
    iarg++;
  }

  // delete old tables, since cannot just change settings

  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);

  if (allocated) {
    memory->destroy(setflag);

    d_table_const.tabindex = d_table->tabindex = typename AT::t_int_2d();
    h_table->tabindex = HAT::t_int_2d();

    d_table_const.cutsq = d_table->cutsq = typename AT::t_float_2d();
    h_table->cutsq = HAT::t_float_2d();
  }
  allocated = 0;

  ntables = 0;
  tables = NULL;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
double PairTableKokkos<Space>::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  tabindex[j][i] = tabindex[i][j];

  if(i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_cutsq[j][i] = m_cutsq[i][j] = tables[tabindex[i][j]].cut*tables[tabindex[i][j]].cut;
  }

  return tables[tabindex[i][j]].cut;
}

/* ----------------------------------------------------------------------
   compute r,e,f vectors from splined values
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairTableKokkos<Space>::compute_table(Table *tb)
{
  update_table = 1;
  PairTable::compute_table(tb);
}

template<ExecutionSpace Space>
void PairTableKokkos<Space>::init_style()
{
  neighbor->request(this,instance_me);
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
    error->all(FLERR,"Cannot use chosen neighbor list style with lj/cut/kk");
  }
}

template<ExecutionSpace Space>
void PairTableKokkos<Space>::cleanup_copy() {
  // WHY needed: this prevents parent copy from deallocating any arrays
  allocated = 0;
  cutsq = NULL;
  eatom = NULL;
  vatom = NULL;
  h_table=NULL; d_table=NULL;
}

namespace LAMMPS_NS {
template class PairTableKokkos<Device>;
template class PairTableKokkos<Host>;

}

