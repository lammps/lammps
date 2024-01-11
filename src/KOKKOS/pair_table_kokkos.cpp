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
   Contributing author: Christian Trott (SNL)
------------------------------------------------------------------------- */

#include "pair_table_kokkos.h"

#include "atom.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairTableKokkos<DeviceType>::PairTableKokkos(LAMMPS *lmp) : PairTable(lmp)
{
  update_table = 0;
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
  h_table = new TableHost();
  d_table = new TableDevice();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairTableKokkos<DeviceType>::~PairTableKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memory->destroy(setflag);
    memoryKK->destroy_kokkos(d_table->cutsq,cutsq);
    memoryKK->destroy_kokkos(d_table->tabindex,tabindex);
  }

  delete h_table;
  h_table = nullptr;
  delete d_table;
  d_table = nullptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairTableKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  if (update_table)
    create_kokkos_tables();
  if (tabstyle == LOOKUP)
    compute_style<LOOKUP>(eflag_in,vflag_in);
  if (tabstyle == LINEAR)
    compute_style<LINEAR>(eflag_in,vflag_in);
  if (tabstyle == SPLINE)
    compute_style<SPLINE>(eflag_in,vflag_in);
  if (tabstyle == BITMAP)
    compute_style<BITMAP>(eflag_in,vflag_in);
}

template<class DeviceType>
template<int TABSTYLE>
void PairTableKokkos<DeviceType>::compute_style(int eflag_in, int vflag_in)
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
  //k_cutsq.template sync<DeviceType>();
  //k_params.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];
  newton_pair = force->newton_pair;
  d_cutsq = d_table->cutsq;
  // loop over neighbors of my atoms

  copymode = 1;

  EV_FLOAT ev;
  if (atom->ntypes > MAX_TYPES_STACKPARAMS) {
    if (neighflag == FULL) {
      PairComputeFunctor<PairTableKokkos<DeviceType>,FULL,false,0,S_TableCompute<DeviceType,TABSTYLE> >
        ff(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,ff,ev);
      else Kokkos::parallel_for(list->inum,ff);
      ff.contribute();
    } else if (neighflag == HALFTHREAD) {
      PairComputeFunctor<PairTableKokkos<DeviceType>,HALFTHREAD,false,0,S_TableCompute<DeviceType,TABSTYLE> >
        ff(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,ff,ev);
      else Kokkos::parallel_for(list->inum,ff);
      ff.contribute();
    } else if (neighflag == HALF) {
      PairComputeFunctor<PairTableKokkos<DeviceType>,HALF,false,0,S_TableCompute<DeviceType,TABSTYLE> >
        f(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,f,ev);
      else Kokkos::parallel_for(list->inum,f);
      f.contribute();
    }
  } else {
    if (neighflag == FULL) {
      PairComputeFunctor<PairTableKokkos<DeviceType>,FULL,true,0,S_TableCompute<DeviceType,TABSTYLE> >
        f(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,f,ev);
      else Kokkos::parallel_for(list->inum,f);
      f.contribute();
    } else if (neighflag == HALFTHREAD) {
      PairComputeFunctor<PairTableKokkos<DeviceType>,HALFTHREAD,true,0,S_TableCompute<DeviceType,TABSTYLE> >
        f(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,f,ev);
      else Kokkos::parallel_for(list->inum,f);
      f.contribute();
    } else if (neighflag == HALF) {
      PairComputeFunctor<PairTableKokkos<DeviceType>,HALF,true,0,S_TableCompute<DeviceType,TABSTYLE> >
        f(this,(NeighListKokkos<DeviceType>*) list);
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

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairTableKokkos<DeviceType>::
compute_fpair(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
  (void) i;
  (void) j;
  union_int_float_t rsq_lookup;
  double fpair;
  const int tidx = d_table_const.tabindex(itype,jtype);
  if (Specialisation::TabStyle == LOOKUP) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    fpair = d_table_const.f(tidx,itable);
  } else if (Specialisation::TabStyle == LINEAR) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    const double fraction = (rsq - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx);
    fpair = d_table_const.f(tidx,itable) + fraction*d_table_const.df(tidx,itable);
  } else if (Specialisation::TabStyle == SPLINE) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    const double b = (rsq - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx);
    const double a = 1.0 - b;
    fpair = a * d_table_const.f(tidx,itable) + b * d_table_const.f(tidx,itable+1) +
      ((a*a*a-a)*d_table_const.f2(tidx,itable) + (b*b*b-b)*d_table_const.f2(tidx,itable+1)) *
      d_table_const.deltasq6(tidx);
  } else {
    rsq_lookup.f = rsq;
    int itable = rsq_lookup.i & d_table_const.nmask(tidx);
    itable >>= d_table_const.nshiftbits(tidx);
    const double fraction = (rsq_lookup.f - d_table_const.rsq(tidx,itable)) * d_table_const.drsq(tidx,itable);
    fpair = d_table_const.f(tidx,itable) + fraction*d_table_const.df(tidx,itable);
  }
  return fpair;
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairTableKokkos<DeviceType>::
compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
  (void) i;
  (void) j;
  double evdwl;
  union_int_float_t rsq_lookup;
  const int tidx = d_table_const.tabindex(itype,jtype);
  if (Specialisation::TabStyle == LOOKUP) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    evdwl = d_table_const.e(tidx,itable);
  } else if (Specialisation::TabStyle == LINEAR) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    const double fraction = (rsq - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx);
    evdwl = d_table_const.e(tidx,itable) + fraction*d_table_const.de(tidx,itable);
  } else if (Specialisation::TabStyle == SPLINE) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    const double b = (rsq - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx);
    const double a = 1.0 - b;
    evdwl = a * d_table_const.e(tidx,itable) + b * d_table_const.e(tidx,itable+1) +
        ((a*a*a-a)*d_table_const.e2(tidx,itable) + (b*b*b-b)*d_table_const.e2(tidx,itable+1)) *
        d_table_const.deltasq6(tidx);
  } else {
    rsq_lookup.f = rsq;
    int itable = rsq_lookup.i & d_table_const.nmask(tidx);
    itable >>= d_table_const.nshiftbits(tidx);
    const double fraction = (rsq_lookup.f - d_table_const.rsq(tidx,itable)) * d_table_const.drsq(tidx,itable);
    evdwl = d_table_const.e(tidx,itable) + fraction*d_table_const.de(tidx,itable);
  }
  return evdwl;
}

template<class DeviceType>
void PairTableKokkos<DeviceType>::create_kokkos_tables()
{
  const int tlm1 = tablength-1;

  memoryKK->create_kokkos(d_table->nshiftbits,h_table->nshiftbits,ntables,"Table::nshiftbits");
  memoryKK->create_kokkos(d_table->nmask,h_table->nmask,ntables,"Table::nmask");
  memoryKK->create_kokkos(d_table->innersq,h_table->innersq,ntables,"Table::innersq");
  memoryKK->create_kokkos(d_table->invdelta,h_table->invdelta,ntables,"Table::invdelta");
  memoryKK->create_kokkos(d_table->deltasq6,h_table->deltasq6,ntables,"Table::deltasq6");

  if (tabstyle == LOOKUP) {
    memoryKK->create_kokkos(d_table->e,h_table->e,ntables,tlm1,"Table::e");
    memoryKK->create_kokkos(d_table->f,h_table->f,ntables,tlm1,"Table::f");
  }

  if (tabstyle == LINEAR) {
    memoryKK->create_kokkos(d_table->rsq,h_table->rsq,ntables,tablength,"Table::rsq");
    memoryKK->create_kokkos(d_table->e,h_table->e,ntables,tablength,"Table::e");
    memoryKK->create_kokkos(d_table->f,h_table->f,ntables,tablength,"Table::f");
    memoryKK->create_kokkos(d_table->de,h_table->de,ntables,tlm1,"Table::de");
    memoryKK->create_kokkos(d_table->df,h_table->df,ntables,tlm1,"Table::df");
  }

  if (tabstyle == SPLINE) {
    memoryKK->create_kokkos(d_table->rsq,h_table->rsq,ntables,tablength,"Table::rsq");
    memoryKK->create_kokkos(d_table->e,h_table->e,ntables,tablength,"Table::e");
    memoryKK->create_kokkos(d_table->f,h_table->f,ntables,tablength,"Table::f");
    memoryKK->create_kokkos(d_table->e2,h_table->e2,ntables,tablength,"Table::e2");
    memoryKK->create_kokkos(d_table->f2,h_table->f2,ntables,tablength,"Table::f2");
  }

  if (tabstyle == BITMAP) {
    int ntable = 1 << tablength;
    memoryKK->create_kokkos(d_table->rsq,h_table->rsq,ntables,ntable,"Table::rsq");
    memoryKK->create_kokkos(d_table->e,h_table->e,ntables,ntable,"Table::e");
    memoryKK->create_kokkos(d_table->f,h_table->f,ntables,ntable,"Table::f");
    memoryKK->create_kokkos(d_table->de,h_table->de,ntables,ntable,"Table::de");
    memoryKK->create_kokkos(d_table->df,h_table->df,ntables,ntable,"Table::df");
    memoryKK->create_kokkos(d_table->drsq,h_table->drsq,ntables,ntable,"Table::drsq");
  }



  for (int i=0; i < ntables; i++) {
    Table* tb = &tables[i];

    h_table->nshiftbits[i] = tb->nshiftbits;
    h_table->nmask[i] = tb->nmask;
    h_table->innersq[i] = tb->innersq;
    h_table->invdelta[i] = tb->invdelta;
    h_table->deltasq6[i] = tb->deltasq6;

    for (int j = 0; j < (int) h_table->rsq.extent(1); j++)
      h_table->rsq(i,j) = tb->rsq[j];
    for (int j = 0; j < (int) h_table->drsq.extent(1); j++)
      h_table->drsq(i,j) = tb->drsq[j];
    for (int j = 0; j < (int) h_table->e.extent(1); j++)
      h_table->e(i,j) = tb->e[j];
    for (int j = 0; j < (int) h_table->de.extent(1); j++)
      h_table->de(i,j) = tb->de[j];
    for (int j = 0; j < (int) h_table->f.extent(1); j++)
      h_table->f(i,j) = tb->f[j];
    for (int j = 0; j < (int) h_table->df.extent(1); j++)
      h_table->df(i,j) = tb->df[j];
    for (int j = 0; j < (int) h_table->e2.extent(1); j++)
      h_table->e2(i,j) = tb->e2[j];
    for (int j = 0; j < (int) h_table->f2.extent(1); j++)
      h_table->f2(i,j) = tb->f2[j];
  }


  Kokkos::deep_copy(d_table->nshiftbits,h_table->nshiftbits);
  d_table_const.nshiftbits = d_table->nshiftbits;
  Kokkos::deep_copy(d_table->nmask,h_table->nmask);
  d_table_const.nmask = d_table->nmask;
  Kokkos::deep_copy(d_table->innersq,h_table->innersq);
  d_table_const.innersq = d_table->innersq;
  Kokkos::deep_copy(d_table->invdelta,h_table->invdelta);
  d_table_const.invdelta = d_table->invdelta;
  Kokkos::deep_copy(d_table->deltasq6,h_table->deltasq6);
  d_table_const.deltasq6 = d_table->deltasq6;

  if (tabstyle == LOOKUP) {
    Kokkos::deep_copy(d_table->e,h_table->e);
    d_table_const.e = d_table->e;
    Kokkos::deep_copy(d_table->f,h_table->f);
    d_table_const.f = d_table->f;
  }

  if (tabstyle == LINEAR) {
    Kokkos::deep_copy(d_table->rsq,h_table->rsq);
    d_table_const.rsq = d_table->rsq;
    Kokkos::deep_copy(d_table->e,h_table->e);
    d_table_const.e = d_table->e;
    Kokkos::deep_copy(d_table->f,h_table->f);
    d_table_const.f = d_table->f;
    Kokkos::deep_copy(d_table->de,h_table->de);
    d_table_const.de = d_table->de;
    Kokkos::deep_copy(d_table->df,h_table->df);
    d_table_const.df = d_table->df;
  }

  if (tabstyle == SPLINE) {
    Kokkos::deep_copy(d_table->rsq,h_table->rsq);
    d_table_const.rsq = d_table->rsq;
    Kokkos::deep_copy(d_table->e,h_table->e);
    d_table_const.e = d_table->e;
    Kokkos::deep_copy(d_table->f,h_table->f);
    d_table_const.f = d_table->f;
    Kokkos::deep_copy(d_table->e2,h_table->e2);
    d_table_const.e2 = d_table->e2;
    Kokkos::deep_copy(d_table->f2,h_table->f2);
    d_table_const.f2 = d_table->f2;
  }

  if (tabstyle == BITMAP) {
    Kokkos::deep_copy(d_table->rsq,h_table->rsq);
    d_table_const.rsq = d_table->rsq;
    Kokkos::deep_copy(d_table->e,h_table->e);
    d_table_const.e = d_table->e;
    Kokkos::deep_copy(d_table->f,h_table->f);
    d_table_const.f = d_table->f;
    Kokkos::deep_copy(d_table->de,h_table->de);
    d_table_const.de = d_table->de;
    Kokkos::deep_copy(d_table->df,h_table->df);
    d_table_const.df = d_table->df;
    Kokkos::deep_copy(d_table->drsq,h_table->drsq);
    d_table_const.drsq = d_table->drsq;
  }

  Kokkos::deep_copy(d_table->cutsq,h_table->cutsq);
  d_table_const.cutsq = d_table->cutsq;
  Kokkos::deep_copy(d_table->tabindex,h_table->tabindex);
  d_table_const.tabindex = d_table->tabindex;

  update_table = 0;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableKokkos<DeviceType>::allocate()
{
  allocated = 1;
  const int nt = atom->ntypes + 1;

  memory->create(setflag,nt,nt,"pair:setflag");
  memoryKK->create_kokkos(d_table->cutsq,h_table->cutsq,cutsq,nt,nt,"pair:cutsq");
  memoryKK->create_kokkos(d_table->tabindex,h_table->tabindex,tabindex,nt,nt,"pair:tabindex");
  d_table_const.cutsq = d_table->cutsq;
  d_table_const.tabindex = d_table->tabindex;

  memset(&setflag[0][0],0,nt*nt*sizeof(int));
  memset(&cutsq[0][0],0,nt*nt*sizeof(double));
  memset(&tabindex[0][0],0,nt*nt*sizeof(int));
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableKokkos<DeviceType>::settings(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal pair_style command");

  // new settings

  if (strcmp(arg[0],"lookup") == 0) tabstyle = LOOKUP;
  else if (strcmp(arg[0],"linear") == 0) tabstyle = LINEAR;
  else if (strcmp(arg[0],"spline") == 0) tabstyle = SPLINE;
  else if (strcmp(arg[0],"bitmap") == 0) tabstyle = BITMAP;
  else error->all(FLERR,"Unknown table style in pair_style command");

  tablength = utils::inumeric(FLERR,arg[1],false,lmp);
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

    d_table_const.tabindex = d_table->tabindex = typename ArrayTypes<DeviceType>::t_int_2d();
    h_table->tabindex = typename ArrayTypes<LMPHostType>::t_int_2d();

    d_table_const.cutsq = d_table->cutsq = typename ArrayTypes<DeviceType>::t_ffloat_2d();
    h_table->cutsq = typename ArrayTypes<LMPHostType>::t_ffloat_2d();
  }
  allocated = 0;

  ntables = 0;
  tables = nullptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairTableKokkos<DeviceType>::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  tabindex[j][i] = tabindex[i][j];

  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_cutsq[j][i] = m_cutsq[i][j] = tables[tabindex[i][j]].cut*tables[tabindex[i][j]].cut;
  }

  return tables[tabindex[i][j]].cut;
}

/* ----------------------------------------------------------------------
   compute r,e,f vectors from splined values
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableKokkos<DeviceType>::compute_table(Table *tb)
{
  update_table = 1;
  PairTable::compute_table(tb);
}

template<class DeviceType>
void PairTableKokkos<DeviceType>::init_style()
{
  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->add_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();
}

namespace LAMMPS_NS {
template class PairTableKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairTableKokkos<LMPHostType>;
#endif

}

