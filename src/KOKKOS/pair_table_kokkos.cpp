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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "pair_table_kokkos.h"
#include "kokkos.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

enum{NONE,RLINEAR,RSQ,BMP};
enum{FULL,HALFTHREAD,HALF};

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairTableKokkos<DeviceType>::PairTableKokkos(LAMMPS *lmp) : Pair(lmp)
{
  update_table = 0;
  atomKK = (AtomKokkos *) atom;
  ntables = 0;
  tables = NULL;
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
/*  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(tabindex);
  }*/
  delete h_table;
  delete d_table;

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairTableKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

template<class DeviceType>
template<int TABSTYLE>
void PairTableKokkos<DeviceType>::compute_style(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL || neighflag == FULLCLUSTER) no_virial_fdotr_compute = 1;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

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

  EV_FLOAT ev;
  if(atom->ntypes > MAX_TYPES_STACKPARAMS) {
    if (neighflag == FULL) {
      PairComputeFunctor<PairTableKokkos<DeviceType>,FULL,false,S_TableCompute<DeviceType,TABSTYLE> >
        ff(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,ff,ev);
      else Kokkos::parallel_for(list->inum,ff);
    } else if (neighflag == HALFTHREAD) {
      PairComputeFunctor<PairTableKokkos<DeviceType>,HALFTHREAD,false,S_TableCompute<DeviceType,TABSTYLE> >
        ff(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,ff,ev);
      else Kokkos::parallel_for(list->inum,ff);
    } else if (neighflag == HALF) {
      PairComputeFunctor<PairTableKokkos<DeviceType>,HALF,false,S_TableCompute<DeviceType,TABSTYLE> >
        f(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,f,ev);
      else Kokkos::parallel_for(list->inum,f);
    } else if (neighflag == N2) {
      PairComputeFunctor<PairTableKokkos<DeviceType>,N2,false,S_TableCompute<DeviceType,TABSTYLE> >
        f(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(nlocal,f,ev);
      else Kokkos::parallel_for(nlocal,f);
    } else if (neighflag == FULLCLUSTER) {
      typedef PairComputeFunctor<PairTableKokkos<DeviceType>,FULLCLUSTER,false,S_TableCompute<DeviceType,TABSTYLE> >
        f_type;
      f_type f(this,(NeighListKokkos<DeviceType>*) list);
      #ifdef KOKKOS_HAVE_CUDA
        const int teamsize = Kokkos::Impl::is_same<typename f_type::device_type, Kokkos::Cuda>::value ? 256 : 1;
      #else
        const int teamsize = 1;
      #endif
      const int nteams = (list->inum*f_type::vectorization::increment+teamsize-1)/teamsize;
      Kokkos::TeamPolicy<DeviceType> config(nteams,teamsize);
      if (eflag || vflag) Kokkos::parallel_reduce(config,f,ev);
      else Kokkos::parallel_for(config,f);
    }
  } else {
    if (neighflag == FULL) {
      PairComputeFunctor<PairTableKokkos<DeviceType>,FULL,true,S_TableCompute<DeviceType,TABSTYLE> >
        f(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,f,ev);
      else Kokkos::parallel_for(list->inum,f);
    } else if (neighflag == HALFTHREAD) {
      PairComputeFunctor<PairTableKokkos<DeviceType>,HALFTHREAD,true,S_TableCompute<DeviceType,TABSTYLE> >
        f(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,f,ev);
      else Kokkos::parallel_for(list->inum,f);
    } else if (neighflag == HALF) {
      PairComputeFunctor<PairTableKokkos<DeviceType>,HALF,true,S_TableCompute<DeviceType,TABSTYLE> >
        f(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,f,ev);
      else Kokkos::parallel_for(list->inum,f);
    } else if (neighflag == N2) {
      PairComputeFunctor<PairTableKokkos<DeviceType>,N2,true,S_TableCompute<DeviceType,TABSTYLE> >
        f(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(nlocal,f,ev);
      else Kokkos::parallel_for(nlocal,f);
    } else if (neighflag == FULLCLUSTER) {
      typedef PairComputeFunctor<PairTableKokkos<DeviceType>,FULLCLUSTER,true,S_TableCompute<DeviceType,TABSTYLE> >
        f_type;
      f_type f(this,(NeighListKokkos<DeviceType>*) list);
      #ifdef KOKKOS_HAVE_CUDA
        const int teamsize = Kokkos::Impl::is_same<typename f_type::device_type, Kokkos::Cuda>::value ? 256 : 1;
      #else
        const int teamsize = 1;
      #endif
      const int nteams = (list->inum*f_type::vectorization::increment+teamsize-1)/teamsize;
      Kokkos::TeamPolicy<DeviceType> config(nteams,teamsize);
      if (eflag || vflag) Kokkos::parallel_reduce(config,f,ev);
      else Kokkos::parallel_for(config,f);
    }
  }
  DeviceType::fence();

  if (eflag) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);
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
  //const Table* const tb = &tables[tabindex[itype][jtype]];

  //if (rsq < d_table_const.innersq(tidx))
  //  error->one(FLERR,"Pair distance < table inner cutoff");

  if (Specialisation::TabStyle == LOOKUP) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    //if (itable >= tlm1)
    //  error->one(FLERR,"Pair distance > table outer cutoff");
    fpair = d_table_const.f(tidx,itable);
  } else if (Specialisation::TabStyle == LINEAR) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    //if (itable >= tlm1)
    //  error->one(FLERR,"Pair distance > table outer cutoff");
    const double fraction = (rsq - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx);
    fpair = d_table_const.f(tidx,itable) + fraction*d_table_const.df(tidx,itable);
  } else if (Specialisation::TabStyle == SPLINE) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    //if (itable >= tlm1)
    //  error->one(FLERR,"Pair distance > table outer cutoff");
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
  //const Table* const tb = &tables[tabindex[itype][jtype]];

  //if (rsq < d_table_const.innersq(tidx))
  //  error->one(FLERR,"Pair distance < table inner cutoff");

  if (Specialisation::TabStyle == LOOKUP) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    //if (itable >= tlm1)
    //  error->one(FLERR,"Pair distance > table outer cutoff");
    evdwl = d_table_const.e(tidx,itable);
  } else if (Specialisation::TabStyle == LINEAR) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    //if (itable >= tlm1)
    //  error->one(FLERR,"Pair distance > table outer cutoff");
    const double fraction = (rsq - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx);
    evdwl = d_table_const.e(tidx,itable) + fraction*d_table_const.de(tidx,itable);
  } else if (Specialisation::TabStyle == SPLINE) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    //if (itable >= tlm1)
    //  error->one(FLERR,"Pair distance > table outer cutoff");
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

  memory->create_kokkos(d_table->nshiftbits,h_table->nshiftbits,ntables,"Table::nshiftbits");
  memory->create_kokkos(d_table->nmask,h_table->nmask,ntables,"Table::nmask");
  memory->create_kokkos(d_table->innersq,h_table->innersq,ntables,"Table::innersq");
  memory->create_kokkos(d_table->invdelta,h_table->invdelta,ntables,"Table::invdelta");
  memory->create_kokkos(d_table->deltasq6,h_table->deltasq6,ntables,"Table::deltasq6");

  if(tabstyle == LOOKUP) {
    memory->create_kokkos(d_table->e,h_table->e,ntables,tlm1,"Table::e");
    memory->create_kokkos(d_table->f,h_table->f,ntables,tlm1,"Table::f");
  }

  if(tabstyle == LINEAR) {
    memory->create_kokkos(d_table->rsq,h_table->rsq,ntables,tablength,"Table::rsq");
    memory->create_kokkos(d_table->e,h_table->e,ntables,tablength,"Table::e");
    memory->create_kokkos(d_table->f,h_table->f,ntables,tablength,"Table::f");
    memory->create_kokkos(d_table->de,h_table->de,ntables,tlm1,"Table::de");
    memory->create_kokkos(d_table->df,h_table->df,ntables,tlm1,"Table::df");
  }

  if(tabstyle == SPLINE) {
    memory->create_kokkos(d_table->rsq,h_table->rsq,ntables,tablength,"Table::rsq");
    memory->create_kokkos(d_table->e,h_table->e,ntables,tablength,"Table::e");
    memory->create_kokkos(d_table->f,h_table->f,ntables,tablength,"Table::f");
    memory->create_kokkos(d_table->e2,h_table->e2,ntables,tablength,"Table::e2");
    memory->create_kokkos(d_table->f2,h_table->f2,ntables,tablength,"Table::f2");
  }

  if(tabstyle == BITMAP) {
    int ntable = 1 << tablength;
    memory->create_kokkos(d_table->rsq,h_table->rsq,ntables,ntable,"Table::rsq");
    memory->create_kokkos(d_table->e,h_table->e,ntables,ntable,"Table::e");
    memory->create_kokkos(d_table->f,h_table->f,ntables,ntable,"Table::f");
    memory->create_kokkos(d_table->de,h_table->de,ntables,ntable,"Table::de");
    memory->create_kokkos(d_table->df,h_table->df,ntables,ntable,"Table::df");
    memory->create_kokkos(d_table->drsq,h_table->drsq,ntables,ntable,"Table::drsq");
  }



  for(int i=0; i < ntables; i++) {
    Table* tb = &tables[i];

    h_table->nshiftbits[i] = tb->nshiftbits;
    h_table->nmask[i] = tb->nmask;
    h_table->innersq[i] = tb->innersq;
    h_table->invdelta[i] = tb->invdelta;
    h_table->deltasq6[i] = tb->deltasq6;

    for(int j = 0; j<h_table->rsq.dimension_1(); j++)
      h_table->rsq(i,j) = tb->rsq[j];
    for(int j = 0; j<h_table->drsq.dimension_1(); j++)
      h_table->drsq(i,j) = tb->drsq[j];
    for(int j = 0; j<h_table->e.dimension_1(); j++)
      h_table->e(i,j) = tb->e[j];
    for(int j = 0; j<h_table->de.dimension_1(); j++)
      h_table->de(i,j) = tb->de[j];
    for(int j = 0; j<h_table->f.dimension_1(); j++)
      h_table->f(i,j) = tb->f[j];
    for(int j = 0; j<h_table->df.dimension_1(); j++)
      h_table->df(i,j) = tb->df[j];
    for(int j = 0; j<h_table->e2.dimension_1(); j++)
      h_table->e2(i,j) = tb->e2[j];
    for(int j = 0; j<h_table->f2.dimension_1(); j++)
      h_table->f2(i,j) = tb->f2[j];
  }

  
  Kokkos::deep_copy(d_table->nshiftbits,h_table->nshiftbits);
  Kokkos::deep_copy(d_table->nmask,h_table->nmask);
  Kokkos::deep_copy(d_table->innersq,h_table->innersq);
  Kokkos::deep_copy(d_table->invdelta,h_table->invdelta);
  Kokkos::deep_copy(d_table->deltasq6,h_table->deltasq6);
  Kokkos::deep_copy(d_table->rsq,h_table->rsq);
  Kokkos::deep_copy(d_table->drsq,h_table->drsq);
  Kokkos::deep_copy(d_table->e,h_table->e);
  Kokkos::deep_copy(d_table->de,h_table->de);
  Kokkos::deep_copy(d_table->f,h_table->f);
  Kokkos::deep_copy(d_table->df,h_table->df);
  Kokkos::deep_copy(d_table->e2,h_table->e2);
  Kokkos::deep_copy(d_table->f2,h_table->f2);
  Kokkos::deep_copy(d_table->tabindex,h_table->tabindex);

  d_table_const.nshiftbits = d_table->nshiftbits;
  d_table_const.nmask = d_table->nmask;
  d_table_const.innersq = d_table->innersq;
  d_table_const.invdelta = d_table->invdelta;
  d_table_const.deltasq6 = d_table->deltasq6;
  d_table_const.rsq = d_table->rsq;
  d_table_const.drsq = d_table->drsq;
  d_table_const.e = d_table->e;
  d_table_const.de = d_table->de;
  d_table_const.f = d_table->f;
  d_table_const.df = d_table->df;
  d_table_const.e2 = d_table->e2;
  d_table_const.f2 = d_table->f2;


  Kokkos::deep_copy(d_table->cutsq,h_table->cutsq);
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
  memory->create_kokkos(d_table->cutsq,h_table->cutsq,cutsq,nt,nt,"pair:cutsq");
  memory->create_kokkos(d_table->tabindex,h_table->tabindex,tabindex,nt,nt,"pair:tabindex");

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
    
    d_table_const.tabindex = d_table->tabindex = typename ArrayTypes<DeviceType>::t_int_2d();
    h_table->tabindex = typename ArrayTypes<LMPHostType>::t_int_2d();

    d_table_const.cutsq = d_table->cutsq = typename ArrayTypes<DeviceType>::t_ffloat_2d();  
    h_table->cutsq = typename ArrayTypes<LMPHostType>::t_ffloat_2d();  
  }
  allocated = 0;

  ntables = 0;
  tables = NULL;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableKokkos<DeviceType>::coeff(int narg, char **arg)
{
  if (narg != 4 && narg != 5) error->all(FLERR,"Illegal pair_coeff command");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  int me;
  MPI_Comm_rank(world,&me);
  tables = (Table *)
    memory->srealloc(tables,(ntables+1)*sizeof(Table),"pair:tables");
  Table *tb = &tables[ntables];
  null_table(tb);
  if (me == 0) read_table(tb,arg[2],arg[3]);
  bcast_table(tb);

  // set table cutoff

  if (narg == 5) tb->cut = force->numeric(FLERR,arg[4]);
  else if (tb->rflag) tb->cut = tb->rhi;
  else tb->cut = tb->rfile[tb->ninput-1];

  // error check on table parameters
  // insure cutoff is within table
  // for BITMAP tables, file values can be in non-ascending order

  if (tb->ninput <= 1) error->one(FLERR,"Invalid pair table length");
  double rlo,rhi;
  if (tb->rflag == 0) {
    rlo = tb->rfile[0];
    rhi = tb->rfile[tb->ninput-1];
  } else {
    rlo = tb->rlo;
    rhi = tb->rhi;
  }
  if (tb->cut <= rlo || tb->cut > rhi)
    error->all(FLERR,"Invalid pair table cutoff");
  if (rlo <= 0.0) error->all(FLERR,"Invalid pair table cutoff");

  // match = 1 if don't need to spline read-in tables
  // this is only the case if r values needed by final tables
  //   exactly match r values read from file
  // for tabstyle SPLINE, always need to build spline tables

  tb->match = 0;
  if (tabstyle == LINEAR && tb->ninput == tablength &&
      tb->rflag == RSQ && tb->rhi == tb->cut) tb->match = 1;
  if (tabstyle == BITMAP && tb->ninput == 1 << tablength &&
      tb->rflag == BMP && tb->rhi == tb->cut) tb->match = 1;
  if (tb->rflag == BMP && tb->match == 0)
    error->all(FLERR,"Bitmapped table in file does not match requested table");

  // spline read-in values and compute r,e,f vectors within table

  if (tb->match == 0) spline_table(tb);
  compute_table(tb);

  // store ptr to table in tabindex

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      tabindex[i][j] = ntables;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Illegal pair_coeff command");
  ntables++;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairTableKokkos<DeviceType>::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  tabindex[j][i] = tabindex[i][j];

  if(i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_cutsq[j][i] = m_cutsq[i][j] = tables[tabindex[i][j]].cut*tables[tabindex[i][j]].cut;
  }

  return tables[tabindex[i][j]].cut;
}

/* ----------------------------------------------------------------------
   read a table section from a tabulated potential file
   only called by proc 0
   this function sets these values in Table:
     ninput,rfile,efile,ffile,rflag,rlo,rhi,fpflag,fplo,fphi,ntablebits
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableKokkos<DeviceType>::read_table(Table *tb, char *file, char *keyword)
{
  char line[MAXLINE];

  // open file

  FILE *fp = force->open_potential(file);
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",file);
    error->one(FLERR,str);
  }

  // loop until section found with matching keyword

  while (1) {
    if (fgets(line,MAXLINE,fp) == NULL)
      error->one(FLERR,"Did not find keyword in table file");
    if (strspn(line," \t\n\r") == strlen(line)) continue;  // blank line
    if (line[0] == '#') continue;                          // comment
    char *word = strtok(line," \t\n\r");
    if (strcmp(word,keyword) == 0) break;           // matching keyword
    fgets(line,MAXLINE,fp);                         // no match, skip section
    param_extract(tb,line);
    fgets(line,MAXLINE,fp);
    for (int i = 0; i < tb->ninput; i++) fgets(line,MAXLINE,fp);
  }

  // read args on 2nd line of section
  // allocate table arrays for file values

  fgets(line,MAXLINE,fp);
  param_extract(tb,line);
  memory->create(tb->rfile,tb->ninput,"pair:rfile");
  memory->create(tb->efile,tb->ninput,"pair:efile");
  memory->create(tb->ffile,tb->ninput,"pair:ffile");

  // setup bitmap parameters for table to read in

  tb->ntablebits = 0;
  int masklo,maskhi,nmask,nshiftbits;
  if (tb->rflag == BMP) {
    while (1 << tb->ntablebits < tb->ninput) tb->ntablebits++;
    if (1 << tb->ntablebits != tb->ninput)
      error->one(FLERR,"Bitmapped table is incorrect length in table file");
    init_bitmap(tb->rlo,tb->rhi,tb->ntablebits,masklo,maskhi,nmask,nshiftbits);
  }

  // read r,e,f table values from file
  // if rflag set, compute r
  // if rflag not set, use r from file

  int itmp;
  double rtmp;
  union_int_float_t rsq_lookup;

  fgets(line,MAXLINE,fp);
  for (int i = 0; i < tb->ninput; i++) {
    fgets(line,MAXLINE,fp);
    sscanf(line,"%d %lg %lg %lg",&itmp,&rtmp,&tb->efile[i],&tb->ffile[i]);

    if (tb->rflag == RLINEAR)
      rtmp = tb->rlo + (tb->rhi - tb->rlo)*i/(tb->ninput-1);
    else if (tb->rflag == RSQ) {
      rtmp = tb->rlo*tb->rlo +
        (tb->rhi*tb->rhi - tb->rlo*tb->rlo)*i/(tb->ninput-1);
      rtmp = sqrt(rtmp);
    } else if (tb->rflag == BMP) {
      rsq_lookup.i = i << nshiftbits;
      rsq_lookup.i |= masklo;
      if (rsq_lookup.f < tb->rlo*tb->rlo) {
        rsq_lookup.i = i << nshiftbits;
        rsq_lookup.i |= maskhi;
      }
      rtmp = sqrtf(rsq_lookup.f);
    }

    tb->rfile[i] = rtmp;
  }

  // close file

  fclose(fp);
}

/* ----------------------------------------------------------------------
   broadcast read-in table info from proc 0 to other procs
   this function communicates these values in Table:
     ninput,rfile,efile,ffile,rflag,rlo,rhi,fpflag,fplo,fphi
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableKokkos<DeviceType>::bcast_table(Table *tb)
{
  MPI_Bcast(&tb->ninput,1,MPI_INT,0,world);

  int me;
  MPI_Comm_rank(world,&me);
  if (me > 0) {
    memory->create(tb->rfile,tb->ninput,"pair:rfile");
    memory->create(tb->efile,tb->ninput,"pair:efile");
    memory->create(tb->ffile,tb->ninput,"pair:ffile");
  }

  MPI_Bcast(tb->rfile,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->efile,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->ffile,tb->ninput,MPI_DOUBLE,0,world);

  MPI_Bcast(&tb->rflag,1,MPI_INT,0,world);
  if (tb->rflag) {
    MPI_Bcast(&tb->rlo,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&tb->rhi,1,MPI_DOUBLE,0,world);
  }
  MPI_Bcast(&tb->fpflag,1,MPI_INT,0,world);
  if (tb->fpflag) {
    MPI_Bcast(&tb->fplo,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&tb->fphi,1,MPI_DOUBLE,0,world);
  }
}

/* ----------------------------------------------------------------------
   build spline representation of e,f over entire range of read-in table
   this function sets these values in Table: e2file,f2file
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableKokkos<DeviceType>::spline_table(Table *tb)
{
  memory->create(tb->e2file,tb->ninput,"pair:e2file");
  memory->create(tb->f2file,tb->ninput,"pair:f2file");

  double ep0 = - tb->ffile[0];
  double epn = - tb->ffile[tb->ninput-1];
  spline(tb->rfile,tb->efile,tb->ninput,ep0,epn,tb->e2file);

  if (tb->fpflag == 0) {
    tb->fplo = (tb->ffile[1] - tb->ffile[0]) / (tb->rfile[1] - tb->rfile[0]);
    tb->fphi = (tb->ffile[tb->ninput-1] - tb->ffile[tb->ninput-2]) /
      (tb->rfile[tb->ninput-1] - tb->rfile[tb->ninput-2]);
  }

  double fp0 = tb->fplo;
  double fpn = tb->fphi;
  spline(tb->rfile,tb->ffile,tb->ninput,fp0,fpn,tb->f2file);
}

/* ----------------------------------------------------------------------
   extract attributes from parameter line in table section
   format of line: N value R/RSQ/BITMAP lo hi FP fplo fphi
   N is required, other params are optional
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableKokkos<DeviceType>::param_extract(Table *tb, char *line)
{
  tb->ninput = 0;
  tb->rflag = NONE;
  tb->fpflag = 0;

  char *word = strtok(line," \t\n\r\f");
  while (word) {
    if (strcmp(word,"N") == 0) {
      word = strtok(NULL," \t\n\r\f");
      tb->ninput = atoi(word);
    } else if (strcmp(word,"R") == 0 || strcmp(word,"RSQ") == 0 ||
               strcmp(word,"BITMAP") == 0) {
      if (strcmp(word,"R") == 0) tb->rflag = RLINEAR;
      else if (strcmp(word,"RSQ") == 0) tb->rflag = RSQ;
      else if (strcmp(word,"BITMAP") == 0) tb->rflag = BMP;
      word = strtok(NULL," \t\n\r\f");
      tb->rlo = atof(word);
      word = strtok(NULL," \t\n\r\f");
      tb->rhi = atof(word);
    } else if (strcmp(word,"FP") == 0) {
      tb->fpflag = 1;
      word = strtok(NULL," \t\n\r\f");
      tb->fplo = atof(word);
      word = strtok(NULL," \t\n\r\f");
      tb->fphi = atof(word);
    } else {
      error->one(FLERR,"Invalid keyword in pair table parameters");
    }
    word = strtok(NULL," \t\n\r\f");
  }

  if (tb->ninput == 0) error->one(FLERR,"Pair table parameters did not set N");
}

/* ----------------------------------------------------------------------
   compute r,e,f vectors from splined values
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableKokkos<DeviceType>::compute_table(Table *tb)
{
  update_table = 1;
  int tlm1 = tablength-1;

  // inner = inner table bound
  // cut = outer table bound
  // delta = table spacing in rsq for N-1 bins

  double inner;
  if (tb->rflag) inner = tb->rlo;
  else inner = tb->rfile[0];
  tb->innersq = inner*inner;
  tb->delta = (tb->cut*tb->cut - tb->innersq) / tlm1;
  tb->invdelta = 1.0/tb->delta;

  // direct lookup tables
  // N-1 evenly spaced bins in rsq from inner to cut
  // e,f = value at midpt of bin
  // e,f are N-1 in length since store 1 value at bin midpt
  // f is converted to f/r when stored in f[i]
  // e,f are never a match to read-in values, always computed via spline interp

  if (tabstyle == LOOKUP) {
    memory->create(tb->e,tlm1,"pair:e");
    memory->create(tb->f,tlm1,"pair:f");

    double r,rsq;
    for (int i = 0; i < tlm1; i++) {
      rsq = tb->innersq + (i+0.5)*tb->delta;
      r = sqrt(rsq);
      tb->e[i] = splint(tb->rfile,tb->efile,tb->e2file,tb->ninput,r);
      tb->f[i] = splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,r)/r;
    }
  }

  // linear tables
  // N-1 evenly spaced bins in rsq from inner to cut
  // rsq,e,f = value at lower edge of bin
  // de,df values = delta from lower edge to upper edge of bin
  // rsq,e,f are N in length so de,df arrays can compute difference
  // f is converted to f/r when stored in f[i]
  // e,f can match read-in values, else compute via spline interp

  if (tabstyle == LINEAR) {
    memory->create(tb->rsq,tablength,"pair:rsq");
    memory->create(tb->e,tablength,"pair:e");
    memory->create(tb->f,tablength,"pair:f");
    memory->create(tb->de,tlm1,"pair:de");
    memory->create(tb->df,tlm1,"pair:df");

    double r,rsq;
    for (int i = 0; i < tablength; i++) {
      rsq = tb->innersq + i*tb->delta;
      r = sqrt(rsq);
      tb->rsq[i] = rsq;
      if (tb->match) {
        tb->e[i] = tb->efile[i];
        tb->f[i] = tb->ffile[i]/r;
      } else {
        tb->e[i] = splint(tb->rfile,tb->efile,tb->e2file,tb->ninput,r);
        tb->f[i] = splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,r)/r;
      }
    }

    for (int i = 0; i < tlm1; i++) {
      tb->de[i] = tb->e[i+1] - tb->e[i];
      tb->df[i] = tb->f[i+1] - tb->f[i];
    }
  }

  // cubic spline tables
  // N-1 evenly spaced bins in rsq from inner to cut
  // rsq,e,f = value at lower edge of bin
  // e2,f2 = spline coefficient for each bin
  // rsq,e,f,e2,f2 are N in length so have N-1 spline bins
  // f is converted to f/r after e is splined
  // e,f can match read-in values, else compute via spline interp

  if (tabstyle == SPLINE) {
    memory->create(tb->rsq,tablength,"pair:rsq");
    memory->create(tb->e,tablength,"pair:e");
    memory->create(tb->f,tablength,"pair:f");
    memory->create(tb->e2,tablength,"pair:e2");
    memory->create(tb->f2,tablength,"pair:f2");

    tb->deltasq6 = tb->delta*tb->delta / 6.0;

    double r,rsq;
    for (int i = 0; i < tablength; i++) {
      rsq = tb->innersq + i*tb->delta;
      r = sqrt(rsq);
      tb->rsq[i] = rsq;
      if (tb->match) {
        tb->e[i] = tb->efile[i];
        tb->f[i] = tb->ffile[i]/r;
      } else {
        tb->e[i] = splint(tb->rfile,tb->efile,tb->e2file,tb->ninput,r);
        tb->f[i] = splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,r);
      }
    }

    // ep0,epn = dh/dg at inner and at cut
    // h(r) = e(r) and g(r) = r^2
    // dh/dg = (de/dr) / 2r = -f/2r

    double ep0 = - tb->f[0] / (2.0 * sqrt(tb->innersq));
    double epn = - tb->f[tlm1] / (2.0 * tb->cut);
    spline(tb->rsq,tb->e,tablength,ep0,epn,tb->e2);

    // fp0,fpn = dh/dg at inner and at cut
    // h(r) = f(r)/r and g(r) = r^2
    // dh/dg = (1/r df/dr - f/r^2) / 2r
    // dh/dg in secant approx = (f(r2)/r2 - f(r1)/r1) / (g(r2) - g(r1))

    double fp0,fpn;
    double secant_factor = 0.1;
    if (tb->fpflag) fp0 = (tb->fplo/sqrt(tb->innersq) - tb->f[0]/tb->innersq) /
      (2.0 * sqrt(tb->innersq));
    else {
      double rsq1 = tb->innersq;
      double rsq2 = rsq1 + secant_factor*tb->delta;
      fp0 = (splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,sqrt(rsq2)) /
             sqrt(rsq2) - tb->f[0] / sqrt(rsq1)) / (secant_factor*tb->delta);
    }

    if (tb->fpflag && tb->cut == tb->rfile[tb->ninput-1]) fpn =
      (tb->fphi/tb->cut - tb->f[tlm1]/(tb->cut*tb->cut)) / (2.0 * tb->cut);
    else {
      double rsq2 = tb->cut * tb->cut;
      double rsq1 = rsq2 - secant_factor*tb->delta;
      fpn = (tb->f[tlm1] / sqrt(rsq2) -
             splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,sqrt(rsq1)) /
             sqrt(rsq1)) / (secant_factor*tb->delta);
    }

    for (int i = 0; i < tablength; i++) tb->f[i] /= sqrt(tb->rsq[i]);
    spline(tb->rsq,tb->f,tablength,fp0,fpn,tb->f2);
  }

  // bitmapped linear tables
  // 2^N bins from inner to cut, spaced in bitmapped manner
  // f is converted to f/r when stored in f[i]
  // e,f can match read-in values, else compute via spline interp

  if (tabstyle == BITMAP) {
    double r;
    union_int_float_t rsq_lookup;
    int masklo,maskhi;

    // linear lookup tables of length ntable = 2^n
    // stored value = value at lower edge of bin

    init_bitmap(inner,tb->cut,tablength,masklo,maskhi,tb->nmask,tb->nshiftbits);
    int ntable = 1 << tablength;
    int ntablem1 = ntable - 1;

    memory->create(tb->rsq,ntable,"pair:rsq");
    memory->create(tb->e,ntable,"pair:e");
    memory->create(tb->f,ntable,"pair:f");
    memory->create(tb->de,ntable,"pair:de");
    memory->create(tb->df,ntable,"pair:df");
    memory->create(tb->drsq,ntable,"pair:drsq");

    union_int_float_t minrsq_lookup;
    minrsq_lookup.i = 0 << tb->nshiftbits;
    minrsq_lookup.i |= maskhi;

    for (int i = 0; i < ntable; i++) {
      rsq_lookup.i = i << tb->nshiftbits;
      rsq_lookup.i |= masklo;
      if (rsq_lookup.f < tb->innersq) {
        rsq_lookup.i = i << tb->nshiftbits;
        rsq_lookup.i |= maskhi;
      }
      r = sqrtf(rsq_lookup.f);
      tb->rsq[i] = rsq_lookup.f;
      if (tb->match) {
        tb->e[i] = tb->efile[i];
        tb->f[i] = tb->ffile[i]/r;
      } else {
        tb->e[i] = splint(tb->rfile,tb->efile,tb->e2file,tb->ninput,r);
        tb->f[i] = splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,r)/r;
      }
      minrsq_lookup.f = MIN(minrsq_lookup.f,rsq_lookup.f);
    }

    tb->innersq = minrsq_lookup.f;

    for (int i = 0; i < ntablem1; i++) {
      tb->de[i] = tb->e[i+1] - tb->e[i];
      tb->df[i] = tb->f[i+1] - tb->f[i];
      tb->drsq[i] = 1.0/(tb->rsq[i+1] - tb->rsq[i]);
    }

    // get the delta values for the last table entries
    // tables are connected periodically between 0 and ntablem1

    tb->de[ntablem1] = tb->e[0] - tb->e[ntablem1];
    tb->df[ntablem1] = tb->f[0] - tb->f[ntablem1];
    tb->drsq[ntablem1] = 1.0/(tb->rsq[0] - tb->rsq[ntablem1]);

    // get the correct delta values at itablemax
    // smallest r is in bin itablemin
    // largest r is in bin itablemax, which is itablemin-1,
    //   or ntablem1 if itablemin=0

    // deltas at itablemax only needed if corresponding rsq < cut*cut
    // if so, compute deltas between rsq and cut*cut
    //   if tb->match, data at cut*cut is unavailable, so we'll take
    //   deltas at itablemax-1 as a good approximation

    double e_tmp,f_tmp;
    int itablemin = minrsq_lookup.i & tb->nmask;
    itablemin >>= tb->nshiftbits;
    int itablemax = itablemin - 1;
    if (itablemin == 0) itablemax = ntablem1;
    int itablemaxm1 = itablemax - 1;
    if (itablemax == 0) itablemaxm1 = ntablem1;
    rsq_lookup.i = itablemax << tb->nshiftbits;
    rsq_lookup.i |= maskhi;
    if (rsq_lookup.f < tb->cut*tb->cut) {
      if (tb->match) {
        tb->de[itablemax] = tb->de[itablemaxm1];
        tb->df[itablemax] = tb->df[itablemaxm1];
        tb->drsq[itablemax] = tb->drsq[itablemaxm1];
      } else {
            rsq_lookup.f = tb->cut*tb->cut;
        r = sqrtf(rsq_lookup.f);
        e_tmp = splint(tb->rfile,tb->efile,tb->e2file,tb->ninput,r);
        f_tmp = splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,r)/r;
        tb->de[itablemax] = e_tmp - tb->e[itablemax];
        tb->df[itablemax] = f_tmp - tb->f[itablemax];
        tb->drsq[itablemax] = 1.0/(rsq_lookup.f - tb->rsq[itablemax]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   set all ptrs in a table to NULL, so can be freed safely
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableKokkos<DeviceType>::null_table(Table *tb)
{
  tb->rfile = tb->efile = tb->ffile = NULL;
  tb->e2file = tb->f2file = NULL;
  tb->rsq = tb->drsq = tb->e = tb->de = NULL;
  tb->f = tb->df = tb->e2 = tb->f2 = NULL;
}

/* ----------------------------------------------------------------------
   free all arrays in a table
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableKokkos<DeviceType>::free_table(Table *tb)
{
  memory->destroy(tb->rfile);
  memory->destroy(tb->efile);
  memory->destroy(tb->ffile);
  memory->destroy(tb->e2file);
  memory->destroy(tb->f2file);

  memory->destroy(tb->rsq);
  memory->destroy(tb->drsq);
  memory->destroy(tb->e);
  memory->destroy(tb->de);
  memory->destroy(tb->f);
  memory->destroy(tb->df);
  memory->destroy(tb->e2);
  memory->destroy(tb->f2);
}

/* ----------------------------------------------------------------------
   spline and splint routines modified from Numerical Recipes
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableKokkos<DeviceType>::spline(double *x, double *y, int n,
                       double yp1, double ypn, double *y2)
{
  int i,k;
  double p,qn,sig,un;
  double *u = new double[n];

  if (yp1 > 0.99e30) y2[0] = u[0] = 0.0;
  else {
    y2[0] = -0.5;
    u[0] = (3.0/(x[1]-x[0])) * ((y[1]-y[0]) / (x[1]-x[0]) - yp1);
  }
  for (i = 1; i < n-1; i++) {
    sig = (x[i]-x[i-1]) / (x[i+1]-x[i-1]);
    p = sig*y2[i-1] + 2.0;
    y2[i] = (sig-1.0) / p;
    u[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]) - (y[i]-y[i-1]) / (x[i]-x[i-1]);
    u[i] = (6.0*u[i] / (x[i+1]-x[i-1]) - sig*u[i-1]) / p;
  }
  if (ypn > 0.99e30) qn = un = 0.0;
  else {
    qn = 0.5;
    un = (3.0/(x[n-1]-x[n-2])) * (ypn - (y[n-1]-y[n-2]) / (x[n-1]-x[n-2]));
  }
  y2[n-1] = (un-qn*u[n-2]) / (qn*y2[n-2] + 1.0);
  for (k = n-2; k >= 0; k--) y2[k] = y2[k]*y2[k+1] + u[k];

  delete [] u;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double PairTableKokkos<DeviceType>::splint(double *xa, double *ya, double *y2a, int n, double x)
{
  int klo,khi,k;
  double h,b,a,y;

  klo = 0;
  khi = n-1;
  while (khi-klo > 1) {
    k = (khi+klo) >> 1;
    if (xa[k] > x) khi = k;
    else klo = k;
  }
  h = xa[khi]-xa[klo];
  a = (xa[khi]-x) / h;
  b = (x-xa[klo]) / h;
  y = a*ya[klo] + b*ya[khi] +
    ((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi]) * (h*h)/6.0;
  return y;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableKokkos<DeviceType>::write_restart(FILE *fp)
{
  write_restart_settings(fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableKokkos<DeviceType>::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableKokkos<DeviceType>::write_restart_settings(FILE *fp)
{
  fwrite(&tabstyle,sizeof(int),1,fp);
  fwrite(&tablength,sizeof(int),1,fp);
  fwrite(&ewaldflag,sizeof(int),1,fp);
  fwrite(&pppmflag,sizeof(int),1,fp);
  fwrite(&msmflag,sizeof(int),1,fp);
  fwrite(&dispersionflag,sizeof(int),1,fp);
  fwrite(&tip4pflag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableKokkos<DeviceType>::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&tabstyle,sizeof(int),1,fp);
    fread(&tablength,sizeof(int),1,fp);
    fread(&ewaldflag,sizeof(int),1,fp);
    fread(&pppmflag,sizeof(int),1,fp);
    fread(&msmflag,sizeof(int),1,fp);
    fread(&dispersionflag,sizeof(int),1,fp);
    fread(&tip4pflag,sizeof(int),1,fp);
  }
  MPI_Bcast(&tabstyle,1,MPI_INT,0,world);
  MPI_Bcast(&tablength,1,MPI_INT,0,world);
  MPI_Bcast(&ewaldflag,1,MPI_INT,0,world);
  MPI_Bcast(&pppmflag,1,MPI_INT,0,world);
  MPI_Bcast(&msmflag,1,MPI_INT,0,world);
  MPI_Bcast(&dispersionflag,1,MPI_INT,0,world);
  MPI_Bcast(&tip4pflag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double PairTableKokkos<DeviceType>::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  int itable;
  double fraction,value,a,b,phi;
  int tlm1 = tablength - 1;

  Table *tb = &tables[tabindex[itype][jtype]];
  if (rsq < tb->innersq) error->one(FLERR,"Pair distance < table inner cutoff");

  if (tabstyle == LOOKUP) {
    itable = static_cast<int> ((rsq-tb->innersq) * tb->invdelta);
    if (itable >= tlm1) error->one(FLERR,"Pair distance > table outer cutoff");
    fforce = factor_lj * tb->f[itable];
  } else if (tabstyle == LINEAR) {
    itable = static_cast<int> ((rsq-tb->innersq) * tb->invdelta);
    if (itable >= tlm1) error->one(FLERR,"Pair distance > table outer cutoff");
    fraction = (rsq - tb->rsq[itable]) * tb->invdelta;
    value = tb->f[itable] + fraction*tb->df[itable];
    fforce = factor_lj * value;
  } else if (tabstyle == SPLINE) {
    itable = static_cast<int> ((rsq-tb->innersq) * tb->invdelta);
    if (itable >= tlm1) error->one(FLERR,"Pair distance > table outer cutoff");
    b = (rsq - tb->rsq[itable]) * tb->invdelta;
    a = 1.0 - b;
    value = a * tb->f[itable] + b * tb->f[itable+1] +
      ((a*a*a-a)*tb->f2[itable] + (b*b*b-b)*tb->f2[itable+1]) *
      tb->deltasq6;
    fforce = factor_lj * value;
  } else {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    itable = rsq_lookup.i & tb->nmask;
    itable >>= tb->nshiftbits;
    fraction = (rsq_lookup.f - tb->rsq[itable]) * tb->drsq[itable];
    value = tb->f[itable] + fraction*tb->df[itable];
    fforce = factor_lj * value;
  }

  if (tabstyle == LOOKUP)
    phi = tb->e[itable];
  else if (tabstyle == LINEAR || tabstyle == BITMAP)
    phi = tb->e[itable] + fraction*tb->de[itable];
  else
    phi = a * tb->e[itable] + b * tb->e[itable+1] +
      ((a*a*a-a)*tb->e2[itable] + (b*b*b-b)*tb->e2[itable+1]) * tb->deltasq6;
  return factor_lj*phi;
}

/* ----------------------------------------------------------------------
   return the Coulomb cutoff for tabled potentials
   called by KSpace solvers which require that all pairwise cutoffs be the same
   loop over all tables not just those indexed by tabindex[i][j] since
     no way to know which tables are active since pair::init() not yet called
------------------------------------------------------------------------- */

template<class DeviceType>
void *PairTableKokkos<DeviceType>::extract(const char *str, int &dim)
{
  if (strcmp(str,"cut_coul") != 0) return NULL;
  if (ntables == 0) error->all(FLERR,"All pair coeffs are not set");

  double cut_coul = tables[0].cut;
  for (int m = 1; m < ntables; m++)
    if (tables[m].cut != cut_coul)
      error->all(FLERR,
                 "Pair table cutoffs must all be equal to use with KSpace");
  dim = 0;
  return &tables[0].cut;
}

template<class DeviceType>
void PairTableKokkos<DeviceType>::init_style()
{
  neighbor->request(this,instance_me);
  neighflag = lmp->kokkos->neighflag;
  int irequest = neighbor->nrequest - 1;

  neighbor->requests[irequest]->
    kokkos_host = Kokkos::Impl::is_same<DeviceType,LMPHostType>::value &&
    !Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;
  neighbor->requests[irequest]->
    kokkos_device = Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;

  if (neighflag == FULL) {
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full_cluster = 0;
  } else if (neighflag == HALF || neighflag == HALFTHREAD) {
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 1;
    neighbor->requests[irequest]->full_cluster = 0;
  } else if (neighflag == N2) {
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full_cluster = 0;
  } else if (neighflag == FULLCLUSTER) {
    neighbor->requests[irequest]->full_cluster = 1;
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
  } else {
    error->all(FLERR,"Cannot use chosen neighbor list style with lj/cut/kk");
  }
}

/*
template <class DeviceType> template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairTableKokkos<DeviceType>::
ev_tally(EV_FLOAT &ev, const int &i, const int &j, const F_FLOAT &fpair,
         const F_FLOAT &delx, const F_FLOAT &dely, const F_FLOAT &delz) const
{
  const int EFLAG = eflag;
  const int NEWTON_PAIR = newton_pair;
  const int VFLAG = vflag_either;

  if (EFLAG) {
    if (eflag_atom) {
      E_FLOAT epairhalf = 0.5 * (ev.evdwl + ev.ecoul);
      if (NEWTON_PAIR || i < nlocal) eatom[i] += epairhalf;
      if (NEWTON_PAIR || j < nlocal) eatom[j] += epairhalf;
    }
  }

  if (VFLAG) {
    const E_FLOAT v0 = delx*delx*fpair;
    const E_FLOAT v1 = dely*dely*fpair;
    const E_FLOAT v2 = delz*delz*fpair;
    const E_FLOAT v3 = delx*dely*fpair;
    const E_FLOAT v4 = delx*delz*fpair;
    const E_FLOAT v5 = dely*delz*fpair;

    if (vflag_global) {
      if (NEIGHFLAG) {
        if (NEWTON_PAIR) {
          ev.v[0] += v0;
          ev.v[1] += v1;
          ev.v[2] += v2;
          ev.v[3] += v3;
          ev.v[4] += v4;
          ev.v[5] += v5;
        } else {
          if (i < nlocal) {
            ev.v[0] += 0.5*v0;
            ev.v[1] += 0.5*v1;
            ev.v[2] += 0.5*v2;
            ev.v[3] += 0.5*v3;
            ev.v[4] += 0.5*v4;
            ev.v[5] += 0.5*v5;
          }
          if (j < nlocal) {
            ev.v[0] += 0.5*v0;
            ev.v[1] += 0.5*v1;
            ev.v[2] += 0.5*v2;
            ev.v[3] += 0.5*v3;
            ev.v[4] += 0.5*v4;
            ev.v[5] += 0.5*v5;
          }
        }
      } else {
        ev.v[0] += 0.5*v0;
        ev.v[1] += 0.5*v1;
        ev.v[2] += 0.5*v2;
        ev.v[3] += 0.5*v3;
        ev.v[4] += 0.5*v4;
        ev.v[5] += 0.5*v5;
      }
    }

    if (vflag_atom) {
      if (NEWTON_PAIR || i < nlocal) {
        d_vatom(i,0) += 0.5*v0;
        d_vatom(i,1) += 0.5*v1;
        d_vatom(i,2) += 0.5*v2;
        d_vatom(i,3) += 0.5*v3;
        d_vatom(i,4) += 0.5*v4;
        d_vatom(i,5) += 0.5*v5;
      }
      if (NEWTON_PAIR || (NEIGHFLAG && j < nlocal)) {
        d_vatom(j,0) += 0.5*v0;
        d_vatom(j,1) += 0.5*v1;
        d_vatom(j,2) += 0.5*v2;
        d_vatom(j,3) += 0.5*v3;
        d_vatom(j,4) += 0.5*v4;
        d_vatom(j,5) += 0.5*v5;
      }
    }
  }
}
*/
template<class DeviceType>
void PairTableKokkos<DeviceType>::cleanup_copy() {
  // WHY needed: this prevents parent copy from deallocating any arrays
  allocated = 0;
  cutsq = NULL;
  eatom = NULL;
  vatom = NULL;
  h_table=NULL; d_table=NULL;
}

template class PairTableKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class PairTableKokkos<LMPHostType>;
#endif

