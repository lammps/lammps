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
   Contributing author: Dan Ibanez (SNL)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "pair_table_rx_kokkos.h"
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
#include "fix.h"

using namespace LAMMPS_NS;

enum{NONE,RLINEAR,RSQ,BMP};

#ifdef DBL_EPSILON
  #define MY_EPSILON (10.0*DBL_EPSILON)
#else
  #define MY_EPSILON (10.0*2.220446049250313e-16)
#endif

#define OneFluidValue (-1)
#define isOneFluid(_site_) ( (_site_) == OneFluidValue )

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairTableRXKokkos<DeviceType>::PairTableRXKokkos(LAMMPS *lmp) : PairTable(lmp)
{
  update_table = 0;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
  h_table = new TableHost();
  d_table = new TableDevice();
  fractionalWeighting = true;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairTableRXKokkos<DeviceType>::~PairTableRXKokkos()
{
  if (copymode) return;
  delete h_table;
  h_table = nullptr;
  delete d_table;
  d_table = nullptr;
  copymode = true; //prevents base class destructor from running
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairTableRXKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
void PairTableRXKokkos<DeviceType>::compute_style(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

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
      PairComputeFunctor<PairTableRXKokkos<DeviceType>,FULL,false,S_TableCompute<DeviceType,TABSTYLE> >
        ff(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,ff,ev);
      else Kokkos::parallel_for(list->inum,ff);
    } else if (neighflag == HALFTHREAD) {
      PairComputeFunctor<PairTableRXKokkos<DeviceType>,HALFTHREAD,false,S_TableCompute<DeviceType,TABSTYLE> >
        ff(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,ff,ev);
      else Kokkos::parallel_for(list->inum,ff);
    } else if (neighflag == HALF) {
      PairComputeFunctor<PairTableRXKokkos<DeviceType>,HALF,false,S_TableCompute<DeviceType,TABSTYLE> >
        f(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,f,ev);
      else Kokkos::parallel_for(list->inum,f);
    } else if (neighflag == N2) {
      PairComputeFunctor<PairTableRXKokkos<DeviceType>,N2,false,S_TableCompute<DeviceType,TABSTYLE> >
        f(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(nlocal,f,ev);
      else Kokkos::parallel_for(nlocal,f);
    }
  } else {
    if (neighflag == FULL) {
      PairComputeFunctor<PairTableRXKokkos<DeviceType>,FULL,true,S_TableCompute<DeviceType,TABSTYLE> >
        f(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,f,ev);
      else Kokkos::parallel_for(list->inum,f);
    } else if (neighflag == HALFTHREAD) {
      PairComputeFunctor<PairTableRXKokkos<DeviceType>,HALFTHREAD,true,S_TableCompute<DeviceType,TABSTYLE> >
        f(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,f,ev);
      else Kokkos::parallel_for(list->inum,f);
    } else if (neighflag == HALF) {
      PairComputeFunctor<PairTableRXKokkos<DeviceType>,HALF,true,S_TableCompute<DeviceType,TABSTYLE> >
        f(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(list->inum,f,ev);
      else Kokkos::parallel_for(list->inum,f);
    } else if (neighflag == N2) {
      PairComputeFunctor<PairTableRXKokkos<DeviceType>,N2,true,S_TableCompute<DeviceType,TABSTYLE> >
        f(this,(NeighListKokkos<DeviceType>*) list);
      if (eflag || vflag) Kokkos::parallel_reduce(nlocal,f,ev);
      else Kokkos::parallel_for(nlocal,f);
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

  if (vflag_fdotr) pair_virial_fdotr_compute(this);
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairTableRXKokkos<DeviceType>::
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
F_FLOAT PairTableRXKokkos<DeviceType>::
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
void PairTableRXKokkos<DeviceType>::create_kokkos_tables()
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
  d_table_const.nshiftbits = d_table->nshiftbits;
  Kokkos::deep_copy(d_table->nmask,h_table->nmask);
  d_table_const.nmask = d_table->nmask;
  Kokkos::deep_copy(d_table->innersq,h_table->innersq);
  d_table_const.innersq = d_table->innersq;
  Kokkos::deep_copy(d_table->invdelta,h_table->invdelta);
  d_table_const.invdelta = d_table->invdelta;
  Kokkos::deep_copy(d_table->deltasq6,h_table->deltasq6);
  d_table_const.deltasq6 = d_table->deltasq6;

  if(tabstyle == LOOKUP) {
    Kokkos::deep_copy(d_table->e,h_table->e);
    d_table_const.e = d_table->e;
    Kokkos::deep_copy(d_table->f,h_table->f);
    d_table_const.f = d_table->f;
  }

  if(tabstyle == LINEAR) {
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

  if(tabstyle == SPLINE) {
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

  if(tabstyle == BITMAP) {
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
void PairTableRXKokkos<DeviceType>::allocate()
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
void PairTableRXKokkos<DeviceType>::settings(int narg, char **arg)
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
    else if (strcmp(arg[iarg],"fractional") == 0) fractionalWeighting = true;
    else if (strcmp(arg[iarg],"molecular") == 0) fractionalWeighting = false;
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
void PairTableRXKokkos<DeviceType>::coeff(int narg, char **arg)
{
  if (narg != 6 && narg != 7) error->all(FLERR,"Illegal pair_coeff command");
  if (!allocated) allocate();

  bool rx_flag = false;
  for (int i = 0; i < modify->nfix; i++)
    if (strncmp(modify->fix[i]->style,"rx",2) == 0) rx_flag = true;
  if (!rx_flag) error->all(FLERR,"PairTableRX requires a fix rx command.");

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  int me;
  MPI_Comm_rank(world,&me);
  tables = (Table *)
    memory->srealloc(tables,(ntables+1)*sizeof(Table),"pair:tables");
  Table *tb = &tables[ntables];
  null_table(tb);
  if (me == 0) read_table(tb,arg[2],arg[3]);
  bcast_table(tb);

  nspecies = atom->nspecies_dpd;
  if(nspecies==0) error->all(FLERR,"There are no rx species specified.");
  int n;
  n = strlen(arg[3]) + 1;
  site1 = new char[n];
  strcpy(site1,arg[4]);

  int ispecies;
  for (ispecies = 0; ispecies < nspecies; ispecies++){
    if (strcmp(site1,&atom->dname[ispecies][0]) == 0) break;
  }
  if (ispecies == nspecies && strcmp(site1,"1fluid") != 0)
    error->all(FLERR,"Site1 name not recognized in pair coefficients");

  n = strlen(arg[4]) + 1;
  site2 = new char[n];
  strcpy(site2,arg[5]);

  for (ispecies = 0; ispecies < nspecies; ispecies++){
    if (strcmp(site2,&atom->dname[ispecies][0]) == 0) break;
  }
  if (ispecies == nspecies && strcmp(site2,"1fluid") != 0)
    error->all(FLERR,"Site2 name not recognized in pair coefficients");

  // set table cutoff

  if (narg == 7) tb->cut = force->numeric(FLERR,arg[6]);
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

  {
     if ( strcmp(site1,"1fluid") == 0 )
       isite1 = OneFluidValue;
     else {
       isite1 = nspecies;

       for (int k = 0; k < nspecies; k++){
         if (strcmp(site1, atom->dname[k]) == 0){
           isite1 = k;
           break;
         }
       }

       if (isite1 == nspecies) error->all(FLERR,"isite1 == nspecies");
     }

     if ( strcmp(site2,"1fluid") == 0 )
       isite2 = OneFluidValue;
     else {
       isite2 = nspecies;

       for (int k = 0; k < nspecies; k++){
         if (strcmp(site2, atom->dname[k]) == 0){
           isite2 = ispecies;
           break;
         }
       }

       if (isite2 == nspecies)
         error->all(FLERR,"isite2 == nspecies");
     }
  }

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairTableRXKokkos<DeviceType>::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  tabindex[j][i] = tabindex[i][j];

  if(i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_cutsq[j][i] = m_cutsq[i][j] = tables[tabindex[i][j]].cut*tables[tabindex[i][j]].cut;
  }

  return tables[tabindex[i][j]].cut;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double PairTableRXKokkos<DeviceType>::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  int itable;
  double fraction,value,a,b,phi;
  int tlm1 = tablength - 1;

  Table *tb = &tables[tabindex[itype][jtype]];
  double mixWtSite1_i, mixWtSite1_j;
  double mixWtSite2_i, mixWtSite2_j;
  double mixWtSite1old_i, mixWtSite1old_j;
  double mixWtSite2old_i, mixWtSite2old_j;

  fraction = 0.0;
  a = 0.0;
  b = 0.0;

  typename ArrayTypes<LMPHostType>::t_float_2d_randomread h_dvector =
    atomKK->k_dvector.view<LMPHostType>();
  getMixingWeights<LMPHostType>(h_dvector,i,mixWtSite1old_i,mixWtSite2old_i,
      mixWtSite1_i,mixWtSite2_i);
  getMixingWeights<LMPHostType>(h_dvector,j,mixWtSite1old_j,mixWtSite2old_j,
      mixWtSite1_j,mixWtSite2_j);

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

  if (isite1 == isite2) fforce = sqrt(mixWtSite1_i*mixWtSite2_j)*fforce;
  else fforce = (sqrt(mixWtSite1_i*mixWtSite2_j) + sqrt(mixWtSite2_i*mixWtSite1_j))*fforce;

  if (tabstyle == LOOKUP)
    phi = tb->e[itable];
  else if (tabstyle == LINEAR || tabstyle == BITMAP)
    phi = tb->e[itable] + fraction*tb->de[itable];
  else
    phi = a * tb->e[itable] + b * tb->e[itable+1] +
      ((a*a*a-a)*tb->e2[itable] + (b*b*b-b)*tb->e2[itable+1]) * tb->deltasq6;

  if (isite1 == isite2) phi = sqrt(mixWtSite1_i*mixWtSite2_j)*phi;
  else phi = (sqrt(mixWtSite1_i*mixWtSite2_j) + sqrt(mixWtSite2_i*mixWtSite1_j))*phi;

  return factor_lj*phi;
}

/* ----------------------------------------------------------------------
   compute r,e,f vectors from splined values
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableRXKokkos<DeviceType>::compute_table(Table *tb)
{
  update_table = 1;
  PairTable::compute_table(tb);
}

template<class DeviceType>
void PairTableRXKokkos<DeviceType>::init_style()
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
  } else if (neighflag == HALF || neighflag == HALFTHREAD) {
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 1;
  } else if (neighflag == N2) {
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 0;
  } else {
    error->all(FLERR,"Cannot use chosen neighbor list style with lj/cut/kk");
  }
}

template<class DeviceType>
void PairTableRXKokkos<DeviceType>::cleanup_copy() {
  // WHY needed: this prevents parent copy from deallocating any arrays
  allocated = 0;
  cutsq = NULL;
  eatom = NULL;
  vatom = NULL;
  h_table=NULL; d_table=NULL;
}

template<class DeviceType>
template<class ExecDevice>
KOKKOS_INLINE_FUNCTION
void PairTableRXKokkos<DeviceType>::getMixingWeights(
    typename ArrayTypes<ExecDevice>::t_float_2d_randomread dvector,
    int id,
    double &mixWtSite1old, double &mixWtSite2old,
    double &mixWtSite1, double &mixWtSite2) {
  double fractionOFAold, fractionOFA;
  double fractionOld1, fraction1;
  double fractionOld2, fraction2;
  double nMoleculesOFAold, nMoleculesOFA;
  double nMoleculesOld1, nMolecules1;
  double nMoleculesOld2, nMolecules2;
  double nTotal, nTotalOld;

  nTotal = 0.0;
  nTotalOld = 0.0;
  for (int ispecies = 0; ispecies < nspecies; ++ispecies){
    nTotal += dvector(ispecies,id);
    nTotalOld += dvector(ispecies+nspecies,id);
  }
  if(nTotal < MY_EPSILON || nTotalOld < MY_EPSILON)
    error->all(FLERR,"The number of molecules in CG particle is less than 10*DBL_EPSILON.");

  if (isOneFluid(isite1) == false){
    nMoleculesOld1 = dvector(isite1+nspecies,id);
    nMolecules1 = dvector(isite1,id);
    fractionOld1 = nMoleculesOld1/nTotalOld;
    fraction1 = nMolecules1/nTotal;
  }
  if (isOneFluid(isite2) == false){
    nMoleculesOld2 = dvector(isite2+nspecies,id);
    nMolecules2 = dvector(isite2,id);
    fractionOld2 = nMoleculesOld2/nTotalOld;
    fraction2 = nMolecules2/nTotal;
  }

  if (isOneFluid(isite1) || isOneFluid(isite2)){
    nMoleculesOFAold  = 0.0;
    nMoleculesOFA  = 0.0;
    fractionOFAold  = 0.0;
    fractionOFA  = 0.0;

    for (int ispecies = 0; ispecies < nspecies; ispecies++){
      if (isite1 == ispecies || isite2 == ispecies) continue;
      nMoleculesOFAold += dvector(ispecies+nspecies,id);
      nMoleculesOFA += dvector(ispecies,id);
      fractionOFAold += dvector(ispecies+nspecies,id)/nTotalOld;
      fractionOFA += dvector(ispecies,id)/nTotal;
    }
    if(isOneFluid(isite1)){
      nMoleculesOld1 = 1.0-(nTotalOld-nMoleculesOFAold);
      nMolecules1 = 1.0-(nTotal-nMoleculesOFA);
      fractionOld1 = fractionOFAold;
      fraction1 = fractionOFA;
    }
    if(isOneFluid(isite2)){
      nMoleculesOld2 = 1.0-(nTotalOld-nMoleculesOFAold);
      nMolecules2 = 1.0-(nTotal-nMoleculesOFA);
      fractionOld2 = fractionOFAold;
      fraction2 = fractionOFA;
    }
  }

  if(fractionalWeighting){
    mixWtSite1old = fractionOld1;
    mixWtSite1 = fraction1;
    mixWtSite2old = fractionOld2;
    mixWtSite2 = fraction2;
  } else {
    mixWtSite1old = nMoleculesOld1;
    mixWtSite1 = nMolecules1;
    mixWtSite2old = nMoleculesOld2;
    mixWtSite2 = nMolecules2;
  }
}

namespace LAMMPS_NS {
template class PairTableRXKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class PairTableRXKokkos<LMPHostType>;
#endif

}

