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

/* ----------------------------------------------------------------------------------------
   Contributing authors:
   Stan Moore (Sandia)

   Please cite the related publications:
   J.D. Moore, B.C. Barnes, S. Izvekov, M. Lisal, M.S. Sellers, D.E. Taylor & J.K. Brennan
   "A coarse-grain force field for RDX: Density dependent and energy conserving"
   The Journal of Chemical Physics, 2016, 144, 104501.
------------------------------------------------------------------------------------------- */

#include <mpi.h>
#include <math.h>
#include "math_const.h"
#include <stdlib.h>
#include <string.h>
#include "pair_multi_lucy_rx_kokkos.h"
#include "atom_kokkos.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "citeme.h"
#include "modify.h"
#include "fix.h"
#include "atom_masks.h"
#include "neigh_request.h"

using namespace LAMMPS_NS;

enum{NONE,RLINEAR,RSQ};

#define MAXLINE 1024

#define oneFluidParameter (-1)
#define isOneFluid(_site) ( (_site) == oneFluidParameter )

static const char cite_pair_multi_lucy_rx[] =
  "pair_style multi/lucy/rx command:\n\n"
  "@Article{Moore16,\n"
  " author = {J.D. Moore, B.C. Barnes, S. Izvekov, M. Lisal, M.S. Sellers, D.E. Taylor and J. K. Brennan},\n"
  " title = {A coarse-grain force field for RDX:  Density dependent and energy conserving},\n"
  " journal = {J. Chem. Phys.},\n"
  " year =    2016,\n"
  " volume =  144\n"
  " pages =   {104501}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairMultiLucyRXKokkos<DeviceType>::PairMultiLucyRXKokkos(LAMMPS *lmp) : PairMultiLucyRX(lmp)
{
  respa_enable = 0;

  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  k_error_flag = DAT::tdual_int_scalar("pair:error_flag");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairMultiLucyRXKokkos<DeviceType>::~PairMultiLucyRXKokkos()
{

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMultiLucyRXKokkos<DeviceType>::init_style()
{
  PairMultiLucyRX::init_style();

  // irequest = neigh request made by parent class

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
    neighbor->requests[irequest]->ghost = 1;
  } else if (neighflag == HALF || neighflag == HALFTHREAD) {
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 1;
    neighbor->requests[irequest]->ghost = 1;
  } else {
    error->all(FLERR,"Cannot use chosen neighbor list style with reax/c/kk");
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMultiLucyRXKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  double evdwl,evdwlOld;

  evdwlOld = 0.0;
  evdwl = 0.0;
  if (neighflag == FULL) no_virial_fdotr_compute = 1;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memory->destroy_kokkos(k_eatom,eatom);
    memory->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.d_view;
  }
  if (vflag_atom) {
    memory->destroy_kokkos(k_vatom,vatom);
    memory->create_kokkos(k_vatom,vatom,maxvatom,6,"pair:vatom");
    d_vatom = k_vatom.d_view;
  }

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  uCG = atomKK->k_uCG.view<DeviceType>();
  uCGnew = atomKK->k_uCGnew.view<DeviceType>();
  dvector = atomKK->k_dvector.view<DeviceType>();
  rho = atomKK->k_rho.view<DeviceType>();

  nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int newton_pair = force->newton_pair;

  {
    const int ntotal = nlocal + nghost;
    d_fractionOld1 = typename AT::t_float_1d("PairMultiLucyRX::fractionOld1",ntotal);
    d_fractionOld2 = typename AT::t_float_1d("PairMultiLucyRX::fractionOld2",ntotal);
    d_fraction1 = typename AT::t_float_1d("PairMultiLucyRX::fraction1",ntotal);
    d_fraction2 = typename AT::t_float_1d("PairMultiLucyRX::fraction2",ntotal);

    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXgetParams>(0,ntotal),*this);
  }

  const int inum = list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  computeLocalDensity();

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  if (evflag) {
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXCompute<HALF,1,1> >(0,inum),*this,ev);
  } else {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXCompute<HALF,1,0> >(0,inum),*this);
  }

  k_error_flag.template modify<DeviceType>();
  k_error_flag.template sync<LMPHostType>();
  if (k_error_flag.h_view() == 1)
    error->one(FLERR,"Density < table inner cutoff");
  else if (k_error_flag.h_view() == 2)
    error->one(FLERR,"Density > table outer cutoff");
  else if (k_error_flag.h_view() == 3)
    error->one(FLERR,"Only LOOKUP and LINEAR table styles have been implemented for pair multi/lucy/rx");

  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

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
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<DeviceType>::operator()(TagPairMultiLucyRXgetParams, const int &i) const {
  getParams(i, d_fractionOld1[i], d_fractionOld2[i], d_fraction1[i], d_fraction2[i]);
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<DeviceType>::operator()(TagPairMultiLucyRXCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int &ii, EV_FLOAT& ev) const {
  int i,j,jj,inum,jnum,itype,jtype,itable;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,evdwlOld,fpair;
  double rsq;

  double fractionOld1_i,fractionOld1_j;
  double fractionOld2_i,fractionOld2_j;
  double fraction1_i;

  double pi = MathConst::MY_PI;
  double A_i, A_j;
  double fraction_i,fraction_j;
  int jtable;

  Table *tb;

  int tlm1 = tablength - 1;

  i = d_ilist[ii];
  xtmp = x(i,0);
  ytmp = x(i,1);
  ztmp = x(i,2);
  itype = type[i];
  jnum = d_numneigh[i];

  double fx_i = 0.0;
  double fy_i = 0.0;
  double fz_i = 0.0;

  fractionOld1_i = d_fractionOld1[i];
  fractionOld2_i = d_fractionOld2[i];
  fraction1_i = d_fraction1[i];

  for (jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;

    delx = xtmp - x(j,0);
    dely = ytmp - x(j,1);
    delz = ztmp - x(j,2);
    rsq = delx*delx + dely*dely + delz*delz;
    jtype = type[j];

    if (rsq < d_cutsq(itype,jtype)) { // optimize
      fpair = 0.0;

      fractionOld1_j = d_fractionOld1[j];
      fractionOld2_j = d_fractionOld2[j];

      tb = &tables[tabindex[itype][jtype]];
      if (rho[i]*rho[i] < tb->innersq || rho[j]*rho[j] < tb->innersq){
        k_error_flag.d_view() = 1;
      }
      if (tabstyle == LOOKUP) {
        itable = static_cast<int> (((rho[i]*rho[i]) - tb->innersq) * tb->invdelta);
        jtable = static_cast<int> (((rho[j]*rho[j]) - tb->innersq) * tb->invdelta);
        if (itable >= tlm1 || jtable >= tlm1){
          k_error_flag.d_view() = 2;
        }
        A_i = tb->f[itable];
        A_j = tb->f[jtable];

        const double rfactor = 1.0-sqrt(rsq/d_cutsq(itype,jtype));
        fpair = 0.5*(A_i + A_j)*(4.0-3.0*rfactor)*rfactor*rfactor*rfactor;
        fpair /= sqrt(rsq);

      } else if (tabstyle == LINEAR) {
        itable = static_cast<int> ((rho[i]*rho[i] - tb->innersq) * tb->invdelta);
        jtable = static_cast<int> (((rho[j]*rho[j]) - tb->innersq) * tb->invdelta);
        if (itable >= tlm1 || jtable >= tlm1){
          k_error_flag.d_view() = 2;
        }
        if(itable<0) itable=0;
        if(itable>=tlm1) itable=tlm1;
        if(jtable<0) jtable=0;
        if(jtable>=tlm1)jtable=tlm1;

        fraction_i = (((rho[i]*rho[i]) - tb->rsq[itable]) * tb->invdelta);
        fraction_j = (((rho[j]*rho[j]) - tb->rsq[jtable]) * tb->invdelta);
        if(itable==0) fraction_i=0.0;
        if(itable==tlm1) fraction_i=0.0;
        if(jtable==0) fraction_j=0.0;
        if(jtable==tlm1) fraction_j=0.0;

        A_i = tb->f[itable] + fraction_i*tb->df[itable];
        A_j = tb->f[jtable] + fraction_j*tb->df[jtable];

        const double rfactor = 1.0-sqrt(rsq/d_cutsq(itype,jtype));
        fpair = 0.5*(A_i + A_j)*(4.0-3.0*rfactor)*rfactor*rfactor*rfactor;
        fpair /= sqrt(rsq);

      } else k_error_flag.d_view() = 3;

      if (isite1 == isite2) fpair = sqrt(fractionOld1_i*fractionOld2_j)*fpair;
      else fpair = (sqrt(fractionOld1_i*fractionOld2_j) + sqrt(fractionOld2_i*fractionOld1_j))*fpair;

      fx_i += delx*fpair;
      fy_i += dely*fpair;
      fz_i += delz*fpair;
      if (NEWTON_PAIR || j < nlocal) {
        f(j,0) -= delx*fpair;
        f(j,1) -= dely*fpair;
        f(j,2) -= delz*fpair;
      }
      //if (evflag) ev_tally(i,j,nlocal,newton_pair,0.0,0.0,fpair,delx,dely,delz);
      if (EVFLAG) this->template ev_tally<NEIGHFLAG,NEWTON_PAIR>(ev,i,j,0.0,fpair,delx,dely,delz);
    }
  }

  f(i,0) += fx_i;
  f(i,1) += fy_i;
  f(i,2) += fz_i;

  tb = &tables[tabindex[itype][itype]];
  itable = static_cast<int> (((rho[i]*rho[i]) - tb->innersq) * tb->invdelta);
  if (tabstyle == LOOKUP) evdwl = tb->e[itable];
  else if (tabstyle == LINEAR){
    if (itable >= tlm1){
      k_error_flag.d_view() = 2;
    }
    if(itable==0) fraction_i=0.0;
    else fraction_i = (((rho[i]*rho[i]) - tb->rsq[itable]) * tb->invdelta);
    evdwl = tb->e[itable] + fraction_i*tb->de[itable];
  } else k_error_flag.d_view() = 3;

  evdwl *=(pi*d_cutsq(itype,itype)*d_cutsq(itype,itype))/84.0;
  evdwlOld = fractionOld1_i*evdwl;
  evdwl = fraction1_i*evdwl;

  uCG[i] += evdwlOld;
  uCGnew[i] += evdwl;

  evdwl = evdwlOld;

  //if (evflag) ev_tally(0,0,nlocal,newton_pair,evdwl,0.0,0.0,0.0,0.0,0.0);
  if (EVFLAG) ev.evdwl += evdwl;
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<DeviceType>::operator()(TagPairMultiLucyRXCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(TagPairMultiLucyRXCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(), ii, ev);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType>
void PairMultiLucyRXKokkos<DeviceType>::coeff(int narg, char **arg)
{
  if (narg != 6 && narg != 7) error->all(FLERR,"Illegal pair_coeff command");

  bool rx_flag = false;
  for (int i = 0; i < modify->nfix; i++)
    if (strncmp(modify->fix[i]->style,"rx",2) == 0) rx_flag = true;
  if (!rx_flag) error->all(FLERR,"PairMultiLucyRXKokkos<DeviceType> requires a fix rx command.");

  if (!allocated) allocate();

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
  int n;
  n = strlen(arg[3]) + 1;
  site1 = new char[n];
  strcpy(site1,arg[4]);

  n = strlen(arg[4]) + 1;
  site2 = new char[n];
  strcpy(site2,arg[5]);

  // set table cutoff

  if (narg == 7) tb->cut = force->numeric(FLERR,arg[6]);
  else if (tb->rflag) tb->cut = tb->rhi;
  else tb->cut = tb->rfile[tb->ninput-1];

  // error check on table parameters
  // insure cutoff is within table

  if (tb->ninput <= 1) error->one(FLERR,"Invalid pair table length");
  if (tb->rflag == 0) {
    rho_0 = tb->rfile[0];
  } else {
    rho_0 = tb->rlo;
  }

  tb->match = 0;
  if (tabstyle == LINEAR && tb->ninput == tablength &&
      tb->rflag == RSQ) tb->match = 1;

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

  // Match site* to isite values.

  if (strcmp(site1, "1fluid") == 0)
     isite1 = oneFluidParameter;
  else {
     isite1 = nspecies;
     for (int ispecies = 0; ispecies < nspecies; ++ispecies)
        if (strcmp(site1, atom->dname[ispecies]) == 0){
           isite1 = ispecies;
           break;
        }

     if (isite1 == nspecies)
        error->all(FLERR,"Pair_multi_lucy_rx site1 is invalid.");
  }

  if (strcmp(site2, "1fluid") == 0)
     isite2 = oneFluidParameter;
  else {
     isite2 = nspecies;
     for (int ispecies = 0; ispecies < nspecies; ++ispecies)
        if (strcmp(site2, atom->dname[ispecies]) == 0){
           isite2 = ispecies;
           break;
        }

     if (isite2 == nspecies)
        error->all(FLERR,"Pair_multi_lucy_rx site2 is invalid.");
  }

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMultiLucyRXKokkos<DeviceType>::computeLocalDensity()
{
  x = atomKK->k_x.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  rho = atomKK->k_rho.view<DeviceType>();
  nlocal = atom->nlocal;

  //sync

  const int inum = list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  const double pi = MathConst::MY_PI;

  const bool newton_pair = force->newton_pair;
  one_type = (atom->ntypes == 1);

  // Special cut-off values for when there's only one type.
  cutsq_type11 = cutsq[1][1];
  rcut_type11 = sqrt(cutsq_type11);
  factor_type11 = 84.0/(5.0*pi*rcut_type11*rcut_type11*rcut_type11);

  // zero out density
  int m = nlocal;
  if (newton_pair) m += atom->nghost;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXZero>(0,m),*this);

// rho = density at each atom
// loop over neighbors of my atoms
  if (newton_pair)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXComputeLocalDensity<HALF,1> >(0,inum),*this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXComputeLocalDensity<HALF,0> >(0,inum),*this);

  if (newton_pair) comm->reverse_comm_pair(this);

  comm->forward_comm_pair(this);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<DeviceType>::operator()(TagPairMultiLucyRXZero, const int &i) const {
  rho[i] = 0.0;
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<DeviceType>::operator()(TagPairMultiLucyRXComputeLocalDensity<NEIGHFLAG,NEWTON_PAIR>, const int &ii) const {
  const int i = d_ilist[ii];

  const double xtmp = x(i,0);
  const double ytmp = x(i,1);
  const double ztmp = x(i,2);

  double rho_i = rho[i];

  const int itype = type[i];
  const int jnum = d_numneigh[i];

  const double pi = MathConst::MY_PI;

  for (int jj = 0; jj < jnum; jj++){
    const int j = (d_neighbors(i,jj) & NEIGHMASK);
    const int jtype = type[j];

    const double delx = xtmp - x(j,0);
    const double dely = ytmp - x(j,1);
    const double delz = ztmp - x(j,2);
    const double rsq = delx*delx + dely*dely + delz*delz;

    if (one_type) {
      if (rsq < cutsq_type11) {
        const double rcut = rcut_type11;
        const double r_over_rcut = sqrt(rsq) / rcut;
        const double tmpFactor = 1.0 - r_over_rcut;
        const double tmpFactor4 = tmpFactor*tmpFactor*tmpFactor*tmpFactor;
        const double factor = factor_type11*(1.0 + 1.5*r_over_rcut)*tmpFactor4;
        rho_i += factor;
        if (NEWTON_PAIR || j < nlocal)
          rho[j] += factor;
      } else if (rsq < d_cutsq(itype,jtype)) {
        const double rcut = sqrt(d_cutsq(itype,jtype));
        const double tmpFactor = 1.0-sqrt(rsq)/rcut;
        const double tmpFactor4 = tmpFactor*tmpFactor*tmpFactor*tmpFactor;
        const double factor = (84.0/(5.0*pi*rcut*rcut*rcut))*(1.0+3.0*sqrt(rsq)/(2.0*rcut))*tmpFactor4;
        rho_i += factor;
        if (NEWTON_PAIR || j < nlocal)
          rho[j] += factor;
      }
    }
  }

  rho[i] = rho_i;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<DeviceType>::getParams(int id, double &fractionOld1, double &fractionOld2, double &fraction1, double &fraction2) const
{
  double fractionOld, fraction;
  double nTotal, nTotalOld;

  nTotal = 0.0;
  nTotalOld = 0.0;
  for (int ispecies = 0; ispecies < nspecies; ispecies++){
    nTotal += dvector(ispecies,id);
    nTotalOld += dvector(ispecies+nspecies,id);
  }

  if (isOneFluid(isite1) == false){
    fractionOld1 = dvector(isite1+nspecies,id)/nTotalOld;
    fraction1 = dvector(isite1,id)/nTotal;
  }
  if (isOneFluid(isite2) == false){
    fractionOld2 = dvector(isite2+nspecies,id)/nTotalOld;
    fraction2 = dvector(isite2,id)/nTotal;
  }

  if (isOneFluid(isite1) || isOneFluid(isite2)){
    fractionOld  = 0.0;
    fraction  = 0.0;

    for (int ispecies = 0; ispecies < nspecies; ispecies++){
      if (isite1 == ispecies || isite2 == ispecies) continue;
      fractionOld += dvector(ispecies+nspecies,id) / nTotalOld;
      fraction += dvector(ispecies,id) / nTotal;
    }
    if (isOneFluid(isite1)){
      fractionOld1 = fractionOld;
      fraction1 = fraction;
    }
    if (isOneFluid(isite2)){
      fractionOld2 = fractionOld;
      fraction2 = fraction;
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairMultiLucyRXKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;
  rho = atomKK->k_rho.view<DeviceType>();

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = rho[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMultiLucyRXKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;
  rho = atomKK->k_rho.view<DeviceType>();

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) rho[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairMultiLucyRXKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  rho = atomKK->k_rho.view<DeviceType>();

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = rho[i];
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMultiLucyRXKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  rho = atomKK->k_rho.view<DeviceType>();

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    rho[j] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                const F_FLOAT &dely, const F_FLOAT &delz) const
{
  const int EFLAG = eflag;
  const int VFLAG = vflag_either;

  // The eatom and vatom arrays are atomic for Half/Thread neighbor style
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_eatom = k_eatom.view<DeviceType>();
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom = k_vatom.view<DeviceType>();

  if (EFLAG) {
    if (eflag_atom) {
      const E_FLOAT epairhalf = 0.5 * epair;
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) v_eatom[i] += epairhalf;
        if (NEWTON_PAIR || j < nlocal) v_eatom[j] += epairhalf;
      } else {
        v_eatom[i] += epairhalf;
      }
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
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) {
          ev.v[0] += 0.5*v0;
          ev.v[1] += 0.5*v1;
          ev.v[2] += 0.5*v2;
          ev.v[3] += 0.5*v3;
          ev.v[4] += 0.5*v4;
          ev.v[5] += 0.5*v5;
        }
        if (NEWTON_PAIR || j < nlocal) {
        ev.v[0] += 0.5*v0;
        ev.v[1] += 0.5*v1;
        ev.v[2] += 0.5*v2;
        ev.v[3] += 0.5*v3;
        ev.v[4] += 0.5*v4;
        ev.v[5] += 0.5*v5;
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
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) {
          v_vatom(i,0) += 0.5*v0;
          v_vatom(i,1) += 0.5*v1;
          v_vatom(i,2) += 0.5*v2;
          v_vatom(i,3) += 0.5*v3;
          v_vatom(i,4) += 0.5*v4;
          v_vatom(i,5) += 0.5*v5;
        }
        if (NEWTON_PAIR || j < nlocal) {
        v_vatom(j,0) += 0.5*v0;
        v_vatom(j,1) += 0.5*v1;
        v_vatom(j,2) += 0.5*v2;
        v_vatom(j,3) += 0.5*v3;
        v_vatom(j,4) += 0.5*v4;
        v_vatom(j,5) += 0.5*v5;
        }
      } else {
        v_vatom(i,0) += 0.5*v0;
        v_vatom(i,1) += 0.5*v1;
        v_vatom(i,2) += 0.5*v2;
        v_vatom(i,3) += 0.5*v3;
        v_vatom(i,4) += 0.5*v4;
        v_vatom(i,5) += 0.5*v5;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairMultiLucyRXKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class PairMultiLucyRXKokkos<LMPHostType>;
#endif
}
