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
   Contributing author: Greg Wagner (SNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
//KK*
#include "meam_kokkos.h"
#include "kokkos.h"
#include "pair_kokkos.h"
//#include "pair_meamc.h"

#include "pair_meam_kokkos.h"
#include "atom_kokkos.h"
//*KK
#include "force.h"
#include "comm.h"
//KK*
//#include "memory.h"
#include "memory_kokkos.h"
//*KK
#include "neighbor.h"
//KK*
//#include "neigh_list.h"
#include "neigh_list_kokkos.h"
//*KK
#include "neigh_request.h"
#include "error.h"
//*KK
#include "atom_masks.h"
//*KK

using namespace LAMMPS_NS;

#if 0
static const int nkeywords = 21;
static const char *keywords[] = {
  "Ec","alpha","rho0","delta","lattce",
  "attrac","repuls","nn2","Cmin","Cmax","rc","delr",
  "augt1","gsmooth_factor","re","ialloy",
  "mixture_ref_t","erose_form","zbl",
  "emb_lin_neg","bkgd_dyn"};
#endif

/* ---------------------------------------------------------------------- */
template<class DeviceType>
PairMEAMKokkos<DeviceType>::PairMEAMKokkos(LAMMPS *lmp) : PairMEAM(lmp)
{
  respa_enable = 0;

  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
  meam_inst_kk = new MEAMKokkos<DeviceType>(memory);
  delete meam_inst;
  meam_inst = meam_inst_kk;
}
/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairMEAMKokkos<DeviceType>::~PairMEAMKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    delete meam_inst_kk;
  }
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
template<class DeviceType>
void PairMEAMKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,6,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);


  // neighbor list info

  NeighListKokkos<DeviceType>* k_halflist = static_cast<NeighListKokkos<DeviceType>*>(listhalf);
  int inum_half = listhalf->inum;
  int* numneigh_half = listhalf->numneigh;
  int* ilist_half = listhalf->ilist;

  d_ilist_half = k_halflist->d_ilist;
  d_numneigh_half = k_halflist->d_numneigh;
  d_neighbors_half = k_halflist->d_neighbors;
  NeighListKokkos<DeviceType>* k_fulllist = static_cast<NeighListKokkos<DeviceType>*>(listfull);
  d_numneigh_full = k_fulllist->d_numneigh;
  d_neighbors_full = k_fulllist->d_neighbors;

  copymode = 1;

  // strip neighbor lists of any special bond flags before using with MEAM
  // necessary before doing neigh_f2c and neigh_c2f conversions each step
  if (neighbor->ago == 0) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMEAMKernelNeighStrip >(0,inum_half),*this);
  }

  // check size of scrfcn based on half neighbor list

  nlocal = atom->nlocal;
  nall = nlocal + atom->nghost;

  int n = 0;
  //for (ii = 0; ii < inum_half; ii++) n += numneigh_half[ilist_half[ii]];
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairMEAMKernelA>(0,inum_half), *this, n);

  meam_inst_kk->meam_dens_setup(atom->nmax, nall, n);

  //double **x = atom->x;
  x = atomKK->k_x.view<DeviceType>();

  //double **f = atom->f;
  f = atomKK->k_f.view<DeviceType>();

  //int *type = atom->type;
  type = atomKK->k_type.view<DeviceType>();

  int ntype = atom->ntypes;

  // 3 stages of MEAM calculation
  // loop over my atoms followed by communication

  int offset = 0;
  int errorflag = 0;
#if 0
  for (ii = 0; ii < inum_half; ii++) {
    i = ilist_half[ii];
    meam_inst->meam_dens_init(i,ntype,type,map,x,
                    numneigh_half[i],firstneigh_half[i],
                    numneigh_full[i],firstneigh_full[i],
                    offset);
    offset += numneigh_half[i];
  }
#endif
  // To do: create the cumulative offset array in host and device
  
  k_offset = DAT::tdual_int_1d("pair:offset",inum_half+1);
  h_offset = k_offset.h_view;
  d_offset = k_offset.template view<DeviceType>();
  ArrayTypes<LMPHostType>::t_int_1d h_ilist;
  ArrayTypes<LMPHostType>::t_int_1d h_numneigh;
  h_ilist = Kokkos::create_mirror_view(k_halflist->d_ilist);
  h_numneigh = Kokkos::create_mirror_view(k_halflist->d_numneigh);
  Kokkos::deep_copy(h_ilist,k_halflist->d_ilist);
  Kokkos::deep_copy(h_numneigh,k_halflist->d_numneigh);

  h_offset[0] = 0;
  for (int ii = 0; ii < inum_half; ii++) {
    int i = h_ilist[ii];
    h_offset[ii+1] = h_offset[ii] + h_numneigh[i]; 
  }
  k_offset.template modify<LMPHostType>();
  k_offset.template sync<DeviceType>();
  meam_inst_kk->meam_dens_init(inum_half,ntype,type,d_map,x,d_numneigh_half,d_numneigh_full,&offset,d_ilist_half,d_neighbors_half, d_neighbors_full, d_offset, neighflag);
  meam_inst_kk->k_rho0.template modify<DeviceType>();
  meam_inst_kk->k_rho0.template sync<LMPHostType>();

  meam_inst_kk->k_arho2b.template modify<DeviceType>();
  meam_inst_kk->k_arho2b.template sync<LMPHostType>();

  meam_inst_kk->k_arho1.template modify<DeviceType>();
  meam_inst_kk->k_arho1.template sync<LMPHostType>();

  meam_inst_kk->k_arho2.template modify<DeviceType>();
  meam_inst_kk->k_arho2.template sync<LMPHostType>();

  meam_inst_kk->k_arho3.template modify<DeviceType>();
  meam_inst_kk->k_arho3.template sync<LMPHostType>();

  meam_inst_kk->k_arho3b.template modify<DeviceType>();
  meam_inst_kk->k_arho3b.template sync<LMPHostType>();

  meam_inst_kk->k_t_ave.template modify<DeviceType>();
  meam_inst_kk->k_t_ave.template sync<LMPHostType>();

  meam_inst_kk->k_tsq_ave.template modify<DeviceType>();
  meam_inst_kk->k_tsq_ave.template sync<LMPHostType>();

  comm->reverse_comm(this);

  meam_inst_kk->k_rho0.template modify<LMPHostType>();
  meam_inst_kk->k_rho0.template sync<DeviceType>();

  meam_inst_kk->k_arho2b.template modify<LMPHostType>();
  meam_inst_kk->k_arho2b.template sync<DeviceType>();

  meam_inst_kk->k_arho1.template modify<LMPHostType>();
  meam_inst_kk->k_arho1.template sync<DeviceType>();

  meam_inst_kk->k_arho2.template modify<LMPHostType>();
  meam_inst_kk->k_arho2.template sync<DeviceType>();

  meam_inst_kk->k_arho3.template modify<LMPHostType>();
  meam_inst_kk->k_arho3.template sync<DeviceType>();

  meam_inst_kk->k_arho3b.template modify<LMPHostType>();
  meam_inst_kk->k_arho3b.template sync<DeviceType>();

  meam_inst_kk->k_t_ave.template modify<LMPHostType>();
  meam_inst_kk->k_t_ave.template sync<DeviceType>();

  meam_inst_kk->k_tsq_ave.template modify<LMPHostType>();
  meam_inst_kk->k_tsq_ave.template sync<DeviceType>();


  meam_inst_kk->meam_dens_final(nlocal,eflag_either,eflag_global,eflag_atom,
                   &eng_vdwl,d_eatom,ntype,type,d_map,errorflag);
  if (errorflag) {
    char str[128];
    sprintf(str,"MEAM library error %d",errorflag);
    error->one(FLERR,str);
  }

  comm->forward_comm(this);

  offset = 0;

  // vptr is first value in vatom if it will be used by meam_force()
  // else vatom may not exist, so pass dummy ptr

#if 0 // To do: is this correct? vflag_atom is used to access vatom
  typename ArrayTypes<DeviceType>::t_virial_array vptr;
  if (vflag_atom) vptr = d_vatom;
  else vptr = NULL;
  for (ii = 0; ii < inum_half; ii++) {
    i = ilist_half[ii];
    meam_inst->meam_force(i,eflag_either,eflag_global,eflag_atom,
                vflag_atom,&eng_vdwl,eatom,ntype,type,map,x,
                numneigh_half[i],firstneigh_half[i],
                numneigh_full[i],firstneigh_full[i],
                offset,f,vptr);
    offset += numneigh_half[i];
  }
#endif
  meam_inst_kk->meam_force(inum_half, eflag_either,eflag_global,eflag_atom,
                vflag_atom,&eng_vdwl,d_eatom,ntype,type,d_map,x,
                d_numneigh_half, d_numneigh_full,f,d_vatom,d_ilist_half, d_offset, d_neighbors_half, d_neighbors_full, neighflag);

  if (vflag_fdotr) pair_virial_fdotr_compute(this);
  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */
template<class DeviceType>
void PairMEAMKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairMEAM::coeff(narg,arg);

  //sync map

  int n = atom->ntypes;
  k_map = DAT::tdual_int_1d("pair:map",n+1);
  HAT::t_int_1d h_map = k_map.h_view;

  for (int i = 1; i <= n; i++)
    h_map[i] = map[i];

  k_map.template modify<LMPHostType>();
  k_map.template sync<DeviceType>();

  d_map = k_map.template view<DeviceType>();

  // To do: need to synchronize phirar variables
   
  meam_inst_kk->meam_setup_done(); 
  
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */
template<class DeviceType>
void PairMEAMKokkos<DeviceType>::init_style()
{

  PairMEAM::init_style();

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);

  // MEAM needs both a full and half neighbor list? Not sure how to get that.
  if (!(neighflag == FULL || neighflag == HALF || neighflag == HALFTHREAD))
    error->all(FLERR, "Cannot use chosen neighbor list style with pair meam/kk");

  request->set_kokkos_host(std::is_same<DeviceType,LMPHostType>::value &&
                           !std::is_same<DeviceType,LMPDeviceType>::value);
  request->set_kokkos_device(std::is_same<DeviceType,LMPDeviceType>::value);
  if (neighflag == FULL) request->enable_full();

}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
int PairMEAMKokkos<DeviceType>::pack_forward_comm_kokkos(int n, DAT::tdual_int_2d k_sendlist, int iswap_in, DAT::tdual_xfloat_1d &buf,
                                int pbc_flag, int *pbc)
{
  d_sendlist = k_sendlist.view<DeviceType>();
  iswap = iswap_in;
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMEAMPackForwardComm>(0,n),*this);
  return n;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairMEAMKokkos<DeviceType>::operator()(TagPairMEAMPackForwardComm, const int &i) const {
  int j = d_sendlist(iswap, i);
  int m = i*38;
  v_buf[m++] = meam_inst_kk->d_rho0[j];
  v_buf[m++] = meam_inst_kk->d_rho1[j];
  v_buf[m++] = meam_inst_kk->d_rho2[j];
    v_buf[m++] = meam_inst_kk->d_rho3[j];
    v_buf[m++] = meam_inst_kk->d_frhop[j];
    v_buf[m++] = meam_inst_kk->d_gamma[j];
    v_buf[m++] = meam_inst_kk->d_dgamma1[j];
    v_buf[m++] = meam_inst_kk->d_dgamma2[j];
    v_buf[m++] = meam_inst_kk->d_dgamma3[j];
    v_buf[m++] = meam_inst_kk->d_arho2b[j];
    v_buf[m++] = meam_inst_kk->d_arho1(j,0);
    v_buf[m++] = meam_inst_kk->d_arho1(j,1);
    v_buf[m++] = meam_inst_kk->d_arho1(j,2);
    v_buf[m++] = meam_inst_kk->d_arho2(j,0);
    v_buf[m++] = meam_inst_kk->d_arho2(j,1);
    v_buf[m++] = meam_inst_kk->d_arho2(j,2);
    v_buf[m++] = meam_inst_kk->d_arho2(j,3);
    v_buf[m++] = meam_inst_kk->d_arho2(j,4);
    v_buf[m++] = meam_inst_kk->d_arho2(j,5);
    for (int k = 0; k < 10; k++) v_buf[m++] = meam_inst_kk->d_arho3(j,k);
    v_buf[m++] = meam_inst_kk->d_arho3b(j,0);
    v_buf[m++] = meam_inst_kk->d_arho3b(j,1);
    v_buf[m++] = meam_inst_kk->d_arho3b(j,2);
    v_buf[m++] = meam_inst_kk->d_t_ave(j,0);
    v_buf[m++] = meam_inst_kk->d_t_ave(j,1);
    v_buf[m++] = meam_inst_kk->d_t_ave(j,2);
    v_buf[m++] = meam_inst_kk->d_tsq_ave(j,0);
    v_buf[m++] = meam_inst_kk->d_tsq_ave(j,1);
    v_buf[m++] = meam_inst_kk->d_tsq_ave(j,2);
} 

/* ---------------------------------------------------------------------- */
template<class DeviceType>
void PairMEAMKokkos<DeviceType>::unpack_forward_comm_kokkos(int n, int first_in, DAT::tdual_xfloat_1d &buf)
{
  first = first_in;
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMEAMUnpackForwardComm>(0,n),*this);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairMEAMKokkos<DeviceType>::operator()(TagPairMEAMUnpackForwardComm, const int &i) const{
   int m = i*38;

    meam_inst_kk->d_rho0[i+first] = v_buf[m++];
    meam_inst_kk->d_rho1[i+first] = v_buf[m++];
    meam_inst_kk->d_rho2[i+first] = v_buf[m++];
    meam_inst_kk->d_rho3[i+first] = v_buf[m++];
    meam_inst_kk->d_frhop[i+first] = v_buf[m++];
    meam_inst_kk->d_gamma[i+first] = v_buf[m++];
    meam_inst_kk->d_dgamma1[i+first] = v_buf[m++];
    meam_inst_kk->d_dgamma2[i+first] = v_buf[m++];
    meam_inst_kk->d_dgamma3[i+first] = v_buf[m++];
    meam_inst_kk->d_arho2b[i+first] = v_buf[m++];
    meam_inst_kk->d_arho1(i+first,0) = v_buf[m++];
    meam_inst_kk->d_arho1(i+first,1) = v_buf[m++];
    meam_inst_kk->d_arho1(i+first,2) = v_buf[m++];
    meam_inst_kk->d_arho2(i+first,0) = v_buf[m++];
    meam_inst_kk->d_arho2(i+first,1) = v_buf[m++];
    meam_inst_kk->d_arho2(i+first,2) = v_buf[m++];
    meam_inst_kk->d_arho2(i+first,3) = v_buf[m++];
    meam_inst_kk->d_arho2(i+first,4) = v_buf[m++];
    meam_inst_kk->d_arho2(i+first,5) = v_buf[m++];
    for (int k = 0; k < 10; k++) meam_inst_kk->d_arho3(i+first,k) = v_buf[m++];
    meam_inst_kk->d_arho3b(i+first,0) = v_buf[m++];
    meam_inst_kk->d_arho3b(i+first,1) = v_buf[m++];
    meam_inst_kk->d_arho3b(i+first,2) = v_buf[m++];
    meam_inst_kk->d_t_ave(i+first,0) = v_buf[m++];
    meam_inst_kk->d_t_ave(i+first,1) = v_buf[m++];
    meam_inst_kk->d_t_ave(i+first,2) = v_buf[m++];
    meam_inst_kk->d_tsq_ave(i+first,0) = v_buf[m++];
    meam_inst_kk->d_tsq_ave(i+first,1) = v_buf[m++];
    meam_inst_kk->d_tsq_ave(i+first,2) = v_buf[m++];
 } 

template<class DeviceType>
int PairMEAMKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf,
                               int pbc_flag, int *pbc)
{
    int i,j,k,m;

    m = 0;

    for (i = 0; i < n; i++) {
       j = list[i];
    buf[m++] = meam_inst_kk->h_rho0[j];
    buf[m++] = meam_inst_kk->h_rho1[j];
    buf[m++] = meam_inst_kk->h_rho2[j];
    buf[m++] = meam_inst_kk->h_rho3[j];
    buf[m++] = meam_inst_kk->h_frhop[j];
    buf[m++] = meam_inst_kk->h_gamma[j];
    buf[m++] = meam_inst_kk->h_dgamma1[j];
    buf[m++] = meam_inst_kk->h_dgamma2[j];
    buf[m++] = meam_inst_kk->h_dgamma3[j];
    buf[m++] = meam_inst_kk->h_arho2b[j];
    buf[m++] = meam_inst_kk->h_arho1(j,0);
    buf[m++] = meam_inst_kk->h_arho1(j,1);
    buf[m++] = meam_inst_kk->h_arho1(j,2);
    buf[m++] = meam_inst_kk->h_arho2(j,0);
    buf[m++] = meam_inst_kk->h_arho2(j,1);
    buf[m++] = meam_inst_kk->h_arho2(j,2);
    buf[m++] = meam_inst_kk->h_arho2(j,3);
    buf[m++] = meam_inst_kk->h_arho2(j,4);
    buf[m++] = meam_inst_kk->h_arho2(j,5);
    for (k = 0; k < 10; k++) buf[m++] = meam_inst_kk->h_arho3(j,k);
    buf[m++] = meam_inst_kk->h_arho3b(j,0);
    buf[m++] = meam_inst_kk->h_arho3b(j,1);
    buf[m++] = meam_inst_kk->h_arho3b(j,2);
    buf[m++] = meam_inst_kk->h_t_ave(j,0);
    buf[m++] = meam_inst_kk->h_t_ave(j,1);
    buf[m++] = meam_inst_kk->h_t_ave(j,2);
    buf[m++] = meam_inst_kk->h_tsq_ave(j,0);
    buf[m++] = meam_inst_kk->h_tsq_ave(j,1);
    buf[m++] = meam_inst_kk->h_tsq_ave(j,2);
  }

  return m;
}

template<class DeviceType>
void PairMEAMKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
    int i,k,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    meam_inst_kk->h_rho0[i] = buf[m++];
    meam_inst_kk->h_rho1[i] = buf[m++];
    meam_inst_kk->h_rho2[i] = buf[m++];
    meam_inst_kk->h_rho3[i] = buf[m++];
    meam_inst_kk->h_frhop[i] = buf[m++];
    meam_inst_kk->h_gamma[i] = buf[m++];
    meam_inst_kk->h_dgamma1[i] = buf[m++];
    meam_inst_kk->h_dgamma2[i] = buf[m++];
    meam_inst_kk->h_dgamma3[i] = buf[m++];
    meam_inst_kk->h_arho2b[i] = buf[m++];
    meam_inst_kk->h_arho1(i,0) = buf[m++];
    meam_inst_kk->h_arho1(i,1) = buf[m++];
    meam_inst_kk->h_arho1(i,2) = buf[m++];
    meam_inst_kk->h_arho2(i,0) = buf[m++];
    meam_inst_kk->h_arho2(i,1) = buf[m++];
    meam_inst_kk->h_arho2(i,2) = buf[m++];
    meam_inst_kk->h_arho2(i,3) = buf[m++];
    meam_inst_kk->h_arho2(i,4) = buf[m++];
    meam_inst_kk->h_arho2(i,5) = buf[m++];
    for (k = 0; k < 10; k++) meam_inst_kk->h_arho3(i,k) = buf[m++];
    meam_inst_kk->h_arho3b(i,0) = buf[m++];
    meam_inst_kk->h_arho3b(i,1) = buf[m++];
    meam_inst_kk->h_arho3b(i,2) = buf[m++];
    meam_inst_kk->h_t_ave(i,0) = buf[m++];
    meam_inst_kk->h_t_ave(i,1) = buf[m++];
    meam_inst_kk->h_t_ave(i,2) = buf[m++];
    meam_inst_kk->h_tsq_ave(i,0) = buf[m++];
    meam_inst_kk->h_tsq_ave(i,1) = buf[m++];
    meam_inst_kk->h_tsq_ave(i,2) = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
int PairMEAMKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  int i,k,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = meam_inst_kk->h_rho0[i];
    buf[m++] = meam_inst_kk->h_arho2b[i];
    buf[m++] = meam_inst_kk->h_arho1(i,0);
    buf[m++] = meam_inst_kk->h_arho1(i,1);
    buf[m++] = meam_inst_kk->h_arho1(i,2);
    buf[m++] = meam_inst_kk->h_arho2(i,0);
    buf[m++] = meam_inst_kk->h_arho2(i,1);
    buf[m++] = meam_inst_kk->h_arho2(i,2);
    buf[m++] = meam_inst_kk->h_arho2(i,3);
    buf[m++] = meam_inst_kk->h_arho2(i,4);
    buf[m++] = meam_inst_kk->h_arho2(i,5);
    for (k = 0; k < 10; k++) buf[m++] = meam_inst_kk->h_arho3(i,k);
    buf[m++] = meam_inst_kk->h_arho3b(i,0);
    buf[m++] = meam_inst_kk->h_arho3b(i,1);
    buf[m++] = meam_inst_kk->h_arho3b(i,2);
    buf[m++] = meam_inst_kk->h_t_ave(i,0);
    buf[m++] = meam_inst_kk->h_t_ave(i,1);
    buf[m++] = meam_inst_kk->h_t_ave(i,2);
    buf[m++] = meam_inst_kk->h_tsq_ave(i,0);
    buf[m++] = meam_inst_kk->h_tsq_ave(i,1);
    buf[m++] = meam_inst_kk->h_tsq_ave(i,2);
  }

  return m;
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
void PairMEAMKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,k,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    meam_inst_kk->h_rho0[j] += buf[m++];
    meam_inst_kk->h_arho2b[j] += buf[m++];
    meam_inst_kk->h_arho1(j,0) += buf[m++];
    meam_inst_kk->h_arho1(j,1) += buf[m++];
    meam_inst_kk->h_arho1(j,2) += buf[m++];
    meam_inst_kk->h_arho2(j,0) += buf[m++];
    meam_inst_kk->h_arho2(j,1) += buf[m++];
    meam_inst_kk->h_arho2(j,2) += buf[m++];
    meam_inst_kk->h_arho2(j,3) += buf[m++];
    meam_inst_kk->h_arho2(j,4) += buf[m++];
    meam_inst_kk->h_arho2(j,5) += buf[m++];
    for (k = 0; k < 10; k++) meam_inst_kk->h_arho3(j,k) += buf[m++];
    meam_inst_kk->h_arho3b(j,0) += buf[m++];
    meam_inst_kk->h_arho3b(j,1) += buf[m++];
    meam_inst_kk->h_arho3b(j,2) += buf[m++];
    meam_inst_kk->h_t_ave(j,0) += buf[m++];
    meam_inst_kk->h_t_ave(j,1) += buf[m++];
    meam_inst_kk->h_t_ave(j,2) += buf[m++];
    meam_inst_kk->h_tsq_ave(j,0) += buf[m++];
    meam_inst_kk->h_tsq_ave(j,1) += buf[m++];
    meam_inst_kk->h_tsq_ave(j,2) += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */
template<class DeviceType>
double PairMEAMKokkos<DeviceType>::memory_usage()
{
  double bytes = 11 * meam_inst_kk->nmax * sizeof(double);
  bytes += (3 + 6 + 10 + 3 + 3 + 3) * meam_inst_kk->nmax * sizeof(double);
  bytes += 3 * meam_inst_kk->maxneigh * sizeof(double);
  return bytes;
}
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairMEAMKokkos<DeviceType>::operator()(TagPairMEAMKernelNeighStrip, const int &ii) const {
/* ----------------------------------------------------------------------
 *    strip special bond flags from neighbor list entries
 *       are not used with MEAM
 *          need to do here so Fortran lib doesn't see them
 *             done once per reneighbor so that neigh_f2c and neigh_c2f don't see them
 *             ------------------------------------------------------------------------- */
    const int i = d_ilist_half[ii];
    const int jnum_half = d_numneigh_half[i];
    const int jnum_full = d_numneigh_full[i];
    for (int jj = 0; jj < jnum_half; jj++) {
        d_neighbors_half(i,jj) &= NEIGHMASK;
    }
    for (int jj = 0; jj < jnum_full; jj++) {
        d_neighbors_full(i,jj) &= NEIGHMASK;
    }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairMEAMKokkos<DeviceType>::operator()(TagPairMEAMKernelA, const int ii, int &n) const {
    const int i = d_ilist_half[ii];
    n += d_numneigh_half[i]; 
}

namespace LAMMPS_NS {
template class PairMEAMKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class PairMEAMKokkos<LMPHostType>;
#endif
}

