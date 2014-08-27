;/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "neighbor_kokkos.h"
#include "atom.h"
#include "pair.h"
#include "neigh_request.h"
#include "memory.h"

using namespace LAMMPS_NS;

enum{NSQ,BIN,MULTI};     // also in neigh_list.cpp

/* ---------------------------------------------------------------------- */

NeighborKokkos::NeighborKokkos(LAMMPS *lmp) : Neighbor(lmp)
{
  atoms_per_bin = 16;

  nlist_host = 0;
  lists_host = NULL;
  pair_build_host = NULL;
  stencil_create_host = NULL;
  nlist_device = 0;
  lists_device = NULL;
  pair_build_device = NULL;
  stencil_create_device = NULL;
}

/* ---------------------------------------------------------------------- */

NeighborKokkos::~NeighborKokkos()
{
  memory->destroy_kokkos(k_cutneighsq,cutneighsq);
  cutneighsq = NULL;

  for (int i = 0; i < nlist_host; i++) delete lists_host[i];
  delete [] lists_host;
  for (int i = 0; i < nlist_device; i++) delete lists_device[i];
  delete [] lists_device;

  delete [] pair_build_device;
  delete [] pair_build_host;

  memory->destroy_kokkos(k_ex1_type,ex1_type);
  memory->destroy_kokkos(k_ex2_type,ex2_type);
  memory->destroy_kokkos(k_ex1_group,ex1_group);
  memory->destroy_kokkos(k_ex2_group,ex2_group);
  memory->destroy_kokkos(k_ex_mol_group,ex_mol_group);
  memory->destroy_kokkos(k_ex1_bit,ex1_bit);
  memory->destroy_kokkos(k_ex2_bit,ex2_bit);
  memory->destroy_kokkos(k_ex_mol_bit,ex_mol_bit);

}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::init()
{
  atomKK = (AtomKokkos *) atom;
  Neighbor::init();
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::init_cutneighsq_kokkos(int n)
{
  memory->create_kokkos(k_cutneighsq,cutneighsq,n+1,n+1,"neigh:cutneighsq");
  k_cutneighsq.modify<LMPHostType>();
}

/* ---------------------------------------------------------------------- */

int NeighborKokkos::init_lists_kokkos()
{
  int i;

  for (i = 0; i < nlist_host; i++) delete lists_host[i];
  delete [] lists_host;
  delete [] pair_build_host;
  delete [] stencil_create_host;
  nlist_host = 0;

  for (i = 0; i < nlist_device; i++) delete lists_device[i];
  delete [] lists_device;
  delete [] pair_build_device;
  delete [] stencil_create_device;
  nlist_device = 0;

  nlist = 0;
  for (i = 0; i < nrequest; i++) {
    if (requests[i]->kokkos_device) nlist_device++;
    else if (requests[i]->kokkos_host) nlist_host++;
    else nlist++;
  }

  lists_host = new NeighListKokkos<LMPHostType>*[nrequest];
  pair_build_host = new PairPtrHost[nrequest];
  stencil_create_host = new StencilPtrHost[nrequest];
  for (i = 0; i < nrequest; i++) {
    lists_host[i] = NULL;
    pair_build_host[i] = NULL;
    stencil_create_host[i] = NULL;
  }

  for (i = 0; i < nrequest; i++) {
    if (!requests[i]->kokkos_host) continue;
    lists_host[i] = new NeighListKokkos<LMPHostType>(lmp);
    lists_host[i]->index = i;
    lists_host[i]->dnum = requests[i]->dnum;
    if (requests[i]->pair) {
      Pair *pair = (Pair *) requests[i]->requestor;
      pair->init_list(requests[i]->id,lists_host[i]);
    }
  }

  lists_device = new NeighListKokkos<LMPDeviceType>*[nrequest];
  pair_build_device = new PairPtrDevice[nrequest];
  stencil_create_device = new StencilPtrDevice[nrequest];
  for (i = 0; i < nrequest; i++) {
    lists_device[i] = NULL;
    pair_build_device[i] = NULL;
    stencil_create_device[i] = NULL;
  }

  for (i = 0; i < nrequest; i++) {
    if (!requests[i]->kokkos_device) continue;
    lists_device[i] = new NeighListKokkos<LMPDeviceType>(lmp);
    lists_device[i]->index = i;
    lists_device[i]->dnum = requests[i]->dnum;
    if (requests[i]->pair) {
      Pair *pair = (Pair *) requests[i]->requestor;
      pair->init_list(requests[i]->id,lists_device[i]);
    }
  }

  // return # of non-Kokkos lists

  return nlist;
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::init_list_flags1_kokkos(int i)
{ 
  if (lists_host[i]) {
    lists_host[i]->buildflag = 1;
    if (pair_build_host[i] == NULL) lists_host[i]->buildflag = 0;
    if (requests[i]->occasional) lists_host[i]->buildflag = 0;
    
    lists_host[i]->growflag = 1;
    if (requests[i]->copy) lists_host[i]->growflag = 0;
    
    lists_host[i]->stencilflag = 1;
    if (style == NSQ) lists_host[i]->stencilflag = 0;
    if (stencil_create[i] == NULL) lists_host[i]->stencilflag = 0;
    
    lists_host[i]->ghostflag = 0;
    if (requests[i]->ghost) lists_host[i]->ghostflag = 1;
    if (requests[i]->ghost && !requests[i]->occasional) anyghostlist = 1;
  }
  
  if (lists_device[i]) {
    lists_device[i]->buildflag = 1;
    if (pair_build_device[i] == NULL) lists_device[i]->buildflag = 0;
    if (requests[i]->occasional) lists_device[i]->buildflag = 0;
    
    lists_device[i]->growflag = 1;
    if (requests[i]->copy) lists_device[i]->growflag = 0;
    
    lists_device[i]->stencilflag = 1;
    if (style == NSQ) lists_device[i]->stencilflag = 0;
    if (stencil_create[i] == NULL) lists_device[i]->stencilflag = 0;
    
    lists_device[i]->ghostflag = 0;
    if (requests[i]->ghost) lists_device[i]->ghostflag = 1;
    if (requests[i]->ghost && !requests[i]->occasional) anyghostlist = 1;
  }
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::init_list_flags2_kokkos(int i)
{ 
  if (lists_host[i]) {
    if (lists_host[i]->buildflag) blist[nblist++] = i;
    if (lists_host[i]->growflag && requests[i]->occasional == 0)
      glist[nglist++] = i;
    if (lists_host[i]->stencilflag && requests[i]->occasional == 0)
      slist[nslist++] = i;
  }

  if (lists_device[i]) {
    if (lists_device[i]->buildflag) blist[nblist++] = i;
    if (lists_device[i]->growflag && requests[i]->occasional == 0)
      glist[nglist++] = i;
    if (lists_device[i]->stencilflag && requests[i]->occasional == 0)
      slist[nslist++] = i;
  }
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::init_list_grow_kokkos(int i)
{
  if (lists_host[i]!=NULL && lists_host[i]->growflag)
    lists_host[i]->grow(maxatom);
  if (lists_device[i]!=NULL && lists_device[i]->growflag)
    lists_device[i]->grow(maxatom);
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::init_ex_type_kokkos(int n)
{
  memory->create_kokkos(k_ex_type,ex_type,n+1,n+1,"neigh:ex_type");
  k_ex_type.modify<LMPHostType>();
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::init_ex_bit_kokkos()
{
  memory->create_kokkos(k_ex1_bit, ex1_bit, nex_group, "neigh:ex1_bit");
  k_ex1_bit.modify<LMPHostType>();
  memory->create_kokkos(k_ex2_bit, ex2_bit, nex_group, "neigh:ex2_bit");
  k_ex2_bit.modify<LMPHostType>();
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::init_ex_mol_bit_kokkos()
{
  memory->create_kokkos(k_ex_mol_bit, ex_mol_bit, nex_mol, "neigh:ex_mol_bit");
  k_ex_mol_bit.modify<LMPHostType>();
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::choose_build(int index, NeighRequest *rq)
{
  if (rq->kokkos_host != 0) {
    PairPtrHost pb = NULL;
    if (rq->full) pb = &NeighborKokkos::full_bin_kokkos<LMPHostType,0>;
    else if (rq->half) pb = &NeighborKokkos::full_bin_kokkos<LMPHostType,1>;
    pair_build_host[index] = pb;
    return;
  }
  if (rq->kokkos_device != 0) {
    PairPtrDevice pb = NULL;
    if (rq->full) {
      if (rq->full_cluster) pb = &NeighborKokkos::full_bin_cluster_kokkos<LMPDeviceType>;
      else pb = &NeighborKokkos::full_bin_kokkos<LMPDeviceType,0>;
    }
    else if (rq->half) pb = &NeighborKokkos::full_bin_kokkos<LMPDeviceType,1>;
    pair_build_device[index] = pb;
    return;
  }

  Neighbor::choose_build(index,rq);
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::build_kokkos(int i)
{
  if (lists_host[blist[i]])
    (this->*pair_build_host[blist[i]])(lists_host[blist[i]]);
  else if (lists_device[blist[i]])
    (this->*pair_build_device[blist[i]])(lists_device[blist[i]]);
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::setup_bins_kokkos(int i)
{
  if (lists_host[slist[i]]) {
    lists_host[slist[i]]->stencil_allocate(smax,style);
    (this->*stencil_create[slist[i]])(lists_host[slist[i]],sx,sy,sz);
  } else if (lists_device[slist[i]]) {
    lists_device[slist[i]]->stencil_allocate(smax,style);
    (this->*stencil_create[slist[i]])(lists_device[slist[i]],sx,sy,sz);
  }

  if (i < nslist-1) return;

  if (maxhead > k_bins.d_view.dimension_0()) {
    k_bins = DAT::tdual_int_2d("Neighbor::d_bins",maxhead,atoms_per_bin);
    k_bincount = DAT::tdual_int_1d("Neighbor::d_bincount",maxhead);
  }
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::modify_ex_type_grow_kokkos(){
  memory->grow_kokkos(k_ex1_type,ex1_type,maxex_type,"neigh:ex1_type");
  k_ex1_type.modify<LMPHostType>();
  memory->grow_kokkos(k_ex2_type,ex2_type,maxex_type,"neigh:ex2_type");
  k_ex2_type.modify<LMPHostType>();
}

/* ---------------------------------------------------------------------- */
void NeighborKokkos::modify_ex_group_grow_kokkos(){
  memory->grow_kokkos(k_ex1_group,ex1_group,maxex_group,"neigh:ex1_group");
  k_ex1_group.modify<LMPHostType>();
  memory->grow_kokkos(k_ex2_group,ex2_group,maxex_group,"neigh:ex2_group");
  k_ex2_group.modify<LMPHostType>();
}

/* ---------------------------------------------------------------------- */
void NeighborKokkos::modify_mol_group_grow_kokkos(){
  memory->grow_kokkos(k_ex_mol_group,ex_mol_group,maxex_mol,"neigh:ex_mol_group");
  k_ex_mol_group.modify<LMPHostType>();
}

// include to trigger instantiation of templated functions

#include "neigh_full_kokkos.h"
