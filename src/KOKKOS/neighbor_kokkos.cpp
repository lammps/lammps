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
#include "update.h"
#include "atom_masks.h"
#include "error.h"
#include "kokkos.h"

using namespace LAMMPS_NS;

enum{NSQ,BIN,MULTI};     // also in neigh_list.cpp

/* ---------------------------------------------------------------------- */

NeighborKokkos::NeighborKokkos(LAMMPS *lmp) : Neighbor(lmp),
  neighbond_host(lmp),neighbond_device(lmp)
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

  device_flag = 0;
}

/* ---------------------------------------------------------------------- */

NeighborKokkos::~NeighborKokkos()
{
  if (!copymode) {
    memory->destroy_kokkos(k_cutneighsq,cutneighsq);
    cutneighsq = NULL;

    for (int i = 0; i < nlist_host; i++) delete lists_host[i];
    delete [] lists_host;
    for (int i = 0; i < nlist_device; i++) delete lists_device[i];
    delete [] lists_device;

    delete [] pair_build_device;
    delete [] pair_build_host;

    memory->destroy_kokkos(k_ex_type,ex_type);
    memory->destroy_kokkos(k_ex1_type,ex1_type);
    memory->destroy_kokkos(k_ex2_type,ex2_type);
    memory->destroy_kokkos(k_ex1_group,ex1_group);
    memory->destroy_kokkos(k_ex2_group,ex2_group);
    memory->destroy_kokkos(k_ex_mol_group,ex_mol_group);
    memory->destroy_kokkos(k_ex1_bit,ex1_bit);
    memory->destroy_kokkos(k_ex2_bit,ex2_bit);
    memory->destroy_kokkos(k_ex_mol_bit,ex_mol_bit);

    memory->destroy_kokkos(k_bondlist,bondlist);
    memory->destroy_kokkos(k_anglelist,anglelist);
    memory->destroy_kokkos(k_dihedrallist,dihedrallist);
    memory->destroy_kokkos(k_improperlist,improperlist);
  }
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

  // 1st time allocation of xhold

  if (dist_check)
      xhold = DAT::tdual_x_array("neigh:xhold",maxhold);

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
    if (rq->ghost) {
      if (rq->full) pb = &NeighborKokkos::full_bin_kokkos<LMPHostType,0,1>;
      else if (rq->half) error->one(FLERR,"Cannot (yet) request ghost atoms with Kokkos half neighbor list");
      pair_build_host[index] = pb;
    } else {
      if (rq->full) pb = &NeighborKokkos::full_bin_kokkos<LMPHostType,0,0>;
      else if (rq->half) pb = &NeighborKokkos::full_bin_kokkos<LMPHostType,1,0>;
      pair_build_host[index] = pb;
    }
    return;
  }
  if (rq->kokkos_device != 0) {
    PairPtrDevice pb = NULL;
    if (rq->ghost) {
      if (rq->full) {
        if (rq->full_cluster) pb = &NeighborKokkos::full_bin_cluster_kokkos<LMPDeviceType>;
        else pb = &NeighborKokkos::full_bin_kokkos<LMPDeviceType,0,1>;
      }
      else if (rq->half) pb = &NeighborKokkos::full_bin_kokkos<LMPDeviceType,1,1>;
    } else {
      if (rq->full) {
        if (rq->full_cluster) pb = &NeighborKokkos::full_bin_cluster_kokkos<LMPDeviceType>;
        else pb = &NeighborKokkos::full_bin_kokkos<LMPDeviceType,0,0>;
      }
      else if (rq->half) pb = &NeighborKokkos::full_bin_kokkos<LMPDeviceType,1,0>;
    }
    pair_build_device[index] = pb;
    return;
  }

  Neighbor::choose_build(index,rq);
}

/* ----------------------------------------------------------------------
   if any atom moved trigger distance (half of neighbor skin) return 1
   shrink trigger distance if box size has changed
   conservative shrink procedure:
     compute distance each of 8 corners of box has moved since last reneighbor
     reduce skin distance by sum of 2 largest of the 8 values
     new trigger = 1/2 of reduced skin distance
   for orthogonal box, only need 2 lo/hi corners
   for triclinic, need all 8 corners since deformations can displace all 8
------------------------------------------------------------------------- */

int NeighborKokkos::check_distance()
{
  if (nlist_device)
    check_distance_kokkos<LMPDeviceType>();
  else
    check_distance_kokkos<LMPHostType>();
}

template<class DeviceType>
int NeighborKokkos::check_distance_kokkos()
{
  typedef DeviceType device_type;

  double delx,dely,delz,rsq;
  double delta,delta1,delta2;

  if (boxcheck) {
    if (triclinic == 0) {
      delx = bboxlo[0] - boxlo_hold[0];
      dely = bboxlo[1] - boxlo_hold[1];
      delz = bboxlo[2] - boxlo_hold[2];
      delta1 = sqrt(delx*delx + dely*dely + delz*delz);
      delx = bboxhi[0] - boxhi_hold[0];
      dely = bboxhi[1] - boxhi_hold[1];
      delz = bboxhi[2] - boxhi_hold[2];
      delta2 = sqrt(delx*delx + dely*dely + delz*delz);
      delta = 0.5 * (skin - (delta1+delta2));
      deltasq = delta*delta;
    } else {
      domain->box_corners();
      delta1 = delta2 = 0.0;
      for (int i = 0; i < 8; i++) {
        delx = corners[i][0] - corners_hold[i][0];
        dely = corners[i][1] - corners_hold[i][1];
        delz = corners[i][2] - corners_hold[i][2];
        delta = sqrt(delx*delx + dely*dely + delz*delz);
        if (delta > delta1) delta1 = delta;
        else if (delta > delta2) delta2 = delta;
      }
      delta = 0.5 * (skin - (delta1+delta2));
      deltasq = delta*delta;
    }
  } else deltasq = triggersq;

  atomKK->sync(ExecutionSpaceFromDevice<DeviceType>::space,X_MASK);
  x = atomKK->k_x;
  xhold.sync<DeviceType>();
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int flag = 0;
  copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagNeighborCheckDistance<DeviceType> >(0,nlocal),*this,flag);
  DeviceType::fence();
  copymode = 0;

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  if (flagall && ago == MAX(every,delay)) ndanger++;
  return flagall;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NeighborKokkos::operator()(TagNeighborCheckDistance<DeviceType>, const int &i, int &flag) const {
  typedef DeviceType device_type;
  const X_FLOAT delx = x.view<DeviceType>()(i,0) - xhold.view<DeviceType>()(i,0);
  const X_FLOAT dely = x.view<DeviceType>()(i,1) - xhold.view<DeviceType>()(i,1);
  const X_FLOAT delz = x.view<DeviceType>()(i,2) - xhold.view<DeviceType>()(i,2);
  const X_FLOAT rsq = delx*delx + dely*dely + delz*delz;
  if (rsq > deltasq) flag = 1;
}

/* ----------------------------------------------------------------------
   build perpetuals neighbor lists
   called at setup and every few timesteps during run or minimization
   topology lists also built if topoflag = 1, USER-CUDA calls with topoflag = 0
------------------------------------------------------------------------- */


void NeighborKokkos::build(int topoflag)
{
  if (nlist_device)
    build_kokkos<LMPDeviceType>(topoflag);
  else
    build_kokkos<LMPHostType>(topoflag);
}

template<class DeviceType>
void NeighborKokkos::build_kokkos(int topoflag)
{
  typedef DeviceType device_type;

  int i;

  ago = 0;
  ncalls++;
  lastcall = update->ntimestep;

  // store current atom positions and box size if needed

  if (dist_check) {
    atomKK->sync(ExecutionSpaceFromDevice<DeviceType>::space,X_MASK);
    x = atomKK->k_x;
    int nlocal = atom->nlocal;
    if (includegroup) nlocal = atom->nfirst;
    int maxhold_kokkos = xhold.view<DeviceType>().dimension_0();
    if (nlocal > maxhold || maxhold_kokkos < maxhold) {
      maxhold = atom->nmax;
      xhold = DAT::tdual_x_array("neigh:xhold",maxhold);
    }
    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagNeighborXhold<DeviceType> >(0,nlocal),*this);
    DeviceType::fence();
    copymode = 0;
    xhold.modify<DeviceType>();
    if (boxcheck) {
      if (triclinic == 0) {
        boxlo_hold[0] = bboxlo[0];
        boxlo_hold[1] = bboxlo[1];
        boxlo_hold[2] = bboxlo[2];
        boxhi_hold[0] = bboxhi[0];
        boxhi_hold[1] = bboxhi[1];
        boxhi_hold[2] = bboxhi[2];
      } else {
        domain->box_corners();
        corners = domain->corners;
        for (i = 0; i < 8; i++) {
          corners_hold[i][0] = corners[i][0];
          corners_hold[i][1] = corners[i][1];
          corners_hold[i][2] = corners[i][2];
        }
      }
    }
  }

  // if any lists store neighbors of ghosts:
  //   invoke grow() if nlocal+nghost exceeds previous list size
  // else only invoke grow() if nlocal exceeds previous list size
  // only for lists with growflag set and which are perpetual (glist)

  if (anyghostlist && atom->nlocal+atom->nghost > maxatom) {
    maxatom = atom->nmax;
    for (i = 0; i < nglist; i++)
      if (lists[glist[i]]) lists[glist[i]]->grow(maxatom);
      else init_list_grow_kokkos(glist[i]);
  } else if (atom->nlocal > maxatom) {
    maxatom = atom->nmax;
    for (i = 0; i < nglist; i++)
      if (lists[glist[i]]) lists[glist[i]]->grow(maxatom);
      else init_list_grow_kokkos(glist[i]);
  }

  // extend atom bin list if necessary

  if (style != NSQ && atom->nmax > maxbin) {
    maxbin = atom->nmax;
    memory->destroy(bins);
    memory->create(bins,maxbin,"bins");
  }

  // check that using special bond flags will not overflow neigh lists

  if (atom->nlocal+atom->nghost > NEIGHMASK)
    error->one(FLERR,"Too many local+ghost atoms for neighbor list");

  // invoke building of pair and molecular topology neighbor lists
  // only for pairwise lists with buildflag set
  // blist is for standard neigh lists, otherwise is a Kokkos list

  for (i = 0; i < nblist; i++) {
    if (lists[blist[i]]) {
      atomKK->sync(Host,ALL_MASK);
      (this->*pair_build[blist[i]])(lists[blist[i]]);
    } else {
      if (lists_host[blist[i]])
        (this->*pair_build_host[blist[i]])(lists_host[blist[i]]);
      else if (lists_device[blist[i]])
        (this->*pair_build_device[blist[i]])(lists_device[blist[i]]);
    }
  }

  if (atom->molecular && topoflag)
    build_topology_kokkos();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NeighborKokkos::operator()(TagNeighborXhold<DeviceType>, const int &i) const {
  typedef DeviceType device_type;
  xhold.view<DeviceType>()(i,0) = x.view<DeviceType>()(i,0);
  xhold.view<DeviceType>()(i,1) = x.view<DeviceType>()(i,1);
  xhold.view<DeviceType>()(i,2) = x.view<DeviceType>()(i,2);
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

  //if (i < nslist-1) return; // this won't work if a non-kokkos neighbor list is last

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

/* ---------------------------------------------------------------------- */

void NeighborKokkos::init_topology_kokkos() {
  if (nlist_device) {
    neighbond_device.init_topology_kk();
  } else {
    neighbond_host.init_topology_kk();
  }
}

/* ----------------------------------------------------------------------
   build all topology neighbor lists every few timesteps
   normally built with pair lists, but USER-CUDA separates them
------------------------------------------------------------------------- */

void NeighborKokkos::build_topology_kokkos() {
  if (nlist_device) {
    neighbond_device.build_topology_kk();

    k_bondlist = neighbond_device.k_bondlist;
    k_anglelist = neighbond_device.k_anglelist;
    k_dihedrallist = neighbond_device.k_dihedrallist;
    k_improperlist = neighbond_device.k_improperlist;

    k_bondlist.modify<LMPDeviceType>();
    k_anglelist.modify<LMPDeviceType>();
    k_dihedrallist.modify<LMPDeviceType>();
    k_improperlist.modify<LMPDeviceType>();
  } else {
    neighbond_host.build_topology_kk();

    k_bondlist = neighbond_host.k_bondlist;
    k_anglelist = neighbond_host.k_anglelist;
    k_dihedrallist = neighbond_host.k_dihedrallist;
    k_improperlist = neighbond_host.k_improperlist;

    k_bondlist.modify<LMPHostType>();
    k_anglelist.modify<LMPHostType>();
    k_dihedrallist.modify<LMPHostType>();
    k_improperlist.modify<LMPHostType>();
  }
}

// include to trigger instantiation of templated functions

#include "neigh_full_kokkos.h"
