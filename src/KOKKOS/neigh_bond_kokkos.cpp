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
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "neigh_bond_kokkos.h"
#include "atom_kokkos.h"
#include "atom_vec.h"
#include "molecule.h"
#include "force.h"
#include "update.h"
#include "domain_kokkos.h"
#include "output.h"
#include "thermo.h"
#include "memory_kokkos.h"
#include "error.h"
#include "modify.h"
#include "fix.h"
#include <cstring>
#include "atom_masks.h"
#include "domain.h"

using namespace LAMMPS_NS;

#define BONDDELTA 10000
#define LB_FACTOR 1.5

/* ---------------------------------------------------------------------- */

template<class DeviceType>
NeighBondKokkos<DeviceType>::NeighBondKokkos(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  k_nlist = DAT::tdual_int_scalar("NeighBond:nlist");
  d_nlist = k_nlist.view<DeviceType>();
  h_nlist = k_nlist.h_view;

  k_fail_flag = DAT::tdual_int_scalar("NeighBond:fail_flag");
  d_fail_flag = k_fail_flag.view<DeviceType>();
  h_fail_flag = k_fail_flag.h_view;

  maxbond = 0;
  maxangle = 0;
  maxdihedral = 0;
  maximproper = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NeighBondKokkos<DeviceType>::init_topology_kk() {

  atomKK = (AtomKokkos *) atom;
  atomKK->sync(Host,BOND_MASK | ANGLE_MASK | DIHEDRAL_MASK | IMPROPER_MASK);

  // topology lists

  // 1st time allocation of topology lists

  if (atom->molecular && atom->nbonds && maxbond == 0) {
    if (nprocs == 1) maxbond = atom->nbonds;
    else maxbond = static_cast<int> (LB_FACTOR * atom->nbonds / nprocs);
    memoryKK->create_kokkos(k_bondlist,neighbor->bondlist,maxbond,3,"neigh:neighbor->bondlist");
  }

  if (atom->molecular && atom->nangles && maxangle == 0) {
    if (nprocs == 1) maxangle = atom->nangles;
    else maxangle = static_cast<int> (LB_FACTOR * atom->nangles / nprocs);
    memoryKK->create_kokkos(k_anglelist,neighbor->anglelist,maxangle,4,"neigh:neighbor->anglelist");
  }

  if (atom->molecular && atom->ndihedrals && maxdihedral == 0) {
    if (nprocs == 1) maxdihedral = atom->ndihedrals;
    else maxdihedral = static_cast<int>
           (LB_FACTOR * atom->ndihedrals / nprocs);
    memoryKK->create_kokkos(k_dihedrallist,neighbor->dihedrallist,maxdihedral,5,"neigh:neighbor->dihedrallist");
  }

  if (atom->molecular && atom->nimpropers && maximproper == 0) {
    if (nprocs == 1) maximproper = atom->nimpropers;
    else maximproper = static_cast<int>
           (LB_FACTOR * atom->nimpropers / nprocs);
    memoryKK->create_kokkos(k_improperlist,neighbor->improperlist,maximproper,5,"neigh:neighbor->improperlist");
  }

  // set flags that determine which topology neighboring routines to use
  // bonds,etc can only be broken for atom->molecular = 1, not 2
  // SHAKE sets bonds and angles negative
  // gcmc sets all bonds, angles, etc negative
  // bond_quartic sets bonds to 0
  // delete_bonds sets all interactions negative

  int i,m;
  int bond_off = 0;
  int angle_off = 0;
  for (i = 0; i < modify->nfix; i++)
    if ((strcmp(modify->fix[i]->style,"shake") == 0)
        || (strcmp(modify->fix[i]->style,"rattle") == 0))
      bond_off = angle_off = 1;
  if (force->bond && force->bond_match("quartic")) bond_off = 1;

  if (atom->avec->bonds_allow && atom->molecular == 1) {
    for (i = 0; i < atom->nlocal; i++) {
      if (bond_off) break;
      for (m = 0; m < atom->num_bond[i]; m++)
        if (atom->bond_type[i][m] <= 0) bond_off = 1;
    }
  }

  if (atom->avec->angles_allow && atom->molecular == 1) {
    for (i = 0; i < atom->nlocal; i++) {
      if (angle_off) break;
      for (m = 0; m < atom->num_angle[i]; m++)
        if (atom->angle_type[i][m] <= 0) angle_off = 1;
    }
  }

  int dihedral_off = 0;
  if (atom->avec->dihedrals_allow && atom->molecular == 1) {
    for (i = 0; i < atom->nlocal; i++) {
      if (dihedral_off) break;
      for (m = 0; m < atom->num_dihedral[i]; m++)
        if (atom->dihedral_type[i][m] <= 0) dihedral_off = 1;
    }
  }

  int improper_off = 0;
  if (atom->avec->impropers_allow && atom->molecular == 1) {
    for (i = 0; i < atom->nlocal; i++) {
      if (improper_off) break;
      for (m = 0; m < atom->num_improper[i]; m++)
        if (atom->improper_type[i][m] <= 0) improper_off = 1;
    }
  }

  for (i = 0; i < modify->nfix; i++)
    if ((strcmp(modify->fix[i]->style,"gcmc") == 0))
      bond_off = angle_off = dihedral_off = improper_off = 1;

  // sync on/off settings across all procs

  int on_or_off = bond_off;
  MPI_Allreduce(&on_or_off,&bond_off,1,MPI_INT,MPI_MAX,world);
  on_or_off = angle_off;
  MPI_Allreduce(&on_or_off,&angle_off,1,MPI_INT,MPI_MAX,world);
  on_or_off = dihedral_off;
  MPI_Allreduce(&on_or_off,&dihedral_off,1,MPI_INT,MPI_MAX,world);
  on_or_off = improper_off;
  MPI_Allreduce(&on_or_off,&improper_off,1,MPI_INT,MPI_MAX,world);

  // set ptrs to topology build functions

  if (atom->molecular == 2) bond_build_kk = &NeighBondKokkos<DeviceType>::bond_template;
  else if (bond_off) bond_build_kk = &NeighBondKokkos<DeviceType>::bond_partial;
  else bond_build_kk = &NeighBondKokkos<DeviceType>::bond_all;

  if (atom->molecular == 2) angle_build_kk = &NeighBondKokkos<DeviceType>::angle_template;
  else if (angle_off) angle_build_kk = &NeighBondKokkos<DeviceType>::angle_partial;
  else angle_build_kk = &NeighBondKokkos<DeviceType>::angle_all;

  if (atom->molecular == 2) dihedral_build_kk = &NeighBondKokkos<DeviceType>::dihedral_template;
  else if (dihedral_off) dihedral_build_kk = &NeighBondKokkos<DeviceType>::dihedral_partial;
  else dihedral_build_kk = &NeighBondKokkos<DeviceType>::dihedral_all;

  if (atom->molecular == 2) improper_build_kk = &NeighBondKokkos<DeviceType>::improper_template;
  else if (improper_off) improper_build_kk = &NeighBondKokkos<DeviceType>::improper_partial;
  else improper_build_kk = &NeighBondKokkos<DeviceType>::improper_all;

  // set topology neighbor list counts to 0
  // in case all are turned off but potential is still defined

  neighbor->nbondlist = neighbor->nanglelist = neighbor->ndihedrallist = neighbor->nimproperlist = 0;
}

/* ----------------------------------------------------------------------
   build all topology neighbor lists every few timesteps
   normally built with pair lists, but USER-CUDA separates them
------------------------------------------------------------------------- */

template<class DeviceType>
void NeighBondKokkos<DeviceType>::build_topology_kk()
{
  atomKK->sync(execution_space, X_MASK | TAG_MASK);
  int nall = atom->nlocal + atom->nghost;

  nlocal = atom->nlocal;
  x = atomKK->k_x.view<DeviceType>();
  tag = atomKK->k_tag.view<DeviceType>();
  newton_bond = force->newton_bond;

  lostbond = output->thermo->lostbond;

  // don't yet have atom_map_kokkos routines, so move data from host to device

  if (atom->map_style != 1)
    error->all(FLERR,"Must use atom map style array with Kokkos");

  int* map_array_host = atom->get_map_array();
  int map_size = atom->get_map_size();
  k_map_array = DAT::tdual_int_1d("NeighBond:map_array",map_size);
  for (int i=0; i<map_size; i++)
    k_map_array.h_view[i] = map_array_host[i];
  k_map_array.template modify<LMPHostType>();
  k_map_array.template sync<DeviceType>();
  map_array = k_map_array.view<DeviceType>();

  int* sametag_host = atomKK->sametag;
  k_sametag = DAT::tdual_int_1d("NeighBond:sametag",nall);
  for (int i=0; i<nall; i++)
    k_sametag.h_view[i] = sametag_host[i];
  k_sametag.template modify<LMPHostType>();
  k_sametag.template sync<DeviceType>();
  sametag = k_sametag.view<DeviceType>();

  if (force->bond) (this->*bond_build_kk)();
  if (force->angle) (this->*angle_build_kk)();
  if (force->dihedral) (this->*dihedral_build_kk)();
  if (force->improper) (this->*improper_build_kk)();
}

/* ---------------------------------------------------------------------- */

// bondlist, anglelist, dihedrallist, improperlist
//   no longer store map_array() of the bond partners
// instead store domain->closest_image() of the bond partners of atom I
// this enables distances between list atoms to be calculated
//   w/out invoking domain->minimium_image(), e.g. in bond->compute()

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NeighBondKokkos<DeviceType>::bond_all()
{
  atomKK->sync(execution_space, BOND_MASK);
  v_bondlist = k_bondlist.view<DeviceType>();
  num_bond = atomKK->k_num_bond.view<DeviceType>();
  bond_atom = atomKK->k_bond_atom.view<DeviceType>();
  bond_type = atomKK->k_bond_type.view<DeviceType>();

  // Cannot grow a Kokkos view in a parallel loop, so
  //  if the capacity of the list is exceeded, count the size
  //  needed, reallocate on the host, and then
  //  repeat the parallel loop again

  do {
    nmissing = 0;

    h_nlist() = 0;
    k_nlist.template modify<LMPHostType>();
    k_nlist.template sync<DeviceType>();

    h_fail_flag() = 0;
    k_fail_flag.template modify<LMPHostType>();
    k_fail_flag.template sync<DeviceType>();

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagNeighBondBondAll>(0,nlocal),*this,nmissing);

    k_nlist.template modify<DeviceType>();
    k_nlist.template sync<LMPHostType>();
    neighbor->nbondlist = h_nlist();

    k_fail_flag.template modify<DeviceType>();
    k_fail_flag.template sync<LMPHostType>();
    if (h_fail_flag()) {
      maxbond = neighbor->nbondlist + BONDDELTA;
      memoryKK->grow_kokkos(k_bondlist,neighbor->bondlist,maxbond,3,"neighbor:neighbor->bondlist");
      v_bondlist = k_bondlist.view<DeviceType>();
    }
  } while (h_fail_flag());

  if (nmissing && lostbond == Thermo::ERROR) {
    char str[128];
    sprintf(str,"Bond atoms missing on proc %d at step " BIGINT_FORMAT,
            me,update->ntimestep);
    error->one(FLERR,str);
  }

  if (neighbor->cluster_check) bond_check();
  if (lostbond == Thermo::IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Bond atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }

  k_bondlist.modify<DeviceType>();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NeighBondKokkos<DeviceType>::operator()(TagNeighBondBondAll, const int &i, int &nmissing) const {
  for (int m = 0; m < num_bond[i]; m++) {
    int atom1 = map_array(bond_atom(i,m));
    if (atom1 == -1) {
      nmissing++;
      if (lostbond == Thermo::ERROR) return;
      continue;
    }
    atom1 = closest_image(i,atom1);
    if (newton_bond || i < atom1) {
      const int nbondlist = Kokkos::atomic_fetch_add(&d_nlist(),1);
      if (nbondlist >= maxbond && !d_fail_flag())
        Kokkos::atomic_fetch_add(&d_fail_flag(),1);
      if (d_fail_flag()) continue;
      v_bondlist(nbondlist,0) = i;
      v_bondlist(nbondlist,1) = atom1;
      v_bondlist(nbondlist,2) = bond_type(i,m);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NeighBondKokkos<DeviceType>::bond_template()
{
  error->all(FLERR,"Cannot (yet) use molecular templates with Kokkos");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NeighBondKokkos<DeviceType>::bond_partial()
{
  atomKK->sync(execution_space, BOND_MASK);
  v_bondlist = k_bondlist.view<DeviceType>();
  num_bond = atomKK->k_num_bond.view<DeviceType>();
  bond_atom = atomKK->k_bond_atom.view<DeviceType>();
  bond_type = atomKK->k_bond_type.view<DeviceType>();

  // Cannot grow a Kokkos view in a parallel loop, so
  //  if the capacity of the list is exceeded, count the size
  //  needed, reallocate on the host, and then
  //  repeat the parallel loop again

  do {
    nmissing = 0;

    h_nlist() = 0;
    k_nlist.template modify<LMPHostType>();
    k_nlist.template sync<DeviceType>();

    h_fail_flag() = 0;
    k_fail_flag.template modify<LMPHostType>();
    k_fail_flag.template sync<DeviceType>();

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagNeighBondBondPartial>(0,nlocal),*this,nmissing);

    k_nlist.template modify<DeviceType>();
    k_nlist.template sync<LMPHostType>();
    neighbor->nbondlist = h_nlist();

    k_fail_flag.template modify<DeviceType>();
    k_fail_flag.template sync<LMPHostType>();
    if (h_fail_flag()) {
      maxbond = neighbor->nbondlist + BONDDELTA;
      memoryKK->grow_kokkos(k_bondlist,neighbor->bondlist,maxbond,3,"neighbor:neighbor->bondlist");
      v_bondlist = k_bondlist.view<DeviceType>();
    }
  } while (h_fail_flag());

  if (nmissing && lostbond == Thermo::ERROR) {
    char str[128];
    sprintf(str,"Bond atoms missing on proc %d at step " BIGINT_FORMAT,
            me,update->ntimestep);
    error->one(FLERR,str);
  }

  if (neighbor->cluster_check) bond_check();
  if (lostbond == Thermo::IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Bond atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }

  k_bondlist.modify<DeviceType>();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NeighBondKokkos<DeviceType>::operator()(TagNeighBondBondPartial, const int &i, int &nmissing) const {
  for (int m = 0; m < num_bond[i]; m++) {
    if (bond_type(i,m) <= 0) continue;
    int atom1 = map_array(bond_atom(i,m));
    if (atom1 == -1) {
      nmissing++;
      if (lostbond == Thermo::ERROR) return;
      continue;
    }
    atom1 = closest_image(i,atom1);
    if (newton_bond || i < atom1) {
      const int nbondlist = Kokkos::atomic_fetch_add(&d_nlist(),1);
      if (nbondlist >= maxbond && !d_fail_flag())
        Kokkos::atomic_fetch_add(&d_fail_flag(),1);
      if (d_fail_flag()) continue;
      v_bondlist(nbondlist,0) = i;
      v_bondlist(nbondlist,1) = atom1;
      v_bondlist(nbondlist,2) = bond_type(i,m);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NeighBondKokkos<DeviceType>::bond_check()
{
  int flag = 0;

  update_domain_variables();
  atomKK->sync(execution_space, X_MASK);
  k_bondlist.sync<DeviceType>();

  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagNeighBondBondCheck>(0,neighbor->nbondlist),*this,flag);

  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) error->all(FLERR,"Bond extent > half of periodic box length");
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NeighBondKokkos<DeviceType>::operator()(TagNeighBondBondCheck, const int &m, int &flag) const {
  const int i = v_bondlist(m,0);
  const int j = v_bondlist(m,1);
  X_FLOAT dxstart,dystart,dzstart;
  X_FLOAT dx,dy,dz;
  dxstart = dx = x(i,0) - x(j,0);
  dystart = dy = x(i,1) - x(j,1);
  dzstart = dz = x(i,2) - x(j,2);
  minimum_image(dx,dy,dz);
  if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NeighBondKokkos<DeviceType>::angle_all()
{
  atomKK->sync(execution_space, ANGLE_MASK);
  v_anglelist = k_anglelist.view<DeviceType>();
  num_angle = atomKK->k_num_angle.view<DeviceType>();
  angle_atom1 = atomKK->k_angle_atom1.view<DeviceType>();
  angle_atom2 = atomKK->k_angle_atom2.view<DeviceType>();
  angle_atom3 = atomKK->k_angle_atom3.view<DeviceType>();
  angle_type = atomKK->k_angle_type.view<DeviceType>();

  // Cannot grow a Kokkos view in a parallel loop, so
  //  if the capacity of the list is exceeded, count the size
  //  needed, reallocate on the host, and then
  //  repeat the parallel loop again

  do {
    nmissing = 0;

    h_nlist() = 0;
    k_nlist.template modify<LMPHostType>();
    k_nlist.template sync<DeviceType>();

    h_fail_flag() = 0;
    k_fail_flag.template modify<LMPHostType>();
    k_fail_flag.template sync<DeviceType>();

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagNeighBondAngleAll>(0,nlocal),*this,nmissing);

    k_nlist.template modify<DeviceType>();
    k_nlist.template sync<LMPHostType>();
    neighbor->nanglelist = h_nlist();

    k_fail_flag.template modify<DeviceType>();
    k_fail_flag.template sync<LMPHostType>();
    if (h_fail_flag()) {
      maxangle = neighbor->nanglelist + BONDDELTA;
      memoryKK->grow_kokkos(k_anglelist,neighbor->anglelist,maxangle,4,"neighbor:neighbor->anglelist");
      v_anglelist = k_anglelist.view<DeviceType>();
    }
  } while (h_fail_flag());

  if (nmissing && lostbond == Thermo::ERROR) {
    char str[128];
    sprintf(str,"Angle atoms missing on proc %d at step " BIGINT_FORMAT,
            me,update->ntimestep);
    error->one(FLERR,str);
  }

  if (neighbor->cluster_check) angle_check();
  if (lostbond == Thermo::IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Angle atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }

  k_anglelist.modify<DeviceType>();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NeighBondKokkos<DeviceType>::operator()(TagNeighBondAngleAll, const int &i, int &nmissing) const {
  for (int m = 0; m < num_angle[i]; m++) {
    int atom1 = map_array(angle_atom1(i,m));
    int atom2 = map_array(angle_atom2(i,m));
    int atom3 = map_array(angle_atom3(i,m));
    if (atom1 == -1 || atom2 == -1 || atom3 == -1) {
      nmissing++;
      if (lostbond == Thermo::ERROR) return;
      continue;
    }
    atom1 = closest_image(i,atom1);
    atom2 = closest_image(i,atom2);
    atom3 = closest_image(i,atom3);
    if (newton_bond || (i <= atom1 && i <= atom2 && i <= atom3)) {
      const int nanglelist = Kokkos::atomic_fetch_add(&d_nlist(),1);
      if (nanglelist >= maxangle && !d_fail_flag())
        Kokkos::atomic_fetch_add(&d_fail_flag(),1);
      if (d_fail_flag()) continue;
      v_anglelist(nanglelist,0) = atom1;
      v_anglelist(nanglelist,1) = atom2;
      v_anglelist(nanglelist,2) = atom3;
      v_anglelist(nanglelist,3) = angle_type(i,m);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NeighBondKokkos<DeviceType>::angle_template()
{
  error->all(FLERR,"Cannot (yet) use molecular templates with Kokkos");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NeighBondKokkos<DeviceType>::angle_partial()
{
  atomKK->sync(execution_space, ANGLE_MASK);
  v_anglelist = k_anglelist.view<DeviceType>();
  num_angle = atomKK->k_num_angle.view<DeviceType>();
  angle_atom1 = atomKK->k_angle_atom1.view<DeviceType>();
  angle_atom2 = atomKK->k_angle_atom2.view<DeviceType>();
  angle_atom3 = atomKK->k_angle_atom3.view<DeviceType>();
  angle_type = atomKK->k_angle_type.view<DeviceType>();

  // Cannot grow a Kokkos view in a parallel loop, so
  //  if the capacity of the list is exceeded, count the size
  //  needed, reallocate on the host, and then
  //  repeat the parallel loop again

  do {
    nmissing = 0;

    h_nlist() = 0;
    k_nlist.template modify<LMPHostType>();
    k_nlist.template sync<DeviceType>();

    h_fail_flag() = 0;
    k_fail_flag.template modify<LMPHostType>();
    k_fail_flag.template sync<DeviceType>();

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagNeighBondAnglePartial>(0,nlocal),*this,nmissing);

    k_nlist.template modify<DeviceType>();
    k_nlist.template sync<LMPHostType>();
    neighbor->nanglelist = h_nlist();

    k_fail_flag.template modify<DeviceType>();
    k_fail_flag.template sync<LMPHostType>();
    if (h_fail_flag()) {
      maxangle = neighbor->nanglelist + BONDDELTA;
      memoryKK->grow_kokkos(k_anglelist,neighbor->anglelist,maxangle,4,"neighbor:neighbor->anglelist");
      v_anglelist = k_anglelist.view<DeviceType>();
    }
  } while (h_fail_flag());

  if (nmissing && lostbond == Thermo::ERROR) {
    char str[128];
    sprintf(str,"Angle atoms missing on proc %d at step " BIGINT_FORMAT,
            me,update->ntimestep);
    error->one(FLERR,str);
  }

  if (neighbor->cluster_check) angle_check();
  if (lostbond == Thermo::IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Angle atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }

  k_anglelist.modify<DeviceType>();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NeighBondKokkos<DeviceType>::operator()(TagNeighBondAnglePartial, const int &i, int &nmissing) const {
  for (int m = 0; m < num_angle[i]; m++) {
    if (angle_type(i,m) <= 0) continue;
    int atom1 = map_array(angle_atom1(i,m));
    int atom2 = map_array(angle_atom2(i,m));
    int atom3 = map_array(angle_atom3(i,m));
    if (atom1 == -1 || atom2 == -1 || atom3 == -1) {
      nmissing++;
      if (lostbond == Thermo::ERROR) return;
      continue;
    }
    atom1 = closest_image(i,atom1);
    atom2 = closest_image(i,atom2);
    atom3 = closest_image(i,atom3);
    if (newton_bond || (i <= atom1 && i <= atom2 && i <= atom3)) {
      const int nanglelist = Kokkos::atomic_fetch_add(&d_nlist(),1);
      if (nanglelist >= maxangle && !d_fail_flag())
        Kokkos::atomic_fetch_add(&d_fail_flag(),1);
      if (d_fail_flag()) continue;
      v_anglelist(nanglelist,0) = atom1;
      v_anglelist(nanglelist,1) = atom2;
      v_anglelist(nanglelist,2) = atom3;
      v_anglelist(nanglelist,3) = angle_type(i,m);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NeighBondKokkos<DeviceType>::angle_check()
{
  int flag = 0;

  // check all 3 distances
  // in case angle potential computes any of them

  update_domain_variables();
  atomKK->sync(execution_space, X_MASK);
  k_anglelist.sync<DeviceType>();

  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagNeighBondAngleCheck>(0,neighbor->nanglelist),*this,flag);

  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) error->all(FLERR,"Angle extent > half of periodic box length");
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NeighBondKokkos<DeviceType>::operator()(TagNeighBondAngleCheck, const int &m, int &flag) const {
  const int i = v_anglelist(m,0);
  const int j = v_anglelist(m,1);
  const int k = v_anglelist(m,2);
  X_FLOAT dxstart,dystart,dzstart;
  X_FLOAT dx,dy,dz;
  dxstart = dx = x(i,0) - x(j,0);
  dystart = dy = x(i,1) - x(j,1);
  dzstart = dz = x(i,2) - x(j,2);
  minimum_image(dx,dy,dz);
  if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
  dxstart = dx = x(i,0) - x(k,0);
  dystart = dy = x(i,1) - x(k,1);
  dzstart = dz = x(i,2) - x(k,2);
  minimum_image(dx,dy,dz);
  if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
  dxstart = dx = x(j,0) - x(k,0);
  dystart = dy = x(j,1) - x(k,1);
  dzstart = dz = x(j,2) - x(k,2);
  minimum_image(dx,dy,dz);
  if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NeighBondKokkos<DeviceType>::dihedral_all()
{
  atomKK->sync(execution_space, DIHEDRAL_MASK);
  v_dihedrallist = k_dihedrallist.view<DeviceType>();
  num_dihedral = atomKK->k_num_dihedral.view<DeviceType>();
  dihedral_atom1 = atomKK->k_dihedral_atom1.view<DeviceType>();
  dihedral_atom2 = atomKK->k_dihedral_atom2.view<DeviceType>();
  dihedral_atom3 = atomKK->k_dihedral_atom3.view<DeviceType>();
  dihedral_atom4 = atomKK->k_dihedral_atom4.view<DeviceType>();
  dihedral_type = atomKK->k_dihedral_type.view<DeviceType>();

  // Cannot grow a Kokkos view in a parallel loop, so
  //  if the capacity of the list is exceeded, count the size
  //  needed, reallocate on the host, and then
  //  repeat the parallel loop again

  do {
    nmissing = 0;

    h_nlist() = 0;
    k_nlist.template modify<LMPHostType>();
    k_nlist.template sync<DeviceType>();

    h_fail_flag() = 0;
    k_fail_flag.template modify<LMPHostType>();
    k_fail_flag.template sync<DeviceType>();

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagNeighBondDihedralAll>(0,nlocal),*this,nmissing);

    k_nlist.template modify<DeviceType>();
    k_nlist.template sync<LMPHostType>();
    neighbor->ndihedrallist = h_nlist();

    k_fail_flag.template modify<DeviceType>();
    k_fail_flag.template sync<LMPHostType>();
    if (h_fail_flag()) {
      maxdihedral = neighbor->ndihedrallist + BONDDELTA;
      memoryKK->grow_kokkos(k_dihedrallist,neighbor->dihedrallist,maxdihedral,5,"neighbor:neighbor->dihedrallist");
      v_dihedrallist = k_dihedrallist.view<DeviceType>();
    }
  } while (h_fail_flag());

  if (nmissing && lostbond == Thermo::ERROR) {
    char str[128];
    sprintf(str,"Dihedral atoms missing on proc %d at step " BIGINT_FORMAT,
            me,update->ntimestep);
    error->one(FLERR,str);
  }

  if (neighbor->cluster_check) dihedral_check(neighbor->ndihedrallist,v_dihedrallist);
  if (lostbond == Thermo::IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Dihedral atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }

  k_dihedrallist.modify<DeviceType>();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NeighBondKokkos<DeviceType>::operator()(TagNeighBondDihedralAll, const int &i, int &nmissing) const {
  for (int m = 0; m < num_dihedral[i]; m++) {
    int atom1 = map_array(dihedral_atom1(i,m));
    int atom2 = map_array(dihedral_atom2(i,m));
    int atom3 = map_array(dihedral_atom3(i,m));
    int atom4 = map_array(dihedral_atom4(i,m));
    if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
      nmissing++;
      if (lostbond == Thermo::ERROR) return;
      continue;
    }
    atom1 = closest_image(i,atom1);
    atom2 = closest_image(i,atom2);
    atom3 = closest_image(i,atom3);
    atom4 = closest_image(i,atom4);
    if (newton_bond ||
        (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4)) {
      const int ndihedrallist = Kokkos::atomic_fetch_add(&d_nlist(),1);
      if (ndihedrallist >= maxdihedral && !d_fail_flag())
        Kokkos::atomic_fetch_add(&d_fail_flag(),1);
      if (d_fail_flag()) continue;
      v_dihedrallist(ndihedrallist,0) = atom1;
      v_dihedrallist(ndihedrallist,1) = atom2;
      v_dihedrallist(ndihedrallist,2) = atom3;
      v_dihedrallist(ndihedrallist,3) = atom4;
      v_dihedrallist(ndihedrallist,4) = dihedral_type(i,m);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NeighBondKokkos<DeviceType>::dihedral_template()
{
  error->all(FLERR,"Cannot (yet) use molecular templates with Kokkos");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NeighBondKokkos<DeviceType>::dihedral_partial()
{
  atomKK->sync(execution_space, DIHEDRAL_MASK);
  v_dihedrallist = k_dihedrallist.view<DeviceType>();
  num_dihedral = atomKK->k_num_dihedral.view<DeviceType>();
  dihedral_atom1 = atomKK->k_dihedral_atom1.view<DeviceType>();
  dihedral_atom2 = atomKK->k_dihedral_atom2.view<DeviceType>();
  dihedral_atom3 = atomKK->k_dihedral_atom3.view<DeviceType>();
  dihedral_atom4 = atomKK->k_dihedral_atom4.view<DeviceType>();
  dihedral_type = atomKK->k_dihedral_type.view<DeviceType>();

  // Cannot grow a Kokkos view in a parallel loop, so
  //  if the capacity of the list is exceeded, count the size
  //  needed, reallocate on the host, and then
  //  repeat the parallel loop again

  do {
    nmissing = 0;

    h_nlist() = 0;
    k_nlist.template modify<LMPHostType>();
    k_nlist.template sync<DeviceType>();

    h_fail_flag() = 0;
    k_fail_flag.template modify<LMPHostType>();
    k_fail_flag.template sync<DeviceType>();

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagNeighBondDihedralPartial>(0,nlocal),*this,nmissing);

    k_nlist.template modify<DeviceType>();
    k_nlist.template sync<LMPHostType>();
    neighbor->ndihedrallist = h_nlist();

    k_fail_flag.template modify<DeviceType>();
    k_fail_flag.template sync<LMPHostType>();
    if (h_fail_flag()) {
      maxdihedral = neighbor->ndihedrallist + BONDDELTA;
      memoryKK->grow_kokkos(k_dihedrallist,neighbor->dihedrallist,maxdihedral,5,"neighbor:neighbor->dihedrallist");
      v_dihedrallist = k_dihedrallist.view<DeviceType>();
    }
  } while (h_fail_flag());

  if (nmissing && lostbond == Thermo::ERROR) {
    char str[128];
    sprintf(str,"Dihedral atoms missing on proc %d at step " BIGINT_FORMAT,
            me,update->ntimestep);
    error->one(FLERR,str);
  }

  if (neighbor->cluster_check) dihedral_check(neighbor->ndihedrallist,v_dihedrallist);
  if (lostbond == Thermo::IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Dihedral atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }

  k_dihedrallist.modify<DeviceType>();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NeighBondKokkos<DeviceType>::operator()(TagNeighBondDihedralPartial, const int &i, int &nmissing) const {
  for (int m = 0; m < num_dihedral[i]; m++) {
    if (dihedral_type(i,m) <= 0) continue;
    int atom1 = map_array(dihedral_atom1(i,m));
    int atom2 = map_array(dihedral_atom2(i,m));
    int atom3 = map_array(dihedral_atom3(i,m));
    int atom4 = map_array(dihedral_atom4(i,m));
    if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
      nmissing++;
      if (lostbond == Thermo::ERROR) return;
      continue;
    }
    atom1 = closest_image(i,atom1);
    atom2 = closest_image(i,atom2);
    atom3 = closest_image(i,atom3);
    atom4 = closest_image(i,atom4);
    if (newton_bond ||
        (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4)) {
      const int ndihedrallist = Kokkos::atomic_fetch_add(&d_nlist(),1);
      if (ndihedrallist >= maxdihedral && !d_fail_flag())
        Kokkos::atomic_fetch_add(&d_fail_flag(),1);
      if (d_fail_flag()) continue;
      v_dihedrallist(ndihedrallist,0) = atom1;
      v_dihedrallist(ndihedrallist,1) = atom2;
      v_dihedrallist(ndihedrallist,2) = atom3;
      v_dihedrallist(ndihedrallist,3) = atom4;
      v_dihedrallist(ndihedrallist,4) = dihedral_type(i,m);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NeighBondKokkos<DeviceType>::dihedral_check(int nlist, typename AT::t_int_2d list_in)
{
  list = list_in;
  int flag = 0;

  // check all 6 distances
  // in case dihedral/improper potential computes any of them

  update_domain_variables();
  atomKK->sync(execution_space, X_MASK);
  k_dihedrallist.sync<DeviceType>();

  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagNeighBondDihedralCheck>(0,nlist),*this,flag);

  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all)
    error->all(FLERR,"Dihedral/improper extent > half of periodic box length");
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NeighBondKokkos<DeviceType>::operator()(TagNeighBondDihedralCheck, const int &m, int &flag) const {
  const int i = list(m,0);
  const int j = list(m,1);
  const int k = list(m,2);
  const int l = list(m,3);
  X_FLOAT dxstart,dystart,dzstart;
  X_FLOAT dx,dy,dz;
  dxstart = dx = x(i,0) - x(j,0);
  dystart = dy = x(i,1) - x(j,1);
  dzstart = dz = x(i,2) - x(j,2);
  minimum_image(dx,dy,dz);
  if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
  dxstart = dx = x(i,0) - x(k,0);
  dystart = dy = x(i,1) - x(k,1);
  dzstart = dz = x(i,2) - x(k,2);
  minimum_image(dx,dy,dz);
  if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
  dxstart = dx = x(i,0) - x(l,0);
  dystart = dy = x(i,1) - x(l,1);
  dzstart = dz = x(i,2) - x(l,2);
  minimum_image(dx,dy,dz);
  if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
  dxstart = dx = x(j,0) - x(k,0);
  dystart = dy = x(j,1) - x(k,1);
  dzstart = dz = x(j,2) - x(k,2);
  minimum_image(dx,dy,dz);
  if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
  dxstart = dx = x(j,0) - x(l,0);
  dystart = dy = x(j,1) - x(l,1);
  dzstart = dz = x(j,2) - x(l,2);
  minimum_image(dx,dy,dz);
  if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
  dxstart = dx = x(k,0) - x(l,0);
  dystart = dy = x(k,1) - x(l,1);
  dzstart = dz = x(k,2) - x(l,2);
  minimum_image(dx,dy,dz);
  if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NeighBondKokkos<DeviceType>::improper_all()
{
  atomKK->sync(execution_space, IMPROPER_MASK);
  v_improperlist = k_improperlist.view<DeviceType>();
  num_improper = atomKK->k_num_improper.view<DeviceType>();
  improper_atom1 = atomKK->k_improper_atom1.view<DeviceType>();
  improper_atom2 = atomKK->k_improper_atom2.view<DeviceType>();
  improper_atom3 = atomKK->k_improper_atom3.view<DeviceType>();
  improper_atom4 = atomKK->k_improper_atom4.view<DeviceType>();
  improper_type = atomKK->k_improper_type.view<DeviceType>();

  // Cannot grow a Kokkos view in a parallel loop, so
  //  if the capacity of the list is exceeded, count the size
  //  needed, reallocate on the host, and then
  //  repeat the parallel loop again

  do {
    nmissing = 0;

    h_nlist() = 0;
    k_nlist.template modify<LMPHostType>();
    k_nlist.template sync<DeviceType>();

    h_fail_flag() = 0;
    k_fail_flag.template modify<LMPHostType>();
    k_fail_flag.template sync<DeviceType>();

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagNeighBondImproperAll>(0,nlocal),*this,nmissing);

    k_nlist.template modify<DeviceType>();
    k_nlist.template sync<LMPHostType>();
    neighbor->nimproperlist = h_nlist();

    k_fail_flag.template modify<DeviceType>();
    k_fail_flag.template sync<LMPHostType>();
    if (h_fail_flag()) {
      maximproper = neighbor->nimproperlist + BONDDELTA;
      memoryKK->grow_kokkos(k_improperlist,neighbor->improperlist,maximproper,5,"neighbor:neighbor->improperlist");
      v_improperlist = k_improperlist.view<DeviceType>();
    }
  } while (h_fail_flag());

  if (nmissing && lostbond == Thermo::ERROR) {
    char str[128];
    sprintf(str,"Improper atoms missing on proc %d at step " BIGINT_FORMAT,
            me,update->ntimestep);
    error->one(FLERR,str);
  }

  if (neighbor->cluster_check) dihedral_check(neighbor->nimproperlist,v_improperlist);
  if (lostbond == Thermo::IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Improper atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }

  k_improperlist.modify<DeviceType>();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NeighBondKokkos<DeviceType>::operator()(TagNeighBondImproperAll, const int &i, int &nmissing) const {
  for (int m = 0; m < num_improper[i]; m++) {
    int atom1 = map_array(improper_atom1(i,m));
    int atom2 = map_array(improper_atom2(i,m));
    int atom3 = map_array(improper_atom3(i,m));
    int atom4 = map_array(improper_atom4(i,m));
    if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
      nmissing++;
      if (lostbond == Thermo::ERROR) return;
      continue;
    }
    atom1 = closest_image(i,atom1);
    atom2 = closest_image(i,atom2);
    atom3 = closest_image(i,atom3);
    atom4 = closest_image(i,atom4);
    if (newton_bond ||
        (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4)) {
      const int nimproperlist = Kokkos::atomic_fetch_add(&d_nlist(),1);
      if (nimproperlist >= maximproper && !d_fail_flag())
        Kokkos::atomic_fetch_add(&d_fail_flag(),1);
      if (d_fail_flag()) continue;
      v_improperlist(nimproperlist,0) = atom1;
      v_improperlist(nimproperlist,1) = atom2;
      v_improperlist(nimproperlist,2) = atom3;
      v_improperlist(nimproperlist,3) = atom4;
      v_improperlist(nimproperlist,4) = improper_type(i,m);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NeighBondKokkos<DeviceType>::improper_template()
{
  error->all(FLERR,"Cannot (yet) use molecular templates with Kokkos");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NeighBondKokkos<DeviceType>::improper_partial()
{
  atomKK->sync(execution_space, IMPROPER_MASK);
  v_improperlist = k_improperlist.view<DeviceType>();
  num_improper = atomKK->k_num_improper.view<DeviceType>();
  improper_atom1 = atomKK->k_improper_atom1.view<DeviceType>();
  improper_atom2 = atomKK->k_improper_atom2.view<DeviceType>();
  improper_atom3 = atomKK->k_improper_atom3.view<DeviceType>();
  improper_atom4 = atomKK->k_improper_atom4.view<DeviceType>();
  improper_type = atomKK->k_improper_type.view<DeviceType>();

  // Cannot grow a Kokkos view in a parallel loop, so
  //  if the capacity of the list is exceeded, count the size
  //  needed, reallocate on the host, and then
  //  repeat the parallel loop again

  do {
    nmissing = 0;

    h_nlist() = 0;
    k_nlist.template modify<LMPHostType>();
    k_nlist.template sync<DeviceType>();

    h_fail_flag() = 0;
    k_fail_flag.template modify<LMPHostType>();
    k_fail_flag.template sync<DeviceType>();

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagNeighBondImproperPartial>(0,nlocal),*this,nmissing);

    k_nlist.template modify<DeviceType>();
    k_nlist.template sync<LMPHostType>();
    neighbor->nimproperlist = h_nlist();

    k_fail_flag.template modify<DeviceType>();
    k_fail_flag.template sync<LMPHostType>();
    if (h_fail_flag()) {
      maximproper = neighbor->nimproperlist + BONDDELTA;
      memoryKK->grow_kokkos(k_improperlist,neighbor->improperlist,maximproper,5,"neighbor:neighbor->improperlist");
      v_improperlist = k_improperlist.view<DeviceType>();
    }
  } while (h_fail_flag());

  if (nmissing && lostbond == Thermo::ERROR) {
    char str[128];
    sprintf(str,"Improper atoms missing on proc %d at step " BIGINT_FORMAT,
            me,update->ntimestep);
    error->one(FLERR,str);
  }

  if (neighbor->cluster_check) dihedral_check(neighbor->nimproperlist,v_improperlist);
  if (lostbond == Thermo::IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Improper atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }

  k_improperlist.modify<DeviceType>();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NeighBondKokkos<DeviceType>::operator()(TagNeighBondImproperPartial, const int &i, int &nmissing) const {
  for (int m = 0; m < num_improper[i]; m++) {
    if (improper_type(i,m) <= 0) continue;
    int atom1 = map_array(improper_atom1(i,m));
    int atom2 = map_array(improper_atom2(i,m));
    int atom3 = map_array(improper_atom3(i,m));
    int atom4 = map_array(improper_atom4(i,m));
    if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
      nmissing++;
      if (lostbond == Thermo::ERROR) return;
      continue;
    }
    atom1 = closest_image(i,atom1);
    atom2 = closest_image(i,atom2);
    atom3 = closest_image(i,atom3);
    atom4 = closest_image(i,atom4);
    if (newton_bond ||
        (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4)) {
      const int nimproperlist = Kokkos::atomic_fetch_add(&d_nlist(),1);
      if (nimproperlist >= maximproper && !d_fail_flag())
        Kokkos::atomic_fetch_add(&d_fail_flag(),1);
      if (d_fail_flag()) continue;
      v_improperlist(nimproperlist,0) = atom1;
      v_improperlist(nimproperlist,1) = atom2;
      v_improperlist(nimproperlist,2) = atom3;
      v_improperlist(nimproperlist,3) = atom4;
      v_improperlist(nimproperlist,4) = improper_type(i,m);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int NeighBondKokkos<DeviceType>::closest_image(const int i, int j) const
{
  if (j < 0) return j;

  const X_FLOAT xi0 = x(i,0);
  const X_FLOAT xi1 = x(i,1);
  const X_FLOAT xi2 = x(i,2);

  int closest = j;
  X_FLOAT delx = xi0 - x(j,0);
  X_FLOAT dely = xi1 - x(j,1);
  X_FLOAT delz = xi2 - x(j,2);
  X_FLOAT rsqmin = delx*delx + dely*dely + delz*delz;
  X_FLOAT rsq;

  while (sametag[j] >= 0) {
    j = sametag[j];
    delx = xi0 - x(j,0);
    dely = xi1 - x(j,1);
    delz = xi2 - x(j,2);
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq < rsqmin) {
      rsqmin = rsq;
      closest = j;
    }
  }
  return closest;
}

/* ----------------------------------------------------------------------
   minimum image convention
   use 1/2 of box size as test
   for triclinic, also add/subtract tilt factors in other dims as needed
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NeighBondKokkos<DeviceType>::minimum_image(X_FLOAT &dx, X_FLOAT &dy, X_FLOAT &dz) const
{
  if (triclinic == 0) {
    if (xperiodic) {
      if (fabs(dx) > xprd_half) {
        if (dx < 0.0) dx += xprd;
        else dx -= xprd;
      }
    }
    if (yperiodic) {
      if (fabs(dy) > yprd_half) {
        if (dy < 0.0) dy += yprd;
        else dy -= yprd;
      }
    }
    if (zperiodic) {
      if (fabs(dz) > zprd_half) {
        if (dz < 0.0) dz += zprd;
        else dz -= zprd;
      }
    }

  } else {
    if (zperiodic) {
      if (fabs(dz) > zprd_half) {
        if (dz < 0.0) {
          dz += zprd;
          dy += yz;
          dx += xz;
        } else {
          dz -= zprd;
          dy -= yz;
          dx -= xz;
        }
      }
    }
    if (yperiodic) {
      if (fabs(dy) > yprd_half) {
        if (dy < 0.0) {
          dy += yprd;
          dx += xy;
        } else {
          dy -= yprd;
          dx -= xy;
        }
      }
    }
    if (xperiodic) {
      if (fabs(dx) > xprd_half) {
        if (dx < 0.0) dx += xprd;
        else dx -= xprd;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NeighBondKokkos<DeviceType>::update_domain_variables()
{
  triclinic = domain->triclinic;
  xperiodic = domain->xperiodic;
  xprd_half = domain->xprd_half;
  xprd = domain->xprd;
  yperiodic = domain->yperiodic;
  yprd_half = domain->yprd_half;
  yprd = domain->yprd;
  zperiodic = domain->zperiodic;
  zprd_half = domain->zprd_half;
  zprd = domain->zprd;
  xy = domain->xy;
  xz = domain->xz;
  yz = domain->yz;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class NeighBondKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class NeighBondKokkos<LMPHostType>;
#endif
}

