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

#include "neigh_request.h"
#include "atom.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NeighRequest::NeighRequest(LAMMPS *lmp) : Pointers(lmp)
{
  // default ID = 0

  id = 0;

  // class user of list: default is pair request
  // only one is set to 1

  pair = 1;
  fix = compute = command = neigh = 0;

  // kind of list: default is half neighbor list
  // only one is set to 1

  half = 1;
  full = 0;
  cac = 0;

  // attribute flags, mutiple can be set to 1
  // default is every reneighboring, not occasional
  // default is use newton_pair setting in force
  // default is no neighbors of ghosts
  // default is use cutoffs, not size of particles
  // default is no associated neighbor history info in FixNeighHistory
  // default is no one-sided sphere/surface interactions (when size = 1)
  // default is neighbors of atoms, not bonds
  // default is no multilevel rRESPA neighbors
  // default is no OpenMP multi-threaded neighbor list build
  // default is no Intel-specific neighbor list build
  // default is no Kokkos neighbor list build
  // default is no Shardlow Splitting Algorithm (SSA) neighbor list build
  // default is no list-specific cutoff
  // default is no storage of auxiliary floating point values

  occasional = 0;
  newton = 0;
  ghost = 0;
  size = 0;
  history = 0;
  granonesided = 0;
  respainner = respamiddle = respaouter = 0;
  bond = 0;
  omp = 0;
  intel = 0;
  kokkos_host = kokkos_device = 0;
  ssa = 0;
  cut = 0;
  cutoff = 0.0;

  // skip info, default is no skipping

  skip = 0;
  iskip = NULL;
  ijskip = NULL;

  // only set when command = 1;

  command_style = NULL;

  // info set by Neighbor class when morphing original requests

  skiplist = -1;
  off2on = 0;
  copy = 0;
  copylist = -1;
  halffull = 0;
  halffulllist = -1;
  unique = 0;

  // internal settings

  index_bin = index_stencil = index_pair = -1;
}

/* ---------------------------------------------------------------------- */

NeighRequest::~NeighRequest()
{
  delete [] iskip;
  memory->destroy(ijskip);
}

/* ----------------------------------------------------------------------
   compare this request to other request
   identical means requestor identity and all params it sets are the same
   do not check other params that Neighbor can change after requests are made
   return 1 if identical, 0 if not
------------------------------------------------------------------------- */

int NeighRequest::identical(NeighRequest *other)
{
  int same = 1;

  // check for match of requestor_instance and instance counter
  // prevents an old fix from being unfix/refix in same memory location
  //   stored in requestor, and thus appearing old, when really new
  // only needed for classes with persistent neigh lists: Pair, Fix, Compute

  if (requestor != other->requestor) same = 0;
  if (requestor_instance != other->requestor_instance) same = 0;
  if (id != other->id) same = 0;

  // only compare settings made by requestors
  // not settings made later by Neighbor class

  if (pair != other->pair) same = 0;
  if (fix != other->fix) same = 0;
  if (compute != other->compute) same = 0;
  if (command != other->command) same = 0;
  if (neigh != other->neigh) same = 0;

  if (half != other->half) same = 0;
  if (full != other->full) same = 0;

  if (occasional != other->occasional) same = 0;
  if (newton != other->newton) same = 0;
  if (ghost != other->ghost) same = 0;
  if (size != other->size) same = 0;
  if (history != other->history) same = 0;
  if (granonesided != other->granonesided) same = 0;
  if (respainner != other->respainner) same = 0;
  if (respamiddle != other->respamiddle) same = 0;
  if (respaouter != other->respaouter) same = 0;
  if (bond != other->bond) same = 0;
  if (omp != other->omp) same = 0;
  if (intel != other->intel) same = 0;
  if (kokkos_host != other->kokkos_host) same = 0;
  if (kokkos_device != other->kokkos_device) same = 0;
  if (ssa != other->ssa) same = 0;
  if (copy != other->copy) same = 0;
  if (cutoff != other->cutoff) same = 0;

  if (skip != other->skip) same = 0;
  if (same && skip && other->skip) same = same_skip(other);

  return same;
}

/* ----------------------------------------------------------------------
   compare same info for two requests that have skip = 1
   return 1 if identical, 0 if not
------------------------------------------------------------------------- */

int NeighRequest::same_skip(NeighRequest *other)
{
  int i,j;

  int ntypes = atom->ntypes;
  int same = 1;

  for (i = 1; i <= ntypes; i++)
    if (iskip[i] != other->iskip[i]) same = 0;
  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      if (ijskip[i][j] != other->ijskip[i][j]) same = 0;

  return same;
}

/* ----------------------------------------------------------------------
   set params in this request to those of other request
   copy same fields that are checked in identical()
   purpose is to allow comparison of new requests to old requests
   skipflag = 1 to copy skip vector/array
------------------------------------------------------------------------- */

void NeighRequest::copy_request(NeighRequest *other, int skipflag)
{
  requestor = other->requestor;
  requestor_instance = other->requestor_instance;
  id = other->id;

  pair = other->pair;
  fix = other->fix;
  compute = other->compute;
  command = other->command;

  half = other->half;
  full = other->full;
  cac = other->cac;
  
  occasional = other->occasional;
  newton = other->newton;
  ghost = other->ghost;
  size = other->size;
  history = other->history;
  granonesided = other->granonesided;
  respainner =  other->respainner;
  respamiddle = other->respamiddle;
  respaouter = other->respaouter;
  bond = other->bond;
  omp = other->omp;
  intel = other->intel;
  kokkos_host = other->kokkos_host;
  kokkos_device = other->kokkos_device;
  ssa = other->ssa;
  cut = other->cut;
  cutoff = other->cutoff;

  iskip = NULL;
  ijskip = NULL;

  if (!skipflag) return;

  int i,j;
  int ntypes = atom->ntypes;

  if (other->iskip) {
    iskip = new int[ntypes+1];
    for (i = 1; i <= ntypes; i++)
      iskip[i] = other->iskip[i];
  }

  if (other->ijskip) {
    memory->create(ijskip,ntypes+1,ntypes+1,"neigh_request:ijskip");
    for (i = 1; i <= ntypes; i++)
      for (j = 1; j <= ntypes; j++)
        ijskip[i][j] = other->ijskip[i][j];
  }
}

