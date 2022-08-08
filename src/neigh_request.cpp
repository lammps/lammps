/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
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
#include "neighbor.h"

using namespace LAMMPS_NS;
using namespace NeighConst;

/* ---------------------------------------------------------------------- */

NeighRequest::NeighRequest(LAMMPS *_lmp) : Pointers(_lmp)
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

  // attribute flags, multiple can be set to 1
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
  iskip = nullptr;
  ijskip = nullptr;

  // only set when command = 1;

  command_style = nullptr;

  // info set by Neighbor class when morphing original requests

  skiplist = -1;
  off2on = 0;
  copy = 0;
  trim = 0;
  copylist = -1;
  halffull = 0;
  halffulllist = -1;
  unique = 0;

  // internal settings

  index_bin = index_stencil = index_pair = -1;
}

/* ---------------------------------------------------------------------- */

NeighRequest::NeighRequest(LAMMPS *_lmp, int idx, void *ptr, int num) : NeighRequest(_lmp)
{
  index = idx;
  requestor = ptr;
  requestor_instance = num;
}

/* ---------------------------------------------------------------------- */

NeighRequest::NeighRequest(NeighRequest *old) : NeighRequest(old->lmp)
{
  copy_request(old, 1);
}

/* ---------------------------------------------------------------------- */

NeighRequest::~NeighRequest()
{
  delete[] iskip;
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
  const int ntypes = atom->ntypes;
  int same = 1;

  for (int i = 1; i <= ntypes; i++)
    if (iskip[i] != other->iskip[i]) same = 0;
  for (int i = 1; i <= ntypes; i++)
    for (int j = 1; j <= ntypes; j++)
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

  occasional = other->occasional;
  newton = other->newton;
  ghost = other->ghost;
  size = other->size;
  history = other->history;
  granonesided = other->granonesided;
  respainner = other->respainner;
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

  iskip = nullptr;
  ijskip = nullptr;

  if (!skipflag) return;

  int i, j;
  int ntp1 = atom->ntypes + 1;

  skip = other->skip;

  if (other->iskip) {
    iskip = new int[ntp1];
    for (i = 1; i < ntp1; i++) iskip[i] = other->iskip[i];
  }

  if (other->ijskip) {
    memory->create(ijskip, ntp1, ntp1, "neigh_request:ijskip");
    for (i = 1; i < ntp1; i++)
      for (j = 1; j < ntp1; j++) ijskip[i][j] = other->ijskip[i][j];
  }
}

/* ---------------------------------------------------------------------- */

void NeighRequest::apply_flags(int flags)
{
  if (flags & REQ_FULL) {
    half = 0;
    full = 1;
  }
  // clang-format off
  if (flags & REQ_GHOST)       { ghost = 1; }
  if (flags & REQ_SIZE)        { size = 1; }
  if (flags & REQ_HISTORY)     { history = 1; }
  if (flags & REQ_NEWTON_ON)   { newton = 1; }
  if (flags & REQ_NEWTON_OFF)  { newton = 2; }
  if (flags & REQ_OCCASIONAL)  { occasional = 1; }
  if (flags & REQ_RESPA_INOUT) { respainner = respaouter = 1; }
  if (flags & REQ_RESPA_ALL)   { respainner = respamiddle = respaouter = 1; }
  if (flags & REQ_SSA)         { ssa = 1; }
  // clang-format on
}

/* ---------------------------------------------------------------------- */

void NeighRequest::set_cutoff(double _cutoff)
{
  cut = 1;
  cutoff = _cutoff;
}

void NeighRequest::set_id(int _id)
{
  id = _id;
}

void NeighRequest::set_kokkos_device(int flag)
{
  kokkos_device = flag;
}

void NeighRequest::set_kokkos_host(int flag)
{
  kokkos_host = flag;
}

void NeighRequest::set_skip(int *_iskip, int **_ijskip)
{
  skip = 1;
  iskip = _iskip;
  ijskip = _ijskip;
}

void NeighRequest::enable_full()
{
  half = 0;
  full = 1;
}

void NeighRequest::enable_ghost()
{
  ghost = 1;
}

void NeighRequest::enable_intel()
{
  intel = 1;
  omp = 0;
}
