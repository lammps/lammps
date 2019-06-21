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

#include <cstring>
#include "comm_tiled_kokkos.h"
#include "comm_brick.h"
#include "atom_kokkos.h"
#include "atom_vec.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "output.h"
#include "dump.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

#define BUFFACTOR 1.5
#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000
#define EPSILON 1.0e-6

#define DELTA_PROCS 16

/* ---------------------------------------------------------------------- */

CommTiledKokkos::CommTiledKokkos(LAMMPS *lmp) : CommTiled(lmp)
{
  comm_style = (const char *) "tiled/kk";
}

/* ---------------------------------------------------------------------- */
//IMPORTANT: we *MUST* pass "*oldcomm" to the Comm initializer here, as
//           the code below *requires* that the (implicit) copy constructor
//           for Comm is run and thus creating a shallow copy of "oldcomm".
//           The call to Comm::copy_arrays() then converts the shallow copy
//           into a deep copy of the class with the new layout.

CommTiledKokkos::CommTiledKokkos(LAMMPS *lmp, Comm *oldcomm) : CommTiled(lmp,oldcomm)
{

}

/* ---------------------------------------------------------------------- */

CommTiledKokkos::~CommTiledKokkos()
{

}

/* ----------------------------------------------------------------------
   forward communication of atom coords every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm(int dummy)
{
  if (comm_x_only) {
    atomKK->sync(Host,X_MASK);
    atomKK->modified(Host,X_MASK);
  } else if (ghost_velocity) {
    atomKK->sync(Host,X_MASK | V_MASK);
    atomKK->modified(Host,X_MASK | V_MASK);
  } else {
    atomKK->sync(Host,ALL_MASK);
    atomKK->modified(Host,ALL_MASK);
  }

  CommTiled::forward_comm(dummy);
}

/* ----------------------------------------------------------------------
   reverse communication of forces on atoms every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm()
{
  if (comm_f_only)
    atomKK->sync(Host,F_MASK);
  else
    atomKK->sync(Host,ALL_MASK);
  CommTiled::reverse_comm();
  if (comm_f_only)
    atomKK->modified(Host,F_MASK);
  else
    atomKK->modified(Host,ALL_MASK);
  atomKK->sync(Device,ALL_MASK);
}

/* ----------------------------------------------------------------------
   exchange: move atoms to correct processors
   atoms exchanged with procs that touch sub-box in each of 3 dims
   send out atoms that have left my box, receive ones entering my box
   atoms will be lost if not inside a touching proc's box
     can happen if atom moves outside of non-periodic bounary
     or if atom moves more than one proc away
   this routine called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before exchange is called
------------------------------------------------------------------------- */

void CommTiledKokkos::exchange()
{
  atomKK->sync(Host,ALL_MASK);
  CommTiled::exchange();
  atomKK->modified(Host,ALL_MASK);
}

/* ----------------------------------------------------------------------
   borders: list nearby atoms to send to neighboring procs at every timestep
   one list is created per swap/proc that will be made
   as list is made, actually do communication
   this does equivalent of a forward_comm(), so don't need to explicitly
     call forward_comm() on reneighboring timestep
   this routine is called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before borders is called
------------------------------------------------------------------------- */

void CommTiledKokkos::borders()
{
  atomKK->sync(Host,ALL_MASK);
  CommTiled::borders();
  atomKK->modified(Host,ALL_MASK);
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Pair
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm_pair(Pair *pair)
{
  CommTiled::forward_comm_pair(pair);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Pair
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm_pair(Pair *pair)
{
  CommTiled::reverse_comm_pair(pair);
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Fix
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_forward from Fix
   size > 0 -> Fix passes max size per atom
   the latter is only useful if Fix does several comm modes,
     some are smaller than max stored in its comm_forward
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm_fix(Fix *fix, int size)
{
  CommTiled::forward_comm_fix(fix,size);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_forward from Fix
   size > 0 -> Fix passes max size per atom
   the latter is only useful if Fix does several comm modes,
     some are smaller than max stored in its comm_forward
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm_fix(Fix *fix, int size)
{
  CommTiled::reverse_comm_fix(fix,size);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix with variable size data
   query fix for all pack sizes to insure buf_send is big enough
   handshake sizes before irregular comm to insure buf_recv is big enough
   NOTE: how to setup one big buf recv with correct offsets ??
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm_fix_variable(Fix *fix)
{
  CommTiled::reverse_comm_fix_variable(fix);
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Compute
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm_compute(Compute *compute)
{
  CommTiled::forward_comm_compute(compute);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Compute
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm_compute(Compute *compute)
{
  CommTiled::reverse_comm_compute(compute);
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Dump
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm_dump(Dump *dump)
{
  CommTiled::forward_comm_dump(dump);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Dump
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm_dump(Dump *dump)
{
  CommTiled::reverse_comm_dump(dump);
}

/* ----------------------------------------------------------------------
   forward communication of Nsize values in per-atom array
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm_array(int nsize, double **array)
{
  CommTiled::forward_comm_array(nsize,array);
}

/* ----------------------------------------------------------------------
   exchange info provided with all 6 stencil neighbors
   NOTE: this method is currently not used
------------------------------------------------------------------------- */

int CommTiledKokkos::exchange_variable(int n, double *inbuf, double *&outbuf)
{
  return CommTiled::exchange_variable(n,inbuf,outbuf);
}
