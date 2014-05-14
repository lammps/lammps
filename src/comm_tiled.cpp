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

#include "lmptype.h"
#include "comm_tiled.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{SINGLE,MULTI};               // same as in Comm

/* ---------------------------------------------------------------------- */

CommTiled::CommTiled(LAMMPS *lmp) : Comm(lmp)
{
  style = 1;
  layout = 0;

  error->all(FLERR,"Comm_style tiled is not yet supported");
}

/* ---------------------------------------------------------------------- */

CommTiled::~CommTiled()
{
}

/* ---------------------------------------------------------------------- */

void CommTiled::init()
{
}

/* ----------------------------------------------------------------------
   setup spatial-decomposition communication patterns
   function of neighbor cutoff(s) & cutghostuser & current box size
   single mode sets slab boundaries (slablo,slabhi) based on max cutoff
   multi mode sets type-dependent slab boundaries (multilo,multihi)
------------------------------------------------------------------------- */

void CommTiled::setup()
{
}

/* ----------------------------------------------------------------------
   forward communication of atom coords every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommTiled::forward_comm(int dummy)
{
}

/* ----------------------------------------------------------------------
   reverse communication of forces on atoms every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommTiled::reverse_comm()
{
}

/* ----------------------------------------------------------------------
   exchange: move atoms to correct processors
   atoms exchanged with all 6 stencil neighbors
   send out atoms that have left my box, receive ones entering my box
   atoms will be lost if not inside some proc's box
     can happen if atom moves outside of non-periodic bounary
     or if atom moves more than one proc away
   this routine called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before exchange is called
------------------------------------------------------------------------- */

void CommTiled::exchange()
{
}

/* ----------------------------------------------------------------------
   borders: list nearby atoms to send to neighboring procs at every timestep
   one list is created for every swap that will be made
   as list is made, actually do swaps
   this does equivalent of a communicate, so don't need to explicitly
     call communicate routine on reneighboring timestep
   this routine is called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before borders is called
------------------------------------------------------------------------- */

void CommTiled::borders()
{
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Pair
------------------------------------------------------------------------- */

void CommTiled::forward_comm_pair(Pair *pair)
{
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Pair
------------------------------------------------------------------------- */

void CommTiled::reverse_comm_pair(Pair *pair)
{
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Fix
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void CommTiled::forward_comm_fix(Fix *fix)
{
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void CommTiled::reverse_comm_fix(Fix *fix)
{
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Fix
   n = total datums for all atoms, allows for variable number/atom
------------------------------------------------------------------------- */

void CommTiled::forward_comm_variable_fix(Fix *fix)
{
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix
   n = total datums for all atoms, allows for variable number/atom
------------------------------------------------------------------------- */

void CommTiled::reverse_comm_variable_fix(Fix *fix)
{
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Compute
------------------------------------------------------------------------- */

void CommTiled::forward_comm_compute(Compute *compute)
{
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Compute
------------------------------------------------------------------------- */

void CommTiled::reverse_comm_compute(Compute *compute)
{
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Dump
------------------------------------------------------------------------- */

void CommTiled::forward_comm_dump(Dump *dump)
{
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Dump
------------------------------------------------------------------------- */

void CommTiled::reverse_comm_dump(Dump *dump)
{
}

/* ----------------------------------------------------------------------
   forward communication of N values in array
------------------------------------------------------------------------- */

void CommTiled::forward_comm_array(int n, double **array)
{
}

/* ----------------------------------------------------------------------
   exchange info provided with all 6 stencil neighbors
------------------------------------------------------------------------- */

int CommTiled::exchange_variable(int n, double *inbuf, double *&outbuf)
{
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint CommTiled::memory_usage()
{
  bigint bytes = 0;
  return bytes;
}
