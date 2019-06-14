/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_NEIGH_REQUEST_H
#define LMP_NEIGH_REQUEST_H

#include "pointers.h"

namespace LAMMPS_NS {

class NeighRequest : protected Pointers {
 public:
  int index;                // index of which neigh request this is
  void *requestor;          // class that made request
  int requestor_instance;   // instance of that class (only Fix, Compute, Pair)
  int id;                   // ID of request as stored by requestor
                            // used to track multiple requests from one class

  // -----------------------------
  // flags set by requesting class for attributes of neighor list they need
  // all must be set appropriately, all have defaults
  // -----------------------------

  // which class style requests the list
  // one flag is 1, others are 0

  int pair;              // pair is set by default
  int fix;
  int compute;
  int command;
  int neigh;

  // half/full setting, determines which neighbors appear in list
  // one flag is 1, other is 0

  int half;              // half neigh list (set by default)
  int full;              // full neigh list
  int CAC;               // CAC style neighbor list
  
  // attribute flags, all are 0 by default

  int occasional;        // how often list is built
                         // 0 if needed every reneighboring during run
                         // 1 if only occasionally needed by a fix, compute, etc

  int newton;            // which owned/ghost pairs appear in list
                         // 0 if use force::newton_pair setting
                         // 1 if override with pair newton on
                         // 2 if override with pair newton off

  int ghost;             // 1 if includes ghost atom neighbors
  int size;              // 1 if pair cutoff set by particle radius
  int history;           // 1 if there is also neigh history info (FixNeighHist)
  int granonesided;      // 1 if one-sided granular list for
                         //   sphere/surf interactions
  int respainner;        // 1 if need a rRESPA inner list
  int respamiddle;       // 1 if need a rRESPA middle list
  int respaouter;        // 1 if need a rRESPA outer list
  int bond;              // 1 if store bond neighbors instead of atom neighs
  int omp;               // set by USER-OMP package
  int intel;             // set by USER-INTEL package
  int kokkos_host;       // set by KOKKOS package
  int kokkos_device;
  int ssa;               // set by USER-DPD package, for Shardlow lists
  int cut;               // 1 if use a non-standard cutoff length
  double cutoff;         // special cutoff distance for this list

  // flags set by pair hybrid

  int skip;              // 1 if this list skips atom types from another list
  int *iskip;            // iskip[i] if atoms of type I are not in list
  int **ijskip;          // ijskip[i][j] if pairs of type I,J are not in list

  // command_style only set if command = 1
  // allows print_pair_info() to access command name

  const char *command_style;

  // -----------------------------
  // flags set by Neighbor class to morph original requests
  // -----------------------------

  int skiplist;          // index of list to skip from
  int off2on;            // 1 if this is newton on list, but skips from off list

  int copy;              // 1 if this list copied from another list
  int copylist;          // index of list to copy from

  int halffull;          // 1 if half list computed from another full list
  int halffulllist;      // index of full list to derive half from

  int unique;            // 1 if this list requires its own
                         // NStencil, Nbin class - because of requestor cutoff

  // -----------------------------
  // internal settings made by Neighbor class
  // -----------------------------

  int index_bin;         // index of NBin class assigned to this request
  int index_stencil;     // index of NStencil class assigned to this request
  int index_pair;        // index of NPair class assigned to this request

  // methods

  NeighRequest(class LAMMPS *);
  ~NeighRequest();
  void archive();
  int identical(NeighRequest *);
  int same_skip(NeighRequest *);
  void copy_request(NeighRequest *, int);
};

}

#endif
