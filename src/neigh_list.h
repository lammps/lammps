/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_NEIGH_LIST_H
#define LMP_NEIGH_LIST_H

#include "pointers.h"    // IWYU pragma: export

namespace LAMMPS_NS {

class NeighList : protected Pointers {
 public:
  enum RequestorType { NONE, PAIR, FIX, COMPUTE };
  void *requestor;                 // object that made request
  RequestorType requestor_type;    // type of requestor

  int index;    // index of which neigh list this is
                // also indexes the request it came from
                // and the npair list of NPair classes

  int bin_method;        // 0 if no binning, else 1-N index into binnames
  int stencil_method;    // 0 if no stencil, else 1-N index into stencilnames
  int pair_method;       // 0 if no pair, else 1-N index into pairnames

  // settings from NeighRequest

  int occasional;     // 0 if build every reneighbor, 1 if not
  int ghost;          // 1 if list stores neighbors of ghosts
  int ssa;            // 1 if list stores Shardlow data
  int history;        // 1 if there is neigh history (FixNeighHist)
  int respaouter;     // 1 if list is a rRespa outer list
  int respamiddle;    // 1 if there is also a rRespa middle list
  int respainner;     // 1 if there is also a rRespa inner list
  int copy;           // 1 if this list is copied from another list
  int trim;           // 1 if this list is trimmed from another list
  int kk2cpu;         // 1 if this list is copied from Kokkos to CPU
  int copymode;       // 1 if this is a Kokkos on-device copy
  int id;             // copied from neighbor list request
  int molskip;        // 1/2 if this is an intra-/inter-molecular skip list

  // data structs to store neighbor pairs I,J and associated values

  int inum;            // # of I atoms neighbors are stored for
  int gnum;            // # of ghost atoms neighbors are stored for
  int *ilist;          // local indices of I atoms
  int *numneigh;       // # of J neighbors for each I atom
  int **firstneigh;    // ptr to 1st J int value of each I atom
  int maxatom;         // size of allocated per-atom arrays

  int pgsize;            // size of each page
  int oneatom;           // max size for one atom
  MyPage<int> *ipage;    // pages of neighbor indices

  // data structs to store rRESPA neighbor pairs I,J and associated values

  int inum_inner;            // # of I atoms neighbors are stored for
  int gnum_inner;            // # of ghost atoms neighbors are stored for
  int *ilist_inner;          // local indices of I atoms
  int *numneigh_inner;       // # of J neighbors for each I atom
  int **firstneigh_inner;    // ptr to 1st J int value of each I atom

  int inum_middle;            // # of I atoms neighbors are stored for
  int gnum_middle;            // # of ghost atoms neighbors are stored for
  int *ilist_middle;          // local indices of I atoms
  int *numneigh_middle;       // # of J neighbors for each I atom
  int **firstneigh_middle;    // ptr to 1st J int value of each I atom

  MyPage<int> *ipage_inner;     // pages of neighbor indices for inner
  MyPage<int> *ipage_middle;    // pages of neighbor indices for middle

  // atom types to skip when building list
  // copied info from corresponding request into realloced vec/array

  int *iskip;      // iskip[i] = 1 if atoms of type I are not in list
  int **ijskip;    // ijskip[i][j] = 1 if pairs of type I,J are not in list

  // settings and pointers for related neighbor lists and fixes

  NeighList *listcopy;    // me = copy list, point to list I copy from
  NeighList *listskip;    // me = skip list, point to list I skip from
  NeighList *listfull;    // me = half list, point to full I derive from

  class Fix *fix_bond;    // fix that stores bond info

  // Kokkos package

  int kokkos;    // 1 if list stores Kokkos data
  ExecutionSpace execution_space;

  // DPD-REACT package and Shardlow Splitting Algorithm (SSA) support

  class NPair *np;    // ptr to NPair instance I depend on

  // methods

  NeighList(class LAMMPS *);
  ~NeighList() override;
  void post_constructor(class NeighRequest *);
  void setup_pages(int, int);    // setup page data structures
  void grow(int, int);           // grow all data structs
  void print_attributes();       // debug routine
  int get_maxlocal() { return maxatom; }
  double memory_usage();
};

}    // namespace LAMMPS_NS

#endif
