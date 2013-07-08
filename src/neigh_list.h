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

#ifndef LMP_NEIGH_LIST_H
#define LMP_NEIGH_LIST_H

#include "pointers.h"
#include "my_page.h"

namespace LAMMPS_NS {

class NeighList : protected Pointers {
 public:
  int index;                       // index of which neigh list it is
                                   // needed when a class invokes it directly
                                   // also indexes the request it came from

  int buildflag;                   // 1 if pair_build invoked every reneigh
  int growflag;                    // 1 if stores atom-based arrays & pages
  int stencilflag;                 // 1 if stores stencil arrays
  int ghostflag;                   // 1 if it stores neighbors of ghosts

  // data structs to store neighbor pairs I,J and associated values

  int inum;                        // # of I atoms neighbors are stored for
  int gnum;                        // # of ghost atoms neighbors are stored for
  int *ilist;                      // local indices of I atoms
  int *numneigh;                   // # of J neighbors for each I atom
  int **firstneigh;                // ptr to 1st J int value of each I atom
  double **firstdouble;            // ptr to 1st J double value of each I atom

  int pgsize;                      // size of each page
  int oneatom;                     // max size for one atom
  int dnum;                        // # of doubles per neighbor, 0 if none
  MyPage<int> *ipage;              // pages of neighbor indices
  MyPage<double> *dpage;           // pages of neighbor doubles, if dnum > 0

  // atom types to skip when building list
  // iskip,ijskip are just ptrs to corresponding request

  int *iskip;         // iskip[i] = 1 if atoms of type I are not in list
  int **ijskip;       // ijskip[i][j] = 1 if pairs of type I,J are not in list

  // settings and pointers for related neighbor lists and fixes

  NeighList *listgranhistory;          // point at history list
  class FixShearHistory *fix_history;  // fix that stores history info

  int respamiddle;              // 1 if this respaouter has middle list
  NeighList *listinner;         // me = respaouter, point to respainner
  NeighList *listmiddle;        // me = respaouter, point to respamiddle
  NeighList *listfull;          // me = half list, point to full I derive from
  NeighList *listcopy;          // me = copy list, point to list I copy from
  NeighList *listskip;          // me = skip list, point to list I skip from

  // stencils of bin indices for neighbor finding

  int maxstencil;                  // max size of stencil
  int nstencil;                    // # of bins in stencil
  int *stencil;                    // list of bin offsets
  int **stencilxyz;                // bin offsets in xyz dims

  int maxstencil_multi;            // max sizes of stencils
  int *nstencil_multi;             // # bins in each type-based multi stencil
  int **stencil_multi;             // list of bin offsets in each stencil
  double **distsq_multi;           // sq distances to bins in each stencil

  class CudaNeighList *cuda_list;  // CUDA neighbor list

  NeighList(class LAMMPS *);
  ~NeighList();
  void setup_pages(int, int, int);      // setup page data structures
  void grow(int);                       // grow maxlocal
  void stencil_allocate(int, int);      // allocate stencil arrays
  void copy_skip_info(int *, int **);   // copy skip info from a neigh request
  void print_attributes();              // debug routine
  int get_maxlocal() {return maxatoms;}
  bigint memory_usage();

 private:
  int maxatoms;                    // size of allocated atom arrays
};

}

#endif
