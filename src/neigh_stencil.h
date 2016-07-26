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

#ifndef LMP_NEIGH_STENCIL_H
#define LMP_NEIGH_STENCIL_H

#include "pointers.h"

namespace LAMMPS_NS {

class NeighStencil : protected Pointers {
 public:
  int style;                       // ID for NeighStencil method this is
  class NeighBin *nb;              // NeighBin instance I depend on

  bigint last_create;              // timesteps for last operations performed
  bigint last_stencil_memory;
  bigint last_copy_bin;

  int nstencil;                    // # of bins in stencil
  int *stencil;                    // list of bin offsets
  int **stencilxyz;                // bin offsets in xyz dims
  int *nstencil_multi;             // # bins in each type-based multi stencil
  int **stencil_multi;             // list of bin offsets in each stencil
  double **distsq_multi;           // sq distances to bins in each stencil

  NeighStencil(class LAMMPS *);
  virtual ~NeighStencil();
  void copy_neighbor_info();
  void create_setup();
  bigint memory_usage();

  virtual void create() = 0;

 protected:

  // data from Neighbor class

  int neighstyle;
  double cutneighmax;
  double cutneighmaxsq;
  double *cuttypesq;

  // data from NeighBin class

  int mbinx,mbiny,mbinz;
  double binsizex,binsizey,binsizez;
  double bininvx,bininvy,bininvz;

  // data common to all NeighStencil variants

  int xyzflag;                     // 1 if stencilxyz is allocated
  int maxstencil;                  // max size of stencil
  int maxstencil_multi;            // max sizes of stencils
  int sx,sy,sz;                    // extent of stencil in each dim

  int dimension;

  // methods for all NeighStencil variants

  void copy_bin_info();                     // copy info from NeighBin class
  double bin_distance(int, int, int);       // distance between bin corners
};

}

#endif

/* ERROR/WARNING messages:

*/
