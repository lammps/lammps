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

#ifndef LMP_NSTENCIL_H
#define LMP_NSTENCIL_H

#include "pointers.h"

namespace LAMMPS_NS {

class NStencil : protected Pointers {
 public:
  int istyle;                      // 1-N index into binnames
  class NBin *nb;                  // ptr to NBin instance I depend on
  bigint last_stencil;             // last timestep stencil was created

  int nstencil;                    // # of bins in stencil
  int *stencil;                    // list of bin offsets
  int **stencilxyz;                // bin offsets in xyz dims
  int *nstencil_multi;             // # bins in each type-based multi stencil
  int **stencil_multi;             // list of bin offsets in each stencil
  double **distsq_multi;           // sq distances to bins in each stencil
  int sx,sy,sz;                    // extent of stencil in each dim
  int post_create_flag;            // 1 if child stencil class implements post create methods

  double cutoff_custom;            // cutoff set by requestor

  NStencil(class LAMMPS *);
  virtual ~NStencil();
  void post_constructor(class NeighRequest *);
  void copy_neighbor_info();
  virtual void create_setup();
  bigint memory_usage();

  virtual void create() = 0;
  virtual void post_create() {}
  virtual void post_create_setup() {}

  inline int get_maxstencil() {return maxstencil;}

 protected:

  // data from Neighbor class

  int neighstyle;
  double cutneighmax;
  double cutneighmaxsq;
  double *cuttypesq;

  // data from NBin class

  int mbinx,mbiny,mbinz;
  double binsizex,binsizey,binsizez;
  double bininvx,bininvy,bininvz;

  // data common to all NStencil variants

  int xyzflag;                     // 1 if stencilxyz is allocated
  int maxstencil;                  // max size of stencil
  int maxstencil_multi;            // max sizes of stencils

  int dimension;

  // methods for all NStencil variants

  void copy_bin_info();                     // copy info from NBin class
  double bin_distance(int, int, int);       // distance between bin corners
};

}

#endif

/* ERROR/WARNING messages:

*/
