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

#ifndef LMP_NSTENCIL_H
#define LMP_NSTENCIL_H

#include "pointers.h"    // IWYU pragma: keep

namespace LAMMPS_NS {

class NStencil : protected Pointers {
 public:
  int istyle;             // 1-N index into binnames
  class NBin *nb;         // ptr to NBin instance I depend on
  bigint last_stencil;    // last timestep stencil was created

  int nstencil;                 // # of bins in stencil
  int *stencil;                 // list of bin offsets
  int **stencilxyz;             // bin offsets in xyz dims
  int *nstencil_multi_old;      // # bins in each type-based old multi stencil
  int **stencil_multi_old;      // list of bin offsets in each stencil
  double **distsq_multi_old;    // sq distances to bins in each stencil
  int **nstencil_multi;         // # bins bins in each igroup-jgroup multi stencil
  int ***stencil_multi;         // list of bin offsets in each multi stencil
  int **maxstencil_multi;       // max stencil size for each multi stencil
  int maxcollections;           // size of multi arrays

  int sx, sy, sz;            // extent of stencil in each dim
  int **stencil_sx_multi;    // analogs for each multi stencil
  int **stencil_sy_multi;
  int **stencil_sz_multi;

  double cutoff_custom;    // cutoff set by requestor

  // Arrays to store options for multi itype-jtype stencils
  bool **flag_half_multi;    // flag creation of a half stencil for icollection-jcollection
  bool **flag_skip_multi;    // skip creation of icollection-jcollection stencils (for newton on)
  bool **flag_same_multi;    // flag same size collection (doesn't always correspond to a half, e.g. newton + tri)
  int **bin_collection_multi;    // what collection to use for bin information

  NStencil(class LAMMPS *);
  ~NStencil() override;
  void post_constructor(class NeighRequest *);
  void copy_neighbor_info();
  virtual void create_setup();
  double memory_usage();

  virtual void create() = 0;

  inline int get_maxstencil() { return maxstencil; }

 protected:
  // data from Neighbor class

  int neighstyle;
  double cutneighmax;
  double cutneighmaxsq;
  double *cuttypesq;
  double **cutneighsq;
  double **cutcollectionsq;
  int ncollections;
  int *collection;

  // data from NBin class

  int mbinx, mbiny, mbinz;
  double binsizex, binsizey, binsizez;
  double bininvx, bininvy, bininvz;

  // data from NBin class for multi

  int *mbinx_multi;
  int *mbiny_multi;
  int *mbinz_multi;
  double *binsizex_multi;
  double *binsizey_multi;
  double *binsizez_multi;
  double *bininvx_multi;
  double *bininvy_multi;
  double *bininvz_multi;

  // Stored bin information for each stencil

  int **stencil_mbinx_multi;
  int **stencil_mbiny_multi;
  int **stencil_mbinz_multi;
  double **stencil_binsizex_multi;
  double **stencil_binsizey_multi;
  double **stencil_binsizez_multi;

  // data common to all NStencil variants

  int xyzflag;                 // 1 if stencilxyz is allocated
  int maxstencil;              // max size of stencil
  int maxstencil_multi_old;    // max sizes of stencils

  int dimension;

  // methods for standard NStencil variants

  void copy_bin_info();                  // copy info from NBin class
  double bin_distance(int, int, int);    // distance between bin corners

  // methods for multi NStencil

  double bin_distance_multi(int, int, int,
                            int);    // distance between bin corners for different collections
  void copy_bin_info_multi();        // copy multi info from NBin class
  virtual void set_stencil_properties() {}    // determine which stencils to build and how
};

}    // namespace LAMMPS_NS

#endif
