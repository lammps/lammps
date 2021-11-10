// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef NPAIR_CLASS
// clang-format off
NPairStyle(skip/intel,
           NPairSkipIntel,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_INTEL);

NPairStyle(skip/ghost/intel,
           NPairSkipIntel,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_GHOST | NP_INTEL);
// clang-format on
#else

#ifndef LMP_NPAIR_SKIP_INTEL_H
#define LMP_NPAIR_SKIP_INTEL_H

#include "fix_intel.h"
#include "npair.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace LAMMPS_NS {

class NPairSkipIntel : public NPair {
 public:
  NPairSkipIntel(class LAMMPS *);
  ~NPairSkipIntel();
  virtual void copy_neighbor_info();
  void build(class NeighList *);

 protected:
  FixIntel *_fix;
  int *_inum_starts, *_inum_counts, *_full_props;

  template <class flt_t, int THREE>
  void build_t(NeighList *, int *numhalf, int *cnumneigh, int *numhalf_skip);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: The 'package intel' command is required for /intel styles

Self explanatory.

*/
