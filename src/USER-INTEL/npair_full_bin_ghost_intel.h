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

/* ----------------------------------------------------------------------
   Contributing authors: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#ifdef NPAIR_CLASS

NPairStyle(full/bin/ghost/intel,
           NPairFullBinGhostIntel,
           NP_FULL | NP_BIN | NP_GHOST | NP_NEWTON | NP_NEWTOFF |
           NP_ORTHO | NP_TRI | NP_INTEL)

#else

#ifndef LMP_NPAIR_FULL_BIN_GHOST_INTEL_H
#define LMP_NPAIR_FULL_BIN_GHOST_INTEL_H

#include "npair_intel.h"

namespace LAMMPS_NS {

class NPairFullBinGhostIntel : public NPairIntel {
 public:
  NPairFullBinGhostIntel(class LAMMPS *);
  ~NPairFullBinGhostIntel() {}
  void build(class NeighList *);
 private:
  template<class flt_t, class acc_t>
  void fbi(NeighList * list, IntelBuffers<flt_t,acc_t> * buffers);
  template<class flt_t, class acc_t, int need_ic>
  void fbi(const int offload, NeighList * list,
           IntelBuffers<flt_t,acc_t> * buffers,
           const int astart, const int aend);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
