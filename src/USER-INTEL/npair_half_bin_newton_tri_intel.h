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

#ifdef NPAIR_CLASS

NPairStyle(half/bin/newton/tri/intel,
           NPairHalfBinNewtonTriIntel,
           NP_HALF | NP_BIN | NP_NEWTON | NP_TRI | NP_INTEL)

#else

#ifndef LMP_NPAIR_HALF_BIN_NEWTON_INTEL_TRI_H
#define LMP_NPAIR_HALF_BIN_NEWTON_INTEL_TRI_H

#include "npair_intel.h"
#include "fix_intel.h"

namespace LAMMPS_NS {

class NPairHalfBinNewtonTriIntel : public NPairIntel {
 public:
  NPairHalfBinNewtonTriIntel(class LAMMPS *);
  ~NPairHalfBinNewtonTriIntel() {}
  void build(class NeighList *);

 private:
  template <class flt_t, class acc_t>
  void hbnti(NeighList *, IntelBuffers<flt_t,acc_t> *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
