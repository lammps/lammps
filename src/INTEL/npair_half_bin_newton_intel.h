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
NPairStyle(half/bin/newton/intel,
           NPairHalfBinNewtonIntel,
           NP_HALF | NP_BIN | NP_NEWTON | NP_ORTHO | NP_INTEL);
// clang-format on
#else

#ifndef LMP_NPAIR_HALF_BIN_NEWTON_INTEL_H
#define LMP_NPAIR_HALF_BIN_NEWTON_INTEL_H

#include "fix_intel.h"
#include "npair_intel.h"

namespace LAMMPS_NS {

class NPairHalfBinNewtonIntel : public NPairIntel {
 public:
  NPairHalfBinNewtonIntel(class LAMMPS *);
  void build(class NeighList *) override;

 private:
  template <class flt_t, class acc_t> void hbni(NeighList *, IntelBuffers<flt_t, acc_t> *);
};

}    // namespace LAMMPS_NS

#endif
#endif
