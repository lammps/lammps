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

/* ----------------------------------------------------------------------
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#ifdef NPAIR_CLASS
// clang-format off
NPairStyle(trim/intel,
           NPairTrimIntel,
           NP_COPY | NP_TRIM | NP_INTEL);
// clang-format on
#else

#ifndef LMP_NPAIR_TRIM_INTEL_H
#define LMP_NPAIR_TRIM_INTEL_H

#include "fix_intel.h"
#include "npair.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace LAMMPS_NS {

class NPairTrimIntel : public NPair {
 public:
  NPairTrimIntel(class LAMMPS *);
  void build(class NeighList *) override;

 protected:
  FixIntel *_fix;

  template <class flt_t, class acc_t> void build_t(NeighList *, IntelBuffers<flt_t, acc_t> *);
};

}    // namespace LAMMPS_NS

#endif
#endif
