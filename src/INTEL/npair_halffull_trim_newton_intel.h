// clang-format off
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

/* ----------------------------------------------------------------------
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#ifdef NPAIR_CLASS
// clang-format off
NPairStyle(halffull/trim/newton/intel,
           NPairHalffullTrimNewtonIntel,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI| NP_TRIM | NP_INTEL);

NPairStyle(halffull/trim/newton/skip/intel,
           NPairHalffullTrimNewtonIntel,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_TRIM | NP_INTEL);
// clang-format on
#else

#ifndef LMP_NPAIR_HALFFULL_TRIM_NEWTON_INTEL_H
#define LMP_NPAIR_HALFFULL_TRIM_NEWTON_INTEL_H

#include "fix_intel.h"
#include "npair.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace LAMMPS_NS {

class NPairHalffullTrimNewtonIntel : public NPair {
 public:
  NPairHalffullTrimNewtonIntel(class LAMMPS *);
  void build(class NeighList *) override;

 protected:
  FixIntel *_fix;

  template <class flt_t, class acc_t> void build_t(NeighList *, IntelBuffers<flt_t, acc_t> *);

  template <class flt_t, class acc_t> void build_t3(NeighList *, int *, IntelBuffers<flt_t, acc_t> *);
};

}    // namespace LAMMPS_NS

#endif
#endif
