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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(zero/intel,PairZeroIntel);
// clang-format on
#else

#ifndef LMP_PAIR_ZERO_INTEL_H
#define LMP_PAIR_ZERO_INTEL_H

#include "fix_intel.h"
#include "pair_zero.h"

namespace LAMMPS_NS {

class PairZeroIntel : public PairZero {
 public:
  PairZeroIntel(class LAMMPS *);
  void compute(int, int) override;
  void init_style() override;

 private:
  FixIntel *fix;

  template <class flt_t, class acc_t>
  void compute(int eflag, int vflag, IntelBuffers<flt_t, acc_t> *buffers);
};
}    // namespace LAMMPS_NS

#endif
#endif
