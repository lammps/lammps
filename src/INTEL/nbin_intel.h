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

#ifdef NBIN_CLASS
// clang-format off
NBinStyle(intel,
          NBinIntel,
          NB_STANDARD | NB_INTEL);
// clang-format on
#else

#ifndef LMP_NBIN_INTEL_H
#define LMP_NBIN_INTEL_H

#include "fix_intel.h"
#include "memory.h"
#include "nbin_standard.h"

namespace LAMMPS_NS {

class NBinIntel : public NBinStandard {
 public:
  NBinIntel(class LAMMPS *);
  ~NBinIntel() override;

  void bin_atoms_setup(int) override;
  void bin_atoms() override;
  int *get_binpacked() { return _binpacked; }

 private:
  FixIntel *_fix;
  int *_atombin, *_binpacked;
  int _precision_mode;
  double memory_usage() override;

  template <class flt_t, class acc_t> void bin_atoms(IntelBuffers<flt_t, acc_t> *);

#ifdef _LMP_INTEL_OFFLOAD
  int _cop, _offload_alloc;
#endif
};

}    // namespace LAMMPS_NS

#endif
#endif
