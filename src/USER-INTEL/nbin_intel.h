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

#ifdef NBIN_CLASS

NBinStyle(intel,
          NBinIntel,
          NB_INTEL)

#else

#ifndef LMP_NBIN_INTEL_H
#define LMP_NBIN_INTEL_H

#include "nbin_standard.h"
#include "fix_intel.h"
#include "memory.h"

namespace LAMMPS_NS {

class NBinIntel : public NBinStandard {
 public:
  NBinIntel(class LAMMPS *);
  ~NBinIntel();
  void bin_atoms_setup(int);
  void bin_atoms();
  int * get_binpacked() { return _binpacked; }

 private:
  FixIntel *_fix;
  int *_atombin, *_binpacked;
  int _precision_mode;
  bigint memory_usage();

  template <class flt_t, class acc_t>
  void bin_atoms(IntelBuffers<flt_t,acc_t> *);

  #ifdef _LMP_INTEL_OFFLOAD
  int _cop, _offload_alloc;
  #endif
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: The 'package intel' command is required for /intel styles

Self-explanatory.

E: Intel package expects no atoms within cutoff of {1e15,1e15,1e15}.

The Intel package can make use of dummy atoms for padding with a large position
that should not be within the cutoff.

*/
