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
   Contributing author: Rodrigo Canales (RWTH Aachen University)
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS

KSpaceStyle(pppm/intel,PPPMIntel)

#else

#ifndef LMP_PPPMINTEL_H
#define LMP_PPPMINTEL_H

#include "pppm.h"
#include "fix_intel.h"

namespace LAMMPS_NS {

class PPPMIntel : public PPPM {
 public:
  PPPMIntel(class LAMMPS *, int, char **);
  virtual ~PPPMIntel();
  virtual void init();
  virtual void compute(int, int);

 protected:
  FixIntel *fix;

  #ifdef _LMP_INTEL_OFFLOAD
  int _use_base;
  #endif

  template<class flt_t, class acc_t>
  void particle_map(IntelBuffers<flt_t,acc_t> *buffers);
  template<class flt_t, class acc_t>
  void make_rho(IntelBuffers<flt_t,acc_t> *buffers);
  template<class flt_t, class acc_t>
  void fieldforce_ik(IntelBuffers<flt_t,acc_t> *buffers);
  template<class flt_t, class acc_t>
  void fieldforce_ad(IntelBuffers<flt_t,acc_t> *buffers);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: PPPM order greater than supported by USER-INTEL

There is a compile time limit on the maximum order for PPPM
in the USER-INTEL package that might be different from LAMMPS

*/
