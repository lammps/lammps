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

#ifdef INTEGRATE_CLASS

IntegrateStyle(verlet/lrt/intel,VerletLRTIntel)

#else

#ifndef LMP_VERLET_LRT_INTEL_H
#define LMP_VERLET_LRT_INTEL_H

#include "verlet.h"
#include "pppm_intel.h"

#ifdef LMP_INTEL_USELRT
  #if defined(LMP_INTEL_LRT11) || defined(__APPLE__)
    #if __cplusplus > 199711L
      #define _LMP_INTEL_LRT_11
      #include <thread>
    #else
      #undef LMP_INTEL_USELRT
      #ifdef LMP_INTEL_LRT11
        #error C++11 support required for LMP_INTEL_LRT11 define
      #endif
    #endif
  #else
    #define _LMP_INTEL_LRT_PTHREAD
    #include <pthread.h>
  #endif
#endif

namespace LAMMPS_NS {

class VerletLRTIntel : public Verlet {
 public:
  VerletLRTIntel(class LAMMPS *, int, char **);
  virtual ~VerletLRTIntel();
  virtual void init();
  virtual void setup(int flag);
  virtual void run(int);

 protected:
  PPPMIntel *_intel_kspace;

  #if defined(_LMP_INTEL_LRT_PTHREAD)
  static void *k_launch_loop(void *context);
  pthread_t _kspace_thread;
  pthread_attr_t _kspace_attr;
  pthread_mutex_t _kmutex;
  pthread_cond_t _kcond;
  int _kspace_ready, _kspace_done, _krun_n;
  #endif
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: LRT otion for Intel package disabled at compile time

This option cannot be used with the Intel package because LAMMPS was not built
with support for it.

W: No fixes defined, atoms won't move

If you are not using a fix like nve, nvt, npt then atom velocities and
coordinates will not be updated during timestepping.

*/
