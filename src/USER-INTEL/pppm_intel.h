/* *- c++ -*- -----------------------------------------------------------
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
   Contributing authors: William McDoniel (RWTH Aachen University)
                         Rodrigo Canales (RWTH Aachen University)
                         Markus Hoehnerbach (RWTH Aachen University)
                         W. Michael Brown (Intel)
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
  PPPMIntel(class LAMMPS *);
  virtual ~PPPMIntel();
  virtual void init();
  virtual void compute(int, int);
  virtual double memory_usage();
  void compute_first(int, int);
  void compute_second(int, int);
  void pack_buffers();

  #ifdef _LMP_INTEL_OFFLOAD
  int use_base();
  #endif

 protected:
  FixIntel *fix;

  int _use_lrt;
  FFT_SCALAR **perthread_density;
  FFT_SCALAR *particle_ekx;
  FFT_SCALAR *particle_eky;
  FFT_SCALAR *particle_ekz;

  int _use_table;
  int rho_points;
  FFT_SCALAR **rho_lookup;
  FFT_SCALAR **drho_lookup;
  FFT_SCALAR half_rho_scale, half_rho_scale_plus;

  #ifdef _LMP_INTEL_OFFLOAD
  int _use_base;
  #endif

  virtual void allocate();

  template<class flt_t, class acc_t>
  void test_function(IntelBuffers<flt_t,acc_t> *buffers);

  void precompute_rho();
  template<class flt_t, class acc_t>
  void particle_map(IntelBuffers<flt_t,acc_t> *buffers);
  template<class flt_t, class acc_t, int use_table>
  void make_rho(IntelBuffers<flt_t,acc_t> *buffers);
  template<class flt_t, class acc_t>
  void make_rho(IntelBuffers<flt_t,acc_t> *buffers) {
    if (_use_table == 1) {
      make_rho<flt_t,acc_t,1>(buffers);
    } else {
      make_rho<flt_t,acc_t,0>(buffers);
    }
  }
  template<class flt_t, class acc_t, int use_table>
  void fieldforce_ik(IntelBuffers<flt_t,acc_t> *buffers);
  template<class flt_t, class acc_t>
  void fieldforce_ik(IntelBuffers<flt_t,acc_t> *buffers) {
    if (_use_table == 1) {
      fieldforce_ik<flt_t, acc_t, 1>(buffers);
    } else {
      fieldforce_ik<flt_t, acc_t, 0>(buffers);
    }
  }
  template<class flt_t, class acc_t, int use_table>
  void fieldforce_ad(IntelBuffers<flt_t,acc_t> *buffers);
  template<class flt_t, class acc_t>
  void fieldforce_ad(IntelBuffers<flt_t,acc_t> *buffers) {
    if (_use_table == 1) {
      fieldforce_ad<flt_t,acc_t,1>(buffers);
    } else {
      fieldforce_ad<flt_t,acc_t,0>(buffers);
    }
  }
  FFT_SCALAR ***create3d_offset(FFT_SCALAR ***&, int, int, int,
                                int, int, int, const char *name);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: PPPM order greater than supported by USER-INTEL

There is a compile time limit on the maximum order for PPPM
in the USER-INTEL package that might be different from LAMMPS

*/
