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

#ifdef KSPACE_CLASS

KSpaceStyle(pppm/tip4p/cg/omp,PPPMTIP4PCGOMP)

#else

#ifndef LMP_PPPM_TIP4P_CG_OMP_H
#define LMP_PPPM_TIP4P_CG_OMP_H

#include "pppm_tip4p_cg.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PPPMTIP4PCGOMP : public PPPMTIP4PCG, public ThrOMP {
 public:
  PPPMTIP4PCGOMP(class LAMMPS *, int, char **);
  virtual ~PPPMTIP4PCGOMP () {};
  virtual void compute(int, int);

 protected:
  virtual void allocate();
  virtual void deallocate();

  virtual void compute_gf_ik();
  virtual void compute_gf_ad();

  virtual void particle_map();
  virtual void make_rho();

  virtual void fieldforce_ik();
  virtual void fieldforce_ad();
  // virtual void fieldforce_peratom();  XXX: need to benchmark first.

 private:
  void compute_rho1d_thr(FFT_SCALAR * const * const, const FFT_SCALAR &,
                         const FFT_SCALAR &, const FFT_SCALAR &);
  void compute_drho1d_thr(FFT_SCALAR * const * const, const FFT_SCALAR &,
			  const FFT_SCALAR &, const FFT_SCALAR &);

  void find_M_thr(const int, int &, int &, dbl3_t &);

//  void slabcorr(int);  // XXX: not (yet) multi-threaded

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Kspace style pppm/tip4p/omp requires newton on

Self-explanatory.

*/
