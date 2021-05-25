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

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(pppm/cg/omp,PPPMCGOMP);
// clang-format on
#else

#ifndef LMP_PPPM_CG_OMP_H
#define LMP_PPPM_CG_OMP_H

#include "pppm_cg.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PPPMCGOMP : public PPPMCG, public ThrOMP {
 public:
  PPPMCGOMP(class LAMMPS *);
  virtual ~PPPMCGOMP();
  virtual void compute(int, int);

 protected:
  virtual void allocate();

  virtual void compute_gf_ik();
  virtual void compute_gf_ad();

  virtual void make_rho();
  virtual void fieldforce_ik();
  virtual void fieldforce_ad();
  virtual void fieldforce_peratom();

 private:
  void compute_rho1d_thr(FFT_SCALAR *const *const, const FFT_SCALAR &, const FFT_SCALAR &,
                         const FFT_SCALAR &);
  void compute_drho1d_thr(FFT_SCALAR *const *const, const FFT_SCALAR &, const FFT_SCALAR &,
                          const FFT_SCALAR &);
};

}    // namespace LAMMPS_NS

#endif
#endif
