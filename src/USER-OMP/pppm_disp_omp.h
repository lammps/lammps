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

KSpaceStyle(pppm/disp/omp,PPPMDispOMP)

#else

#ifndef LMP_PPPM_DISP_OMP_H
#define LMP_PPPM_DISP_OMP_H

#include "pppm_disp.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

  class PPPMDispOMP : public PPPMDisp, public ThrOMP {
 public:
  PPPMDispOMP(class LAMMPS *, int, char **);
  virtual ~PPPMDispOMP () {};
  virtual void compute(int, int);

 protected:
  virtual void allocate();
  virtual void deallocate();

  virtual void compute_gf();
  virtual void compute_gf_6();

  virtual void particle_map(double,double,double,
                            double,int**,int,int,
                            int,int,int,int,int,int);
                                

  virtual void fieldforce_c_ik();
  virtual void fieldforce_c_ad();
  virtual void fieldforce_c_peratom();
  virtual void fieldforce_g_ik();
  virtual void fieldforce_g_ad();
  virtual void fieldforce_g_peratom();
  virtual void fieldforce_a_ik();
  virtual void fieldforce_a_ad();
  virtual void fieldforce_a_peratom();

  virtual void make_rho_c();
  virtual void make_rho_g();
  virtual void make_rho_a();

  void compute_rho1d_thr(FFT_SCALAR * const * const, const FFT_SCALAR &,
                         const FFT_SCALAR &, const FFT_SCALAR &,
                         const int, FFT_SCALAR * const * const);
  void compute_drho1d_thr(FFT_SCALAR * const * const, const FFT_SCALAR &,
			  const FFT_SCALAR &, const FFT_SCALAR &,
                          const int, FFT_SCALAR * const * const);
//  void compute_rho_coeff();
//  void slabcorr(int);

};

}

#endif
#endif
