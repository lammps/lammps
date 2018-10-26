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

KSpaceStyle(pppm/disp/tip4p/omp,PPPMDispTIP4POMP)

#else

#ifndef LMP_PPPM_DISP_TIP4P_OMP_H
#define LMP_PPPM_DISP_TIP4P_OMP_H

#include "pppm_disp_tip4p.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

  class PPPMDispTIP4POMP : public PPPMDispTIP4P, public ThrOMP {
 public:
  PPPMDispTIP4POMP(class LAMMPS *);
  virtual ~PPPMDispTIP4POMP ();

 protected:
  virtual void allocate();

  virtual void compute_gf();
  virtual void compute_gf_6();

  virtual void compute(int,int);

  virtual void particle_map(double, double, double,
                            double, int **, int, int,
                            int, int, int, int, int, int);
  virtual void particle_map_c(double, double, double,
                              double, int **, int, int,
                              int, int, int, int, int, int);
  virtual void make_rho_c();  // XXX: not (yet) multi-threaded
  virtual void make_rho_g();
  virtual void make_rho_a();

  virtual void fieldforce_c_ik();
  virtual void fieldforce_c_ad();
  // virtual void fieldforce_peratom();  XXX: need to benchmark first.
  virtual void fieldforce_g_ik();
  virtual void fieldforce_g_ad();
  virtual void fieldforce_g_peratom();
  virtual void fieldforce_a_ik();
  virtual void fieldforce_a_ad();
  virtual void fieldforce_a_peratom();

  private:
  void compute_rho1d_thr(FFT_SCALAR * const * const, const FFT_SCALAR &,
                         const FFT_SCALAR &, const FFT_SCALAR &,
                         const int, FFT_SCALAR * const * const);
  void compute_drho1d_thr(FFT_SCALAR * const * const, const FFT_SCALAR &,
                          const FFT_SCALAR &, const FFT_SCALAR &,
                          const int, FFT_SCALAR * const * const);
  virtual void find_M_thr(int, int &, int &, dbl3_t &);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Kspace style pppm/tip4p/omp requires newton on

Self-explanatory.

*/
