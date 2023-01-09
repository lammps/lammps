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

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(pppm/disp/omp,PPPMDispOMP);
// clang-format on
#else

#ifndef LMP_PPPM_DISP_OMP_H
#define LMP_PPPM_DISP_OMP_H

#include "pppm_disp.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PPPMDispOMP : public PPPMDisp, public ThrOMP {
 public:
  PPPMDispOMP(class LAMMPS *);
  ~PPPMDispOMP() override;
  void compute(int, int) override;

 protected:
  void allocate() override;

  virtual void compute_gf();
  virtual void compute_gf_6();

  void particle_map(double, double, double, double, int **, int, int, int, int, int, int, int,
                    int) override;

  void fieldforce_c_ik() override;
  void fieldforce_c_ad() override;
  void fieldforce_c_peratom() override;
  void fieldforce_g_ik() override;
  void fieldforce_g_ad() override;
  void fieldforce_g_peratom() override;
  void fieldforce_a_ik() override;
  void fieldforce_a_ad() override;
  void fieldforce_a_peratom() override;

  void make_rho_c() override;
  void make_rho_g() override;
  void make_rho_a() override;

  void compute_rho1d_thr(FFT_SCALAR *const *const, const FFT_SCALAR &, const FFT_SCALAR &,
                         const FFT_SCALAR &, const int, FFT_SCALAR *const *const);
  void compute_drho1d_thr(FFT_SCALAR *const *const, const FFT_SCALAR &, const FFT_SCALAR &,
                          const FFT_SCALAR &, const int, FFT_SCALAR *const *const);
  //  void compute_rho_coeff();
  //  void slabcorr(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
