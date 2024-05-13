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
KSpaceStyle(pppm/omp,PPPMOMP);
// clang-format on
#else

#ifndef LMP_PPPM_OMP_H
#define LMP_PPPM_OMP_H

#include "pppm.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PPPMOMP : public PPPM, public ThrOMP {
 public:
  PPPMOMP(class LAMMPS *);
  ~PPPMOMP() override;
  void compute(int, int) override;

 protected:
  void allocate() override;

  void compute_gf_ik() override;
  void compute_gf_ad() override;

  void make_rho() override;
  void fieldforce_ik() override;
  void fieldforce_ad() override;
  void fieldforce_peratom() override;

 private:
  void compute_rho1d_thr(FFT_SCALAR *const *const, const FFT_SCALAR &, const FFT_SCALAR &,
                         const FFT_SCALAR &);
  void compute_drho1d_thr(FFT_SCALAR *const *const, const FFT_SCALAR &, const FFT_SCALAR &,
                          const FFT_SCALAR &);
  //  void slabcorr(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
