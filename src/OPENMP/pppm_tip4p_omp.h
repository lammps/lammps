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
KSpaceStyle(pppm/tip4p/omp,PPPMTIP4POMP);
// clang-format on
#else

#ifndef LMP_PPPM_TIP4P_OMP_H
#define LMP_PPPM_TIP4P_OMP_H

#include "pppm_tip4p.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PPPMTIP4POMP : public PPPMTIP4P, public ThrOMP {
 public:
  PPPMTIP4POMP(class LAMMPS *);
  ~PPPMTIP4POMP() override;
  void compute(int, int) override;

 protected:
  void allocate() override;

  void compute_gf_ik() override;
  void compute_gf_ad() override;

  void particle_map() override;
  void make_rho() override;    // XXX: not (yet) multi-threaded

  void fieldforce_ik() override;
  void fieldforce_ad() override;
  // virtual void fieldforce_peratom();  XXX: need to benchmark first.

 private:
  void compute_rho1d_thr(FFT_SCALAR *const *const, const FFT_SCALAR &, const FFT_SCALAR &,
                         const FFT_SCALAR &);
  void compute_drho1d_thr(FFT_SCALAR *const *const, const FFT_SCALAR &, const FFT_SCALAR &,
                          const FFT_SCALAR &);

  void find_M_thr(const int, int &, int &, dbl3_t &);

  //  void slabcorr(int);  // XXX: not (yet) multi-threaded
};

}    // namespace LAMMPS_NS

#endif
#endif
