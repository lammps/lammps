/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(pppm/disp/dielectric,PPPMDispDielectric);
// clang-format on
#else

#ifndef LMP_PPPM_DISP_DIELECTRIC_H
#define LMP_PPPM_DISP_DIELECTRIC_H

#include "pppm_disp.h"

namespace LAMMPS_NS {

class PPPMDispDielectric : public PPPMDisp {
 public:
  PPPMDispDielectric(class LAMMPS *);
  virtual ~PPPMDispDielectric();
  virtual double memory_usage();
  virtual void compute(int, int);
  void qsum_qsq();
  void slabcorr(int);

  double **efield;
  double *phi;
  int potflag;    // 1/0 if per-atom electrostatic potential phi is needed

 protected:
  virtual void fieldforce_c_ik();
  virtual void fieldforce_c_ad();
  virtual void fieldforce_c_peratom();

  class AtomVecDielectric *avec;
  int mu_flag;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
