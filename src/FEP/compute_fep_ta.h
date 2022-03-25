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

/* ----------------------------------------------------------------------
   Contributing author: Shifeng Ke (Zhejiang University)
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(fep/ta,ComputeFEPTA);
// clang-format on
#else

#ifndef COMPUTE_FEP_TA_H
#define COMPUTE_FEP_TA_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeFEPTA : public Compute {
 public:
  ComputeFEPTA(class LAMMPS *, int, char **); // compute ID groupID fep/ta temp xy/xz/yz scale_factor
  ~ComputeFEPTA() override;
  void init() override;
  void compute_vector() override;

 private:
  int pairflag;
  int tailflag;
  int fepinitflag;
  int eflag, vflag;
  double temp_fep;
  double scale_factor;
  int tan_axis1, tan_axis2, norm_axis;

  double boxlo_orig[3], boxhi_orig[3];
  double area_orig;

  int nmax;
  double **x_orig;
  double **f_orig;
  double eng_vdwl_orig, eng_coul_orig;
  double pvirial_orig[6];
  double *peatom_orig, **pvatom_orig;
  double energy_orig;
  double kvirial_orig[6];
  double *keatom_orig, **kvatom_orig;

  class Fix *fixgpu;

  double compute_epair();
  void change_box();
  void backup_box();
  void restore_box();
  void allocate_storage();
  void deallocate_storage();
  void backup_xfev();
  void restore_xfev();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Energy was not tallied on needed timestep

You are using a thermo keyword that requires potentials to
have tallied energy, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

*/
