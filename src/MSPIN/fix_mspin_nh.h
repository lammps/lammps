/* ----------------------------------------------------------------------
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
   Contributing author: Akhlak Mahmood

   Contact:
     Department of Materials Science and Engineering,
     North Carolina State University,
     Raleigh, NC, USA

     amahmoo3@ncsu.edu; mahmoodakhlak@gmail.com
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(rigid/nh/mspin, FixMspinNH)
// clang-format on
#else

#ifndef LMP_FIX_MSPIN_NH_H
#define LMP_FIX_MSPIN_NH_H

#include "fix_rigid_nh.h"

namespace LAMMPS_NS {

class FixMspinNH : public FixRigidNH {
 public:
  FixMspinNH(class LAMMPS *, int, char **);
  virtual ~FixMspinNH();
  virtual int setmask();
  virtual void init();
  virtual void setup(int);
  virtual void post_force(int);
  virtual void compute_zeeman();
  virtual void compute_dipolar();
  virtual double compute_scalar();

  double extract_zeeman_pe();
  double extract_dipolar_pe();
  double extract_distance(int, int);

 protected:
  void calculate_dipoles(int);

  double **mu;
  double **dq;
  double *qm;
  int *qmcount;
  int qm_icustom;

  int nsum;    // total number of rigid atoms

  int zeeman_flag;
  int dipolar_flag;
  int uniform_field;
  double dipole_cutoff;

  double alpha;    // dipole interaction scaling factor
  double beta;     // zeeman+dipolar scaling factor

  double qb2f;
  double mu_0;    // force/Ampere^2 in Real
  double bxdx, bxdy, bxdz, bydx, bydy, bydz, bzdx, bzdy, bzdz;

  double zeeman_pe, dipolar_pe;
};
}    // namespace LAMMPS_NS

#endif
#endif
