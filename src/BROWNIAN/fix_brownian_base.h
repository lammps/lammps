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

#ifndef LMP_FIX_BROWNIAN_BASE_H
#define LMP_FIX_BROWNIAN_BASE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBrownianBase : public Fix {
 public:
  FixBrownianBase(class LAMMPS *, int, char **);
  ~FixBrownianBase() override;
  void init() override;
  int setmask() override;
  void reset_dt() override;

 protected:
  int seed;                  // RNG seed
  double dt, sqrtdt;         // time step interval and its sqrt
  int gamma_t_flag;          // 0/1 if isotropic translational damping is unset/set
  int gamma_r_flag;          // 0/1 if isotropic rotational damping is unset/set
  int gamma_t_eigen_flag;    // 0/1 if anisotropic translational damping is unset/set
  int gamma_r_eigen_flag;    // 0/1 if anisotropic rotational damping is unset/set
  int rot_temp_flag;         // 0/1 if rotational temperature is unset/set
  int planar_rot_flag;       // 0/1 if rotation is constrained to 2D (xy) plane

  double gamma_t, gamma_r;    // translational and rotational (isotropic) damping params
  double *gamma_t_inv;        // anisotropic damping parameter eigenvalues
  double *gamma_r_inv;
  double *gamma_t_invsqrt;
  double *gamma_r_invsqrt;

  int dipole_flag;        // set if dipole is used for asphere
  double *dipole_body;    // direction dipole is slaved to in body frame

  int noise_flag;             // 0/1 for noise off/on
  int gaussian_noise_flag;    // 0/1 for uniform/gaussian noise

  double temp;        // temperature
  double rot_temp;    // temperature
  double g1, g2;      // prefactors in time stepping

  class RanMars *rng;
};

}    // namespace LAMMPS_NS
#endif
