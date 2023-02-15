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

#ifdef FIX_CLASS
// clang-format off
FixStyle(rheo,FixRHEO)
// clang-format on
#else

#ifndef LMP_FIX_RHEO_H
#define LMP_FIX_RHEO_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRHEO : public Fix {
 public:
  FixRHEO(class LAMMPS *, int, char **);
  virtual ~FixRHEO();
  int setmask();
  virtual void post_constructor();
  virtual void init();
  virtual void setup_pre_force(int);
  virtual void pre_force(int);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  void reset_dt();

  int kernel_style;
  int thermal_flag;
  int rhosum_flag;
  int shift_flag;
  int solid_flag;

  int thermal_fix_defined;
  int viscosity_fix_defined;
  int pressure_fix_defined;

  int *status, *surface;
  double *conductivity, *viscosity, *pressure;
  double **f_pressure;

  class ComputeRHEOGrad *compute_grad;
  class ComputeRHEOKernel *compute_kernel;
  class ComputeRHEOInterface *compute_interface;
  class ComputeRHEORhoSum *compute_rhosum;
  class ComputeRHEOVShift *compute_vshift;

  enum {QUINTIC, CRK0, CRK1, CRK2};
  enum {LINEAR, CUBIC, TAITWATER};

  enum {
    // Phase status
    STATUS_FLUID = 1 << 0,
    STATUS_REACTIVE = 1 << 1,
    STATUS_SOLID = 1 << 2,
    STATUS_FREEZING = 1 << 3

    // Temporary status options - reset in preforce
    STATUS_SHIFT = 1 << 4,
    STATUS_NO_FORCE = 1 << 5,

    // Surface status
    STATUS_BULK = 1 << 6,
    STATUS_LAYER = 1 << 7,
    STATUS_SURFACE = 1 << 8,
    STATUS_SPLASH = 1 << 9,
  };

 protected:
  double cut, rho0, csq;
  int zmin_kernel, rhosum_zmin;

  double dtv, dtf;
};

}    // namespace LAMMPS_NS

#endif
#endif
