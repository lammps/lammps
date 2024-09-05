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
  ~FixRHEO() override;
  int setmask() override;
  void post_constructor() override;
  void init() override;
  void setup_pre_force(int) override;
  void setup(int) override;
  void pre_force(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void reset_dt() override;

  // Model parameters
  double cut;
  double *rho0, *csq;
  int self_mass_flag;
  int zmin_kernel, zmin_surface, zmin_splash;
  int kernel_style, surface_style;
  double divr_surface;

  // Accessory fixes/computes
  int thermal_flag;
  int rhosum_flag;
  int shift_flag;
  int interface_flag;
  int surface_flag;
  int oxidation_flag;

  int viscosity_fix_defined;
  int pressure_fix_defined;
  int thermal_fix_defined;
  int oxidation_fix_defined;

  class ComputeRHEOGrad *compute_grad;
  class ComputeRHEOKernel *compute_kernel;
  class ComputeRHEOInterface *compute_interface;
  class ComputeRHEOSurface *compute_surface;
  class ComputeRHEORhoSum *compute_rhosum;
  class ComputeRHEOVShift *compute_vshift;

 protected:
  double dtv, dtf;
};

namespace RHEO_NS {

  enum { QUINTIC, WENDLANDC4, RK0, RK1, RK2 };
  enum { COORDINATION, DIVR };

  // Status variables
  enum Status {
    // Phase status
    STATUS_SOLID = 1 << 0,
    // Gap for future phase: STATUS_ = 1 << 1,

    // Surface status
    STATUS_BULK = 1 << 2,
    STATUS_LAYER = 1 << 3,
    STATUS_SURFACE = 1 << 4,
    STATUS_SPLASH = 1 << 5,

    // Temporary status options - reset in preforce
    STATUS_NO_SHIFT = 1 << 6,
    STATUS_NO_INTEGRATION = 1 << 7,
    STATUS_FREEZING = 1 << 8,
    STATUS_MELTING = 1 << 9
  };

  // Masks and their inverses
  enum {
    PHASECHECK = 0x00000003,      // 00000000000000000000000000000011
    SURFACECHECK = 0x0000003C,    // 00000000000000000000000000111100
    OPTIONSMASK = 0xFFFFFC3F,     // 11111111111111111111110000111111
    SURFACEMASK = 0xFFFFFFC3,     // 11111111111111111111111111000011
    PHASEMASK = 0xFFFFFFFC        // 11111111111111111111111111111100
  };
}    // namespace RHEO_NS
}    // namespace LAMMPS_NS

#endif
#endif
