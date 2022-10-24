/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    *** Smooth Mach Dynamics ***
 *
 * This file is part of the MACHDYN package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
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
FixStyle(smd/move_tri_surf,FixSMDMoveTriSurf);
// clang-format on
#else

#ifndef LMP_FIX_SMD_INTEGRATE_TRIANGULAR_SURFACE_H
#define LMP_FIX_SMD_INTEGRATE_TRIANGULAR_SURFACE_H

#include "fix.h"
#include <Eigen/Eigen>

namespace LAMMPS_NS {

class FixSMDMoveTriSurf : public Fix {
 public:
  FixSMDMoveTriSurf(class LAMMPS *, int, char **);

  int setmask() override;
  void init() override;
  void initial_integrate(int) override;
  void reset_dt() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 protected:
  double dtv;
  bool linearFlag, wiggleFlag, rotateFlag;
  double vx, vy, vz;
  Eigen::Vector3d rotation_axis, origin;
  double rotation_period;
  Eigen::Matrix3d u_cross, uxu;
  double wiggle_travel, wiggle_max_travel, wiggle_direction;
};

}    // namespace LAMMPS_NS

#endif
#endif
