/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef TEST_DEM_COMMON_H
#define TEST_DEM_COMMON_H

#include <string>

// Shared trajectory-comparison body for all DEM test drivers.  Builds a
// LAMMPS instance from the global test_config (newton on or off), runs the
// configured 'run_segments' and compares per-atom positions, velocities,
// torques and angular velocities -- plus any enabled analytic model --
// against the recorded reference.  The per-driver TEST() fixtures only set
// up the package/skip guards and then call this with the newton flag.
void run_dem_trajectory_test(bool newton, const std::string &label);

#endif
