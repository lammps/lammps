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

#ifndef TEST_BPM_COMMON_H
#define TEST_BPM_COMMON_H

#include <string>

// Shared trajectory-comparison body for the BPM pair+bond test drivers.  Builds
// a LAMMPS instance from the global test_config, applies the configured
// pair_style/pair_coeff and bond_style/bond_coeff, runs the 'run_segments' and
// compares per-atom positions, velocities and forces -- plus any enabled
// analytic model -- against the recorded reference.  The harness is generic over
// BPM (the styles come from the YAML), so the same driver serves bpm/peri,
// bpm/spring, etc., and is extensible to styles layered on BPM such as RHEO.
// BPM bond styles mandate 'newton bond off'; the per-driver TEST() fixtures pass
// that as the newton flag.
void run_bpm_trajectory_test(bool newton, const std::string &label);

#endif
