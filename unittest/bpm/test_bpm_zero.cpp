/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// Verification tests for bond_style bpm/zero (the force-free BPM bond style used
// for topology/breaking debugging).  This reuses the shared test_bpm_common
// driver unchanged -- only the YAML names the styles -- which doubles as a check
// that the harness is generic over the BPM package and not specific to bpm/peri.
// The reference systems live in tests/bpmzero-*.yaml.

#include "test_bpm_common.h"
#include "test_main.h"

#include "info.h"

#include "gtest/gtest.h"

using namespace LAMMPS_NS;

// BPM bond styles mandate 'newton bond off'.  The pair newton flag is varied:
// newton_on runs 'newton on off', newton_off runs 'newton off off'.

TEST(BpmZero, newton_on)
{
    if (!Info::has_package("BPM")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    run_bpm_trajectory_test(true, "newton on");
}

TEST(BpmZero, newton_off)
{
    if (!Info::has_package("BPM")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    run_bpm_trajectory_test(false, "newton off");
}
