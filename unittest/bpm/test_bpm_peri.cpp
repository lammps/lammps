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

// Verification tests for the BPM peridynamics styles (pair_style bpm/peri and
// bond_style bpm/peri, models pmb/lps/ves/eps).  All test logic is shared via
// test_bpm_common; the reference systems live in tests/bpmperi-*.yaml.  The
// driver exercises the pair and bond styles together, comparing the per-atom
// trajectory and any enabled analytic model against the recorded reference.

#include "test_bpm_common.h"
#include "test_main.h"

#include "info.h"

#include "gtest/gtest.h"

using namespace LAMMPS_NS;

// BPM bond styles mandate 'newton bond off'.  The pair newton flag is varied:
// newton_on runs 'newton on off', newton_off runs 'newton off off'.

TEST(BpmPeri, newton_on)
{
    if (!Info::has_package("BPM")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    run_bpm_trajectory_test(true, "newton on");
}

TEST(BpmPeri, newton_off)
{
    if (!Info::has_package("BPM")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    run_bpm_trajectory_test(false, "newton off");
}
