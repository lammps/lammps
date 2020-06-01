/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "fmt/format.h"
#include <string>

using namespace LAMMPS_NS;
using ::testing::Eq;

// this tests a subset of {fmt} that is most relevant to LAMMPS

TEST(Fmtlib, insert_string) {
    const char word[] = "word";
    auto text = fmt::format("word {}",word);
    ASSERT_THAT(text, Eq("word word"));
}

TEST(Fmtlib, insert_int) {
    const int word = 333;
    auto text = fmt::format("word {}",word);
    ASSERT_THAT(text, Eq("word 333"));
}

TEST(Fmtlib, insert_neg_int) {
    const int word = -333;
    auto text = fmt::format("word {}",word);
    ASSERT_THAT(text, Eq("word -333"));
}

TEST(Fmtlib, insert_double) {
    const double word = 1.5;
    auto text = fmt::format("word {}",word);
    ASSERT_THAT(text, Eq("word 1.5"));
}

TEST(Fmtlib, insert_neg_double) {
    const double word = -1.5;
    auto text = fmt::format("word {}",word);
    ASSERT_THAT(text, Eq("word 1.5"));
}
