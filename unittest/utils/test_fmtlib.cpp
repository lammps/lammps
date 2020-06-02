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
#include <exception>

using namespace LAMMPS_NS;
using ::testing::Eq;

// this tests a subset of {fmt} that is most relevant to LAMMPS

TEST(Fmtlib, insert_string) {
    const char val[] = "word";
    auto text = fmt::format("word {}",val);
    ASSERT_THAT(text, Eq("word word"));
}

TEST(Fmtlib, insert_int) {
    const int val = 333;
    auto text = fmt::format("word {}",val);
    ASSERT_THAT(text, Eq("word 333"));
}

TEST(Fmtlib, insert_neg_int) {
    const int val = -333;
    auto text = fmt::format("word {}",val);
    ASSERT_THAT(text, Eq("word -333"));
}

TEST(Fmtlib, insert_bigint) {
    if (sizeof(bigint) == 4) GTEST_SKIP();
    const bigint val = 9945234592L;
    auto text = fmt::format("word {}",val);
    ASSERT_THAT(text, Eq("word 9945234592"));
}

TEST(Fmtlib, insert_neg_bigint) {
    if (sizeof(bigint) == 4) GTEST_SKIP();
    const bigint val = -9945234592L;
    auto text = fmt::format("word {}",val);
    ASSERT_THAT(text, Eq("word -9945234592"));
}

TEST(Fmtlib, insert_tagint) {
    if (sizeof(tagint) == 4) GTEST_SKIP();
    const tagint val = 9945234592L;
    auto text = fmt::format("word {}",val);
    ASSERT_THAT(text, Eq("word 9945234592"));
}

TEST(Fmtlib, insert_neg_tagint) {
    if (sizeof(tagint) == 4) GTEST_SKIP();
    const tagint val = -9945234592L;
    auto text = fmt::format("word {}",val);
    ASSERT_THAT(text, Eq("word -9945234592"));
}

TEST(Fmtlib, insert_imageint) {
    if (sizeof(imageint) == 4) GTEST_SKIP();
    const imageint val = 9945234592L;
    auto text = fmt::format("word {}",val);
    ASSERT_THAT(text, Eq("word 9945234592"));
}

TEST(Fmtlib, insert_neg_imageint) {
    if (sizeof(imageint) == 4) GTEST_SKIP();
    const imageint val = -9945234592L;
    auto text = fmt::format("word {}",val);
    ASSERT_THAT(text, Eq("word -9945234592"));
}

TEST(Fmtlib, insert_double) {
    const double val = 1.5;
    auto text = fmt::format("word {}",val);
    ASSERT_THAT(text, Eq("word 1.5"));
}

TEST(Fmtlib, insert_neg_double) {
    const double val = -1.5;
    auto text = fmt::format("word {}",val);
    ASSERT_THAT(text, Eq("word -1.5"));
}

TEST(Fmtlib, int_for_double) {
    const double val = -1.5;
    ASSERT_THROW(fmt::format("word {:d}",val),std::exception);
}

TEST(Fmtlib, double_for_int) {
    const int val = 15;
    ASSERT_THROW(fmt::format("word {:g}",val),std::exception);
}
