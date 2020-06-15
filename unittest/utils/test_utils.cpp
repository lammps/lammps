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

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "utils.h"
#include <string>
#include <cerrno>
#include <cstdio>
#include <cstdlib>

using namespace LAMMPS_NS;
using ::testing::Eq;

TEST(Utils, trim_comment) {
    auto trimmed = utils::trim_comment("some text # comment");
    ASSERT_THAT(trimmed, Eq("some text "));
}

TEST(Utils, count_words) {
    ASSERT_EQ(utils::count_words("some text # comment"), 4);
}

TEST(Utils, count_words_non_default) {
    ASSERT_EQ(utils::count_words("some text # comment", " #"), 3);
}

TEST(Utils, trim_and_count_words) {
    ASSERT_EQ(utils::trim_and_count_words("some text # comment"), 2);
}

TEST(Utils, count_words_with_extra_spaces) {
    ASSERT_EQ(utils::count_words("   some text # comment   "), 4);
}

TEST(Utils, valid_integer1) {
    ASSERT_TRUE(utils::is_integer("10"));
}

TEST(Utils, valid_integer2) {
    ASSERT_TRUE(utils::is_integer("-10"));
}

TEST(Utils, valid_integer3) {
    ASSERT_TRUE(utils::is_integer("+10"));
}

TEST(Utils, valid_double1) {
    ASSERT_TRUE(utils::is_double("10.0"));
}

TEST(Utils, valid_double2) {
    ASSERT_TRUE(utils::is_double("1."));
}

TEST(Utils, valid_double3) {
    ASSERT_TRUE(utils::is_double(".0"));
}

TEST(Utils, valid_double4) {
    ASSERT_TRUE(utils::is_double("-10.0"));
}

TEST(Utils, valid_double5) {
    ASSERT_TRUE(utils::is_double("-1."));
}

TEST(Utils, valid_double6) {
    ASSERT_TRUE(utils::is_double("-.0"));
}

TEST(Utils, valid_double7) {
    ASSERT_TRUE(utils::is_double("+10.0"));
}

TEST(Utils, valid_double8) {
    ASSERT_TRUE(utils::is_double("+1."));
}

TEST(Utils, valid_double9) {
    ASSERT_TRUE(utils::is_double("+.0"));
}

TEST(Utils, empty_not_an_integer) {
    ASSERT_FALSE(utils::is_integer(""));
}

TEST(Utils, empty_not_a_double) {
    ASSERT_FALSE(utils::is_double(""));
}

TEST(Utils, text_not_an_integer) {
    ASSERT_FALSE(utils::is_integer("one"));
}

TEST(Utils, text_not_a_double) {
    ASSERT_FALSE(utils::is_double("half"));
}

TEST(Utils, double_not_an_integer1) {
    ASSERT_FALSE(utils::is_integer("10.0"));
}

TEST(Utils, double_not_an_integer2) {
    ASSERT_FALSE(utils::is_integer(".0"));
}

TEST(Utils, double_not_an_integer3) {
    ASSERT_FALSE(utils::is_integer("1."));
}

TEST(Utils, integer_is_double1) {
    ASSERT_TRUE(utils::is_double("10"));
}

TEST(Utils, integer_is_double2) {
    ASSERT_TRUE(utils::is_double("-10"));
}

TEST(Utils, is_double_with_exponential) {
    ASSERT_TRUE(utils::is_double("+1e02"));
}

TEST(Utils, is_double_with_neg_exponential) {
    ASSERT_TRUE(utils::is_double("1.0e-22"));
}

TEST(Utils, is_double_with_pos_exponential) {
    ASSERT_TRUE(utils::is_double(".1e+22"));
}

TEST(Utils, signed_double_and_exponential) {
    ASSERT_TRUE(utils::is_double("-10E-22"));
}

TEST(Utils, is_double_with_d_exponential) {
    ASSERT_FALSE(utils::is_double("10d22"));
}

TEST(Utils, is_double_with_neg_d_exponential) {
    ASSERT_FALSE(utils::is_double("10d-22"));
}

TEST(Utils, signed_double_and_d_exponential) {
    ASSERT_FALSE(utils::is_double("-10D-22"));
}

TEST(Utils, strmatch_beg) {
    ASSERT_TRUE(utils::strmatch("rigid/small/omp","^rigid"));
}

TEST(Utils, strmatch_mid1) {
    ASSERT_TRUE(utils::strmatch("rigid/small/omp","small"));
}

TEST(Utils, strmatch_mid2) {
    ASSERT_TRUE(utils::strmatch("rigid/small/omp","omp"));
}

TEST(Utils, strmatch_end) {
    ASSERT_TRUE(utils::strmatch("rigid/small/omp","/omp$"));
}

TEST(Utils, no_strmatch_beg) {
    ASSERT_FALSE(utils::strmatch("rigid/small/omp","^small"));
}

TEST(Utils, no_strmatch_mid) {
    ASSERT_FALSE(utils::strmatch("rigid/small/omp","none"));
}

TEST(Utils, no_strmatch_end) {
    ASSERT_FALSE(utils::strmatch("rigid/small/omp","/opt$"));
}

TEST(Utils, strmatch_whole_line) {
    ASSERT_TRUE(utils::strmatch("ITEM: UNITS\n","^\\s*ITEM: UNITS\\s*$"));
}

TEST(Utils, no_strmatch_whole_line) {
    ASSERT_FALSE(utils::strmatch("ITEM: UNITS\n","^\\s*ITEM: UNIT\\s*$"));
}

TEST(Utils, strmatch_integer_in_line) {
    ASSERT_TRUE(utils::strmatch(" 5  angles\n","^\\s*\\d+\\s+angles\\s"));
}

TEST(Utils, strmatch_float_in_line) {
    ASSERT_TRUE(utils::strmatch(" 5.0  angles\n","^\\s*\\f+\\s+angles\\s"));
}

TEST(Utils, strmatch_int_as_float_in_line) {
    ASSERT_TRUE(utils::strmatch(" 5  angles\n","^\\s*\\f+\\s+angles\\s"));
}

TEST(Utils, strmatch_char_range) {
    ASSERT_TRUE(utils::strmatch("rigid","^[ip-s]+gid"));
}

TEST(Utils, strmatch_opt_range) {
    ASSERT_TRUE(utils::strmatch("rigid","^[0-9]*[p-s]igid"));
}

TEST(Utils, path_join) {
#if defined(_WIN32)
    ASSERT_THAT(utils::path_join("c:\\parent\\folder", "filename"), Eq("c:\\parent\\folder\\filename"));
#else
    ASSERT_THAT(utils::path_join("/parent/folder", "filename"), Eq("/parent/folder/filename"));
#endif
}

TEST(Utils, path_basename) {
#if defined(_WIN32)
    ASSERT_THAT(utils::path_basename("c:\\parent\\folder\\filename"), Eq("filename"));
#else
    ASSERT_THAT(utils::path_basename("/parent/folder/filename"), Eq("filename"));
#endif
}

TEST(Utils, getsyserror) {
#if defined(__linux__)
    errno = ENOENT;
    std::string errmesg = utils::getsyserror();
    ASSERT_THAT(errmesg, Eq("No such file or directory"));
#else
    GTEST_SKIP();
#endif
}

TEST(Utils, potential_file) {
    FILE *fp;
    fp = fopen("ctest.txt","w");
    ASSERT_NE(fp,nullptr);
    fputs("# DATE: 2020-02-20 CONTRIBUTOR: Nessuno\n",fp);
    fclose(fp);

    EXPECT_TRUE(utils::file_is_readable("ctest.txt"));
    EXPECT_FALSE(utils::file_is_readable("no_such_file.txt"));

    EXPECT_THAT(utils::get_potential_file_path("ctest.txt"),Eq("ctest.txt"));
    const char *folder = getenv("LAMMPS_POTENTIALS");
    if (folder != nullptr) {
      std::string path=utils::path_join(folder,"Cu_u3.eam");
      EXPECT_THAT(utils::get_potential_file_path("Cu_u3.eam"),Eq(path));
    }

    EXPECT_THAT(utils::get_potential_date("ctest.txt","Test"),Eq("2020-02-20"));

    remove("ctest.txt");
}
