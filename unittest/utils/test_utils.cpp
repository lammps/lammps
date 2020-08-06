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

#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace LAMMPS_NS;
using ::testing::EndsWith;
using ::testing::Eq;
using ::testing::StrEq;

TEST(Utils, trim)
{
    auto trimmed = utils::trim("\t some text");
    ASSERT_THAT(trimmed, StrEq("some text"));

    trimmed = utils::trim("some text \r\n");
    ASSERT_THAT(trimmed, StrEq("some text"));

    trimmed = utils::trim("\v some text \f");
    ASSERT_THAT(trimmed, StrEq("some text"));

    trimmed = utils::trim("   some\t text    ");
    ASSERT_THAT(trimmed, StrEq("some\t text"));

    trimmed = utils::trim("  \t\n   ");
    ASSERT_THAT(trimmed, StrEq(""));
}

TEST(Utils, trim_comment)
{
    auto trimmed = utils::trim_comment("some text # comment");
    ASSERT_THAT(trimmed, StrEq("some text "));
}

TEST(Utils, count_words)
{
    ASSERT_EQ(utils::count_words("some text # comment"), 4);
}

TEST(Utils, count_words_string)
{
    ASSERT_EQ(utils::count_words(std::string("some text # comment")), 4);
}

TEST(Utils, count_words_non_default)
{
    ASSERT_EQ(utils::count_words("some text # comment", " #"), 3);
}

TEST(Utils, trim_and_count_words)
{
    ASSERT_EQ(utils::trim_and_count_words("some text # comment"), 2);
}

TEST(Utils, count_words_with_extra_spaces)
{
    ASSERT_EQ(utils::count_words("   some text # comment   "), 4);
}

TEST(Utils, split_words_simple)
{
    std::vector<std::string> list = utils::split_words("one two three");
    ASSERT_EQ(list.size(), 3);
    ASSERT_THAT(list[0], StrEq("one"));
    ASSERT_THAT(list[1], StrEq("two"));
    ASSERT_THAT(list[2], StrEq("three"));
}

TEST(Utils, split_words_quoted)
{
    std::vector<std::string> list = utils::split_words("one 'two' \"three\"");
    ASSERT_EQ(list.size(), 3);
    ASSERT_THAT(list[0], StrEq("one"));
    ASSERT_THAT(list[1], StrEq("two"));
    ASSERT_THAT(list[2], StrEq("three"));
}

TEST(Utils, split_words_escaped)
{
    std::vector<std::string> list = utils::split_words("1\\' '\"two\"' 3\\\"");
    ASSERT_EQ(list.size(), 3);
    ASSERT_THAT(list[0], StrEq("1\\'"));
    ASSERT_THAT(list[1], StrEq("\"two\""));
    ASSERT_THAT(list[2], StrEq("3\\\""));
}

TEST(Utils, split_words_quote_in_quoted)
{
    std::vector<std::string> list = utils::split_words("one 't\\'wo' \"th\\\"ree\"");
    ASSERT_EQ(list.size(), 3);
    ASSERT_THAT(list[0], StrEq("one"));
    ASSERT_THAT(list[1], StrEq("t\\'wo"));
    ASSERT_THAT(list[2], StrEq("th\\\"ree"));
}

TEST(Utils, valid_integer1)
{
    ASSERT_TRUE(utils::is_integer("10"));
}

TEST(Utils, valid_integer2)
{
    ASSERT_TRUE(utils::is_integer("-10"));
}

TEST(Utils, valid_integer3)
{
    ASSERT_TRUE(utils::is_integer("+10"));
}

TEST(Utils, valid_double1)
{
    ASSERT_TRUE(utils::is_double("10.0"));
}

TEST(Utils, valid_double2)
{
    ASSERT_TRUE(utils::is_double("1."));
}

TEST(Utils, valid_double3)
{
    ASSERT_TRUE(utils::is_double(".0"));
}

TEST(Utils, valid_double4)
{
    ASSERT_TRUE(utils::is_double("-10.0"));
}

TEST(Utils, valid_double5)
{
    ASSERT_TRUE(utils::is_double("-1."));
}

TEST(Utils, valid_double6)
{
    ASSERT_TRUE(utils::is_double("-.0"));
}

TEST(Utils, valid_double7)
{
    ASSERT_TRUE(utils::is_double("+10.0"));
}

TEST(Utils, valid_double8)
{
    ASSERT_TRUE(utils::is_double("+1."));
}

TEST(Utils, valid_double9)
{
    ASSERT_TRUE(utils::is_double("+.0"));
}

TEST(Utils, empty_not_an_integer)
{
    ASSERT_FALSE(utils::is_integer(""));
}

TEST(Utils, empty_not_a_double)
{
    ASSERT_FALSE(utils::is_double(""));
}

TEST(Utils, text_not_an_integer)
{
    ASSERT_FALSE(utils::is_integer("one"));
}

TEST(Utils, text_not_a_double)
{
    ASSERT_FALSE(utils::is_double("half"));
}

TEST(Utils, double_not_an_integer1)
{
    ASSERT_FALSE(utils::is_integer("10.0"));
}

TEST(Utils, double_not_an_integer2)
{
    ASSERT_FALSE(utils::is_integer(".0"));
}

TEST(Utils, double_not_an_integer3)
{
    ASSERT_FALSE(utils::is_integer("1."));
}

TEST(Utils, integer_is_double1)
{
    ASSERT_TRUE(utils::is_double("10"));
}

TEST(Utils, integer_is_double2)
{
    ASSERT_TRUE(utils::is_double("-10"));
}

TEST(Utils, is_double_with_exponential)
{
    ASSERT_TRUE(utils::is_double("+1e02"));
}

TEST(Utils, is_double_with_neg_exponential)
{
    ASSERT_TRUE(utils::is_double("1.0e-22"));
}

TEST(Utils, is_double_with_pos_exponential)
{
    ASSERT_TRUE(utils::is_double(".1e+22"));
}

TEST(Utils, signed_double_and_exponential)
{
    ASSERT_TRUE(utils::is_double("-10E-22"));
}

TEST(Utils, is_double_with_d_exponential)
{
    ASSERT_FALSE(utils::is_double("10d22"));
}

TEST(Utils, is_double_with_neg_d_exponential)
{
    ASSERT_FALSE(utils::is_double("10d-22"));
}

TEST(Utils, signed_double_and_d_exponential)
{
    ASSERT_FALSE(utils::is_double("-10D-22"));
}

TEST(Utils, strmatch_beg)
{
    ASSERT_TRUE(utils::strmatch("rigid/small/omp", "^rigid"));
}

TEST(Utils, strmatch_mid1)
{
    ASSERT_TRUE(utils::strmatch("rigid/small/omp", "small"));
}

TEST(Utils, strmatch_mid2)
{
    ASSERT_TRUE(utils::strmatch("rigid/small/omp", "omp"));
}

TEST(Utils, strmatch_end)
{
    ASSERT_TRUE(utils::strmatch("rigid/small/omp", "/omp$"));
}

TEST(Utils, no_strmatch_beg)
{
    ASSERT_FALSE(utils::strmatch("rigid/small/omp", "^small"));
}

TEST(Utils, no_strmatch_mid)
{
    ASSERT_FALSE(utils::strmatch("rigid/small/omp", "none"));
}

TEST(Utils, no_strmatch_end)
{
    ASSERT_FALSE(utils::strmatch("rigid/small/omp", "/opt$"));
}

TEST(Utils, strmatch_whole_line)
{
    ASSERT_TRUE(utils::strmatch("ITEM: UNITS\n", "^\\s*ITEM: UNITS\\s*$"));
}

TEST(Utils, no_strmatch_whole_line)
{
    ASSERT_FALSE(utils::strmatch("ITEM: UNITS\n", "^\\s*ITEM: UNIT\\s*$"));
}

TEST(Utils, strmatch_integer_in_line)
{
    ASSERT_TRUE(utils::strmatch(" 5  angles\n", "^\\s*\\d+\\s+angles\\s"));
}

TEST(Utils, strmatch_float_in_line)
{
    ASSERT_TRUE(utils::strmatch(" 5.0  angles\n", "^\\s*\\f+\\s+angles\\s"));
}

TEST(Utils, strmatch_int_as_float_in_line)
{
    ASSERT_TRUE(utils::strmatch(" 5  angles\n", "^\\s*\\f+\\s+angles\\s"));
}

TEST(Utils, strmatch_char_range)
{
    ASSERT_TRUE(utils::strmatch("rigid", "^[ip-s]+gid"));
}

TEST(Utils, strmatch_notchar_range)
{
    ASSERT_TRUE(utils::strmatch("rigid", "^[^a-g]+gid"));
}

TEST(Utils, strmatch_backslash)
{
    ASSERT_TRUE(utils::strmatch("\\rigid", "^\\W\\w+gid"));
}

TEST(Utils, strmatch_opt_range)
{
    ASSERT_TRUE(utils::strmatch("rigid", "^[0-9]*[\\Wp-s]igid"));
}

TEST(Utils, strmatch_opt_char)
{
    ASSERT_TRUE(utils::strmatch("rigid", "^r?igid"));
    ASSERT_TRUE(utils::strmatch("igid", "^r?igid"));
}

TEST(Utils, strmatch_dot)
{
    ASSERT_TRUE(utils::strmatch("rigid", ".igid"));
    ASSERT_TRUE(utils::strmatch("Rigid", ".igid"));
}

TEST(Utils, strmatch_digit_nondigit)
{
    ASSERT_TRUE(utils::strmatch(" 5  angles\n", "^\\s*\\d+\\s+\\D+\\s"));
}

TEST(Utils, strmatch_integer_noninteger)
{
    ASSERT_TRUE(utils::strmatch(" -5  angles\n", "^\\s*\\i+\\s+\\I+\\s"));
}

TEST(Utils, strmatch_float_nonfloat)
{
    ASSERT_TRUE(utils::strmatch(" 5.0  angls\n", "^\\s*\\f+\\s+\\F+\\s"));
}

TEST(Utils, strmatch_whitespace_nonwhitespace)
{
    ASSERT_TRUE(utils::strmatch(" 5.0  angles\n", "^\\s*\\S+\\s+\\S+\\s"));
}

TEST(Utils, guesspath)
{
    char buf[256];
    FILE *fp = fopen("test_guesspath.txt", "w");
#if defined(__linux__)
    const char *path = utils::guesspath(buf, sizeof(buf), fp);
    ASSERT_THAT(path, EndsWith("test_guesspath.txt"));
#else
    const char *path = utils::guesspath(buf, sizeof(buf), fp);
    ASSERT_THAT(path, EndsWith("(unknown)"));
#endif
    fclose(fp);
}

TEST(Utils, path_join)
{
#if defined(_WIN32)
    ASSERT_THAT(utils::path_join("c:\\parent\\folder", "filename"),
                Eq("c:\\parent\\folder\\filename"));
#else
    ASSERT_THAT(utils::path_join("/parent/folder", "filename"), Eq("/parent/folder/filename"));
#endif
}

TEST(Utils, path_basename)
{
#if defined(_WIN32)
    ASSERT_THAT(utils::path_basename("c:\\parent\\folder\\filename"), Eq("filename"));
#else
    ASSERT_THAT(utils::path_basename("/parent/folder/filename"), Eq("filename"));
#endif
}

TEST(Utils, getsyserror)
{
#if defined(__linux__)
    errno               = ENOENT;
    std::string errmesg = utils::getsyserror();
    ASSERT_THAT(errmesg, Eq("No such file or directory"));
#else
    GTEST_SKIP();
#endif
}

TEST(Utils, potential_file)
{
    FILE *fp;
    fp = fopen("ctest1.txt", "w");
    ASSERT_NE(fp, nullptr);
    fputs("# DATE: 2020-02-20 UNITS: real CONTRIBUTOR: Nessuno\n", fp);
    fclose(fp);
    fp = fopen("ctest2.txt", "w");
    ASSERT_NE(fp, nullptr);
    fputs("# CONTRIBUTOR: Pippo\n", fp);
    fclose(fp);

    ASSERT_TRUE(utils::file_is_readable("ctest1.txt"));
    ASSERT_TRUE(utils::file_is_readable("ctest2.txt"));
    ASSERT_FALSE(utils::file_is_readable("no_such_file.txt"));

    ASSERT_THAT(utils::get_potential_file_path("ctest1.txt"), Eq("ctest1.txt"));
    ASSERT_THAT(utils::get_potential_file_path("no_such_file.txt"), Eq(""));

    const char *folder = getenv("LAMMPS_POTENTIALS");
    if (folder != nullptr) {
        std::string path = utils::path_join(folder, "Cu_u3.eam");
        EXPECT_THAT(utils::get_potential_file_path("Cu_u3.eam"), Eq(path));
        EXPECT_THAT(utils::get_potential_units(path, "EAM"), Eq("metal"));
    }

    ASSERT_THAT(utils::get_potential_date("ctest1.txt", "Test"), Eq("2020-02-20"));
    ASSERT_THAT(utils::get_potential_units("ctest1.txt", "Test"), Eq("real"));
    ASSERT_THAT(utils::get_potential_date("ctest2.txt", "Test"), Eq(""));
    ASSERT_THAT(utils::get_potential_units("ctest2.txt", "Test"), Eq(""));

    remove("ctest1.txt");
    remove("ctest2.txt");
}

TEST(Utils, unit_conversion)
{
    double factor;
    int flag;

    flag = utils::get_supported_conversions(utils::UNKNOWN);
    ASSERT_EQ(flag, utils::NOCONVERT);
    flag = utils::get_supported_conversions(utils::ENERGY);
    ASSERT_EQ(flag, utils::METAL2REAL | utils::REAL2METAL);

    factor = utils::get_conversion_factor(utils::UNKNOWN, (1 << 30) - 1);
    ASSERT_DOUBLE_EQ(factor, 0.0);
    factor = utils::get_conversion_factor(utils::UNKNOWN, utils::NOCONVERT);
    ASSERT_DOUBLE_EQ(factor, 0.0);
    factor = utils::get_conversion_factor(utils::ENERGY, utils::NOCONVERT);
    ASSERT_DOUBLE_EQ(factor, 1.0);
    factor = utils::get_conversion_factor(utils::ENERGY, utils::METAL2REAL);
    ASSERT_DOUBLE_EQ(factor, 23.060549);
    factor = utils::get_conversion_factor(utils::ENERGY, utils::REAL2METAL);
    ASSERT_DOUBLE_EQ(factor, 1.0 / 23.060549);
}
