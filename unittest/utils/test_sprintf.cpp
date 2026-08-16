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

#include "lmptype.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>

using namespace LAMMPS_NS;
using ::testing::Eq;
using ::testing::StartsWith;

// this tests the printf() conversions supported by utils::sprintf() that
// are most relevant to LAMMPS, i.e. for dump or thermo style output

TEST(Sprintf, plain_text)
{
    auto text = utils::sprintf("some plain text");
    ASSERT_THAT(text, Eq("some plain text"));
}

TEST(Sprintf, empty_format)
{
    auto text = utils::sprintf("");
    ASSERT_THAT(text, Eq(""));
}

TEST(Sprintf, percent_literal)
{
    auto text = utils::sprintf("100%% done");
    ASSERT_THAT(text, Eq("100% done"));
}

TEST(Sprintf, percent_literal_with_args)
{
    auto text = utils::sprintf("%d%% of %d", 50, 200);
    ASSERT_THAT(text, Eq("50% of 200"));
}

TEST(Sprintf, insert_int)
{
    constexpr int val = 333;
    auto text         = utils::sprintf("word %d", val);
    ASSERT_THAT(text, Eq("word 333"));
    text = utils::sprintf("word %i", val);
    ASSERT_THAT(text, Eq("word 333"));
}

TEST(Sprintf, insert_neg_int)
{
    constexpr int val = -333;
    auto text         = utils::sprintf("word %d", val);
    ASSERT_THAT(text, Eq("word -333"));
}

TEST(Sprintf, insert_int_extremes)
{
    auto text = utils::sprintf("%d %d", std::numeric_limits<int>::max(),
                               std::numeric_limits<int>::min());
    ASSERT_THAT(text, Eq("2147483647 -2147483648"));
}

TEST(Sprintf, int_flags_and_width)
{
    ASSERT_THAT(utils::sprintf("%8d", 333), Eq("     333"));
    ASSERT_THAT(utils::sprintf("%-8d|", 333), Eq("333     |"));
    ASSERT_THAT(utils::sprintf("%08d", 333), Eq("00000333"));
    ASSERT_THAT(utils::sprintf("%+d", 333), Eq("+333"));
    ASSERT_THAT(utils::sprintf("% d", 333), Eq(" 333"));
    ASSERT_THAT(utils::sprintf("%2d", 333), Eq("333"));
}

TEST(Sprintf, insert_unsigned)
{
    constexpr unsigned int val = 4294967295U;
    auto text                  = utils::sprintf("word %u", val);
    ASSERT_THAT(text, Eq("word 4294967295"));
}

TEST(Sprintf, insert_octal_hex)
{
    ASSERT_THAT(utils::sprintf("%o", 8), Eq("10"));
    ASSERT_THAT(utils::sprintf("%#o", 8), Eq("010"));
    ASSERT_THAT(utils::sprintf("%x", 255), Eq("ff"));
    ASSERT_THAT(utils::sprintf("%X", 255), Eq("FF"));
    ASSERT_THAT(utils::sprintf("%#x", 255), Eq("0xff"));
    ASSERT_THAT(utils::sprintf("%08x", 48879), Eq("0000beef"));
}

TEST(Sprintf, insert_char)
{
    ASSERT_THAT(utils::sprintf("%c", 'A'), Eq("A"));
    ASSERT_THAT(utils::sprintf("%c%c%c", 'x', 'y', 'z'), Eq("xyz"));
    ASSERT_THAT(utils::sprintf("%c", 65), Eq("A"));
}

TEST(Sprintf, insert_cstring)
{
    constexpr char val[] = "word";
    auto text            = utils::sprintf("word %s", val);
    ASSERT_THAT(text, Eq("word word"));
    const char *ptr = val;
    text            = utils::sprintf("word %s", ptr);
    ASSERT_THAT(text, Eq("word word"));
}

TEST(Sprintf, insert_stdstring)
{
    const std::string val = "word";
    auto text             = utils::sprintf("word %s", val);
    ASSERT_THAT(text, Eq("word word"));
}

TEST(Sprintf, string_flags_and_width)
{
    const std::string val = "word";
    ASSERT_THAT(utils::sprintf("%10s", val), Eq("      word"));
    ASSERT_THAT(utils::sprintf("%-10s|", val), Eq("word      |"));
    ASSERT_THAT(utils::sprintf("%.2s", val), Eq("wo"));
}

TEST(Sprintf, insert_long_long)
{
    constexpr long long val = 9945234592LL;
    auto text               = utils::sprintf("word %lld", val);
    ASSERT_THAT(text, Eq("word 9945234592"));
}

TEST(Sprintf, insert_size_t)
{
    auto text = utils::sprintf("word %zu", sizeof(double));
    ASSERT_THAT(text, Eq("word 8"));
}

TEST(Sprintf, insert_bigint)
{
#if defined(LAMMPS_BIGBIG) || defined(LAMMPS_SMALLBIG)
    constexpr bigint val = 9945234592L;
    auto text            = utils::sprintf("word " BIGINT_FORMAT, val);
    ASSERT_THAT(text, Eq("word 9945234592"));
#else
    GTEST_SKIP();
#endif
}

TEST(Sprintf, insert_neg_bigint)
{
#if defined(LAMMPS_BIGBIG) || defined(LAMMPS_SMALLBIG)
    constexpr bigint val = -9945234592L;
    auto text            = utils::sprintf("word " BIGINT_FORMAT, val);
    ASSERT_THAT(text, Eq("word -9945234592"));
#else
    GTEST_SKIP();
#endif
}

TEST(Sprintf, insert_tagint)
{
#if defined(LAMMPS_BIGBIG)
    constexpr tagint val = 9945234592L;
    auto text            = utils::sprintf("word " TAGINT_FORMAT, val);
    ASSERT_THAT(text, Eq("word 9945234592"));
#else
    GTEST_SKIP();
#endif
}

TEST(Sprintf, insert_neg_tagint)
{
#if defined(LAMMPS_BIGBIG)
    constexpr tagint val = -9945234592L;
    auto text            = utils::sprintf("word " TAGINT_FORMAT, val);
    ASSERT_THAT(text, Eq("word -9945234592"));
#else
    GTEST_SKIP();
#endif
}

TEST(Sprintf, insert_double)
{
    constexpr double val = 1.5;
    auto text            = utils::sprintf("word %f", val);
    ASSERT_THAT(text, Eq("word 1.500000"));
}

TEST(Sprintf, insert_neg_double)
{
    constexpr double val = -1.5;
    auto text            = utils::sprintf("word %f", val);
    ASSERT_THAT(text, Eq("word -1.500000"));
}

TEST(Sprintf, insert_float)
{
    constexpr float val = 1.5f;
    auto text           = utils::sprintf("word %.1f", val);
    ASSERT_THAT(text, Eq("word 1.5"));
}

TEST(Sprintf, double_flags_and_width)
{
    ASSERT_THAT(utils::sprintf("%.2f", 1.5), Eq("1.50"));
    ASSERT_THAT(utils::sprintf("%8.3f", 1.5), Eq("   1.500"));
    ASSERT_THAT(utils::sprintf("%-8.3f|", 1.5), Eq("1.500   |"));
    ASSERT_THAT(utils::sprintf("%08.3f", 1.5), Eq("0001.500"));
    ASSERT_THAT(utils::sprintf("%+.2f", 1.5), Eq("+1.50"));
    ASSERT_THAT(utils::sprintf("% .2f", 1.5), Eq(" 1.50"));
}

TEST(Sprintf, insert_exponential)
{
    ASSERT_THAT(utils::sprintf("%e", 1.5), Eq("1.500000e+00"));
    ASSERT_THAT(utils::sprintf("%.2E", -1.5), Eq("-1.50E+00"));
    ASSERT_THAT(utils::sprintf("%15.8e", 0.1), Eq(" 1.00000000e-01"));
}

TEST(Sprintf, insert_general)
{
    ASSERT_THAT(utils::sprintf("%g", 0.00001), Eq("1e-05"));
    ASSERT_THAT(utils::sprintf("%G", 0.00001), Eq("1E-05"));
    ASSERT_THAT(utils::sprintf("%g", 100000.0), Eq("100000"));
    ASSERT_THAT(utils::sprintf("%g", 1000000.0), Eq("1e+06"));
    // default thermo floating point format
    ASSERT_THAT(utils::sprintf("%12.8g", 0.1), Eq("         0.1"));
    ASSERT_THAT(utils::sprintf("%12.8g", 1.0 / 3.0), Eq("  0.33333333"));
    ASSERT_THAT(utils::sprintf("%-12.8g|", 1.0 / 3.0), Eq("0.33333333  |"));
}

TEST(Sprintf, insert_hexfloat)
{
    // the exact hexadecimal float representation is implementation specific,
    // so we convert the output back and compare the numbers instead
    constexpr double val = 1.5;
    auto text            = utils::sprintf("%a", val);
    ASSERT_EQ(std::strtod(text.c_str(), nullptr), val);
}

TEST(Sprintf, insert_inf_nan)
{
    constexpr double inf = std::numeric_limits<double>::infinity();
    ASSERT_THAT(utils::sprintf("%f", inf), Eq("inf"));
    ASSERT_THAT(utils::sprintf("%f", -inf), Eq("-inf"));
    ASSERT_THAT(utils::sprintf("%g", std::numeric_limits<double>::quiet_NaN()),
                StartsWith("nan"));
}

TEST(Sprintf, insert_pointer)
{
    // pointer representation is implementation specific: compare with snprintf()
    int val = 0;
    char buf[32];
    snprintf(buf, sizeof(buf), "%p", (void *) &val);
    ASSERT_THAT(utils::sprintf("%p", (void *) &val), Eq(buf));
}

TEST(Sprintf, dynamic_width_precision)
{
    ASSERT_THAT(utils::sprintf("%*d", 8, 42), Eq("      42"));
    ASSERT_THAT(utils::sprintf("%-*d|", 8, 42), Eq("42      |"));
    ASSERT_THAT(utils::sprintf("%.*f", 3, 1.5), Eq("1.500"));
    ASSERT_THAT(utils::sprintf("%*.*f", 10, 3, 1.5), Eq("     1.500"));
}

TEST(Sprintf, multiple_args)
{
    const std::string prop = "PotEng";
    auto text = utils::sprintf("Step %d: %s = %-12.8g temp: %8.3f", 100, prop, -1.5, 300.0);
    ASSERT_THAT(text, Eq("Step 100: PotEng = -1.5         temp:  300.000"));
}

TEST(Sprintf, consecutive_conversions)
{
    auto text = utils::sprintf("%d%d%s", 1, 2, "3");
    ASSERT_THAT(text, Eq("123"));
}

TEST(Sprintf, buffer_boundaries)
{
    // exercise field widths around the size of the internal fixed size buffer
    for (int width : {510, 511, 512, 513, 1024}) {
        auto text = utils::sprintf("%*d", width, 1);
        ASSERT_EQ((int) text.size(), width);
        ASSERT_EQ(text.back(), '1');
        ASSERT_EQ(text.front(), ' ');
    }
}

TEST(Sprintf, oversize_output)
{
    auto text = utils::sprintf("%600d", 42);
    ASSERT_EQ((int) text.size(), 600);
    ASSERT_THAT(text.substr(598), Eq("42"));
    ASSERT_EQ(text.front(), ' ');

    const std::string big(2000, 'x');
    text = utils::sprintf("<%s>", big);
    ASSERT_EQ((int) text.size(), 2002);
    ASSERT_THAT(text, Eq('<' + big + '>'));
}
