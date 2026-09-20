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

// these test the validation of user provided format strings that are
// passed to utils::sprintf() and other printf() style functions

using utils::FmtArg;

TEST(CheckFormat, accept_float)
{
    for (const auto *fmt : {"%g", "%e", "%f", "%G", "%20.15g", "%-15.8e", "%+8.3f", "%08.4f",
                            "% .6g", "%#g", "%.0f", "%La", "%lf", "[%g]", "x = %g", "%g\n"})
        ASSERT_THAT(utils::check_format(fmt, FmtArg::FLOAT), Eq("")) << "format: " << fmt;
}

TEST(CheckFormat, accept_integer)
{
    for (const auto *fmt : {"%d", "%i", "%8d", "%-8d", "%08d", "%+d", "%.5d", "%x", "%X", "%o",
                            "%u", "%ld", "%lld", "%zd", "id=%d", "%d\n"})
        ASSERT_THAT(utils::check_format(fmt, FmtArg::INTEGER), Eq("")) << "format: " << fmt;
}

TEST(CheckFormat, accept_bigint)
{
    // integer conversions are interchangeable, the width is fixed up later
    for (const auto *fmt : {"%d", "%ld", "%lld", "%20d"})
        ASSERT_THAT(utils::check_format(fmt, FmtArg::BIGINT), Eq("")) << "format: " << fmt;
}

TEST(CheckFormat, accept_string)
{
    for (const auto *fmt : {"%s", "%-8s", "%.4s", "<%s>"})
        ASSERT_THAT(utils::check_format(fmt, FmtArg::STRING), Eq("")) << "format: " << fmt;
}

TEST(CheckFormat, percent_literal_is_not_a_conversion)
{
    ASSERT_THAT(utils::check_format("%g%%", FmtArg::FLOAT), Eq(""));
    ASSERT_THAT(utils::check_format("100%% done", std::vector<FmtArg>{}), Eq(""));
    ASSERT_THAT(utils::check_format("%%%d", FmtArg::INTEGER), Eq(""));
}

TEST(CheckFormat, literal_text_is_valid)
{
    // a format string need not use all values, and need not have any
    // conversion at all: printf() ignores the surplus arguments
    ASSERT_THAT(utils::check_format("plain text", std::vector<FmtArg>{}), Eq(""));
    ASSERT_THAT(utils::check_format("", std::vector<FmtArg>{}), Eq(""));
    ASSERT_THAT(utils::check_format("plain text", FmtArg::FLOAT), Eq(""));
    ASSERT_THAT(utils::check_format("x", FmtArg::INTEGER), Eq(""));
    ASSERT_THAT(utils::check_format("", FmtArg::BIGINT), Eq(""));
    ASSERT_THAT(utils::check_format("%g", {FmtArg::FLOAT, FmtArg::FLOAT}), Eq(""));
    ASSERT_THAT(utils::check_format("just text", {FmtArg::FLOAT, FmtArg::STRING}), Eq(""));
}

TEST(CheckFormat, accept_char_conversion)
{
    // %c consumes an int and is harmless for one
    ASSERT_THAT(utils::check_format("%c", FmtArg::INTEGER), Eq(""));
    ASSERT_THAT(utils::check_format("%-3c", FmtArg::INTEGER), Eq(""));
    // but there is no valid length modifier to make it consume a bigint
    ASSERT_THAT(utils::check_format("%c", FmtArg::BIGINT),
                Eq("conversion 1 of '%c' cannot be used for large integer values"));
    ASSERT_THAT(utils::adjust_format("%c", FmtArg::INTEGER), Eq("%c"));
}

TEST(CheckFormat, reject_type_mismatch)
{
    // this is the case that silently produced garbage before
    ASSERT_THAT(utils::check_format("%d", FmtArg::FLOAT), StartsWith("conversion 1"));
    ASSERT_THAT(utils::check_format("%s", FmtArg::FLOAT), StartsWith("conversion 1"));
    ASSERT_THAT(utils::check_format("%x", FmtArg::FLOAT), StartsWith("conversion 1"));
    ASSERT_THAT(utils::check_format("%g", FmtArg::INTEGER), StartsWith("conversion 1"));
    ASSERT_THAT(utils::check_format("%f", FmtArg::BIGINT), StartsWith("conversion 1"));
    ASSERT_THAT(utils::check_format("%s", FmtArg::INTEGER), StartsWith("conversion 1"));
    ASSERT_THAT(utils::check_format("%d", FmtArg::STRING), StartsWith("conversion 1"));
}

TEST(CheckFormat, mismatch_message_is_informative)
{
    ASSERT_THAT(utils::check_format("%d", FmtArg::FLOAT),
                Eq("conversion 1 of '%d' formats integer values, "
                   "but floating-point values are provided"));
}

TEST(CheckFormat, reject_surplus_conversions)
{
    // these would consume arguments that are never passed
    ASSERT_THAT(utils::check_format("%g %g", FmtArg::FLOAT),
                Eq("'%g %g' has 2 conversion(s) but only 1 value(s) are provided"));
    // this slipped past the regular expression check for immediate variables
    ASSERT_THAT(utils::check_format("%.3f%d", FmtArg::FLOAT),
                StartsWith("'%.3f%d' has 2 conversion(s)"));
    ASSERT_THAT(utils::check_format("%g", std::vector<FmtArg>{}),
                StartsWith("'%g' has 1 conversion(s)"));
}

TEST(CheckFormat, reject_malformed)
{
    ASSERT_THAT(utils::check_format("%", FmtArg::FLOAT), StartsWith("incomplete conversion"));
    ASSERT_THAT(utils::check_format("%g %", {FmtArg::FLOAT}), StartsWith("incomplete conversion"));
    ASSERT_THAT(utils::check_format("%12.4", FmtArg::FLOAT), StartsWith("incomplete conversion"));
    ASSERT_THAT(utils::check_format("%y", FmtArg::FLOAT), StartsWith("unsupported conversion"));
    // %n writes through a pointer argument and must never be accepted
    ASSERT_THAT(utils::check_format("%n", FmtArg::FLOAT), StartsWith("unsupported conversion"));
    // '*' consumes an extra argument
    ASSERT_THAT(utils::check_format("%*g", FmtArg::FLOAT), StartsWith("variable field width"));
    ASSERT_THAT(utils::check_format("%.*g", FmtArg::FLOAT), StartsWith("variable precision"));
    // a length modifier on a string conversion changes the argument type
    ASSERT_THAT(utils::check_format("%ls", FmtArg::STRING), StartsWith("unsupported length"));
}

TEST(CheckFormat, multiple_conversions)
{
    // whole line format for dump atom without image flags
    const std::vector<FmtArg> line = {FmtArg::BIGINT, FmtArg::INTEGER, FmtArg::FLOAT,
                                      FmtArg::FLOAT, FmtArg::FLOAT};
    ASSERT_THAT(utils::check_format("%d %d %g %g %g", line), Eq(""));
    ASSERT_THAT(utils::check_format("%ld %d %20.15g %20.15g %20.15g", line), Eq(""));
    // printing only the leading values is allowed
    ASSERT_THAT(utils::check_format("%d %d %g %g", line), Eq(""));
    ASSERT_THAT(utils::check_format("%d %d %g %g %g %g", line),
                StartsWith("'%d %d %g %g %g %g' has 6 conversion(s)"));
    ASSERT_THAT(utils::check_format("%d %g %g %g %d", line), StartsWith("conversion 2"));
    // dump xyz passes the element name as a string
    ASSERT_THAT(utils::check_format("%s %g %g %g",
                                    {FmtArg::STRING, FmtArg::FLOAT, FmtArg::FLOAT, FmtArg::FLOAT}),
                Eq(""));
}

TEST(CheckFormat, conversions_are_classified_in_order)
{
    // each conversion must be matched against the value in its own position
    ASSERT_THAT(utils::check_format("%d %g %s", {FmtArg::INTEGER, FmtArg::FLOAT, FmtArg::STRING}),
                Eq(""));
    ASSERT_THAT(utils::check_format("%d %g %s", {FmtArg::INTEGER, FmtArg::STRING, FmtArg::FLOAT}),
                StartsWith("conversion 2"));
    ASSERT_THAT(utils::check_format("%s %g %d", {FmtArg::INTEGER, FmtArg::FLOAT, FmtArg::STRING}),
                StartsWith("conversion 1"));
    // the c conversion belongs to the integer group
    ASSERT_THAT(utils::check_format("id %d = %c", {FmtArg::INTEGER, FmtArg::INTEGER}), Eq(""));
    // a %% sequence is literal text and consumes no value
    ASSERT_THAT(utils::check_format("100%% done", std::vector<FmtArg>{}), Eq(""));
    ASSERT_THAT(utils::check_format("%d%% of %d", {FmtArg::INTEGER, FmtArg::INTEGER}), Eq(""));
}

TEST(AdjustFormat, bigint_length_modifier)
{
    // the resulting format must print a bigint value correctly
    constexpr bigint big = 123456789012345;
    for (const auto *fmt : {"%d", "%8d", "id=%d", "%-10d", "%ld", "%.5d", "%d\n"}) {
        auto adjusted = utils::adjust_format(fmt, FmtArg::BIGINT);
        auto text     = utils::sprintf(adjusted, big);
        ASSERT_THAT(text, ::testing::HasSubstr("123456789012345"))
            << "format: " << fmt << " adjusted: " << adjusted;
    }
}

TEST(AdjustFormat, uses_platform_length_modifier)
{
    // int64_t is "long int" on LP64 and "long long int" on LLP64, so the
    // length modifier must come from PRId64 rather than being hardcoded
    ASSERT_THAT(utils::adjust_format("%d", FmtArg::BIGINT), Eq(std::string("%") + PRId64));
    ASSERT_THAT(utils::adjust_format("%d", FmtArg::BIGINT), Eq(BIGINT_FORMAT));
    // and the result must actually round-trip a value that needs all 64 bits
    ASSERT_THAT(utils::sprintf(utils::adjust_format("%d", FmtArg::BIGINT), MAXBIGINT),
                Eq(std::to_string(MAXBIGINT)));
}

TEST(AdjustFormat, preserves_surrounding_text)
{
    // find('d') used to mangle these into "illd=%d" and "%llld"
    ASSERT_THAT(utils::adjust_format("id=%d", FmtArg::BIGINT),
                Eq(std::string("id=") + BIGINT_FORMAT));
    ASSERT_THAT(utils::adjust_format("%ld", FmtArg::BIGINT), Eq(BIGINT_FORMAT));
    ASSERT_THAT(utils::adjust_format("%d %d", {FmtArg::BIGINT, FmtArg::BIGINT}),
                Eq(std::string(BIGINT_FORMAT) + " " + BIGINT_FORMAT));
}

TEST(AdjustFormat, int_strips_length_modifier)
{
    ASSERT_THAT(utils::adjust_format("%ld", FmtArg::INTEGER), Eq("%d"));
    ASSERT_THAT(utils::adjust_format("%8lld", FmtArg::INTEGER), Eq("%8d"));
    ASSERT_THAT(utils::sprintf(utils::adjust_format("%5d", FmtArg::INTEGER), 42), Eq("   42"));
}

TEST(AdjustFormat, leaves_other_types_alone)
{
    ASSERT_THAT(utils::adjust_format("%20.15g", FmtArg::FLOAT), Eq("%20.15g"));
    ASSERT_THAT(utils::adjust_format("%-8s", FmtArg::STRING), Eq("%-8s"));
    // a format that does not match is returned unchanged
    ASSERT_THAT(utils::adjust_format("%g", FmtArg::BIGINT), Eq("%g"));
    ASSERT_THAT(utils::adjust_format("bogus %", FmtArg::BIGINT), Eq("bogus %"));
}

TEST(AdjustFormat, mixed_line_format)
{
    // dump atom with image flags: tagint, int, 3 doubles, 3 ints
    const std::vector<FmtArg> line = {FmtArg::BIGINT, FmtArg::INTEGER, FmtArg::FLOAT,
                                      FmtArg::FLOAT,  FmtArg::FLOAT,   FmtArg::INTEGER,
                                      FmtArg::INTEGER, FmtArg::INTEGER};
    auto adjusted = utils::adjust_format("%d %d %g %g %g %d %d %d", line);
    ASSERT_THAT(adjusted, Eq(std::string(BIGINT_FORMAT) + " %d %g %g %g %d %d %d"));
    // a format that uses only the leading values is adjusted the same way
    adjusted = utils::adjust_format("id=%d type=%d", line);
    ASSERT_THAT(adjusted, Eq(std::string("id=") + BIGINT_FORMAT + " type=%d"));
}
