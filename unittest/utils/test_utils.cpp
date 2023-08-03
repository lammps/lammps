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
#include "pointers.h"
#include "tokenizer.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cerrno>
#include <string>
#include <vector>

using namespace LAMMPS_NS;
using ::testing::EndsWith;
using ::testing::Eq;
using ::testing::StrEq;

#if !defined(FLERR)
#define FLERR __FILE__, __LINE__
#endif

TEST(Utils, strdup)
{
    std::string original("some_text");
    const char *copy = utils::strdup(original);
    ASSERT_THAT(original, StrEq(copy));
    ASSERT_NE(copy, original.c_str());

    const char *copy2 = utils::strdup(copy);
    ASSERT_THAT(original, StrEq(copy2));
    ASSERT_NE(copy, copy2);

    delete[] copy;
    delete[] copy2;
}

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

TEST(Utils, casemod)
{
    ASSERT_THAT(utils::lowercase("Gba35%*zAKgRvr"), StrEq("gba35%*zakgrvr"));
    ASSERT_THAT(utils::lowercase("A BC DEFG"), StrEq("a bc defg"));
    ASSERT_THAT(utils::uppercase("Gba35%*zAKgRvr"), StrEq("GBA35%*ZAKGRVR"));
    ASSERT_THAT(utils::uppercase("a bc defg"), StrEq("A BC DEFG"));
}

TEST(Utils, trim_comment)
{
    auto trimmed = utils::trim_comment("some text # comment");
    ASSERT_THAT(trimmed, StrEq("some text "));
}

TEST(Utils, star_subst)
{
    std::string starred = "beforeafter";
    std::string subst   = utils::star_subst(starred, 1234, 0);
    ASSERT_THAT(subst, StrEq("beforeafter"));

    starred = "before*after";
    subst   = utils::star_subst(starred, 1234, 6);
    ASSERT_THAT(subst, StrEq("before001234after"));

    starred = "before*";
    subst   = utils::star_subst(starred, 1234, 0);
    ASSERT_THAT(subst, StrEq("before1234"));

    starred = "*after";
    subst   = utils::star_subst(starred, 1234, 2);
    ASSERT_THAT(subst, StrEq("1234after"));
}

TEST(Utils, has_utf8)
{
    const char ascii_string[] = " -2";
    const char utf8_string[]  = " −2";
    ASSERT_FALSE(utils::has_utf8(ascii_string));
    ASSERT_TRUE(utils::has_utf8(utf8_string));
}

TEST(Utils, utf8_subst)
{
    const char ascii_string[] = " -2";
    const char utf8_string[]  = " −2";
    auto ascii                = utils::utf8_subst(ascii_string);
    auto utf8                 = utils::utf8_subst(utf8_string);
    ASSERT_TRUE(ascii == utf8);
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

TEST(Utils, join_words)
{
    std::vector<std::string> words = {"one", "two", "three"};
    auto combined                  = utils::join_words(words, " ");
    ASSERT_THAT(combined, StrEq("one two three"));
    combined = utils::join_words(words, "");
    ASSERT_THAT(combined, StrEq("onetwothree"));
    words[1] = "two ";
    combined = utils::join_words(words, "__");
    ASSERT_THAT(combined, StrEq("one__two __three"));
    words.resize(1);
    combined = utils::join_words(words, "/");
    ASSERT_THAT(combined, StrEq("one"));
    words.emplace_back("");
    combined = utils::join_words(words, "1");
    ASSERT_THAT(combined, StrEq("one1"));
}

TEST(Utils, split_words_simple)
{
    auto list = utils::split_words("one two three");
    ASSERT_EQ(list.size(), 3);
    ASSERT_THAT(list[0], StrEq("one"));
    ASSERT_THAT(list[1], StrEq("two"));
    ASSERT_THAT(list[2], StrEq("three"));
}

TEST(Utils, split_words_leading_whitespace)
{
    auto list = utils::split_words("  one two three");
    ASSERT_EQ(list.size(), 3);
    ASSERT_THAT(list[0], StrEq("one"));
    ASSERT_THAT(list[1], StrEq("two"));
    ASSERT_THAT(list[2], StrEq("three"));
}

TEST(Utils, split_words_trailing_whitespace)
{
    auto list = utils::split_words("one two three  ");
    ASSERT_EQ(list.size(), 3);
    ASSERT_THAT(list[0], StrEq("one"));
    ASSERT_THAT(list[1], StrEq("two"));
    ASSERT_THAT(list[2], StrEq("three"));
}

TEST(Utils, split_words_heredoc)
{
    auto list = utils::split_words("one two three \"\"\"");
    ASSERT_EQ(list.size(), 4);
    ASSERT_THAT(list[0], StrEq("one"));
    ASSERT_THAT(list[1], StrEq("two"));
    ASSERT_THAT(list[2], StrEq("three"));
    ASSERT_THAT(list[3], StrEq("\"\"\""));
}

TEST(Utils, split_words_heredoc_whitespace)
{
    auto list = utils::split_words("one two three \"\"\"   ");
    ASSERT_EQ(list.size(), 4);
    ASSERT_THAT(list[0], StrEq("one"));
    ASSERT_THAT(list[1], StrEq("two"));
    ASSERT_THAT(list[2], StrEq("three"));
    ASSERT_THAT(list[3], StrEq("\"\"\""));
}

TEST(Utils, split_words_quoted)
{
    auto list = utils::split_words("one 'two' \"three\"");
    ASSERT_EQ(list.size(), 3);
    ASSERT_THAT(list[0], StrEq("one"));
    ASSERT_THAT(list[1], StrEq("two"));
    ASSERT_THAT(list[2], StrEq("three"));
}

TEST(Utils, split_words_escaped)
{
    auto list = utils::split_words("1\\' '\"two\"' 3\\\"");
    ASSERT_EQ(list.size(), 3);
    ASSERT_THAT(list[0], StrEq("1\\'"));
    ASSERT_THAT(list[1], StrEq("\"two\""));
    ASSERT_THAT(list[2], StrEq("3\\\""));
}

TEST(Utils, split_words_quote_in_quoted)
{
    auto list = utils::split_words("one 't\\'wo' \"th\\\"ree\"");
    ASSERT_EQ(list.size(), 3);
    ASSERT_THAT(list[0], StrEq("one"));
    ASSERT_THAT(list[1], StrEq("t\\'wo"));
    ASSERT_THAT(list[2], StrEq("th\\\"ree"));
}

TEST(Utils, split_lines)
{
    auto list = utils::split_lines(" line 1\nline 2 \n line 3 \n");
    ASSERT_EQ(list.size(), 3);
    ASSERT_THAT(list[0], StrEq(" line 1"));
    ASSERT_THAT(list[1], StrEq("line 2 "));
    ASSERT_THAT(list[2], StrEq(" line 3 "));
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

TEST(Utils, valid_double10)
{
    ASSERT_TRUE(utils::is_double("-0.15"));
}

TEST(Utils, valid_double11)
{
    ASSERT_TRUE(utils::is_double("-27.5"));
}

TEST(Utils, valid_double12)
{
    ASSERT_TRUE(utils::is_double("+0.15"));
}

TEST(Utils, valid_double13)
{
    ASSERT_TRUE(utils::is_double("+27.5"));
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

TEST(Utils, minus_not_an_integer1)
{
    ASSERT_FALSE(utils::is_integer("1-"));
}

TEST(Utils, plus_not_an_integer1)
{
    ASSERT_FALSE(utils::is_integer("1+"));
}

TEST(Utils, minus_not_an_integer2)
{
    ASSERT_FALSE(utils::is_integer("--1"));
}

TEST(Utils, plus_not_an_integer2)
{
    ASSERT_FALSE(utils::is_integer("++1"));
}

TEST(Utils, plusminus_not_an_integer1)
{
    ASSERT_FALSE(utils::is_integer("-+1"));
}

TEST(Utils, plusminus_not_an_integer2)
{
    ASSERT_FALSE(utils::is_integer("+-1"));
}

TEST(Utils, minus_not_a_double1)
{
    ASSERT_FALSE(utils::is_double("1-"));
}

TEST(Utils, plus_not_a_double1)
{
    ASSERT_FALSE(utils::is_double("1+"));
}

TEST(Utils, minus_not_a_double2)
{
    ASSERT_FALSE(utils::is_double("--1"));
}

TEST(Utils, plus_not_a_double2)
{
    ASSERT_FALSE(utils::is_double("++1"));
}

TEST(Utils, plusminus_not_a_double1)
{
    ASSERT_FALSE(utils::is_double("+-1"));
}

TEST(Utils, plusminus_not_a_double2)
{
    ASSERT_FALSE(utils::is_double("-+1"));
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

TEST(Utils, signed_double_and_broken_exponential)
{
    ASSERT_FALSE(utils::is_double("-10e10-2"));
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

TEST(Utils, valid_id1)
{
    ASSERT_TRUE(utils::is_id("abc"));
}

TEST(Utils, valid_id2)
{
    ASSERT_TRUE(utils::is_id("123"));
}

TEST(Utils, valid_id3)
{
    ASSERT_TRUE(utils::is_id("abc123"));
}

TEST(Utils, valid_id4)
{
    ASSERT_TRUE(utils::is_id("abc_123"));
}

TEST(Utils, valid_id5)
{
    ASSERT_TRUE(utils::is_id("123_abc"));
}

TEST(Utils, valid_id6)
{
    ASSERT_TRUE(utils::is_id("_123"));
}

TEST(Utils, valid_id7)
{
    ASSERT_TRUE(utils::is_id("___"));
}

TEST(Utils, empty_id)
{
    ASSERT_FALSE(utils::is_id(""));
}

TEST(Utils, invalid_id1)
{
    ASSERT_FALSE(utils::is_id("+abc"));
}

TEST(Utils, invalid_id2)
{
    ASSERT_FALSE(utils::is_id("a[1]"));
}

TEST(Utils, invalid_id3)
{
    ASSERT_FALSE(utils::is_id("b(c)"));
}

TEST(Utils, invalid_id4)
{
    ASSERT_FALSE(utils::is_id("a$12"));
}

TEST(Utils, valid_numeric)
{
    ASSERT_EQ(utils::is_type("1"), 0);
    ASSERT_EQ(utils::is_type("21"), 0);
    ASSERT_EQ(utils::is_type("05"), 0);
    ASSERT_EQ(utils::is_type("1*"), 0);
    ASSERT_EQ(utils::is_type("*2"), 0);
    ASSERT_EQ(utils::is_type("1*4"), 0);
}

TEST(Utils, invalid_numeric)
{
    ASSERT_EQ(utils::is_type("1*2*"), -1);
    ASSERT_EQ(utils::is_type("**2"), -1);
    ASSERT_EQ(utils::is_type("*4*"), -1);
    ASSERT_EQ(utils::is_type("30**"), -1);
}

TEST(Utils, valid_label)
{
    ASSERT_EQ(utils::is_type("A"), 1);
    ASSERT_EQ(utils::is_type("c1"), 1);
    ASSERT_EQ(utils::is_type("o1_"), 1);
    ASSERT_EQ(utils::is_type("C1'"), 1);
    ASSERT_EQ(utils::is_type("N2\"-C1'"), 1);
    ASSERT_EQ(utils::is_type("[N2\"][C1']"), 1);
    ASSERT_EQ(utils::is_type("@X2=&X1"), 1);
    ASSERT_EQ(utils::is_type("|Na|Cl|H2O|"), 1);
    ASSERT_EQ(utils::is_type("CA(1)/CB(1)"), 1);
    ASSERT_EQ(utils::is_type("A-B"), 1); // ASCII
    ASSERT_EQ(utils::is_type("A−B"), 1); // UTF-8
}

TEST(Utils, invalid_label)
{
    ASSERT_EQ(utils::is_type("1A"), -1);
    ASSERT_EQ(utils::is_type("#c"), -1);
    ASSERT_EQ(utils::is_type("*B"), -1);
    ASSERT_EQ(utils::is_type(" B"), -1);
    ASSERT_EQ(utils::is_type("A "), -1);
    ASSERT_EQ(utils::is_type("A B"), -1);
    ASSERT_EQ(utils::is_type("\tB"), -1);
    ASSERT_EQ(utils::is_type("C\n"), -1);
    ASSERT_EQ(utils::is_type("d\r"), -1);
    ASSERT_EQ(utils::is_type(""), -1);
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
    ASSERT_TRUE(utils::strmatch("c_name", "^[cfvid]2?_name"));
    ASSERT_TRUE(utils::strmatch("f_name", "^[cfvid]2?_name"));
    ASSERT_TRUE(utils::strmatch("v_name", "^[cfvid]2?_name"));
    ASSERT_TRUE(utils::strmatch("i_name", "^[cfvid]2?_name"));
    ASSERT_TRUE(utils::strmatch("d_name", "^[cfvid]2?_name"));
    ASSERT_TRUE(utils::strmatch("i2_name", "^[cfvid]2?_name"));
    ASSERT_TRUE(utils::strmatch("d2_name", "^[cfvid]2?_name"));
    ASSERT_FALSE(utils::strmatch("d2name", "^[cfvid]2?_name"));
    ASSERT_FALSE(utils::strmatch("i1_name", "^[cfvid]2?_name"));
    ASSERT_FALSE(utils::strmatch("V_name", "^[cfvid]2?_name"));
    ASSERT_FALSE(utils::strmatch("x_name", "^[cfvid]2?_name"));
}

TEST(Utils, strmatch_yaml_suffix)
{
    ASSERT_TRUE(utils::strmatch("test.yaml", "\\.[yY][aA]?[mM][lL]$"));
    ASSERT_TRUE(utils::strmatch("test.yml", "\\.[yY][aA]?[mM][lL]$"));
    ASSERT_TRUE(utils::strmatch("TEST.YAML", "\\.[yY][aA]?[mM][lL]$"));
    ASSERT_TRUE(utils::strmatch("TEST.YML", "\\.[yY][aA]?[mM][lL]$"));
    ASSERT_FALSE(utils::strmatch("test.yamlx", "\\.[yY][aA]?[mM][lL]$"));
    ASSERT_FALSE(utils::strmatch("test.ymlx", "\\.[yY][aA]?[mM][lL]$"));
    ASSERT_FALSE(utils::strmatch("TEST.YAMLX", "\\.[yY][aA]?[mM][lL]$"));
    ASSERT_FALSE(utils::strmatch("TEST.YMLX", "\\.[yY][aA]?[mM][lL]$"));
    ASSERT_FALSE(utils::strmatch("testyaml", "\\.[yY][aA]?[mM][lL]$"));
    ASSERT_FALSE(utils::strmatch("testyml", "\\.[yY][aA]?[mM][lL]$"));
    ASSERT_FALSE(utils::strmatch("TESTYAML", "\\.[yY][aA]?[mM][lL]$"));
    ASSERT_FALSE(utils::strmatch("TESTYML", "\\.[yY][aA]?[mM][lL]$"));
    ASSERT_FALSE(utils::strmatch("yaml.test", "\\.[yY][aA]?[mM][lL]$"));
    ASSERT_FALSE(utils::strmatch("yml.test", "\\.[yY][aA]?[mM][lL]$"));
    ASSERT_FALSE(utils::strmatch("YAML.TEST", "\\.[yY][aA]?[mM][lL]$"));
    ASSERT_FALSE(utils::strmatch("YML.TEST", "\\.[yY][aA]?[mM][lL]$"));
    ASSERT_FALSE(utils::strmatch("test", "\\.[yY][aA]?[mM][lL]$"));
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

TEST(Utils, strmatch_range)
{
    ASSERT_TRUE(utils::strmatch("*11", "^\\d*\\*\\d*$"));
    ASSERT_TRUE(utils::strmatch("2*11", "^\\d*\\*\\d*$"));
    ASSERT_TRUE(utils::strmatch("5*", "^\\d*\\*\\d*$"));
    ASSERT_TRUE(utils::strmatch("*", "^\\d*\\*\\d*$"));
    ASSERT_FALSE(utils::strmatch("x5*", "^\\d*\\*\\d*$"));
    ASSERT_FALSE(utils::strmatch("x*", "^\\d*\\*\\d*$"));
    ASSERT_FALSE(utils::strmatch("*a", "^\\d*\\*\\d*$"));
    ASSERT_FALSE(utils::strmatch("1*2d", "^\\d*\\*\\d*$"));
}

TEST(Utils, strfind_beg)
{
    ASSERT_THAT(utils::strfind("rigid/small/omp", "^rigid"), StrEq("rigid"));
}

TEST(Utils, strfind_mid1)
{
    ASSERT_THAT(utils::strfind("rigid/small/omp", ".small."), StrEq("/small/"));
}

TEST(Utils, strfind_mid2)
{
    ASSERT_THAT(utils::strfind("rigid/small/ompXXX", "omp"), StrEq("omp"));
}

TEST(Utils, strfind_end)
{
    ASSERT_THAT(utils::strfind("rigid/small/omp", "/omp$"), StrEq("/omp"));
}

TEST(Utils, no_strfind_beg)
{
    ASSERT_THAT(utils::strfind("rigid/small/omp", "^small"), StrEq(""));
}

TEST(Utils, no_strfind_mid)
{
    ASSERT_THAT(utils::strfind("rigid/small/omp", "none"), StrEq(""));
}

TEST(Utils, no_strfind_end)
{
    ASSERT_THAT(utils::strfind("rigid/small/omp", "/opt$"), StrEq(""));
}

TEST(Utils, strfind_whole_line)
{
    ASSERT_THAT(utils::strfind("ITEM: UNITS\n", "^\\s*ITEM: UNITS\\s*$"), StrEq("ITEM: UNITS\n"));
}

TEST(Utils, no_strfind_whole_line)
{
    ASSERT_THAT(utils::strfind("ITEM: UNITS\n", "^\\s*ITEM: UNIT\\s*$"), StrEq(""));
}

TEST(Utils, strfind_char_range)
{
    ASSERT_THAT(utils::strfind("rigidXXX", "^[ip-s]+gid"), StrEq("rigid"));
}

TEST(Utils, strfind_notchar_range)
{
    ASSERT_THAT(utils::strfind("rigidYYY", "^[^a-g]+gid"), StrEq("rigid"));
}

TEST(Utils, strfind_backslash)
{
    ASSERT_THAT(utils::strfind("\\rigidZZZ", "^\\W\\w+gid"), StrEq("\\rigid"));
}

TEST(Utils, strfind_opt_range)
{
    ASSERT_THAT(utils::strfind("rigidAAA", "^[0-9]*[\\Wp-s]igid"), StrEq("rigid"));
}

TEST(Utils, strfind_opt_char)
{
    ASSERT_THAT(utils::strfind("rigid111", "^r?igid"), StrEq("rigid"));
    ASSERT_THAT(utils::strfind("igid222", "^r?igid"), StrEq("igid"));
}

TEST(Utils, strfind_dot)
{
    ASSERT_THAT(utils::strfind("AAArigidBBB", ".igid"), StrEq("rigid"));
    ASSERT_THAT(utils::strfind("000Rigid111", ".igid"), StrEq("Rigid"));
}

TEST(Utils, strfind_kim)
{
    ASSERT_THAT(
        utils::strfind("n3409jfse MO_004835508849_000 aslfjiaf", "[MS][MO]_\\d\\d\\d+_\\d\\d\\d"),
        StrEq("MO_004835508849_000"));
    ASSERT_THAT(utils::strfind("VanDuinChakraborty_2003_CHNO__SM_107643900657_000",
                               "[MS][MO]_\\d\\d\\d+_\\d\\d\\d"),
                StrEq("SM_107643900657_000"));
}

TEST(Utils, bounds_case1)
{
    int nlo, nhi;

    nlo = nhi = -1;
    utils::bounds(FLERR, "9", 0, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 9);
    ASSERT_EQ(nhi, 9);
    utils::bounds(FLERR, "1", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 1);
    ASSERT_EQ(nhi, 1);
    utils::bounds(FLERR, "1x", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
    utils::bounds(FLERR, "-1", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
    utils::bounds(FLERR, "+1", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
    utils::bounds(FLERR, "1:3", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
}

TEST(Utils, bounds_case2)
{
    int nlo, nhi;

    nlo = nhi = -1;
    utils::bounds(FLERR, "*", 0, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 0);
    ASSERT_EQ(nhi, 10);
    utils::bounds(FLERR, "*", -10, 5, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -10);
    ASSERT_EQ(nhi, 5);
    utils::bounds(FLERR, "?", -10, 5, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
}

TEST(Utils, bounds_case3)
{
    int nlo, nhi;

    nlo = nhi = -1;
    utils::bounds(FLERR, "2*", 0, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 2);
    ASSERT_EQ(nhi, 10);
    utils::bounds(FLERR, "3*", -10, 5, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 3);
    ASSERT_EQ(nhi, 5);
    utils::bounds(FLERR, "3*:2", -10, 5, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
}

TEST(Utils, boundsbig_case1)
{
    bigint nlo, nhi;

    nlo = nhi = -1;
    utils::bounds(FLERR, "9", 0, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 9);
    ASSERT_EQ(nhi, 9);
    utils::bounds(FLERR, "1", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 1);
    ASSERT_EQ(nhi, 1);
    utils::bounds(FLERR, "1x", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
    utils::bounds(FLERR, "-1", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
    utils::bounds(FLERR, "+1", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
    utils::bounds(FLERR, "1:3", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
}

TEST(Utils, boundsbig_case2)
{
    bigint nlo, nhi;

    nlo = nhi = -1;
    utils::bounds(FLERR, "*", 0, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 0);
    ASSERT_EQ(nhi, 10);
    utils::bounds(FLERR, "*", -10, 5, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -10);
    ASSERT_EQ(nhi, 5);
    utils::bounds(FLERR, "?", -10, 5, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
}

TEST(Utils, boundsbig_case3)
{
    bigint nlo, nhi;

    nlo = nhi = -1;
    utils::bounds(FLERR, "2*", 0, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 2);
    ASSERT_EQ(nhi, 10);
    utils::bounds(FLERR, "3*", -10, 5, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 3);
    ASSERT_EQ(nhi, 5);
    utils::bounds(FLERR, "3*:2", -10, 5, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
}

TEST(Utils, parse_grid_id)
{
    auto words = utils::parse_grid_id(FLERR, "c_1:full:density", nullptr);
    ASSERT_THAT(words[0], StrEq("c_1"));
    ASSERT_THAT(words[1], StrEq("full"));
    ASSERT_THAT(words[2], StrEq("density"));

    words = utils::parse_grid_id(FLERR, "c_1:full:density[1]", nullptr);
    ASSERT_THAT(words[0], StrEq("c_1"));
    ASSERT_THAT(words[1], StrEq("full"));
    ASSERT_THAT(words[2], StrEq("density[1]"));

    words = utils::parse_grid_id(FLERR, "c_1:full:density[*]", nullptr);
    ASSERT_THAT(words[0], StrEq("c_1"));
    ASSERT_THAT(words[1], StrEq("full"));
    ASSERT_THAT(words[2], StrEq("density[*]"));

    words = utils::parse_grid_id(FLERR, "c_1_full_density", nullptr);
    ASSERT_THAT(words[0], StrEq(""));
    ASSERT_THAT(words[1], StrEq(""));
    ASSERT_THAT(words[0], StrEq(""));

    words = utils::parse_grid_id(FLERR, "c_1:full:", nullptr);
    ASSERT_THAT(words[0], StrEq(""));
    ASSERT_THAT(words[1], StrEq(""));
    ASSERT_THAT(words[0], StrEq(""));

    words = utils::parse_grid_id(FLERR, ":full:density", nullptr);
    ASSERT_THAT(words[0], StrEq(""));
    ASSERT_THAT(words[1], StrEq(""));
    ASSERT_THAT(words[0], StrEq(""));

    words = utils::parse_grid_id(FLERR, "c_1:full", nullptr);
    ASSERT_THAT(words[0], StrEq(""));
    ASSERT_THAT(words[1], StrEq(""));
    ASSERT_THAT(words[0], StrEq(""));
}

TEST(Utils, errorurl)
{
    auto errmesg = utils::errorurl(10);
    ASSERT_THAT(errmesg, Eq("\nFor more information see https://docs.lammps.org/err0010"));
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

    ASSERT_TRUE(platform::file_is_readable("ctest1.txt"));
    ASSERT_TRUE(platform::file_is_readable("ctest2.txt"));
    ASSERT_FALSE(platform::file_is_readable("no_such_file.txt"));

    ASSERT_THAT(utils::get_potential_file_path("ctest1.txt"), Eq("ctest1.txt"));
    ASSERT_THAT(utils::get_potential_file_path("no_such_file.txt"), Eq(""));

    const char *folder = getenv("LAMMPS_POTENTIALS");
    if (folder != nullptr) {
        std::string path = platform::path_join(folder, "Cu_u3.eam");
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

TEST(Utils, timespec2seconds_off)
{
    ASSERT_DOUBLE_EQ(utils::timespec2seconds("off"), -1.0);
}

TEST(Utils, timespec2seconds_ss)
{
    ASSERT_DOUBLE_EQ(utils::timespec2seconds("45"), 45.0);
}

TEST(Utils, timespec2seconds_mmss)
{
    ASSERT_DOUBLE_EQ(utils::timespec2seconds("10:45"), 645.0);
}

TEST(Utils, timespec2seconds_hhmmss)
{
    ASSERT_DOUBLE_EQ(utils::timespec2seconds("2:10:45"), 7845.0);
}

TEST(Utils, timespec2seconds_invalid)
{
    ASSERT_DOUBLE_EQ(utils::timespec2seconds("2:aa:45"), -1.0);
}

TEST(Utils, date2num)
{
    ASSERT_EQ(utils::date2num("1Jan05"), 20050101);
    ASSERT_EQ(utils::date2num("10Feb2005"), 20050210);
    ASSERT_EQ(utils::date2num("02Mar10"), 20100302);
    ASSERT_EQ(utils::date2num(" 5Apr1900"), 19000405);
    ASSERT_EQ(utils::date2num("10May22 "), 20220510);
    ASSERT_EQ(utils::date2num("1 Jun 05"), 20050601);
    ASSERT_EQ(utils::date2num("10 Jul 2005"), 20050710);
    ASSERT_EQ(utils::date2num("02 Aug 10"), 20100802);
    ASSERT_EQ(utils::date2num("  5  September  99"), 20990905);
    ASSERT_EQ(utils::date2num("10October22 "), 20221010);
    ASSERT_EQ(utils::date2num("30November 02"), 20021130);
    ASSERT_EQ(utils::date2num("31December100"), 1001231);
}

TEST(Utils, current_date)
{
    auto vals = ValueTokenizer(utils::current_date(), "-");
    int year  = vals.next_int();
    int month = vals.next_int();
    int day   = vals.next_int();
    ASSERT_GT(year, 2020);
    ASSERT_GE(month, 1);
    ASSERT_GE(day, 1);
    ASSERT_LE(month, 12);
    ASSERT_LE(day, 31);
}

TEST(Utils, binary_search)
{
    double data[] = {-2.0, -1.8, -1.0, -1.0, -1.0, -0.5, -0.2, 0.0, 0.1, 0.1,
                     0.2,  0.3,  0.5,  0.5,  0.6,  0.7,  1.0,  1.2, 1.5, 2.0};
    const int n   = sizeof(data) / sizeof(double);
    ASSERT_EQ(utils::binary_search(-5.0, n, data), 0);
    ASSERT_EQ(utils::binary_search(-2.0, n, data), 0);
    ASSERT_EQ(utils::binary_search(-1.9, n, data), 0);
    ASSERT_EQ(utils::binary_search(-1.0, n, data), 4);
    ASSERT_EQ(utils::binary_search(0.0, n, data), 7);
    ASSERT_EQ(utils::binary_search(0.1, n, data), 9);
    ASSERT_EQ(utils::binary_search(0.4, n, data), 11);
    ASSERT_EQ(utils::binary_search(1.1, n, data), 16);
    ASSERT_EQ(utils::binary_search(1.5, n, data), 18);
    ASSERT_EQ(utils::binary_search(2.0, n, data), 19);
    ASSERT_EQ(utils::binary_search(2.5, n, data), 19);
}

static int compare(int a, int b, void *)
{
    if (a < b)
        return -1;
    else if (a > b)
        return 1;
    else
        return 0;
}

TEST(Utils, merge_sort)
{
    int data[] = {773, 405, 490, 830, 632, 96,  428, 728, 912, 840, 878, 745, 213, 219, 249, 380,
                  894, 758, 575, 690, 61,  849, 19,  577, 338, 569, 898, 873, 448, 940, 431, 780,
                  472, 289, 65,  491, 641, 37,  367, 33,  407, 854, 594, 611, 845, 136, 107, 592,
                  275, 865, 158, 626, 399, 703, 686, 734, 188, 559, 781, 558, 737, 281, 638, 664,
                  533, 529, 62,  969, 595, 661, 837, 463, 624, 568, 615, 936, 206, 637, 91,  694,
                  214, 872, 468, 66,  775, 949, 486, 576, 255, 961, 480, 138, 177, 509, 333, 705,
                  10,  375, 321, 952, 210, 111, 475, 268, 708, 864, 244, 121, 988, 540, 942, 682,
                  750, 473, 478, 714, 955, 911, 482, 384, 144, 757, 697, 791, 420, 605, 447, 320};

    const int num = sizeof(data) / sizeof(int);
    utils::merge_sort(data, num, nullptr, &compare);
    bool sorted = true;
    for (int i = 1; i < num; ++i)
        if (data[i - 1] > data[i]) sorted = false;
    ASSERT_TRUE(sorted);
}
