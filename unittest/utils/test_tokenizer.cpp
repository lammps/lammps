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

#include "fmt/format.h"
#include "lmptype.h"
#include "tokenizer.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdint>

using namespace LAMMPS_NS;
using ::testing::Eq;

TEST(Tokenizer, empty_string)
{
    Tokenizer t("", " ");
    ASSERT_EQ(t.count(), 0);
}

TEST(Tokenizer, whitespace_only)
{
    Tokenizer t("    ", " ");
    ASSERT_EQ(t.count(), 0);
}

TEST(Tokenizer, single_word)
{
    Tokenizer t("test", " ");
    ASSERT_EQ(t.count(), 1);
}

TEST(Tokenizer, two_words)
{
    Tokenizer t("test word", " ");
    ASSERT_EQ(t.count(), 2);
}

TEST(Tokenizer, skip)
{
    Tokenizer t("test word", " ");
    ASSERT_TRUE(t.has_next());
    t.skip();
    ASSERT_TRUE(t.has_next());
    t.skip(1);
    ASSERT_FALSE(t.has_next());
    ASSERT_EQ(t.count(), 2);
    ASSERT_THROW(t.skip(), TokenizerException);
    try {
        t.skip();
    } catch (TokenizerException &e) {
        ASSERT_STREQ(e.what(), "No more tokens");
    }
}

TEST(Tokenizer, prefix_separators)
{
    Tokenizer t("  test word", " ");
    ASSERT_EQ(t.count(), 2);
}

TEST(Tokenizer, postfix_separators)
{
    Tokenizer t("test word   ", " ");
    ASSERT_EQ(t.count(), 2);
}

TEST(Tokenizer, iterate_words)
{
    Tokenizer t("  test word   ", " ");
    ASSERT_THAT(t.next(), Eq("test"));
    ASSERT_THAT(t.next(), Eq("word"));
    ASSERT_EQ(t.count(), 2);
}

TEST(Tokenizer, copy_constructor)
{
    Tokenizer t("  test word   ", " ");
    ASSERT_THAT(t.next(), Eq("test"));
    ASSERT_THAT(t.next(), Eq("word"));
    ASSERT_EQ(t.count(), 2);
    Tokenizer u(t);
    ASSERT_THAT(u.next(), Eq("test"));
    ASSERT_THAT(u.next(), Eq("word"));
    ASSERT_EQ(u.count(), 2);
}

TEST(Tokenizer, move_constructor)
{
    Tokenizer t("test new word   ", " ");
    Tokenizer u = std::move(t);
    ASSERT_THAT(u.next(), Eq("test"));
    ASSERT_THAT(u.next(), Eq("new"));
    ASSERT_THAT(u.next(), Eq("word"));
    ASSERT_EQ(u.count(), 3);
}

TEST(Tokenizer, copy_assignment)
{
    Tokenizer t("  test word   ", " ");
    Tokenizer u("  test2 word2 other2 ", " ");
    ASSERT_THAT(t.next(), Eq("test"));
    ASSERT_THAT(t.next(), Eq("word"));
    ASSERT_EQ(t.count(), 2);
    Tokenizer v = u;
    u           = t;
    ASSERT_THAT(u.next(), Eq("test"));
    ASSERT_THAT(u.next(), Eq("word"));
    ASSERT_EQ(u.count(), 2);

    ASSERT_THAT(v.next(), Eq("test2"));
    ASSERT_THAT(v.next(), Eq("word2"));
    ASSERT_THAT(v.next(), Eq("other2"));
    ASSERT_EQ(v.count(), 3);
}

TEST(Tokenizer, move_assignment)
{
    Tokenizer t("  test word   ", " ");
    ASSERT_THAT(t.next(), Eq("test"));
    ASSERT_THAT(t.next(), Eq("word"));
    ASSERT_EQ(t.count(), 2);
    t = Tokenizer("test new word   ", " ");
    ASSERT_THAT(t.next(), Eq("test"));
    ASSERT_THAT(t.next(), Eq("new"));
    ASSERT_THAT(t.next(), Eq("word"));
    ASSERT_EQ(t.count(), 3);
}

TEST(Tokenizer, no_separator_path)
{
    Tokenizer t("one", ":");
    ASSERT_EQ(t.has_next(), true);
    ASSERT_EQ(t.count(), 1);
    ASSERT_THAT(t.next(), Eq("one"));
    ASSERT_EQ(t.has_next(), false);
}

TEST(Tokenizer, unix_paths)
{
    Tokenizer t(":one:two:three:", ":");
    ASSERT_EQ(t.count(), 3);
    ASSERT_THAT(t.next(), Eq("one"));
    ASSERT_THAT(t.next(), Eq("two"));
    ASSERT_EQ(t.has_next(), true);
    ASSERT_THAT(t.next(), Eq("three"));
    ASSERT_EQ(t.has_next(), false);
}

TEST(Tokenizer, windows_paths)
{
    Tokenizer t("c:\\one;\\two\\three;d:four;", ";");
    ASSERT_EQ(t.count(), 3);
    ASSERT_THAT(t.next(), Eq("c:\\one"));
    ASSERT_EQ(t.has_next(), true);
    ASSERT_THAT(t.next(), Eq("\\two\\three"));
    ASSERT_THAT(t.next(), Eq("d:four"));
    ASSERT_EQ(t.has_next(), false);
}

TEST(Tokenizer, default_separators)
{
    Tokenizer t(" \r\n test \t word \f");
    ASSERT_THAT(t.next(), Eq("test"));
    ASSERT_THAT(t.next(), Eq("word"));
    ASSERT_EQ(t.count(), 2);
}

TEST(Tokenizer, contains)
{
    Tokenizer values("test word");
    ASSERT_TRUE(values.contains("test"));
    ASSERT_TRUE(values.contains("word"));
    values = Tokenizer("Triangles");
    ASSERT_TRUE(values.contains("angles"));
    ASSERT_TRUE(values.contains("Triangles"));
}

TEST(Tokenizer, not_contains)
{
    Tokenizer values("test word");
    ASSERT_FALSE(values.contains("test2"));
}

TEST(Tokenizer, matches)
{
    Tokenizer values("test word");
    ASSERT_TRUE(values.matches("test"));
    ASSERT_TRUE(values.matches("^test"));
    ASSERT_TRUE(values.matches("word"));
    ASSERT_TRUE(values.matches("word$"));
    ASSERT_TRUE(values.matches("^\\s*\\S+\\s+word"));
    values = Tokenizer("Triangles");
    ASSERT_TRUE(values.matches("^\\s*Triangles\\s*$"));
    values = Tokenizer("\t20\tatoms");
    ASSERT_TRUE(values.matches("^\\s*\\d+\\s+atoms\\s*$"));
}

TEST(Tokenizer, not_matches)
{
    Tokenizer values("test word");
    ASSERT_FALSE(values.matches("test2"));
    ASSERT_FALSE(values.matches("^word"));
    ASSERT_FALSE(values.matches("^ "));
    ASSERT_FALSE(values.matches(" $"));
    values = Tokenizer("Triangles");
    ASSERT_FALSE(values.matches("^\\s*\\S+\\s+angles"));
    values = Tokenizer("\t0x20\tatoms");
    ASSERT_FALSE(values.matches("^\\s*\\d+\\s+atoms\\s*$"));
}

TEST(Tokenizer, as_vector1)
{
    Tokenizer t(" \r\n test \t word \f");
    std::vector<std::string> list = t.as_vector();
    ASSERT_THAT(list[0], Eq("test"));
    ASSERT_THAT(list[1], Eq("word"));
}

TEST(Tokenizer, as_vector2)
{
    auto list = Tokenizer("a\\b\\c", "\\").as_vector();
    ASSERT_THAT(list[0], Eq("a"));
    ASSERT_THAT(list[1], Eq("b"));
    ASSERT_THAT(list[2], Eq("c"));
    ASSERT_EQ(list.size(), 3);
}

TEST(Tokenizer, as_vector3)
{
    auto list = Tokenizer("a\\", "\\").as_vector();
    ASSERT_THAT(list[0], Eq("a"));
    ASSERT_EQ(list.size(), 1);
}

TEST(Tokenizer, as_vector4)
{
    auto list = Tokenizer("\\a", "\\").as_vector();
    ASSERT_THAT(list[0], Eq("a"));
    ASSERT_EQ(list.size(), 1);
}

TEST(ValueTokenizer, empty_string)
{
    ValueTokenizer values("");
    ASSERT_FALSE(values.has_next());
}

TEST(ValueTokenizer, two_words)
{
    ValueTokenizer t("test word", " ");
    ASSERT_THAT(t.next_string(), Eq("test"));
    ASSERT_THAT(t.next_string(), Eq("word"));
    ASSERT_THROW(t.next_string(), TokenizerException);
}

TEST(ValueTokenizer, skip)
{
    ValueTokenizer t("test word", " ");
    ASSERT_TRUE(t.has_next());
    t.skip();
    ASSERT_TRUE(t.has_next());
    t.skip(1);
    ASSERT_FALSE(t.has_next());
    ASSERT_EQ(t.count(), 2);
    ASSERT_THROW(t.skip(), TokenizerException);
    try {
        t.skip();
    } catch (TokenizerException &e) {
        ASSERT_STREQ(e.what(), "No more tokens");
    }
}

TEST(ValueTokenizer, copy_constructor)
{
    ValueTokenizer t("  test word   ", " ");
    ASSERT_THAT(t.next_string(), Eq("test"));
    ASSERT_THAT(t.next_string(), Eq("word"));
    ASSERT_EQ(t.count(), 2);
    ValueTokenizer u(t);
    ASSERT_THAT(u.next_string(), Eq("test"));
    ASSERT_THAT(u.next_string(), Eq("word"));
    ASSERT_EQ(u.count(), 2);
}

TEST(ValueTokenizer, move_constructor)
{
    ValueTokenizer t("  test new word   ", " ");
    ValueTokenizer u = std::move(t);
    ASSERT_THAT(u.next_string(), Eq("test"));
    ASSERT_THAT(u.next_string(), Eq("new"));
    ASSERT_THAT(u.next_string(), Eq("word"));
    ASSERT_EQ(u.count(), 3);
}

TEST(ValueTokenizer, copy_assignment)
{
    ValueTokenizer t("  test word   ", " ");
    ValueTokenizer u("  test2 word2 other2 ", " ");
    ASSERT_THAT(t.next_string(), Eq("test"));
    ASSERT_THAT(t.next_string(), Eq("word"));
    ASSERT_EQ(t.count(), 2);
    ValueTokenizer v = u;
    u                = t;
    ASSERT_THAT(u.next_string(), Eq("test"));
    ASSERT_THAT(u.next_string(), Eq("word"));
    ASSERT_EQ(u.count(), 2);

    ASSERT_THAT(v.next_string(), Eq("test2"));
    ASSERT_THAT(v.next_string(), Eq("word2"));
    ASSERT_THAT(v.next_string(), Eq("other2"));
    ASSERT_EQ(v.count(), 3);
}

TEST(ValueTokenizer, move_assignment)
{
    ValueTokenizer t("  test word   ", " ");
    ASSERT_THAT(t.next_string(), Eq("test"));
    ASSERT_THAT(t.next_string(), Eq("word"));
    ASSERT_EQ(t.count(), 2);
    t = ValueTokenizer("test new word   ", " ");
    ASSERT_THAT(t.next_string(), Eq("test"));
    ASSERT_THAT(t.next_string(), Eq("new"));
    ASSERT_THAT(t.next_string(), Eq("word"));
    ASSERT_EQ(t.count(), 3);
}

TEST(ValueTokenizer, bad_integer)
{
    ValueTokenizer values("f10 f11 f12 0xff 109951162777 "
                          "36893488147419103232 36893488147419103232");
    ASSERT_THROW(values.next_int(), InvalidIntegerException);
    ASSERT_THROW(values.next_bigint(), InvalidIntegerException);
    ASSERT_THROW(values.next_tagint(), InvalidIntegerException);
    ASSERT_THROW(values.next_int(), InvalidIntegerException);
    ASSERT_THROW(values.next_int(), InvalidIntegerException);
    ASSERT_THROW(values.next_tagint(), InvalidIntegerException);
    ASSERT_THROW(values.next_bigint(), InvalidIntegerException);
}

TEST(ValueTokenizer, bad_double)
{
    ValueTokenizer values("1a.0 --2.0 2.4d3 -1e20000 1.0e-1.0");
    ASSERT_THROW(values.next_double(), InvalidFloatException);
    ASSERT_THROW(values.next_double(), InvalidFloatException);
    ASSERT_THROW(values.next_double(), InvalidFloatException);
    ASSERT_THROW(values.next_double(), InvalidFloatException);
    ASSERT_THROW(values.next_double(), InvalidFloatException);
}

TEST(ValueTokenizer, valid_int)
{
    ValueTokenizer values(fmt::format("10 {} {}", -MAXSMALLINT - 1, MAXSMALLINT));
    ASSERT_EQ(values.next_int(), 10);
    ASSERT_EQ(values.next_int(), -MAXSMALLINT - 1);
    ASSERT_EQ(values.next_int(), MAXSMALLINT);
}

TEST(ValueTokenizer, valid_tagint)
{
    ValueTokenizer values(
        fmt::format("42 {} {} {} {}", -MAXSMALLINT - 1, MAXSMALLINT, -MAXTAGINT - 1, MAXTAGINT));
    ASSERT_EQ(values.next_tagint(), 42);
    ASSERT_EQ(values.next_tagint(), -MAXSMALLINT - 1);
    ASSERT_EQ(values.next_tagint(), MAXSMALLINT);
    ASSERT_EQ(values.next_tagint(), -MAXTAGINT - 1);
    ASSERT_EQ(values.next_tagint(), MAXTAGINT);
}

TEST(ValueTokenizer, valid_bigint)
{
    ValueTokenizer values(
        fmt::format("42 {} {} {} {}", -MAXSMALLINT - 1, MAXSMALLINT, -MAXBIGINT - 1, MAXBIGINT));
    ASSERT_EQ(values.next_bigint(), 42);
    ASSERT_EQ(values.next_bigint(), -MAXSMALLINT - 1);
    ASSERT_EQ(values.next_bigint(), MAXSMALLINT);
    ASSERT_EQ(values.next_bigint(), -MAXBIGINT - 1);
    ASSERT_EQ(values.next_bigint(), MAXBIGINT);
}

TEST(ValueTokenizer, valid_double)
{
    ValueTokenizer values("3.14 -0.00002 .1 0xff " + std::to_string(MAXBIGINT));
    ASSERT_DOUBLE_EQ(values.next_double(), 3.14);
    ASSERT_DOUBLE_EQ(values.next_double(), -0.00002);
    ASSERT_DOUBLE_EQ(values.next_double(), 0.1);
    ASSERT_DOUBLE_EQ(values.next_double(), 255);
    ASSERT_DOUBLE_EQ(values.next_double(), MAXBIGINT);
}

TEST(ValueTokenizer, valid_double_with_exponential)
{
    ValueTokenizer values(fmt::format("3.14e22 {} {}", DBL_MAX, DBL_MIN));
    ASSERT_DOUBLE_EQ(values.next_double(), 3.14e22);
    ASSERT_DOUBLE_EQ(values.next_double(), DBL_MAX);
    ASSERT_DOUBLE_EQ(values.next_double(), DBL_MIN);
}

TEST(ValueTokenizer, contains)
{
    ValueTokenizer values("test word");
    ASSERT_TRUE(values.contains("test"));
    ASSERT_TRUE(values.contains("word"));
    values = ValueTokenizer("Triangles");
    ASSERT_TRUE(values.contains("angles"));
    ASSERT_TRUE(values.contains("Triangles"));
}

TEST(ValueTokenizer, not_contains)
{
    ValueTokenizer values("test word");
    ASSERT_FALSE(values.contains("test2"));
}

TEST(ValueTokenizer, matches)
{
    ValueTokenizer values("test word");
    ASSERT_TRUE(values.matches("test"));
    ASSERT_TRUE(values.matches("^test"));
    ASSERT_TRUE(values.matches("word"));
    ASSERT_TRUE(values.matches("word$"));
    ASSERT_TRUE(values.matches("^\\s*\\S+\\s+word"));
    values = ValueTokenizer("Triangles");
    ASSERT_TRUE(values.matches("^\\s*Triangles\\s*$"));
    values = ValueTokenizer("\t20\tatoms");
    ASSERT_TRUE(values.matches("^\\s*\\d+\\s+atoms\\s*$"));
}

TEST(ValueTokenizer, not_matches)
{
    ValueTokenizer values("test word");
    ASSERT_FALSE(values.matches("test2"));
    ASSERT_FALSE(values.matches("^word"));
    ASSERT_FALSE(values.matches("^ "));
    ASSERT_FALSE(values.matches(" $"));
    values = ValueTokenizer("Triangles");
    ASSERT_FALSE(values.matches("^\\s*\\S+\\s+angles"));
    values = ValueTokenizer("\t0x20\tatoms");
    ASSERT_FALSE(values.matches("^\\s*\\d+\\s+atoms\\s*$"));
}

TEST(ValueTokenizer, missing_int)
{
    ValueTokenizer values("10");
    ASSERT_EQ(values.next_int(), 10);
    ASSERT_THROW(values.next_int(), TokenizerException);
}

TEST(ValueTokenizer, missing_tagint)
{
    ValueTokenizer values("42");
    ASSERT_EQ(values.next_tagint(), 42);
    ASSERT_THROW(values.next_tagint(), TokenizerException);
}

TEST(ValueTokenizer, missing_bigint)
{
    ValueTokenizer values("42");
    ASSERT_EQ(values.next_bigint(), 42);
    ASSERT_THROW(values.next_bigint(), TokenizerException);
}

TEST(ValueTokenizer, missing_double)
{
    ValueTokenizer values("3.14");
    ASSERT_DOUBLE_EQ(values.next_double(), 3.14);
    ASSERT_THROW(values.next_double(), TokenizerException);
}
