#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "utils.h"
#include <string>

using namespace LAMMPS_NS;
using ::testing::Eq;

TEST(Utils, trim_comment) {
    auto trimmed = utils::trim_comment("some text # comment");
    ASSERT_THAT(trimmed, Eq("some text "));
}

TEST(Utils, count_words) {
    ASSERT_EQ(utils::count_words("some text # comment"), 2);
}

TEST(Utils, valid_integer) {
    ASSERT_TRUE(utils::is_integer("10"));
}

TEST(Utils, valid_double) {
    ASSERT_TRUE(utils::is_double("10.0"));
}

TEST(Utils, empty_not_an_integer) {
    ASSERT_FALSE(utils::is_integer(""));
}

TEST(Utils, empty_not_a_double) {
    ASSERT_FALSE(utils::is_double(""));
}

TEST(Utils, double_not_an_integer) {
    ASSERT_FALSE(utils::is_integer("10.0"));
}

TEST(Utils, integer_is_double) {
    ASSERT_TRUE(utils::is_double("10"));
}

TEST(Utils, is_double_with_exponential) {
    ASSERT_TRUE(utils::is_double("10e22"));
}

TEST(Utils, is_double_with_neg_exponential) {
    ASSERT_TRUE(utils::is_double("10e-22"));
}

TEST(Utils, signed_double_and_exponential) {
    ASSERT_TRUE(utils::is_double("-10E-22"));
}
