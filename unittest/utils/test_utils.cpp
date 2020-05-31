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
