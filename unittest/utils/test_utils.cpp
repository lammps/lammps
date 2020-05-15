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
