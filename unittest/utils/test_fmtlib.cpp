#include "lmptype.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "fmt/format.h"
#include <string>

using namespace LAMMPS_NS;
using ::testing::Eq;

// this tests a subset of {fmt} that is most relevant to LAMMPS

TEST(Fmtlib, insert_string) {
    const char word[] = "word";
    auto text = fmt::format("word {}",word);
    ASSERT_THAT(text, Eq("word word"));
}

TEST(Fmtlib, insert_int) {
    const int word = 333;
    auto text = fmt::format("word {}",word);
    ASSERT_THAT(text, Eq("word 333"));
}

TEST(Fmtlib, insert_neg_int) {
    const int word = -333;
    auto text = fmt::format("word {}",word);
    ASSERT_THAT(text, Eq("word -333"));
}

TEST(Fmtlib, insert_double) {
    const double word = 1.5;
    auto text = fmt::format("word {}",word);
    ASSERT_THAT(text, Eq("word 1.5"));
}

TEST(Fmtlib, insert_neg_double) {
    const double word = -1.5;
    auto text = fmt::format("word {}",word);
    ASSERT_THAT(text, Eq("word 1.5"));
}
