#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "tokenizer.h"

using namespace LAMMPS_NS;
using ::testing::Eq;

TEST(Tokenizer, empty_string) {
    Tokenizer t("", " ");
    ASSERT_EQ(t.count(), 0);
}

TEST(Tokenizer, whitespace_only) {
    Tokenizer t("    ", " ");
    ASSERT_EQ(t.count(), 0);
}

TEST(Tokenizer, single_word) {
    Tokenizer t("test", " ");
    ASSERT_EQ(t.count(), 1);
}

TEST(Tokenizer, two_words) {
    Tokenizer t("test word", " ");
    ASSERT_EQ(t.count(), 2);
}

TEST(Tokenizer, prefix_seperators) {
    Tokenizer t("  test word", " ");
    ASSERT_EQ(t.count(), 2);
}

TEST(Tokenizer, postfix_seperators) {
    Tokenizer t("test word   ", " ");
    ASSERT_EQ(t.count(), 2);
}

TEST(Tokenizer, iterate_words) {
    Tokenizer t("  test word   ", " ");
    ASSERT_THAT(t[0], Eq("test"));
    ASSERT_THAT(t[1], Eq("word"));
    ASSERT_EQ(t.count(), 2);
}

TEST(Tokenizer, default_seperators) {
    Tokenizer t(" \r\n test \t word \f");
    ASSERT_THAT(t[0], Eq("test"));
    ASSERT_THAT(t[1], Eq("word"));
    ASSERT_EQ(t.count(), 2);
}

TEST(Tokenizer, for_loop) {
    Tokenizer t(" \r\n test \t word \f");
    std::vector<std::string> list;

    for(auto word : t) {
        list.push_back(word);
    }
    ASSERT_THAT(list[0], Eq("test"));
    ASSERT_THAT(list[1], Eq("word"));
}