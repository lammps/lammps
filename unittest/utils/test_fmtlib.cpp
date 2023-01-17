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
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <exception>

using namespace LAMMPS_NS;
using ::testing::Eq;

// this tests a subset of {fmt} that is most relevant to LAMMPS

TEST(FmtLib, insert_string)
{
    const char val[] = "word";
    auto text        = fmt::format("word {}", val);
    ASSERT_THAT(text, Eq("word word"));
}

TEST(FmtLib, insert_int)
{
    const int val = 333;
    auto text     = fmt::format("word {}", val);
    ASSERT_THAT(text, Eq("word 333"));
}

TEST(FmtLib, insert_neg_int)
{
    const int val = -333;
    auto text     = fmt::format("word {}", val);
    ASSERT_THAT(text, Eq("word -333"));
}

TEST(FmtLib, insert_bigint)
{
#if defined(LAMMPS_BIGBIG) || defined(LAMMPS_SMALLBIG)
    const bigint val = 9945234592L;
    auto text        = fmt::format("word {}", val);
    ASSERT_THAT(text, Eq("word 9945234592"));
#else
    GTEST_SKIP();
#endif
}

TEST(FmtLib, insert_neg_bigint)
{
#if defined(LAMMPS_BIGBIG) || defined(LAMMPS_SMALLBIG)
    const bigint val = -9945234592L;
    auto text        = fmt::format("word {}", val);
    ASSERT_THAT(text, Eq("word -9945234592"));
#else
    GTEST_SKIP();
#endif
}

TEST(FmtLib, insert_tagint)
{
#if defined(LAMMPS_BIGBIG)
    const tagint val = 9945234592L;
    auto text        = fmt::format("word {}", val);
    ASSERT_THAT(text, Eq("word 9945234592"));
#else
    GTEST_SKIP();
#endif
}

TEST(FmtLib, insert_neg_tagint)
{
#if defined(LAMMPS_BIGBIG)
    const tagint val = -9945234592L;
    auto text        = fmt::format("word {}", val);
    ASSERT_THAT(text, Eq("word -9945234592"));
#else
    GTEST_SKIP();
#endif
}

TEST(FmtLib, insert_imageint)
{
#if defined(LAMMPS_BIGBIG)
    const imageint val = 9945234592L;
    auto text          = fmt::format("word {}", val);
    ASSERT_THAT(text, Eq("word 9945234592"));
#else
    GTEST_SKIP();
#endif
}

TEST(FmtLib, insert_neg_imageint)
{
#if defined(LAMMPS_BIGBIG)
    const imageint val = -9945234592L;
    auto text          = fmt::format("word {}", val);
    ASSERT_THAT(text, Eq("word -9945234592"));
#else
    GTEST_SKIP();
#endif
}

TEST(FmtLib, insert_double)
{
    const double val = 1.5;
    auto text        = fmt::format("word {}", val);
    ASSERT_THAT(text, Eq("word 1.5"));
}

TEST(FmtLib, insert_neg_double)
{
    const double val = -1.5;
    auto text        = fmt::format("word {}", val);
    ASSERT_THAT(text, Eq("word -1.5"));
}

TEST(FmtLib, int_for_double)
{
    const double val = -1.5;
    ASSERT_THROW(fmt::format("word {:d}", val), std::exception);
}

TEST(FmtLib, double_for_int)
{
    const int val = 15;
    ASSERT_THROW(fmt::format("word {:g}", val), std::exception);
}
