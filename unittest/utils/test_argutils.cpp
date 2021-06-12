/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "arg_info.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <string>

using namespace LAMMPS_NS;
using ::testing::StrEq;

TEST(ArgInfo, plain)
{
    ArgInfo arg("text");
    ASSERT_EQ(arg.get_dim(), 0);
    ASSERT_EQ(arg.get_type(), ArgInfo::NONE);
    ASSERT_EQ(arg.get_index1(), 0);
    ASSERT_EQ(arg.get_index2(), -1);
    ASSERT_THAT(arg.get_name(), StrEq("text"));
}

TEST(ArgInfo, copy_name)
{
    char *name = nullptr;
    ArgInfo arg("text");
    ASSERT_THAT(arg.get_name(), StrEq("text"));
    name = arg.copy_name();
    ASSERT_THAT(name, StrEq("text"));
    delete[] name;
}

TEST(ArgInfo, compute0)
{
    ArgInfo arg("c_text");
    ASSERT_EQ(arg.get_dim(), 0);
    ASSERT_EQ(arg.get_type(), ArgInfo::COMPUTE);
    ASSERT_EQ(arg.get_index1(), 0);
    ASSERT_EQ(arg.get_index2(), -1);
    ASSERT_THAT(arg.get_name(), StrEq("text"));
}

TEST(ArgInfo, compute1)
{
    ArgInfo arg("c_1[5]", ArgInfo::COMPUTE);
    ASSERT_EQ(arg.get_dim(), 1);
    ASSERT_EQ(arg.get_type(), ArgInfo::COMPUTE);
    ASSERT_EQ(arg.get_index1(), 5);
    ASSERT_EQ(arg.get_index2(), -1);
    ASSERT_THAT(arg.get_name(), StrEq("1"));
}

TEST(ArgInfo, compute2)
{
    ArgInfo arg("c_text[02][05]", ArgInfo::COMPUTE | ArgInfo::FIX);
    ASSERT_EQ(arg.get_dim(), 2);
    ASSERT_EQ(arg.get_type(), ArgInfo::COMPUTE);
    ASSERT_EQ(arg.get_index1(), 2);
    ASSERT_EQ(arg.get_index2(), 5);
    ASSERT_THAT(arg.get_name(), StrEq("text"));
}

TEST(ArgInfo, fix0)
{
    ArgInfo arg("f_2");
    ASSERT_EQ(arg.get_dim(), 0);
    ASSERT_EQ(arg.get_type(), ArgInfo::FIX);
    ASSERT_EQ(arg.get_index1(), 0);
    ASSERT_EQ(arg.get_index2(), -1);
    ASSERT_THAT(arg.get_name(), StrEq("2"));
}

TEST(ArgInfo, fix1)
{
    ArgInfo arg("f_text[5]", ArgInfo::FIX | ArgInfo::VARIABLE);
    ASSERT_EQ(arg.get_dim(), 1);
    ASSERT_EQ(arg.get_type(), ArgInfo::FIX);
    ASSERT_EQ(arg.get_index1(), 5);
    ASSERT_EQ(arg.get_index2(), -1);
    ASSERT_THAT(arg.get_name(), StrEq("text"));
}

TEST(ArgInfo, fix2)
{
    ArgInfo arg("f_text[02][05]", ArgInfo::FIX);
    ASSERT_EQ(arg.get_dim(), 2);
    ASSERT_EQ(arg.get_type(), ArgInfo::FIX);
    ASSERT_EQ(arg.get_index1(), 2);
    ASSERT_EQ(arg.get_index2(), 5);
    ASSERT_THAT(arg.get_name(), StrEq("text"));
}

TEST(ArgInfo, variable0)
{
    ArgInfo arg("v_text");
    ASSERT_EQ(arg.get_dim(), 0);
    ASSERT_EQ(arg.get_type(), ArgInfo::VARIABLE);
    ASSERT_EQ(arg.get_index1(), 0);
    ASSERT_EQ(arg.get_index2(), -1);
    ASSERT_THAT(arg.get_name(), StrEq("text"));
}

TEST(ArgInfo, variable1)
{
    ArgInfo arg("v_text_1[5]", ArgInfo::VARIABLE);
    ASSERT_EQ(arg.get_dim(), 1);
    ASSERT_EQ(arg.get_type(), ArgInfo::VARIABLE);
    ASSERT_EQ(arg.get_index1(), 5);
    ASSERT_EQ(arg.get_index2(), -1);
    ASSERT_THAT(arg.get_name(), StrEq("text_1"));
}

TEST(ArgInfo, variable2)
{
    ArgInfo arg("v_x[02][05]");
    ASSERT_EQ(arg.get_dim(), 2);
    ASSERT_EQ(arg.get_type(), ArgInfo::VARIABLE);
    ASSERT_EQ(arg.get_index1(), 2);
    ASSERT_EQ(arg.get_index2(), 5);
    ASSERT_THAT(arg.get_name(), StrEq("x"));
}

TEST(ArgInfo, dname0)
{
    ArgInfo arg("d_text", ArgInfo::DNAME);
    ASSERT_EQ(arg.get_dim(), 0);
    ASSERT_EQ(arg.get_type(), ArgInfo::DNAME);
    ASSERT_EQ(arg.get_index1(), 0);
    ASSERT_EQ(arg.get_index2(), -1);
    ASSERT_THAT(arg.get_name(), StrEq("text"));
}

TEST(ArgInfo, iname0)
{
    ArgInfo arg("i_text", ArgInfo::INAME);
    ASSERT_EQ(arg.get_dim(), 0);
    ASSERT_EQ(arg.get_type(), ArgInfo::INAME);
    ASSERT_EQ(arg.get_index1(), 0);
    ASSERT_EQ(arg.get_index2(), -1);
    ASSERT_THAT(arg.get_name(), StrEq("text"));
}

TEST(ArgInfo, unsupported1)
{
    ArgInfo arg("v_text[02][05]", ArgInfo::COMPUTE | ArgInfo::FIX);
    ASSERT_EQ(arg.get_type(), ArgInfo::NONE);
}

TEST(ArgInfo, unsupported2)
{
    ArgInfo arg("d_text");
    ASSERT_EQ(arg.get_type(), ArgInfo::NONE);
}

TEST(ArgInfo, unsupported3)
{
    ArgInfo arg("i_text");
    ASSERT_EQ(arg.get_type(), ArgInfo::NONE);
}

TEST(ArgInfo, no_bracket1)
{
    ArgInfo arg("v_text[2");
    ASSERT_EQ(arg.get_type(), ArgInfo::UNKNOWN);
}

TEST(ArgInfo, no_bracket2)
{
    ArgInfo arg("v_text[2][1");
    ASSERT_EQ(arg.get_type(), ArgInfo::UNKNOWN);
}

TEST(ArgInfo, no_bracket3)
{
    ArgInfo arg("v_text[2[1]");
    ASSERT_EQ(arg.get_type(), ArgInfo::UNKNOWN);
}

TEST(ArgInfo, none)
{
    ArgInfo arg("x_text");
    ASSERT_EQ(arg.get_type(), ArgInfo::NONE);
}

TEST(ArgInfo, bad_idx1)
{
    ArgInfo arg("c_1[a]");
    ASSERT_EQ(arg.get_type(), ArgInfo::UNKNOWN);
}

TEST(ArgInfo, bad_idx2)
{
    ArgInfo arg("c_1[1][b]");
    ASSERT_EQ(arg.get_type(), ArgInfo::UNKNOWN);
}
