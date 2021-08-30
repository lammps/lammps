/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "../testing/core.h"
#include "atom.h"
#include "input.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstring>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {
using ::testing::Eq;

class ImageFlagsTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "ImageFlagsTest";
        LAMMPSTest::SetUp();
        ASSERT_NE(lmp, nullptr);
        BEGIN_HIDE_OUTPUT();
        command("units real");
        command("dimension 3");
        command("region box block -2 2 -2 2 -2 2");
        command("create_box 1 box");
        command("create_atoms 1 single  0.0  0.0 0.0    units box");
        command("create_atoms 1 single  1.9 -1.9 1.9999 units box");
        command("pair_style zero 2.0");
        command("pair_coeff * *");
        command("mass * 1.0");
        command("set atom 1 image -1 2 3");
        command("set atom 2 image -2 1 -1");
        command("write_data test_image_flags.data");
        END_HIDE_OUTPUT();
    }

    void TearDown() override
    {
        LAMMPSTest::TearDown();
        remove("test_image_flags.data");
    }
};

TEST_F(ImageFlagsTest, change_box)
{
    auto image = lmp->atom->image;
    int imx    = (image[0] & IMGMASK) - IMGMAX;
    int imy    = (image[0] >> IMGBITS & IMGMASK) - IMGMAX;
    int imz    = (image[0] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx, -1);
    ASSERT_EQ(imy, 2);
    ASSERT_EQ(imz, 3);

    imx = (image[1] & IMGMASK) - IMGMAX;
    imy = (image[1] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[1] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx, -2);
    ASSERT_EQ(imy, 1);
    ASSERT_EQ(imz, -1);

    BEGIN_HIDE_OUTPUT();
    command("change_box all boundary f p p");
    END_HIDE_OUTPUT();

    image = lmp->atom->image;
    imx   = (image[0] & IMGMASK) - IMGMAX;
    imy   = (image[0] >> IMGBITS & IMGMASK) - IMGMAX;
    imz   = (image[0] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx, 0);
    ASSERT_EQ(imy, 2);
    ASSERT_EQ(imz, 3);

    imx = (image[1] & IMGMASK) - IMGMAX;
    imy = (image[1] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[1] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx, 0);
    ASSERT_EQ(imy, 1);
    ASSERT_EQ(imz, -1);

    BEGIN_HIDE_OUTPUT();
    command("change_box all boundary f s p");
    END_HIDE_OUTPUT();

    image = lmp->atom->image;
    imx   = (image[0] & IMGMASK) - IMGMAX;
    imy   = (image[0] >> IMGBITS & IMGMASK) - IMGMAX;
    imz   = (image[0] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx, 0);
    ASSERT_EQ(imy, 0);
    ASSERT_EQ(imz, 3);

    imx = (image[1] & IMGMASK) - IMGMAX;
    imy = (image[1] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[1] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx, 0);
    ASSERT_EQ(imy, 0);
    ASSERT_EQ(imz, -1);

    BEGIN_HIDE_OUTPUT();
    command("change_box all boundary p p m");
    END_HIDE_OUTPUT();

    image = lmp->atom->image;
    imx   = (image[0] & IMGMASK) - IMGMAX;
    imy   = (image[0] >> IMGBITS & IMGMASK) - IMGMAX;
    imz   = (image[0] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx, 0);
    ASSERT_EQ(imy, 0);
    ASSERT_EQ(imz, 0);

    imx = (image[1] & IMGMASK) - IMGMAX;
    imy = (image[1] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[1] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx, 0);
    ASSERT_EQ(imy, 0);
    ASSERT_EQ(imz, 0);
}

TEST_F(ImageFlagsTest, read_data)
{
    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("dimension 3");
    command("boundary p p p");
    command("pair_style zero 2.0");
    command("read_data test_image_flags.data");
    END_HIDE_OUTPUT();

    auto image = lmp->atom->image;
    int imx    = (image[0] & IMGMASK) - IMGMAX;
    int imy    = (image[0] >> IMGBITS & IMGMASK) - IMGMAX;
    int imz    = (image[0] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx, -1);
    ASSERT_EQ(imy, 2);
    ASSERT_EQ(imz, 3);

    imx = (image[1] & IMGMASK) - IMGMAX;
    imy = (image[1] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[1] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx, -2);
    ASSERT_EQ(imy, 1);
    ASSERT_EQ(imz, -1);

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("dimension 3");
    command("boundary f p p");
    command("pair_style zero 2.0");
    command("read_data test_image_flags.data");
    END_HIDE_OUTPUT();

    image = lmp->atom->image;
    imx   = (image[0] & IMGMASK) - IMGMAX;
    imy   = (image[0] >> IMGBITS & IMGMASK) - IMGMAX;
    imz   = (image[0] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx, 0);
    ASSERT_EQ(imy, 2);
    ASSERT_EQ(imz, 3);

    imx = (image[1] & IMGMASK) - IMGMAX;
    imy = (image[1] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[1] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx, 0);
    ASSERT_EQ(imy, 1);
    ASSERT_EQ(imz, -1);

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("dimension 3");
    command("boundary p s p");
    command("pair_style zero 2.0");
    command("read_data test_image_flags.data");
    END_HIDE_OUTPUT();

    image = lmp->atom->image;
    imx   = (image[0] & IMGMASK) - IMGMAX;
    imy   = (image[0] >> IMGBITS & IMGMASK) - IMGMAX;
    imz   = (image[0] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx, -1);
    ASSERT_EQ(imy, 0);
    ASSERT_EQ(imz, 3);

    imx = (image[1] & IMGMASK) - IMGMAX;
    imy = (image[1] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[1] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx, -2);
    ASSERT_EQ(imy, 0);
    ASSERT_EQ(imz, -1);

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("dimension 3");
    command("boundary p p m");
    command("pair_style zero 2.0");
    command("read_data test_image_flags.data");
    END_HIDE_OUTPUT();

    image = lmp->atom->image;
    imx   = (image[0] & IMGMASK) - IMGMAX;
    imy   = (image[0] >> IMGBITS & IMGMASK) - IMGMAX;
    imz   = (image[0] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx, -1);
    ASSERT_EQ(imy, 2);
    ASSERT_EQ(imz, 0);

    imx = (image[1] & IMGMASK) - IMGMAX;
    imy = (image[1] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[1] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx, -2);
    ASSERT_EQ(imy, 1);
    ASSERT_EQ(imz, 0);
}

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }
    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
