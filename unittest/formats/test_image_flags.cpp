/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

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


class ImageFlagsTest : public ::testing::Test {
protected:
    LAMMPS *lmp;

    void SetUp() override
    {
        const char *args[] = {"ImageFlagsTest", "-log", "none", "-echo", "screen", "-nocite"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        if (!verbose) ::testing::internal::GetCapturedStdout();
        ASSERT_NE(lmp, nullptr);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp->input->one("units real");
        lmp->input->one("dimension 3");
        lmp->input->one("region box block -2 2 -2 2 -2 2");
        lmp->input->one("create_box 1 box");
        lmp->input->one("create_atoms 1 single  0.0  0.0 0.0    units box");
        lmp->input->one("create_atoms 1 single  1.9 -1.9 1.9999 units box");
        lmp->input->one("pair_style zero 2.0");
        lmp->input->one("pair_coeff * *");
        lmp->input->one("mass * 1.0");
        lmp->input->one("set atom 1 image -1 2 3");
        lmp->input->one("set atom 2 image -2 1 -1");
        lmp->input->one("write_data test_image_flags.data");
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete lmp;
        if (!verbose) ::testing::internal::GetCapturedStdout();
        remove("test_image_flags.data");
    }
};

TEST_F(ImageFlagsTest, change_box)
{
    auto image = lmp->atom->image;
    int imx = (image[0] & IMGMASK) - IMGMAX;
    int imy = (image[0] >> IMGBITS & IMGMASK) - IMGMAX;
    int imz = (image[0] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx,-1);
    ASSERT_EQ(imy,2);
    ASSERT_EQ(imz,3);
        
    imx = (image[1] & IMGMASK) - IMGMAX;
    imy = (image[1] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[1] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx,-2);
    ASSERT_EQ(imy,1);
    ASSERT_EQ(imz,-1);
        
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("change_box all boundary f p p");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    image = lmp->atom->image;
    imx = (image[0] & IMGMASK) - IMGMAX;
    imy = (image[0] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[0] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx,0);
    ASSERT_EQ(imy,2);
    ASSERT_EQ(imz,3);

    imx = (image[1] & IMGMASK) - IMGMAX;
    imy = (image[1] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[1] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx,0);
    ASSERT_EQ(imy,1);
    ASSERT_EQ(imz,-1);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("change_box all boundary f s p");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    image = lmp->atom->image;
    imx = (image[0] & IMGMASK) - IMGMAX;
    imy = (image[0] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[0] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx,0);
    ASSERT_EQ(imy,0);
    ASSERT_EQ(imz,3);

    imx = (image[1] & IMGMASK) - IMGMAX;
    imy = (image[1] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[1] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx,0);
    ASSERT_EQ(imy,0);
    ASSERT_EQ(imz,-1);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("change_box all boundary p p m");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    image = lmp->atom->image;
    imx = (image[0] & IMGMASK) - IMGMAX;
    imy = (image[0] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[0] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx,0);
    ASSERT_EQ(imy,0);
    ASSERT_EQ(imz,0);

    imx = (image[1] & IMGMASK) - IMGMAX;
    imy = (image[1] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[1] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx,0);
    ASSERT_EQ(imy,0);
    ASSERT_EQ(imz,0);
}

TEST_F(ImageFlagsTest, read_data)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("dimension 3");
    lmp->input->one("boundary p p p");
    lmp->input->one("pair_style zero 2.0");
    lmp->input->one("read_data test_image_flags.data");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    
    auto image = lmp->atom->image;
    int imx = (image[0] & IMGMASK) - IMGMAX;
    int imy = (image[0] >> IMGBITS & IMGMASK) - IMGMAX;
    int imz = (image[0] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx,-1);
    ASSERT_EQ(imy,2);
    ASSERT_EQ(imz,3);
        
    imx = (image[1] & IMGMASK) - IMGMAX;
    imy = (image[1] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[1] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx,-2);
    ASSERT_EQ(imy,1);
    ASSERT_EQ(imz,-1);
        
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("dimension 3");
    lmp->input->one("boundary f p p");
    lmp->input->one("pair_style zero 2.0");
    lmp->input->one("read_data test_image_flags.data");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    image = lmp->atom->image;
    imx = (image[0] & IMGMASK) - IMGMAX;
    imy = (image[0] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[0] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx,0);
    ASSERT_EQ(imy,2);
    ASSERT_EQ(imz,3);

    imx = (image[1] & IMGMASK) - IMGMAX;
    imy = (image[1] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[1] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx,0);
    ASSERT_EQ(imy,1);
    ASSERT_EQ(imz,-1);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("dimension 3");
    lmp->input->one("boundary p s p");
    lmp->input->one("pair_style zero 2.0");
    lmp->input->one("read_data test_image_flags.data");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    image = lmp->atom->image;
    imx = (image[0] & IMGMASK) - IMGMAX;
    imy = (image[0] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[0] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx,-1);
    ASSERT_EQ(imy,0);
    ASSERT_EQ(imz,3);

    imx = (image[1] & IMGMASK) - IMGMAX;
    imy = (image[1] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[1] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx,-2);
    ASSERT_EQ(imy,0);
    ASSERT_EQ(imz,-1);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("dimension 3");
    lmp->input->one("boundary p p m");
    lmp->input->one("pair_style zero 2.0");
    lmp->input->one("read_data test_image_flags.data");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    image = lmp->atom->image;
    imx = (image[0] & IMGMASK) - IMGMAX;
    imy = (image[0] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[0] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx,-1);
    ASSERT_EQ(imy,2);
    ASSERT_EQ(imz,0);

    imx = (image[1] & IMGMASK) - IMGMAX;
    imy = (image[1] >> IMGBITS & IMGMASK) - IMGMAX;
    imz = (image[1] >> IMG2BITS) - IMGMAX;

    ASSERT_EQ(imx,-2);
    ASSERT_EQ(imy,1);
    ASSERT_EQ(imz,0);
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
