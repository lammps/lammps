// unit tests for checking and changing simulation properties through the library interface

#include "library.h"
#include "lammps.h"
#include <string>

#include "gtest/gtest.h"

const char *demo_input[] = {
      "region       box block 0 $x 0 2 0 2",
      "create_box 1 box",
      "create_atoms 1 single 1.0 1.0 ${zpos}" };
const char *cont_input[] = {
    "create_atoms 1 single &",
    "0.2 0.1 0.1" };

class LAMMPS_properties : public ::testing::Test
{
protected:
    void *lmp;
    LAMMPS_properties() {};
    ~LAMMPS_properties() override {};

    void SetUp() override {
        const char *args[] = {"LAMMPS_test", "-log", "none",
                              "-echo", "screen", "-nocite" };
        char **argv = (char **)args;
        int argc = sizeof(args)/sizeof(char *);

        ::testing::internal::CaptureStdout();
        lmp = lammps_open_no_mpi(argc, argv, NULL);
        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_STREQ(output.substr(0,8).c_str(), "LAMMPS (");
    }
    void TearDown() override {
        ::testing::internal::CaptureStdout();
        lammps_close(lmp);
        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_STREQ(output.substr(0,16).c_str(), "Total wall time:");
        lmp = nullptr;
    }
};

TEST_F(LAMMPS_properties, box) {
};
