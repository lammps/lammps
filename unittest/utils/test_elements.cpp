
#include "elements.h"
#include "utils.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {
using ::testing::Eq;
using ::testing::StrEq;

class ElementsTest : public LAMMPSTest {
protected:
    Error *error;

    void SetUp() override
    {
        testbinary = "ElementsTest";
        LAMMPSTest::SetUp();
        error = lmp->error;
    }

};


/* ------------------------------------------------------------------ */

TEST_F(ElementsTest, symbol)
{
  ASSERT_THAT(LAMMPS_NS::elements::symbol(0,error), StrEq("X"));
  ASSERT_THAT(LAMMPS_NS::elements::symbol(1,error), StrEq("H"));
  ASSERT_THAT(LAMMPS_NS::elements::symbol(74,error), StrEq("W"));
  ASSERT_THAT(LAMMPS_NS::elements::symbol(96,error), StrEq("Cm"));

  TEST_FAILURE("atomic_number -1 out of range \\(0-96\\)",
    LAMMPS_NS::elements::symbol(-1,error););

  TEST_FAILURE("atomic_number 99 out of range \\(0-96\\)",
    elements::symbol(99,error););
}

TEST_F(ElementsTest, name)
{
  ASSERT_THAT(LAMMPS_NS::elements::name(0,error), StrEq("X"));
  ASSERT_THAT(LAMMPS_NS::elements::name(1,error), StrEq("Hydrogen"));
  ASSERT_THAT(LAMMPS_NS::elements::name(81,error), StrEq("Thallium"));
  ASSERT_THAT(LAMMPS_NS::elements::name(96,error), StrEq("Curium"));

  TEST_FAILURE("atomic_number -1 out of range \\(0-96\\)",
    LAMMPS_NS::elements::name(-1,error););

  TEST_FAILURE("atomic_number 99 out of range \\(0-96\\)",
    elements::name(99,error););
}

TEST_F(ElementsTest, cpkHexColor)
{
  ASSERT_THAT(LAMMPS_NS::elements::cpkHexColor(0,error), StrEq("000000"));
  ASSERT_THAT(LAMMPS_NS::elements::cpkHexColor(1,error), StrEq("FFFFFF"));
  ASSERT_THAT(LAMMPS_NS::elements::cpkHexColor(6,error), StrEq("909090"));
  ASSERT_THAT(LAMMPS_NS::elements::cpkHexColor(96,error), StrEq("785CE3"));

  TEST_FAILURE("atomic_number -1 out of range \\(0-96\\)",
    LAMMPS_NS::elements::cpkHexColor(-1,error););

  TEST_FAILURE("atomic_number 99 out of range \\(0-96\\)",
    elements::cpkHexColor(99,error););
}

/* ------------------------------------------------------------------ */

TEST_F(ElementsTest, atomic_mass)
{
  ASSERT_EQ(LAMMPS_NS::elements::atomic_mass(0,error), 0.0);
  ASSERT_EQ(LAMMPS_NS::elements::atomic_mass(1,error), 1.008);
  ASSERT_EQ(LAMMPS_NS::elements::atomic_mass(8,error), 15.999);
  ASSERT_EQ(LAMMPS_NS::elements::atomic_mass(96,error), 247.07035);

  TEST_FAILURE("atomic_number -1 out of range \\(0-96\\)",
    LAMMPS_NS::elements::atomic_mass(-1,error););

  TEST_FAILURE("atomic_number 99 out of range \\(0-96\\)",
    elements::atomic_mass(99,error););
}

TEST_F(ElementsTest, vdw_radius)
{
  ASSERT_EQ(LAMMPS_NS::elements::vdw_radius(0,error), 0.0);
  ASSERT_EQ(LAMMPS_NS::elements::vdw_radius(1,error), 120.0);
  ASSERT_EQ(LAMMPS_NS::elements::vdw_radius(11,error), 227.0);
  ASSERT_EQ(LAMMPS_NS::elements::vdw_radius(96,error), 245.0);

  TEST_FAILURE("atomic_number -1 out of range \\(0-96\\)",
    LAMMPS_NS::elements::vdw_radius(-1,error););

  TEST_FAILURE("atomic_number 99 out of range \\(0-96\\)",
    elements::vdw_radius(99,error););
}

TEST_F(ElementsTest, covalent_radius)
{
  ASSERT_EQ(LAMMPS_NS::elements::covalent_radius(0,error), 0.0);
  ASSERT_EQ(LAMMPS_NS::elements::covalent_radius(1,error), 31.0);
  ASSERT_EQ(LAMMPS_NS::elements::covalent_radius(7,error), 71.0);
  ASSERT_EQ(LAMMPS_NS::elements::covalent_radius(96,error), 169.0);

  TEST_FAILURE("atomic_number -1 out of range \\(0-96\\)",
    LAMMPS_NS::elements::covalent_radius(-1,error););

  TEST_FAILURE("atomic_number 99 out of range \\(0-96\\)",
    elements::covalent_radius(99,error););
}

/* ------------------------------------------------------------------ */

TEST_F(ElementsTest, atomic_number_with_symbol)
{
  ASSERT_EQ(LAMMPS_NS::elements::atomic_number_with_symbol("X",error), 0);
  ASSERT_EQ(LAMMPS_NS::elements::atomic_number_with_symbol("H",error), 1);
  ASSERT_EQ(LAMMPS_NS::elements::atomic_number_with_symbol("Ca",error), 20);

  TEST_FAILURE("symbol FOO not found",
    elements::atomic_number_with_symbol("FOO",error););
}

TEST_F(ElementsTest, atomic_number_with_closest_mass)
{
  ASSERT_EQ(LAMMPS_NS::elements::atomic_number_with_closest_mass(1.008,error), 1);
  ASSERT_EQ(LAMMPS_NS::elements::atomic_number_with_closest_mass(15.5,error), 8);
  ASSERT_EQ(LAMMPS_NS::elements::atomic_number_with_closest_mass(196.96657,error), 79);

  TEST_FAILURE("atomic mass -1.1 is negative, must be >=0",
    elements::atomic_number_with_closest_mass(-1.1,error););

  TEST_FAILURE("atomic mass 999.9 higher than heaviest element Curium \\(247.07035\\) available", elements::atomic_number_with_closest_mass(999.9,error););
}
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
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
