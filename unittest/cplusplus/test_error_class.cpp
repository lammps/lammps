// unit tests for public member functions of the Error class

#include "error.h"
#include "output.h"
#include "thermo.h"
#include "exceptions.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstring>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {
using ::testing::ContainsRegex;
using ::testing::StrEq;

class ErrorTest : public LAMMPSTest {
protected:
    Error *error;
    Thermo *thermo;

    void SetUp() override
    {
        testbinary = "ErrorClass";
        LAMMPSTest::SetUp();
        error  = lmp->error;
        thermo = lmp->output->thermo;
    }
};

TEST_F(ErrorTest, warning)
{
    // standard warning
    auto output = CAPTURE_OUTPUT([&] {
        error->warning(FLERR, "one warning");
    });
    ASSERT_THAT(output, ContainsRegex("WARNING: one warning .*test_error_class.cpp:.*"));
    ASSERT_THAT(error->get_maxwarn(), 100);

    // warnings disabled
    HIDE_OUTPUT([&] {
        command("thermo_modify warn ignore");
    });
    output = CAPTURE_OUTPUT([&] {
        error->warning(FLERR, "one warning");
    });
    ASSERT_THAT(error->get_maxwarn(), -1);

    BEGIN_HIDE_OUTPUT();
    command("thermo_modify warn 2");
    error->warning(FLERR, "one warning");
    error->warning(FLERR, "one warning");
    error->warning(FLERR, "one warning");
    END_HIDE_OUTPUT();
    ASSERT_THAT(error->get_maxwarn(), 2);
    ASSERT_THAT(error->get_numwarn(), 5);

    output = CAPTURE_OUTPUT([&] {
        thermo->lost_check();
    });
    ASSERT_THAT(output, ContainsRegex("WARNING: Too many warnings: 5 vs 2. All future.*"));

    output = CAPTURE_OUTPUT([&] {
        error->warning(FLERR, "one warning");
    });

    ASSERT_EQ(output, "");

    BEGIN_HIDE_OUTPUT();
    command("thermo_modify warn reset");
    thermo->lost_check();
    error->warning(FLERR, "one warning");
    END_HIDE_OUTPUT();
    ASSERT_THAT(error->get_maxwarn(), 2);
    ASSERT_THAT(error->get_numwarn(), 1);

    output = CAPTURE_OUTPUT([&] {
        error->warning(FLERR, "one warning");
    });
    ASSERT_THAT(output, ContainsRegex("WARNING: one warning.*"));

    BEGIN_HIDE_OUTPUT();
    command("thermo_modify warn default");
    thermo->lost_check();
    error->warning(FLERR, "one warning");
    END_HIDE_OUTPUT();
    ASSERT_THAT(error->get_maxwarn(), 100);
    ASSERT_THAT(error->get_numwarn(), 1);

    BEGIN_HIDE_OUTPUT();
    command("thermo_modify warn always");
    thermo->lost_check();
    error->warning(FLERR, "one warning");
    END_HIDE_OUTPUT();
    ASSERT_THAT(error->get_maxwarn(), 0);
    ASSERT_THAT(error->get_numwarn(), 2);
};

TEST_F(ErrorTest, one)
{
    bool caught = false;
    try {
        BEGIN_HIDE_OUTPUT();
        error->one(FLERR, "one error");
    } catch (LAMMPSAbortException &e) {
        END_HIDE_OUTPUT();
        EXPECT_THAT(e.what(), ContainsRegex("ERROR on proc 0: one error.*test_error_class.cpp"));
        caught = true;
    } catch (std::exception &e) {
        END_HIDE_OUTPUT();
        GTEST_FAIL() << "Incorrect exception\n";
    }
    ASSERT_TRUE(caught);
};

TEST_F(ErrorTest, all)
{
    bool caught = false;
    try {
        BEGIN_HIDE_OUTPUT();
        error->all(FLERR, "one error");
    } catch (LAMMPSException &e) {
        END_HIDE_OUTPUT();
        EXPECT_THAT(e.what(), ContainsRegex("ERROR: one error.*test_error_class.cpp"));
        caught = true;
    } catch (std::exception &e) {
        END_HIDE_OUTPUT();
        GTEST_FAIL() << "Incorrect exception\n";
    }
    ASSERT_TRUE(caught);
};

TEST_F(ErrorTest, errorpointer)
{
    lmp->input->line = utils::strdup("some command line with no blanks");
    lmp->input->parse();
    EXPECT_THAT(utils::point_to_error(lmp->input, Error::NOPOINTER),
                StrEq("Last input line: some command line with no blanks\n"));
    EXPECT_THAT(utils::point_to_error(lmp->input, -1),
                StrEq("Last input line: some command line with no blanks \n"
                      "                 ^^^^\n"));
    EXPECT_THAT(utils::point_to_error(lmp->input, 0),
                StrEq("Last input line: some command line with no blanks \n"
                      "                      ^^^^^^^\n"));
    EXPECT_THAT(utils::point_to_error(lmp->input, 1),
                StrEq("Last input line: some command line with no blanks \n"
                      "                              ^^^^\n"));
    EXPECT_THAT(utils::point_to_error(lmp->input, 2),
                StrEq("Last input line: some command line with no blanks \n"
                      "                                   ^^^^\n"));
    EXPECT_THAT(utils::point_to_error(lmp->input, 3),
                StrEq("Last input line: some command line with no blanks \n"
                      "                                        ^^\n"));
    EXPECT_THAT(utils::point_to_error(lmp->input, 4),
                StrEq("Last input line: some command line with no blanks \n"
                      "                                           ^^^^^^\n"));

    delete[] lmp->input->line;
    lmp->input->line = utils::strdup("some command line 'with some blanks'");
    lmp->input->parse();
    EXPECT_THAT(utils::point_to_error(nullptr, 0), StrEq(""));
    EXPECT_THAT(utils::point_to_error(lmp->input, Error::NOPOINTER),
                StrEq("Last input line: some command line 'with some blanks'\n"));
    EXPECT_THAT(utils::point_to_error(lmp->input, -1),
                StrEq("Last input line: some command line 'with some blanks'\n"
                      "--> parsed line: some command line \"with some blanks\" \n"
                      "                 ^^^^\n"));
    EXPECT_THAT(utils::point_to_error(lmp->input, 0),
                StrEq("Last input line: some command line 'with some blanks'\n"
                      "--> parsed line: some command line \"with some blanks\" \n"
                      "                      ^^^^^^^\n"));
    EXPECT_THAT(utils::point_to_error(lmp->input, 1),
                StrEq("Last input line: some command line 'with some blanks'\n"
                      "--> parsed line: some command line \"with some blanks\" \n"
                      "                              ^^^^\n"));
    EXPECT_THAT(utils::point_to_error(lmp->input, 2),
                StrEq("Last input line: some command line 'with some blanks'\n"
                      "--> parsed line: some command line \"with some blanks\" \n"
                      "                                   ^^^^^^^^^^^^^^^^^^\n"));
    delete[] lmp->input->line;
    lmp->input->line = nullptr;
};

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
