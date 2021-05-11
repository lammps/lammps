// unit tests for issuing command to a LAMMPS instance through the Input class

#include "error.h"
#include "info.h"
#include "lammps.h"
#include "output.h"
#include "thermo.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <string>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

using ::testing::MatchesRegex;
using utils::split_words;

class Error_class : public LAMMPSTest {
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

TEST_F(Error_class, message)
{
    auto output = CAPTURE_OUTPUT([&] {
        error->message(FLERR, "one message");
    });
    ASSERT_THAT(output, MatchesRegex("one message .*test_error_class.cpp:.*"));
};

TEST_F(Error_class, warning)
{
    // standard warning
    auto output = CAPTURE_OUTPUT([&] {
        error->warning(FLERR, "one warning");
    });
    ASSERT_THAT(output, MatchesRegex("WARNING: one warning .*test_error_class.cpp:.*"));
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
    ASSERT_THAT(output, MatchesRegex("WARNING: Too many warnings: 5 vs 2. All future.*"));

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
    ASSERT_THAT(output, MatchesRegex("WARNING: one warning.*"));

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

TEST_F(Error_class, one)
{
    TEST_FAILURE("ERROR on proc 0: one error.*test_error_class.cpp:.*",
                 error->one(FLERR, "one error"););
};

TEST_F(Error_class, all)
{
    TEST_FAILURE("ERROR: one error.*test_error_class.cpp:.*", error->all(FLERR, "one error"););
};

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (Info::get_mpi_vendor() == "Open MPI" && !LAMMPS_NS::Info::has_exceptions())
        std::cout << "Warning: using OpenMPI without exceptions. "
                     "Death tests will be skipped\n";

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
