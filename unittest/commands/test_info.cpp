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

// unit tests for the info command and the query API of the Info class

#include "lammps.h"

#include "info.h"
#include "input.h"
#include "library.h"
#include "platform.h"
#include "utils.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {
using ::testing::AnyOf;
using ::testing::ContainsRegex;
using ::testing::HasSubstr;
using ::testing::Not;

// gtest's ContainsRegex() uses different regular expression engines on different
// platforms (POSIX ERE vs. its own limited fallback on Windows) whose common
// subset is too small.  Content checks with patterns therefore use the bundled
// (and thus platform-independent) LAMMPS regex implementation instead.
#define ASSERT_MATCH(text, pattern) \
    ASSERT_TRUE(utils::strmatch(text, pattern)) << "no match for pattern: " << (pattern)

// small LJ system with one of each kind of object the info command reports on.
// only core functionality is used, so the test runs with any package selection.

class InfoTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "InfoTest";
        LAMMPSTest::SetUp();
    }

    void InitSystem() override
    {
        HIDE_OUTPUT([&] {
            command("units lj");
            command("atom_style atomic");
            command("lattice fcc 0.8442");
            command("region box block 0 2 0 2 0 2");
            command("create_box 2 box");
            command("create_atoms 1 box");
            command("mass * 1.0");
            command("set type 1 type/fraction 2 0.5 12345");
            command("velocity all create 1.0 4928459");
            command("pair_style lj/cut 2.5");
            command("pair_coeff * * 1.0 1.0");
            command("region r1 sphere 0.0 0.0 0.0 1.0");
            command("region r2 plane 0.0 0.0 0.0 0.0 0.0 1.0 side out");
            command("group g1 region r1");
            command("compute ke all ke");
            command("fix nve all nve");
            command("variable eq equal 2.0*3.0");
            command("variable str string hello");
            command("variable idx index one two");
            command("variable atm atom x*2.0");
            command("variable vec vector [1.0,2.0,3.0]");
            command("dump d1 all atom 100 info_test.dump");
            command("dump d2 all atom 10 info_test2.dump");
            command("dump_modify d2 every v_eq");
        });
    }

    void TearDown() override
    {
        LAMMPSTest::TearDown();
        platform::unlink("info_test.dump");
        platform::unlink("info_test2.dump");
    }

    std::string run_info(const std::string &args)
    {
        return CAPTURE_OUTPUT([&] {
            command("info " + args);
        });
    }
};

TEST_F(InfoTest, All)
{
    auto output = run_info("all");
    ASSERT_THAT(output, HasSubstr("Info-Info-Info"));
    ASSERT_THAT(output, HasSubstr("Printed on "));
    ASSERT_THAT(output, HasSubstr("LAMMPS version: "));
    ASSERT_THAT(output, HasSubstr("OS information: "));
    ASSERT_THAT(output, HasSubstr("sizeof(tagint):"));
    ASSERT_THAT(output, HasSubstr("Compiler: "));
    ASSERT_THAT(output, HasSubstr("Active compile time flags:"));
    ASSERT_THAT(output, HasSubstr("Installed packages:"));
    ASSERT_THAT(output, HasSubstr("Accelerator configuration:"));
    ASSERT_THAT(output, HasSubstr("Memory allocation information (MPI rank 0):"));
    ASSERT_THAT(output, HasSubstr("Communication information:"));
    ASSERT_THAT(output, HasSubstr("FFT information:"));
    ASSERT_THAT(output, HasSubstr("System information:"));
    ASSERT_THAT(output, HasSubstr("Coeff status information:"));
    ASSERT_THAT(output, HasSubstr("Group information:"));
    ASSERT_THAT(output, HasSubstr("Region information:"));
    ASSERT_THAT(output, HasSubstr("Compute information:"));
    ASSERT_THAT(output, HasSubstr("Dump information:"));
    ASSERT_THAT(output, HasSubstr("Fix information:"));
    ASSERT_THAT(output, HasSubstr("Variable information:"));
    ASSERT_THAT(output, HasSubstr("Total time information (MPI rank 0):"));
    // "all" includes the styles listing
    ASSERT_THAT(output, HasSubstr("Styles information"));
}

TEST_F(InfoTest, System)
{
    auto output = run_info("system");
    ASSERT_THAT(output, HasSubstr("System information:"));
    ASSERT_THAT(output, HasSubstr("Units         = lj"));
    ASSERT_THAT(output, HasSubstr("Atom style    = atomic"));
    ASSERT_THAT(output, HasSubstr("Atom map      = "));
    ASSERT_MATCH(output, "Atoms += +32, +types = +2, +style = lj/cut");
    ASSERT_THAT(output, HasSubstr("Atoms with atom IDs"));
    ASSERT_THAT(output, HasSubstr("Atoms with per-type masses"));
    ASSERT_THAT(output, Not(HasSubstr("Atoms with per-atom charges")));
    ASSERT_THAT(output, HasSubstr("Kspace style = none"));
    ASSERT_THAT(output, HasSubstr("Dimensions = 3"));
    ASSERT_THAT(output, HasSubstr("Orthogonal box = "));
    ASSERT_THAT(output, HasSubstr("Boundaries = p,p p,p p,p"));
    ASSERT_THAT(output, HasSubstr("xlo, xhi = 0, "));
    ASSERT_THAT(output, HasSubstr("Current timestep number = 0"));
    ASSERT_THAT(output, HasSubstr("Current timestep size = 0.005"));
    ASSERT_THAT(output, HasSubstr("Current simulation time = 0"));
    ASSERT_THAT(output, Not(HasSubstr("Group information:")));

    // triclinic box reports the tilt factors
    HIDE_OUTPUT([&] {
        command("undump d1");
        command("undump d2");
        command("change_box all triclinic");
    });
    output = run_info("system");
    ASSERT_THAT(output, HasSubstr("Triclinic box = "));
    ASSERT_THAT(output, HasSubstr("Xy, xz, yz = 0, 0, 0"));
}

TEST_F(InfoTest, Objects)
{
    auto output = run_info("groups regions computes fixes dumps variables");
    ASSERT_THAT(output, HasSubstr("Group information:"));
    ASSERT_MATCH(output, "Group\\[ 0\\]: +all +\\(static\\)");
    ASSERT_MATCH(output, "Group\\[ 1\\]: +g1 +\\(static\\)");
    ASSERT_THAT(output, HasSubstr("Region information:"));
    ASSERT_MATCH(output, "Region\\[ +[0-9]\\]: +box, +style = block, +side = in");
    ASSERT_MATCH(output, "Region\\[ +[0-9]\\]: +r1, +style = sphere, +side = in");
    ASSERT_MATCH(output, "Region\\[ +[0-9]\\]: +r2, +style = plane, +side = out");
    ASSERT_THAT(output, HasSubstr("Region[  0]:"));
    ASSERT_THAT(output, HasSubstr("Region[  1]:"));
    ASSERT_THAT(output, HasSubstr("Region[  2]:"));
    ASSERT_THAT(output, HasSubstr("   Boundary:  lo "));
    ASSERT_THAT(output, HasSubstr("   No Boundary"));
    ASSERT_THAT(output, HasSubstr("Compute information:"));
    ASSERT_MATCH(output, "Compute\\[ +[0-9]\\]: +ke, +style = ke[^,]*, +group = all");
    ASSERT_THAT(output, HasSubstr("Dump information:"));
    ASSERT_MATCH(output, "Dump\\[ +0\\]: +d1, +file = info_test.dump, +style = atom, "
                                      "+group = all, +every = 100");
    ASSERT_MATCH(output, "Dump\\[ +1\\]: +d2, +file = info_test2.dump, +style = atom, "
                                      "+group = all, +every = eq");
    ASSERT_THAT(output, HasSubstr("Fix information:"));
    ASSERT_MATCH(output, "Fix\\[ +0\\]: +nve, +style = nve[^,]*, +group = all");
    ASSERT_THAT(output, HasSubstr("Variable information:"));
    ASSERT_MATCH(output, "Variable\\[ +[0-9]+\\]: +eq, +style = equal,");
    ASSERT_MATCH(output, "Variable\\[ +[0-9]+\\]: +str, +style = string,");
    ASSERT_MATCH(output, "Variable\\[ +[0-9]+\\]: +idx, +style = index,");
    ASSERT_MATCH(output, "Variable\\[ +[0-9]+\\]: +atm, +style = atom,");
    ASSERT_MATCH(output, "Variable\\[ +[0-9]+\\]: +vec, +style = vector,");
    ASSERT_THAT(output, Not(HasSubstr("System information:")));
}

TEST_F(InfoTest, Misc)
{
    auto output = run_info("config");
    ASSERT_THAT(output, HasSubstr("LAMMPS version: "));
    ASSERT_THAT(output, HasSubstr("sizeof(bigint):   64-bit"));
    ASSERT_THAT(output, HasSubstr("C++ standard: "));
    ASSERT_THAT(output, AnyOf(HasSubstr("-DLAMMPS_SMALLBIG"), HasSubstr("-DLAMMPS_BIGBIG")));
    ASSERT_THAT(output, HasSubstr("Installed packages:"));

    output = run_info("communication");
    ASSERT_THAT(output, HasSubstr("Communication information:"));
    ASSERT_THAT(output, HasSubstr("MPI library level: "));
    ASSERT_THAT(output, HasSubstr("Comm style = brick"));
    ASSERT_THAT(output, HasSubstr("Communication mode = single"));
    ASSERT_THAT(output, HasSubstr("Communication cutoff = "));
    ASSERT_THAT(output, HasSubstr("Processor grid = "));

    // the KOKKOS package supports only "bin" neighbor lists
    if (!lmp->suffix_enable) {
        HIDE_OUTPUT([&] {
            command("neighbor 0.3 multi");
            command("comm_modify mode multi");
            command("run 0 post no");
        });
        output = run_info("comm");
        ASSERT_THAT(output, HasSubstr("Communication mode = multi"));
        ASSERT_THAT(output, HasSubstr("Communication cutoff for collection 1 = "));
    }

    output = run_info("coeffs");
    ASSERT_THAT(output, HasSubstr("Coeff status information:"));
    ASSERT_THAT(output, HasSubstr("Pair coeffs"));
    ASSERT_THAT(output, Not(HasSubstr("Bond coeffs")));

    output = run_info("time memory");
    ASSERT_THAT(output, HasSubstr("Total time information (MPI rank 0):"));
    ASSERT_THAT(output, HasSubstr("CPU time:"));
    ASSERT_THAT(output, HasSubstr("Wall time:"));
    ASSERT_THAT(output, HasSubstr("Memory allocation information (MPI rank 0):"));
    ASSERT_THAT(output, HasSubstr("Total dynamically allocated memory:"));

    output = run_info("accelerator fft");
    ASSERT_THAT(output, HasSubstr("Accelerator configuration:"));
    ASSERT_THAT(output, HasSubstr("FFT information:"));
    ASSERT_THAT(output, HasSubstr("FFT precision  = "));

    // unknown flags are reported and ignored, the rest is still printed
    output = run_info("xxx system");
    ASSERT_THAT(output, HasSubstr("WARNING: Ignoring unknown or incorrect info command flag: xxx"));
    ASSERT_THAT(output, HasSubstr("System information:"));

    // no flags at all still prints the frame
    output = run_info("");
    ASSERT_THAT(output, HasSubstr("Info-Info-Info"));
    ASSERT_THAT(output, Not(HasSubstr("information:")));
}

TEST_F(InfoTest, Styles)
{
    auto output = run_info("styles");
    ASSERT_THAT(output, HasSubstr("Styles information"));
    ASSERT_THAT(output, HasSubstr("Atom styles:"));
    ASSERT_THAT(output, HasSubstr("Integrate styles:"));
    ASSERT_THAT(output, HasSubstr("Minimize styles:"));
    ASSERT_THAT(output, HasSubstr("Pair styles:"));
    ASSERT_THAT(output, HasSubstr("Bond styles:"));
    ASSERT_THAT(output, HasSubstr("Angle styles:"));
    ASSERT_THAT(output, HasSubstr("Dihedral styles:"));
    ASSERT_THAT(output, HasSubstr("Improper styles:"));
    ASSERT_THAT(output, HasSubstr("KSpace styles:"));
    ASSERT_THAT(output, HasSubstr("Fix styles:"));
    ASSERT_THAT(output, HasSubstr("Compute styles:"));
    ASSERT_THAT(output, HasSubstr("Region styles:"));
    ASSERT_THAT(output, HasSubstr("Dump styles:"));
    ASSERT_THAT(output, HasSubstr("Command styles"));

    // "styles all" and an unknown category select everything, too
    ASSERT_THAT(run_info("styles all"), HasSubstr("Dump styles:"));
    ASSERT_THAT(run_info("styles xxx"), HasSubstr("Dump styles:"));

    output = run_info("styles atom");
    ASSERT_THAT(output, HasSubstr("Atom styles:"));
    ASSERT_MATCH(output, "[ \n]atomic[ \n]");
    ASSERT_THAT(output, Not(HasSubstr("Pair styles:")));

    output = run_info("styles pair");
    ASSERT_THAT(output, HasSubstr("Pair styles:"));
    ASSERT_MATCH(output, "[ \n]lj/cut[ \n]");
    ASSERT_THAT(output, Not(HasSubstr("Atom styles:")));

    ASSERT_MATCH(run_info("styles integrate"), "Integrate styles:[^:]*[ \n]verlet[ \n]");
    ASSERT_MATCH(run_info("styles minimize"), "Minimize styles:[^:]*[ \n]cg[ \n]");
    ASSERT_MATCH(run_info("styles bond"), "Bond styles:[^:]*[ \n]zero[ \n]");
    ASSERT_MATCH(run_info("styles angle"), "Angle styles:[^:]*[ \n]zero[ \n]");
    ASSERT_MATCH(run_info("styles dihedral"), "Dihedral styles:[^:]*[ \n]zero[ \n]");
    ASSERT_MATCH(run_info("styles improper"), "Improper styles:[^:]*[ \n]zero[ \n]");
    ASSERT_THAT(run_info("styles kspace"), HasSubstr("KSpace styles:"));
    ASSERT_MATCH(run_info("styles fix"), "Fix styles:[^:]*[ \n]nve[ \n]");
    ASSERT_MATCH(run_info("styles compute"), "Compute styles:[^:]*[ \n]temp[ \n]");
    ASSERT_MATCH(run_info("styles region"), "Region styles:[^:]*[ \n]block[ \n]");
    ASSERT_MATCH(run_info("styles dump"), "Dump styles:[^:]*[ \n]atom[ \n]");
    output = run_info("styles command");
    ASSERT_THAT(output, HasSubstr("Command styles (add-on input script commands):"));
    ASSERT_MATCH(output, "[ \n]info[ \n]");
    ASSERT_THAT(output, Not(HasSubstr("Dump styles:")));
}

TEST_F(InfoTest, Output)
{
    // without a log file, "out log" produces nothing
    auto output = run_info("out log system");
    ASSERT_THAT(output, Not(HasSubstr("Info-Info-Info")));

    output = run_info("out screen system");
    ASSERT_THAT(output, HasSubstr("System information:"));

    // write to a file, then append to it
    output = run_info("out overwrite info_test.txt system");
    ASSERT_THAT(output, Not(HasSubstr("Info-Info-Info")));
    std::ifstream in("info_test.txt");
    ASSERT_TRUE(in.good());
    std::stringstream buffer;
    buffer << in.rdbuf();
    in.close();
    auto text = buffer.str();
    ASSERT_THAT(text, HasSubstr("Info-Info-Info"));
    ASSERT_THAT(text, HasSubstr("System information:"));
    ASSERT_THAT(text, Not(HasSubstr("Group information:")));

    run_info("out append info_test.txt groups");
    in.open("info_test.txt");
    ASSERT_TRUE(in.good());
    buffer.str("");
    buffer << in.rdbuf();
    in.close();
    text = buffer.str();
    ASSERT_THAT(text, HasSubstr("System information:"));
    ASSERT_THAT(text, HasSubstr("Group information:"));
    // two complete info blocks => four frame lines
    const std::string frame = "\nInfo-Info-Info-Info-Info-Info-Info-Info-Info-Info-Info\n";
    std::size_t nframe = 0, pos = 0;
    while ((pos = text.find(frame, pos)) != std::string::npos) {
        ++nframe;
        pos += frame.size();
    }
    ASSERT_EQ(nframe, 4);

    run_info("out overwrite info_test.txt computes");
    in.open("info_test.txt");
    buffer.str("");
    buffer << in.rdbuf();
    in.close();
    text = buffer.str();
    ASSERT_THAT(text, Not(HasSubstr("System information:")));
    ASSERT_THAT(text, HasSubstr("Compute information:"));
    platform::unlink("info_test.txt");
}

TEST_F(InfoTest, QueryAPI)
{
    ASSERT_TRUE(info->is_defined("compute", "ke"));
    ASSERT_FALSE(info->is_defined("compute", "xxx"));
    ASSERT_TRUE(info->is_defined("dump", "d1"));
    ASSERT_FALSE(info->is_defined("dump", "xxx"));
    ASSERT_TRUE(info->is_defined("fix", "nve"));
    ASSERT_FALSE(info->is_defined("fix", "xxx"));
    ASSERT_TRUE(info->is_defined("group", "g1"));
    ASSERT_FALSE(info->is_defined("group", "xxx"));
    ASSERT_TRUE(info->is_defined("region", "r1"));
    ASSERT_FALSE(info->is_defined("region", "xxx"));
    ASSERT_TRUE(info->is_defined("variable", "eq"));
    ASSERT_FALSE(info->is_defined("variable", "xxx"));
    ASSERT_FALSE(info->is_defined(nullptr, "xxx"));
    ASSERT_FALSE(info->is_defined("variable", nullptr));
    TEST_FAILURE(".*ERROR: Unknown category for info is_defined\\(\\): xxx.*",
                 info->is_defined("xxx", "xxx"););

    ASSERT_TRUE(info->is_active("newton", "pair"));
    ASSERT_TRUE(info->is_active("newton", "bond"));
    ASSERT_TRUE(info->is_active("newton", "any"));
    ASSERT_FALSE(info->is_active("package", "gpu"));
    ASSERT_FALSE(info->is_active("package", "intel"));
    ASSERT_FALSE(info->is_active("package", "omp"));
    ASSERT_TRUE(info->is_active("pair", "single"));
    // the KOKKOS version of the pair styles does not support r-RESPA
    if (!lmp->suffix_enable) ASSERT_TRUE(info->is_active("pair", "respa"));
    ASSERT_FALSE(info->is_active("pair", "manybody"));
    ASSERT_FALSE(info->is_active("pair", "tail"));
    ASSERT_FALSE(info->is_active("pair", "shift"));
    ASSERT_TRUE(info->is_active("comm_style", "brick"));
    ASSERT_FALSE(info->is_active("comm_style", "tiled"));
    ASSERT_TRUE(info->is_active("min_style", "cg"));
    ASSERT_TRUE(info->is_active("run_style", "verlet"));
    ASSERT_TRUE(info->is_active("atom_style", "atomic"));
    ASSERT_FALSE(info->is_active("atom_style", "full"));
    ASSERT_TRUE(info->is_active("pair_style", "lj/cut"));
    ASSERT_TRUE(info->is_active("bond_style", "none"));
    ASSERT_TRUE(info->is_active("angle_style", "none"));
    ASSERT_TRUE(info->is_active("dihedral_style", "none"));
    ASSERT_TRUE(info->is_active("improper_style", "none"));
    ASSERT_TRUE(info->is_active("kspace_style", "none"));
    ASSERT_FALSE(info->is_active(nullptr, "none"));
    TEST_FAILURE(".*ERROR: Unknown category for info is_active\\(\\): xxx.*",
                 info->is_active("xxx", "xxx"););
    TEST_FAILURE(".*ERROR: Unknown name for info package category: xxx.*",
                 info->is_active("package", "xxx"););
    TEST_FAILURE(".*ERROR: Unknown name for info newton category: xxx.*",
                 info->is_active("newton", "xxx"););
    TEST_FAILURE(".*ERROR: Unknown name for info pair category: xxx.*",
                 info->is_active("pair", "xxx"););

    HIDE_OUTPUT([&] {
        command("pair_modify tail yes");
        command("newton off on");
        command("run 0 post no");
    });
    ASSERT_TRUE(info->is_active("pair", "tail"));
    ASSERT_FALSE(info->is_active("pair", "shift"));
    ASSERT_FALSE(info->is_active("newton", "pair"));
    ASSERT_TRUE(info->is_active("newton", "bond"));
    ASSERT_TRUE(info->is_active("newton", "any"));
    HIDE_OUTPUT([&] {
        command("pair_modify tail no shift yes");
        command("run 0 post no");
    });
    ASSERT_FALSE(info->is_active("pair", "tail"));
    ASSERT_TRUE(info->is_active("pair", "shift"));

    ASSERT_TRUE(info->is_available("command", "info"));
    ASSERT_TRUE(info->is_available("pair", "lj/cut"));
    // an unknown style name in a style category is reported as an error
    TEST_FAILURE(".*ERROR: Unknown category for info is_available\\(\\): pair.*",
                 info->is_available("pair", "xxx"););
    ASSERT_TRUE(info->is_available("feature", "exceptions"));
    ASSERT_EQ(info->is_available("feature", "gzip"), Info::has_gzip_support());
    ASSERT_EQ(info->is_available("feature", "png"), Info::has_png_support());
    ASSERT_EQ(info->is_available("feature", "jpeg"), Info::has_jpeg_support());
    ASSERT_EQ(info->is_available("feature", "ffmpeg"), Info::has_ffmpeg_support());
    ASSERT_EQ(info->is_available("feature", "curl"), Info::has_curl_support());
    ASSERT_EQ(info->is_available("feature", "fft_single"), Info::has_fft_single_support());
    ASSERT_FALSE(info->is_available("feature", "xxx"));
    ASSERT_FALSE(info->is_available(nullptr, "xxx"));
    TEST_FAILURE(".*ERROR: Unknown category for info is_available\\(\\): xxx.*",
                 info->is_available("xxx", "xxx"););

    ASSERT_TRUE(info->has_style("atom", "atomic"));
    ASSERT_FALSE(info->has_style("xxx", "atomic"));
    auto styles = info->get_available_styles("atom");
    ASSERT_NE(std::find(styles.begin(), styles.end(), "atomic"), styles.end());
    ASSERT_TRUE(info->get_available_styles("xxx").empty());

    int num    = 0;
    auto names = info->get_variable_names(num);
    ASSERT_EQ(num, (int)names.size());
    ASSERT_NE(std::find(names.begin(), names.end(), "eq"), names.end());
    ASSERT_NE(std::find(names.begin(), names.end(), "vec"), names.end());
    for (int i = 0; i < num; ++i) {
        auto text = info->get_variable_info(i);
        ASSERT_THAT(text, HasSubstr("Variable["));
        ASSERT_THAT(text, HasSubstr("style = "));
    }
    ASSERT_THAT(info->get_variable_info(num + 10), HasSubstr("(unknown)"));
    ASSERT_THAT(info->get_variable_info(-1), HasSubstr("(unknown)"));
}

// the info command on a fresh instance without a simulation box

class InfoEmptyTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "InfoEmptyTest";
        LAMMPSTest::SetUp();
    }
};

TEST_F(InfoEmptyTest, NoBox)
{
    auto output = CAPTURE_OUTPUT([&] {
        command("info system groups regions computes fixes dumps variables coeffs");
    });
    ASSERT_THAT(output, HasSubstr("System information:"));
    ASSERT_THAT(output, HasSubstr("Box has not yet been created"));
    ASSERT_THAT(output, HasSubstr("Kspace style = none"));
    ASSERT_MATCH(output, "Group\\[ 0\\]: +all +\\(static\\)");
    ASSERT_THAT(output, HasSubstr("Region information:"));
    ASSERT_THAT(output, HasSubstr("Compute information:"));
    ASSERT_THAT(output, HasSubstr("Dump information:"));
    ASSERT_THAT(output, HasSubstr("Fix information:"));
    ASSERT_THAT(output, HasSubstr("Variable information:"));
    ASSERT_THAT(output, Not(HasSubstr("Coeff status information:")));
    ASSERT_FALSE(info->is_defined("group", "g1"));
    ASSERT_TRUE(info->is_active("pair_style", "none"));
    ASSERT_FALSE(info->is_active("pair", "single"));
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

    // finalize the KOKKOS package explicitly: otherwise Kokkos is torn down by
    // static destructors at program exit, leading to segfaults in some cases
    // same workaround as the force-style and FFT3d test drivers

    lammps_kokkos_finalize();

    MPI_Finalize();
    return rv;
}
