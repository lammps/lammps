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

#include "MANYBODY/pair_comb.h"
#include "MANYBODY/pair_comb3.h"
#include "MANYBODY/pair_gw.h"
#include "MANYBODY/pair_gw_zbl.h"
#include "MANYBODY/pair_nb3b_harmonic.h"
#include "MANYBODY/pair_sw.h"
#include "MANYBODY/pair_tersoff.h"
#include "MANYBODY/pair_tersoff_mod.h"
#include "MANYBODY/pair_tersoff_mod_c.h"
#include "MANYBODY/pair_tersoff_table.h"
#include "MANYBODY/pair_tersoff_zbl.h"
#include "MANYBODY/pair_vashishta.h"
#include "info.h"
#include "input.h"
#include "potential_file_reader.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "../testing/core.h"

#include <cstring>
#include <iostream>
#include <mpi.h>
#include <vector>

using namespace LAMMPS_NS;
using utils::split_words;

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

#if __cplusplus < 201703L
constexpr int LAMMPS_NS::PairSW::NPARAMS_PER_LINE;
constexpr int LAMMPS_NS::PairComb::NPARAMS_PER_LINE;
constexpr int LAMMPS_NS::PairComb3::NPARAMS_PER_LINE;
constexpr int LAMMPS_NS::PairTersoff::NPARAMS_PER_LINE;
constexpr int LAMMPS_NS::PairTersoffMOD::NPARAMS_PER_LINE;
constexpr int LAMMPS_NS::PairTersoffMODC::NPARAMS_PER_LINE;
constexpr int LAMMPS_NS::PairTersoffZBL::NPARAMS_PER_LINE;
constexpr int LAMMPS_NS::PairGW::NPARAMS_PER_LINE;
constexpr int LAMMPS_NS::PairGWZBL::NPARAMS_PER_LINE;
constexpr int LAMMPS_NS::PairNb3bHarmonic::NPARAMS_PER_LINE;
constexpr int LAMMPS_NS::PairVashishta::NPARAMS_PER_LINE;
constexpr int LAMMPS_NS::PairTersoffTable::NPARAMS_PER_LINE;
#endif

class PotentialFileReaderTest : public LAMMPSTest {};

// open for native units
TEST_F(PotentialFileReaderTest, Sw_native)
{
    BEGIN_HIDE_OUTPUT();
    command("units metal");
    PotentialFileReader reader(lmp, "Si.sw", "Stillinger-Weber");
    END_HIDE_OUTPUT();

    auto *line = reader.next_line(PairSW::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairSW::NPARAMS_PER_LINE);
}

// open with supported conversion enabled
TEST_F(PotentialFileReaderTest, Sw_conv)
{
    BEGIN_HIDE_OUTPUT();
    command("units real");
    PotentialFileReader reader(lmp, "Si.sw", "Stillinger-Weber", utils::METAL2REAL);
    END_HIDE_OUTPUT();

    auto *line = reader.next_line(PairSW::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairSW::NPARAMS_PER_LINE);
}

// open without conversion enabled
TEST_F(PotentialFileReaderTest, Sw_noconv)
{
    BEGIN_HIDE_OUTPUT();
    command("units real");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*ERROR on proc.*potential.*requires metal units but real.*",
                 PotentialFileReader reader(lmp, "Si.sw", "Stillinger-Weber", utils::REAL2METAL););
}

TEST_F(PotentialFileReaderTest, Comb)
{
    BEGIN_HIDE_OUTPUT();
    command("units metal");
    PotentialFileReader reader(lmp, "ffield.comb", "COMB");
    END_HIDE_OUTPUT();

    auto *line = reader.next_line(PairComb::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairComb::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, Comb3)
{
    BEGIN_HIDE_OUTPUT();
    command("units metal");
    PotentialFileReader reader(lmp, "ffield.comb3", "COMB3");
    END_HIDE_OUTPUT();

    auto *line = reader.next_line(PairComb3::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairComb3::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, Tersoff)
{
    BEGIN_HIDE_OUTPUT();
    command("units metal");
    PotentialFileReader reader(lmp, "Si.tersoff", "Tersoff");
    END_HIDE_OUTPUT();

    auto *line = reader.next_line(PairTersoff::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairTersoff::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, TersoffMod)
{
    BEGIN_HIDE_OUTPUT();
    command("units metal");
    PotentialFileReader reader(lmp, "Si.tersoff.mod", "Tersoff/Mod");
    END_HIDE_OUTPUT();

    auto *line = reader.next_line(PairTersoffMOD::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairTersoffMOD::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, TersoffModC)
{
    BEGIN_HIDE_OUTPUT();
    command("units metal");
    PotentialFileReader reader(lmp, "Si.tersoff.modc", "Tersoff/ModC");
    END_HIDE_OUTPUT();

    auto *line = reader.next_line(PairTersoffMODC::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairTersoffMODC::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, TersoffTable)
{
    BEGIN_HIDE_OUTPUT();
    command("units metal");
    PotentialFileReader reader(lmp, "Si.tersoff", "TersoffTable");
    END_HIDE_OUTPUT();

    auto *line = reader.next_line(PairTersoffTable::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairTersoffTable::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, TersoffZBL)
{
    BEGIN_HIDE_OUTPUT();
    command("units metal");
    PotentialFileReader reader(lmp, "SiC.tersoff.zbl", "Tersoff/ZBL");
    END_HIDE_OUTPUT();

    auto *line = reader.next_line(PairTersoffZBL::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairTersoffZBL::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, GW)
{
    BEGIN_HIDE_OUTPUT();
    command("units metal");
    PotentialFileReader reader(lmp, "SiC.gw", "GW");
    END_HIDE_OUTPUT();

    auto *line = reader.next_line(PairGW::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairGW::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, GWZBL)
{
    BEGIN_HIDE_OUTPUT();
    command("units metal");
    PotentialFileReader reader(lmp, "SiC.gw.zbl", "GW/ZBL");
    END_HIDE_OUTPUT();

    auto *line = reader.next_line(PairGWZBL::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairGWZBL::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, Nb3bHarmonic)
{
    BEGIN_HIDE_OUTPUT();
    command("units real");
    PotentialFileReader reader(lmp, "MOH.nb3b.harmonic", "NB3B Harmonic");
    END_HIDE_OUTPUT();

    auto *line = reader.next_line(PairNb3bHarmonic::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairNb3bHarmonic::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, Vashishta)
{
    BEGIN_HIDE_OUTPUT();
    command("units metal");
    PotentialFileReader reader(lmp, "SiC.vashishta", "Vashishta");
    END_HIDE_OUTPUT();

    auto *line = reader.next_line(PairVashishta::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairVashishta::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, UnitConvert)
{
    PotentialFileReader *reader;
    int unit_convert, flag;

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    reader = new PotentialFileReader(lmp, "Si.sw", "Stillinger-Weber");
    END_HIDE_OUTPUT();

    unit_convert = reader->get_unit_convert();
    ASSERT_EQ(unit_convert, 0);
    delete reader;

    BEGIN_HIDE_OUTPUT();
    flag   = utils::get_supported_conversions(utils::UNKNOWN);
    reader = new PotentialFileReader(lmp, "Si.sw", "Stillinger-Weber", flag);
    END_HIDE_OUTPUT();

    unit_convert = reader->get_unit_convert();
    ASSERT_EQ(unit_convert, 0);
    delete reader;

    BEGIN_HIDE_OUTPUT();
    flag   = utils::get_supported_conversions(utils::ENERGY);
    reader = new PotentialFileReader(lmp, "Si.sw", "Stillinger-Weber", flag);
    END_HIDE_OUTPUT();

    unit_convert = reader->get_unit_convert();
    ASSERT_EQ(unit_convert, 0);
    delete reader;

    BEGIN_HIDE_OUTPUT();
    flag = utils::get_supported_conversions(utils::ENERGY);
    command("units real");
    reader = new PotentialFileReader(lmp, "Si.sw", "Stillinger-Weber", flag);
    END_HIDE_OUTPUT();

    unit_convert = reader->get_unit_convert();
    ASSERT_EQ(unit_convert, utils::METAL2REAL);
    delete reader;
}

TEST_F(PotentialFileReaderTest, convenience_functions)
{
    FILE *fp = fopen("potential_reader.file", "w");
    fmt::print(fp, "123\n");
    fmt::print(fp, "-123.45\n");
    fmt::print(fp, "{}\n", DBL_MIN);
    fmt::print(fp, "{}\n", DBL_MAX);
    fmt::print(fp, "{}\n", -DBL_MIN);
    fmt::print(fp, "{}\n", -DBL_MAX);
    fmt::print(fp, "{}\n", MAXSMALLINT);
    fmt::print(fp, "{}\n", MAXTAGINT);
    fmt::print(fp, "{}\n", MAXBIGINT);
    fmt::print(fp, "fooBAR\n");
    fclose(fp);

    BEGIN_HIDE_OUTPUT();
    command("units real");
    PotentialFileReader reader(lmp, "potential_reader.file", "test");
    END_HIDE_OUTPUT();

    utils::logmesg(lmp,"*** {}\n", reader.next_string());
    utils::logmesg(lmp,"*** {}\n", reader.next_string());
    utils::logmesg(lmp,"*** {}\n", reader.next_string());
    utils::logmesg(lmp,"*** {}\n", reader.next_string());
    utils::logmesg(lmp,"*** {}\n", reader.next_string());
    utils::logmesg(lmp,"*** {}\n", reader.next_string());
    utils::logmesg(lmp,"*** {}\n", reader.next_string());
    utils::logmesg(lmp,"*** {}\n", reader.next_string());
    utils::logmesg(lmp,"*** {}\n", reader.next_string());
    utils::logmesg(lmp,"*** {}\n", reader.next_string());
    reader.rewind();

    ASSERT_DOUBLE_EQ( reader.next_double(), 123 );
    ASSERT_DOUBLE_EQ( reader.next_double(), -123.45 );
    ASSERT_DOUBLE_EQ( reader.next_double(), DBL_MIN );
    ASSERT_DOUBLE_EQ( reader.next_double(), DBL_MAX );
    ASSERT_DOUBLE_EQ( reader.next_double(), -DBL_MIN );
    ASSERT_DOUBLE_EQ( reader.next_double(), -DBL_MAX );
    ASSERT_EQ( reader.next_int(), MAXSMALLINT );
    ASSERT_EQ( reader.next_tagint(), MAXTAGINT );
    ASSERT_EQ( reader.next_bigint(), MAXBIGINT );
    ASSERT_STREQ( reader.next_string().c_str(), "fooBAR" );
    platform::unlink("potential_reader.file");
}

class OpenPotentialTest : public LAMMPSTest {};

// open for native units
TEST_F(OpenPotentialTest, Sw_native)
{
    int convert_flag = utils::get_supported_conversions(utils::ENERGY);
    BEGIN_CAPTURE_OUTPUT();
    command("units metal");
    FILE *fp    = utils::open_potential("Si.sw", lmp, &convert_flag);
    auto text   = END_CAPTURE_OUTPUT();
    double conv = utils::get_conversion_factor(utils::ENERGY, convert_flag);

    ASSERT_NE(fp, nullptr);
    ASSERT_DOUBLE_EQ(conv, 1.0);
    fclose(fp);
}

// open with supported conversion enabled
TEST_F(OpenPotentialTest, Sw_conv)
{
    int convert_flag = utils::get_supported_conversions(utils::ENERGY);
    ASSERT_EQ(convert_flag, utils::METAL2REAL | utils::REAL2METAL);
    BEGIN_CAPTURE_OUTPUT();
    command("units real");
    FILE *fp    = utils::open_potential("Si.sw", lmp, &convert_flag);
    auto text   = END_CAPTURE_OUTPUT();
    double conv = utils::get_conversion_factor(utils::ENERGY, convert_flag);

    ASSERT_NE(fp, nullptr);
    ASSERT_EQ(convert_flag, utils::METAL2REAL);
    ASSERT_DOUBLE_EQ(conv, 23.060549);
    fclose(fp);
}

// open with conversion disabled
TEST_F(OpenPotentialTest, Sw_noconv)
{
    BEGIN_HIDE_OUTPUT();
    command("units real");
    END_HIDE_OUTPUT();
    TEST_FAILURE(".*Potential.*requires metal units but real.*",
                 utils::open_potential("Si.sw", lmp, nullptr););
    BEGIN_HIDE_OUTPUT();
    command("units lj");
    END_HIDE_OUTPUT();
    int convert_flag = utils::get_supported_conversions(utils::UNKNOWN);
    ASSERT_EQ(convert_flag, utils::NOCONVERT);
}

// open non-existing potential
TEST_F(OpenPotentialTest, No_file)
{
    int convert_flag = utils::get_supported_conversions(utils::ENERGY);
    BEGIN_HIDE_OUTPUT();
    command("units metal");
    FILE *fp = utils::open_potential("Unknown.sw", lmp, &convert_flag);
    END_HIDE_OUTPUT();
    ASSERT_EQ(fp, nullptr);
}

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
