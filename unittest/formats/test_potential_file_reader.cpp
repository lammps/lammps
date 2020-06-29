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

#include "MANYBODY/pair_comb.h"
#include "MANYBODY/pair_comb3.h"
#include "MANYBODY/pair_eim.h"
#include "MANYBODY/pair_gw.h"
#include "MANYBODY/pair_gw_zbl.h"
#include "MANYBODY/pair_nb3b_harmonic.h"
#include "MANYBODY/pair_sw.h"
#include "MANYBODY/pair_tersoff.h"
#include "MANYBODY/pair_tersoff_mod.h"
#include "MANYBODY/pair_tersoff_mod_c.h"
#include "MANYBODY/pair_tersoff_zbl.h"
#include "MANYBODY/pair_vashishta.h"
#include "USER-MISC/pair_tersoff_table.h"
#include "input.h"
#include "lammps.h"
#include "potential_file_reader.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstring>
#include <mpi.h>

using namespace LAMMPS_NS;
using utils::split_words;

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

const int LAMMPS_NS::PairSW::NPARAMS_PER_LINE;
const int LAMMPS_NS::PairComb::NPARAMS_PER_LINE;
const int LAMMPS_NS::PairComb3::NPARAMS_PER_LINE;
const int LAMMPS_NS::PairTersoff::NPARAMS_PER_LINE;
const int LAMMPS_NS::PairTersoffMOD::NPARAMS_PER_LINE;
const int LAMMPS_NS::PairTersoffMODC::NPARAMS_PER_LINE;
const int LAMMPS_NS::PairTersoffZBL::NPARAMS_PER_LINE;
const int LAMMPS_NS::PairGW::NPARAMS_PER_LINE;
const int LAMMPS_NS::PairGWZBL::NPARAMS_PER_LINE;
const int LAMMPS_NS::PairNb3bHarmonic::NPARAMS_PER_LINE;
const int LAMMPS_NS::PairVashishta::NPARAMS_PER_LINE;
const int LAMMPS_NS::PairTersoffTable::NPARAMS_PER_LINE;

class PotentialFileReaderTest : public ::testing::Test {
protected:
    LAMMPS *lmp;

    void SetUp() override
    {
        const char *args[] = {
            "PotentialFileReaderTest", "-log", "none", "-echo", "screen", "-nocite"};
        char **argv = (char **)args;
        int argc    = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete lmp;
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }
};

TEST_F(PotentialFileReaderTest, Sw)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    PotentialFileReader reader(lmp, "Si.sw", "Stillinger-Weber");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairSW::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairSW::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, Comb)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    PotentialFileReader reader(lmp, "ffield.comb", "COMB");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairComb::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairComb::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, Comb3)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    PotentialFileReader reader(lmp, "ffield.comb3", "COMB3");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairComb3::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairComb3::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, Tersoff)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    PotentialFileReader reader(lmp, "Si.tersoff", "Tersoff");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairTersoff::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairTersoff::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, TersoffMod)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    PotentialFileReader reader(lmp, "Si.tersoff.mod", "Tersoff/Mod");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairTersoffMOD::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairTersoffMOD::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, TersoffModC)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    PotentialFileReader reader(lmp, "Si.tersoff.modc", "Tersoff/ModC");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairTersoffMODC::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairTersoffMODC::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, TersoffTable)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    PotentialFileReader reader(lmp, "Si.tersoff", "TersoffTable");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairTersoffTable::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairTersoffTable::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, TersoffZBL)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    PotentialFileReader reader(lmp, "SiC.tersoff.zbl", "Tersoff/ZBL");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairTersoffZBL::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairTersoffZBL::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, GW)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    PotentialFileReader reader(lmp, "SiC.gw", "GW");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairGW::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairGW::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, GWZBL)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    PotentialFileReader reader(lmp, "SiC.gw.zbl", "GW/ZBL");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairGWZBL::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairGWZBL::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, Nb3bHarmonic)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units real");
    PotentialFileReader reader(lmp, "MOH.nb3b.harmonic", "NB3B Harmonic");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairNb3bHarmonic::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairNb3bHarmonic::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, Vashishta)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    PotentialFileReader reader(lmp, "SiC.vashishta", "Vashishta");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairVashishta::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairVashishta::NPARAMS_PER_LINE);
}

TEST_F(PotentialFileReaderTest, UnitConvert)
{
    PotentialFileReader *reader;
    int unit_convert, flag;

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    reader = new PotentialFileReader(lmp, "Si.sw", "Stillinger-Weber");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    unit_convert = reader->get_unit_convert();
    ASSERT_EQ(unit_convert, 0);
    delete reader;

    if (!verbose) ::testing::internal::CaptureStdout();
    flag   = utils::get_supported_conversions(utils::UNKNOWN);
    reader = new PotentialFileReader(lmp, "Si.sw", "Stillinger-Weber", flag);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    unit_convert = reader->get_unit_convert();
    ASSERT_EQ(unit_convert, 0);
    delete reader;

    if (!verbose) ::testing::internal::CaptureStdout();
    flag   = utils::get_supported_conversions(utils::ENERGY);
    reader = new PotentialFileReader(lmp, "Si.sw", "Stillinger-Weber", flag);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    unit_convert = reader->get_unit_convert();
    ASSERT_EQ(unit_convert, 0);
    delete reader;

    if (!verbose) ::testing::internal::CaptureStdout();
    flag   = utils::get_supported_conversions(utils::ENERGY);
    lmp->input->one("units real");
    reader = new PotentialFileReader(lmp, "Si.sw", "Stillinger-Weber", flag);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    unit_convert = reader->get_unit_convert();
    ASSERT_EQ(unit_convert, utils::METAL2REAL);
    delete reader;
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
