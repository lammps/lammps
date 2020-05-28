#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "potential_file_reader.h"
#include "lammps.h"
#include "utils.h"
#include "MANYBODY/pair_sw.h"
#include "MANYBODY/pair_comb.h"
#include "MANYBODY/pair_comb3.h"
#include "MANYBODY/pair_tersoff.h"
#include "MANYBODY/pair_tersoff_mod.h"
#include "MANYBODY/pair_tersoff_mod_c.h"
#include "MANYBODY/pair_tersoff_zbl.h"
#include "MANYBODY/pair_gw.h"
#include "MANYBODY/pair_gw_zbl.h"
#include "MANYBODY/pair_nb3b_harmonic.h"

#include <mpi.h>

using namespace LAMMPS_NS;

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

class PotenialFileReaderTest : public ::testing::Test {
protected:
    LAMMPS * lmp;

    void SetUp() override {
        const char *args[] = {"PotentialFileReaderTest", "-log", "none", "-echo", "screen", "-nocite" };
        char **argv = (char **)args;
        int argc = sizeof(args)/sizeof(char *);
        ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        ::testing::internal::GetCapturedStdout();
    }

    void TearDown() override {
        ::testing::internal::CaptureStdout();
        delete lmp;
        ::testing::internal::GetCapturedStdout();
    }
};

TEST_F(PotenialFileReaderTest, Si) {
    ::testing::internal::CaptureStdout();
    PotentialFileReader reader(lmp, "Si.sw", "Stillinger-Weber");
    ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairSW::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairSW::NPARAMS_PER_LINE);
}

TEST_F(PotenialFileReaderTest, Comb) {
    ::testing::internal::CaptureStdout();
    PotentialFileReader reader(lmp, "ffield.comb", "COMB");
    ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairComb::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairComb::NPARAMS_PER_LINE);
}

TEST_F(PotenialFileReaderTest, Comb3) {
    ::testing::internal::CaptureStdout();
    PotentialFileReader reader(lmp, "ffield.comb3", "COMB3");
    ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairComb3::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairComb3::NPARAMS_PER_LINE);
}

TEST_F(PotenialFileReaderTest, Tersoff) {
    ::testing::internal::CaptureStdout();
    PotentialFileReader reader(lmp, "Si.tersoff", "Tersoff");
    ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairTersoff::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairTersoff::NPARAMS_PER_LINE);
}

TEST_F(PotenialFileReaderTest, TersoffMod) {
    ::testing::internal::CaptureStdout();
    PotentialFileReader reader(lmp, "Si.tersoff.mod", "Tersoff/Mod");
    ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairTersoffMOD::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairTersoffMOD::NPARAMS_PER_LINE);
}

TEST_F(PotenialFileReaderTest, TersoffModC) {
    ::testing::internal::CaptureStdout();
    PotentialFileReader reader(lmp, "Si.tersoff.modc", "Tersoff/ModC");
    ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairTersoffMODC::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairTersoffMODC::NPARAMS_PER_LINE);
}

TEST_F(PotenialFileReaderTest, TersoffZBL) {
    ::testing::internal::CaptureStdout();
    PotentialFileReader reader(lmp, "SiC.tersoff.zbl", "Tersoff/ZBL");
    ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairTersoffZBL::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairTersoffZBL::NPARAMS_PER_LINE);
}

TEST_F(PotenialFileReaderTest, GW) {
    ::testing::internal::CaptureStdout();
    PotentialFileReader reader(lmp, "SiC.gw", "GW");
    ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairGW::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairGW::NPARAMS_PER_LINE);
}

TEST_F(PotenialFileReaderTest, GWZBL) {
    ::testing::internal::CaptureStdout();
    PotentialFileReader reader(lmp, "SiC.gw.zbl", "GW/ZBL");
    ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairGWZBL::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairGWZBL::NPARAMS_PER_LINE);
}

TEST_F(PotenialFileReaderTest, Nb3bHarmonic) {
    ::testing::internal::CaptureStdout();
    PotentialFileReader reader(lmp, "MOH.nb3b.harmonic", "NB3B Harmonic");
    ::testing::internal::GetCapturedStdout();

    auto line = reader.next_line(PairNb3bHarmonic::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairNb3bHarmonic::NPARAMS_PER_LINE);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);
    return RUN_ALL_TESTS();
}
