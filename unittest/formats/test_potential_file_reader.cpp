#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "potential_file_reader.h"
#include "lammps.h"
#include "utils.h"
#include "MANYBODY/pair_sw.h"
#include "MANYBODY/pair_comb.h"
#include "MANYBODY/pair_comb3.h"
#include "MANYBODY/pair_tersoff.h"

#include <mpi.h>

using namespace LAMMPS_NS;

const int LAMMPS_NS::PairSW::Param::NPARAMS_PER_LINE;
const int LAMMPS_NS::PairComb::Param::NPARAMS_PER_LINE;
const int LAMMPS_NS::PairComb3::Param::NPARAMS_PER_LINE;
const int LAMMPS_NS::PairTersoff::Param::NPARAMS_PER_LINE;

class PotenialFileReaderTest : public ::testing::Test {
protected:
    LAMMPS * lmp;

    void SetUp() override {
        const char *args[] = {"PotentialFileReaderTest", "-log", "none", "-echo", "screen", "-nocite" };
        char **argv = (char **)args;
        int argc = sizeof(args)/sizeof(char *);
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
    }

    void TearDown() override {
        delete lmp;
    }
};

TEST_F(PotenialFileReaderTest, Si) {
    PotentialFileReader reader(lmp, "Si.sw", "Stillinger-Weber");
    auto line = reader.next_line(PairSW::Param::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairSW::Param::NPARAMS_PER_LINE);
}

TEST_F(PotenialFileReaderTest, Comb) {
    PotentialFileReader reader(lmp, "ffield.comb", "COMB");
    auto line = reader.next_line(PairComb::Param::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairComb::Param::NPARAMS_PER_LINE);
}

TEST_F(PotenialFileReaderTest, Comb3) {
    PotentialFileReader reader(lmp, "ffield.comb3", "COMB3");
    auto line = reader.next_line(PairComb3::Param::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairComb3::Param::NPARAMS_PER_LINE);
}

TEST_F(PotenialFileReaderTest, Tersoff) {
    PotentialFileReader reader(lmp, "Si.tersoff", "Tersoff");
    auto line = reader.next_line(PairTersoff::Param::NPARAMS_PER_LINE);
    ASSERT_EQ(utils::count_words(line), PairTersoff::Param::NPARAMS_PER_LINE);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);
    return RUN_ALL_TESTS();
}
