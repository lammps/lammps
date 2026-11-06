// SPDX-License-Identifier: GPL-2.0-or-later OR LGPL-2.1-or-later

#include "input.h"
#include "lammps.h"
#include "library.h"
#include "lmptype.h"

#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <mpi.h>
#include <string>
#include <vector>

using namespace LAMMPS_NS;

namespace {

double force_component(LAMMPS &lmp, tagint id, int component)
{
    const int nlocal = lammps_extract_setting(&lmp, "nlocal");
    const auto *ids = static_cast<tagint *>(lammps_extract_atom(&lmp, "id"));
    const auto *forces = static_cast<double **>(lammps_extract_atom(&lmp, "f"));
    double local = 0.0;
    for (int i = 0; i < nlocal; ++i)
        if (ids[i] == id) local = forces[i][component];
    double total = 0.0;
    MPI_Allreduce(&local, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return total;
}

class DispersionD3SpecialBonds : public ::testing::TestWithParam<const char *> {};

TEST_P(DispersionD3SpecialBonds, ForcesMatchEnergyGradient)
{
    for (const std::string newton : {"on", "off"}) {
        SCOPED_TRACE(newton);
        LAMMPS::argv args = {"DispersionD3SpecialBonds", "-log", "none", "-screen", "none", "-nocite"};
        LAMMPS lmp(args, MPI_COMM_WORLD);
        for (const auto &command : std::vector<std::string>{
                 "units metal", "atom_style bond", "atom_modify map array", "boundary f f f",
                 "newton " + newton, "region box block -5 15 -5 15 -5 15",
                 "create_box 1 box bond/types 1 extra/bond/per/atom 2 extra/special/per/atom 3",
                 "mass 1 28.0855", "create_atoms 1 single 0 0 0",
                 "create_atoms 1 single 2.3 0.1 0", "create_atoms 1 single 4.6 0.2 0.4",
                 "create_atoms 1 single 2.0 2.4 0.3", "create_bonds single/bond 1 1 2",
                 "create_bonds single/bond 1 2 3", "bond_style zero", "bond_coeff 1",
                 "special_bonds lj " + std::string(GetParam()), "pair_style dispersion/d3 bj pbe 8.0 8.0",
                 "pair_coeff * * Si", "neighbor 0.5 bin", "neigh_modify every 1 delay 0 check no",
                 "thermo_style custom step pe"}) {
            lmp.input->one(command);
        }

        constexpr double delta = 1.0e-5;
        for (int id = 1; id <= 4; ++id) {
            SCOPED_TRACE(id);
            lmp.input->one("group probe id " + std::to_string(id));
            for (int d = 0; d < 3; ++d) {
                SCOPED_TRACE(d);
                lmp.input->one("run 0 post no");
                const double analytic = force_component(lmp, id, d);
                auto displace = [&](double step) {
                    std::array<double, 3> move = {0.0, 0.0, 0.0};
                    move[d] = step;
                    lmp.input->one("displace_atoms probe move " + std::to_string(move[0]) + " " +
                                   std::to_string(move[1]) + " " + std::to_string(move[2]) + " units box");
                    lmp.input->one("run 0 post no");
                };
                displace(delta);
                const double plus = lammps_get_thermo(&lmp, "pe");
                displace(-2.0 * delta);
                const double minus = lammps_get_thermo(&lmp, "pe");
                displace(delta);
                ASSERT_TRUE(std::isfinite(analytic));
                ASSERT_TRUE(std::isfinite(plus));
                ASSERT_TRUE(std::isfinite(minus));
                EXPECT_NEAR(analytic, -(plus - minus) / (2.0 * delta), 1.0e-7);
            }
            lmp.input->one("group probe delete");
        }
    }
}

INSTANTIATE_TEST_SUITE_P(Weights, DispersionD3SpecialBonds,
                        ::testing::Values("1.0 1.0 1.0", "0.5 0.8 1.0", "0.0 0.0 1.0"));

} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    lammps_kokkos_finalize();
    MPI_Finalize();
    return result;
}
