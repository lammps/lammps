/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   This file is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

// Unit tests for the dense direct-Ewald operator used to update the periodic
// QM-image potential during every xTB SCC iteration.

#include "qmmm_xtb_ewald.h"
#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <vector>

using namespace LAMMPS_NS;

namespace {

QMMMXTBEwald::CellMatrix restricted_cell(double xprd, double yprd, double zprd, double xy = 0.0,
                                         double xz = 0.0, double yz = 0.0)
{
    return {xprd, xy, xz, 0.0, yprd, yz, 0.0, 0.0, zprd};
}

TEST(QMMMXTBEwald, ResponseIsSymmetricAndTranslationInvariant)
{
    QMMMXTBEwald ewald;
    ewald.setup(restricted_cell(18.0, 19.0, 20.0), 0.55, {10, 10, 10}, 100);

    const std::vector<double> x    = {4.2, 5.1, 6.3, 9.4, 8.7, 7.6, 13.1, 3.8, 11.2};
    std::vector<double> translated = x;
    for (std::size_t i = 0; i < translated.size(); i += 3) {
        translated[i] += 2.75;
        translated[i + 1] -= 1.25;
        translated[i + 2] += 0.80;
    }

    std::vector<double> response, shifted_response;
    ewald.response(x, response);
    ewald.response(translated, shifted_response);
    ASSERT_EQ(response.size(), 9U);
    ASSERT_EQ(shifted_response.size(), response.size());

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(response[3 * i + j], response[3 * j + i], 1.0e-14);
    }
    for (std::size_t i = 0; i < response.size(); ++i)
        EXPECT_NEAR(response[i], shifted_response[i], 2.0e-14);
}

TEST(QMMMXTBEwald, AnalyticForceMatchesEnergyFiniteDifference)
{
    QMMMXTBEwald ewald;
    ewald.setup(restricted_cell(18.0, 19.0, 20.0), 0.55, {12, 12, 12}, 144);

    std::vector<double> x       = {4.2, 5.1, 6.3, 9.4, 8.7, 7.6, 13.1, 3.8, 11.2};
    const std::vector<double> q = {-0.35, 0.20, 0.15};
    std::vector<double> force(x.size(), 0.0);
    ewald.add_force(x, q, force);

    constexpr double step = 1.0e-5;
    for (std::size_t coordinate = 0; coordinate < x.size(); ++coordinate) {
        x[coordinate] += step;
        const double energy_plus = ewald.energy(x, q);
        x[coordinate] -= 2.0 * step;
        const double energy_minus = ewald.energy(x, q);
        x[coordinate] += step;
        const double finite_difference_force = -(energy_plus - energy_minus) / (2.0 * step);
        EXPECT_NEAR(force[coordinate], finite_difference_force, 2.0e-10);
    }

    for (int dim = 0; dim < 3; ++dim)
        EXPECT_NEAR(force[dim] + force[3 + dim] + force[6 + dim], 0.0, 2.0e-14);
}

TEST(QMMMXTBEwald, TriclinicAnalyticForceMatchesEnergyFiniteDifference)
{
    QMMMXTBEwald ewald;
    ewald.setup(restricted_cell(18.0, 19.0, 20.0, 3.2, -1.4, 2.1), 0.55, {12, 12, 12}, 144);

    std::vector<double> x       = {4.2, 5.1, 6.3, 9.4, 8.7, 7.6, 13.1, 3.8, 11.2};
    const std::vector<double> q = {-0.35, 0.20, 0.15};
    std::vector<double> force(x.size(), 0.0);
    ewald.add_force(x, q, force);

    constexpr double step = 1.0e-5;
    for (std::size_t coordinate = 0; coordinate < x.size(); ++coordinate) {
        x[coordinate] += step;
        const double energy_plus = ewald.energy(x, q);
        x[coordinate] -= 2.0 * step;
        const double energy_minus = ewald.energy(x, q);
        x[coordinate] += step;
        const double finite_difference_force = -(energy_plus - energy_minus) / (2.0 * step);
        EXPECT_NEAR(force[coordinate], finite_difference_force, 2.0e-10);
    }

    for (int dim = 0; dim < 3; ++dim)
        EXPECT_NEAR(force[dim] + force[3 + dim] + force[6 + dim], 0.0, 2.0e-14);
}

TEST(QMMMXTBEwald, TriclinicResponseMatchesReciprocalLatticeReference)
{
    QMMMXTBEwald ewald;
    ewald.setup(restricted_cell(18.0, 19.0, 20.0, 3.2, -1.4, 2.1), 0.55, {10, 10, 10}, 100);

    const std::vector<double> x = {4.2, 5.1, 6.3, 9.4, 8.7, 7.6, 13.1, 3.8, 11.2};
    std::vector<double> response;
    ewald.response(x, response);

    // Independent values from the reciprocal-lattice convention
    // k = 2*pi*H^-T*n with the same half-space integer-vector enumeration.
    const std::vector<double> reference = {
        -0.14934592489933957, -0.13448377698551009, -0.11316607730419709,
        -0.13448377698551009, -0.14934592489933957, -0.13659403366664649,
        -0.11316607730419709, -0.13659403366664649, -0.14934592489933957};
    ASSERT_EQ(response.size(), reference.size());
    for (std::size_t i = 0; i < response.size(); ++i)
        EXPECT_NEAR(response[i], reference[i], 2.0e-14);
}

TEST(QMMMXTBEwald, ErfCorrectionForceMatchesEnergyFiniteDifference)
{
    std::array<double, 3> xi       = {4.2, 5.1, 6.3};
    const std::array<double, 3> xj = {9.4, 8.7, 7.6};
    constexpr double qi            = -0.35;
    constexpr double qj            = 0.20;
    constexpr double alpha         = 0.55;

    double fi[3] = {0.0, 0.0, 0.0};
    double fj[3] = {0.0, 0.0, 0.0};
    QMMMXTBEwald::add_erf_pair_force(xi.data(), xj.data(), qi, qj, alpha, fi, fj);

    auto correction_energy = [&](const std::array<double, 3> &x) {
        const double dx = x[0] - xj[0];
        const double dy = x[1] - xj[1];
        const double dz = x[2] - xj[2];
        const double r  = std::sqrt(dx * dx + dy * dy + dz * dz);
        return -qi * qj * std::erf(alpha * r) / r;
    };

    constexpr double step = 1.0e-5;
    for (int dim = 0; dim < 3; ++dim) {
        xi[dim] += step;
        const double energy_plus = correction_energy(xi);
        xi[dim] -= 2.0 * step;
        const double energy_minus = correction_energy(xi);
        xi[dim] += step;
        const double finite_difference_force = -(energy_plus - energy_minus) / (2.0 * step);
        EXPECT_NEAR(fi[dim], finite_difference_force, 1.0e-12);
        EXPECT_NEAR(fi[dim] + fj[dim], 0.0, 1.0e-15);
    }
}

} // namespace
