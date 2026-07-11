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

#include "../../src/ASPHERE/math_extra_superellipsoids.h"
#include "../../src/math_extra.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cmath>
#include <limits>
#include <vector>
// TODO: consider making a fixture with several setup functions?

static constexpr double EPSILON      = 1e-4;
static constexpr double SOLV_EPSILON = std::numeric_limits<double>::epsilon() * 100;

TEST(HandwrittenSolver, invertible)
{
    double A[16] = {4, 2, 1, 3, 0, 5, 2, 1, 1, 0, 3, 2, 2, 1, 0, 4};

    double b[4] = {23.0, 20.0, 18.0, 20.0};

    double expected_solution[4] = {1.0, 2.0, 3.0, 4.0};

    bool success = MathExtraSuperellipsoids::solve_4x4_robust_unrolled(A, b);

    ASSERT_TRUE(success) << "The solver falsely flagged an invertible matrix as singular.";

    for (int i = 0; i < 4; ++i) {
        ASSERT_NEAR(b[i], expected_solution[i], SOLV_EPSILON) << "Failed at index " << i;
    }
}

TEST(ContactPointAndNormal, sphere)
{
    // First grain
    double xci[3]    = {1.0, 5.246, 3.123};
    double ri        = 2.5;
    double shapei[3] = {ri, ri, ri};
    double Ri[3][3]  = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    double blocki[2] = {2.0, 2.0};
    int flagi        = 0;

    // Second grains
    double xcj[3]    = {2.0, -1.562, 4.607};
    double rj        = 1.25;
    double shapej[3] = {rj, rj, rj};
    double Rj[3][3]  = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    double blockj[2] = {2.0, 2.0};
    int flagj        = 0;

    // Analytical solution
    double X0_analytical[4]  = {rj * xci[0] / (ri + rj) + ri * xcj[0] / (ri + rj),
                                rj * xci[1] / (ri + rj) + ri * xcj[1] / (ri + rj),
                                rj * xci[2] / (ri + rj) + ri * xcj[2] / (ri + rj), rj / ri};
    double nij_analytical[3] = {xcj[0] - xci[0], xcj[1] - xci[1], xcj[2] - xci[2]};
    MathExtra::norm3(nij_analytical);

    int method = MathExtraSuperellipsoids::FORMULATION_ALGEBRAIC;

    // Contact detection
    double X0[4] = {0.0, 0.0, 0.0, 0.0}, nij[3];
    MathExtraSuperellipsoids::determine_contact_point(xci, Ri, shapei, blocki, flagi, xcj, Rj,
                                                      shapej, blockj, flagj, X0, nij, method);

    ASSERT_NEAR(X0[0], X0_analytical[0], EPSILON);
    ASSERT_NEAR(X0[1], X0_analytical[1], EPSILON);
    ASSERT_NEAR(X0[2], X0_analytical[2], EPSILON);
    ASSERT_NEAR(X0[3], X0_analytical[3], EPSILON);

    ASSERT_NEAR(nij[0], nij_analytical[0], EPSILON);
    ASSERT_NEAR(nij[1], nij_analytical[1], EPSILON);
    ASSERT_NEAR(nij[2], nij_analytical[2], EPSILON);

    // Rotational invariance
    double anglei   = 0.456;
    double axisi[3] = {1, 2, 3};
    MathExtra::norm3(axisi);
    double quati[4] = {std::cos(anglei), std::sin(anglei) * axisi[0], std::sin(anglei) * axisi[1],
                       std::sin(anglei) * axisi[2]};
    MathExtra::quat_to_mat(quati, Ri);

    double anglej   = 0.123;
    double axisj[3] = {-1, 2, 1};
    MathExtra::norm3(axisj);
    double quatj[4] = {std::cos(anglej), std::sin(anglej) * axisj[0], std::sin(anglej) * axisj[1],
                       std::sin(anglej) * axisj[2]};
    MathExtra::quat_to_mat(quatj, Rj);

    X0[0] = X0[1] = X0[2] = X0[3] = 0.0;
    MathExtraSuperellipsoids::determine_contact_point(xci, Ri, shapei, blocki, flagi, xcj, Rj,
                                                      shapej, blockj, flagj, X0, nij, method);

    ASSERT_NEAR(X0[0], X0_analytical[0], EPSILON) << "Method: " << method;
    ASSERT_NEAR(X0[1], X0_analytical[1], EPSILON) << "Method: " << method;
    ASSERT_NEAR(X0[2], X0_analytical[2], EPSILON) << "Method: " << method;
    ASSERT_NEAR(X0[3], X0_analytical[3], EPSILON) << "Method: " << method;

    ASSERT_NEAR(nij[0], nij_analytical[0], EPSILON);
    ASSERT_NEAR(nij[1], nij_analytical[1], EPSILON);
    ASSERT_NEAR(nij[2], nij_analytical[2], EPSILON);
}

TEST(ContactPointAndNormal, supersphere_mono)
{
    double r        = 3.456;
    double xci[3]   = {-2 * r, 0.0, 0.0};
    double xcj[3]   = {2 * r, 0.0, 0.0};
    double shape[3] = {r, r, r};
    double R[3][3]  = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

    std::vector<double> blocks = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0};
    int method                 = MathExtraSuperellipsoids::FORMULATION_ALGEBRAIC;

    // Analytical solution
    double X0_analytical[4]  = {0.0, 0.0, 0.0, 1.0};
    double nij_analytical[3] = {1.0, 0.0, 0.0};

    for (auto n : blocks) {
        double block[2] = {n, n};
        int flag        = (n < 2.01) ? 0 : 1;

        // Contact detection
        // Some starting point away from (0,0,0). Possibly bad initial guess so test is demanding
        double X0[4] = {r, -r, 2 * r, 0.0}, nij[3];

        int status = MathExtraSuperellipsoids::determine_contact_point(
            xci, R, shape, block, flag, xcj, R, shape, block, flag, X0, nij, method);

        ASSERT_NEAR(X0[0], X0_analytical[0], EPSILON)
            << "Method: " << method << " | n: " << n << " | status: " << status << " | X0: ["
            << X0[0] << ", " << X0[1] << ", " << X0[2] << ", " << X0[3] << "]";
        ASSERT_NEAR(X0[0], X0_analytical[0], EPSILON) << "Method: " << method;
        ASSERT_NEAR(X0[1], X0_analytical[1], EPSILON) << "Method: " << method;
        ASSERT_NEAR(X0[2], X0_analytical[2], EPSILON) << "Method: " << method;
        ASSERT_NEAR(X0[3], X0_analytical[3], EPSILON) << "Method: " << method;

        ASSERT_NEAR(nij[0], nij_analytical[0], EPSILON);
        ASSERT_NEAR(nij[1], nij_analytical[1], EPSILON);
        ASSERT_NEAR(nij[2], nij_analytical[2], EPSILON);
    }
}

TEST(ContactPointAndNormal, sphere_geometric)
{
    // First grain
    double ri        = 2.5;
    double rj        = 1.25;
    double overlap   = -0.5;
    double xci[3]    = {-(ri - overlap / 2.0), 0.0, 0.0};
    double shapei[3] = {ri, ri, ri};
    double Ri[3][3]  = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    double blocki[2] = {2.0, 2.0};
    int flagi        = 0;

    // Second grains
    double xcj[3] = {rj - overlap / 2.0, 0.0, 0.0};

    double shapej[3] = {rj, rj, rj};
    double Rj[3][3]  = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    double blockj[2] = {2.0, 2.0};
    int flagj        = 0;

    // Analytical solution
    double X0_analytical[4]  = {0.0, 0.0, 0.0, 1.0};
    double nij_analytical[3] = {xcj[0] - xci[0], xcj[1] - xci[1], xcj[2] - xci[2]};
    MathExtra::norm3(nij_analytical);

    int method = MathExtraSuperellipsoids::FORMULATION_GEOMETRIC;

    // Contact detection
    double X0[4] = {.1, .1, .1, 1.0}, nij[3];
    MathExtraSuperellipsoids::determine_contact_point(xci, Ri, shapei, blocki, flagi, xcj, Rj,
                                                      shapej, blockj, flagj, X0, nij, method);

    ASSERT_NEAR(X0[0], X0_analytical[0], EPSILON);
    ASSERT_NEAR(X0[1], X0_analytical[1], EPSILON);
    ASSERT_NEAR(X0[2], X0_analytical[2], EPSILON);
    ASSERT_NEAR(X0[3], X0_analytical[3], EPSILON);

    ASSERT_NEAR(nij[0], nij_analytical[0], EPSILON);
    ASSERT_NEAR(nij[1], nij_analytical[1], EPSILON);
    ASSERT_NEAR(nij[2], nij_analytical[2], EPSILON);

    // Rotational invariance
    double anglei   = 0.456;
    double axisi[3] = {1, 2, 3};
    MathExtra::norm3(axisi);
    double quati[4] = {std::cos(anglei), std::sin(anglei) * axisi[0], std::sin(anglei) * axisi[1],
                       std::sin(anglei) * axisi[2]};
    MathExtra::quat_to_mat(quati, Ri);

    double anglej   = 0.123;
    double axisj[3] = {-1, 2, 1};
    MathExtra::norm3(axisj);
    double quatj[4] = {std::cos(anglej), std::sin(anglej) * axisj[0], std::sin(anglej) * axisj[1],
                       std::sin(anglej) * axisj[2]};
    MathExtra::quat_to_mat(quatj, Rj);

    X0[0] = X0[1] = X0[2] = X0[3] = 0.0;
    MathExtraSuperellipsoids::determine_contact_point(xci, Ri, shapei, blocki, flagi, xcj, Rj,
                                                      shapej, blockj, flagj, X0, nij, method);

    ASSERT_NEAR(X0[0], X0_analytical[0], EPSILON) << "Method: " << method;
    ASSERT_NEAR(X0[1], X0_analytical[1], EPSILON) << "Method: " << method;
    ASSERT_NEAR(X0[2], X0_analytical[2], EPSILON) << "Method: " << method;
    ASSERT_NEAR(X0[3], X0_analytical[3], EPSILON) << "Method: " << method;

    ASSERT_NEAR(nij[0], nij_analytical[0], EPSILON);
    ASSERT_NEAR(nij[1], nij_analytical[1], EPSILON);
    ASSERT_NEAR(nij[2], nij_analytical[2], EPSILON);
}

TEST(ContactPointAndNormal, supersphere_poly_geometric)
{
    double r1      = 3.456;
    double r2      = 2.0 * r1; // Polydisperse: radius_2 = 3 * radius_1
    double overlap = r1 / 10.0;
    double xci[3]  = {-(r1 - overlap / 2.0), 0.0, 0.0};
    double xcj[3]  = {r2 - overlap / 2.0, 0.0, 0.0};

    double shapei[3] = {r1, r1, r1};
    double shapej[3] = {r2, r2, r2};

    // Identity Rotation
    double R[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

    std::vector<double> blocks = {
        2.0, 3.0, 4.0, 5.0, 6.0,
        7.0, 8.0, 9.0, 10.0}; // test would no converge for higher n if not starting along the line
                              // connecting the centers
    int method = MathExtraSuperellipsoids::FORMULATION_GEOMETRIC;

    double nij_analytical[3] = {1.0, 0.0, 0.0};
    double X0_analytical[4]  = {0.0, 0.0, 0.0, 1.0};

    for (auto n : blocks) {
        double block[2] = {n, n};
        int flag        = (n < 2.01) ? 0 : 1;

        // Initial Guess: Offset from 0 to test convergence
        double X0[4] = {overlap, overlap, overlap, 1.0 / 2.0}, nij[3];
        int status   = MathExtraSuperellipsoids::determine_contact_point(
            xci, R, shapei, block, flag, xcj, R, shapej, block, flag, X0, nij, method);

        ASSERT_NEAR(X0[0], X0_analytical[0], EPSILON)
            << "Method: " << method << " | n: " << n << " | status: " << status << " | X0: ["
            << X0[0] << ", " << X0[1] << ", " << X0[2] << ", " << X0[3] << "]";

        ASSERT_EQ(status, 0) << "Failed to converge/detect contact for n=" << n;

        ASSERT_NEAR(X0[0], X0_analytical[0], EPSILON) << "Position X failed for n=" << n;
        ASSERT_NEAR(X0[1], X0_analytical[1], EPSILON) << "Position Y failed for n=" << n;
        ASSERT_NEAR(X0[2], X0_analytical[2], EPSILON) << "Position Z failed for n=" << n;
        ASSERT_NEAR(X0[3], X0_analytical[3], EPSILON) << "Lagrange Multiplier failed for n=" << n;

        ASSERT_NEAR(nij[0], nij_analytical[0], EPSILON) << "Normal X failed for n=" << n;
        ASSERT_NEAR(nij[1], nij_analytical[1], EPSILON) << "Normal Y failed for n=" << n;
        ASSERT_NEAR(nij[2], nij_analytical[2], EPSILON) << "Normal Z failed for n=" << n;
    }
}

// TODO: supersphere_mono with grains overlapping
// TODO: supersphere_poly with grains overlapping
// TODO: more
// for polydisperse solution should be at the radii ratio
