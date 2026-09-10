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

#include "math_const.h"
#include "math_extra.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathConst::MY_PI;
using MathConst::MY_PI2;

// =========================================================================
// 3-vector operations
// =========================================================================

TEST(MathExtraVec3, copy2)
{
    double v[2]   = {1.5, -2.3};
    double ans[2] = {0.0, 0.0};
    MathExtra::copy2(v, ans);
    EXPECT_DOUBLE_EQ(ans[0], 1.5);
    EXPECT_DOUBLE_EQ(ans[1], -2.3);
}

TEST(MathExtraVec3, copy3)
{
    double v[3]   = {1.0, 2.0, 3.0};
    double ans[3] = {0.0, 0.0, 0.0};
    MathExtra::copy3(v, ans);
    EXPECT_DOUBLE_EQ(ans[0], 1.0);
    EXPECT_DOUBLE_EQ(ans[1], 2.0);
    EXPECT_DOUBLE_EQ(ans[2], 3.0);
}

TEST(MathExtraVec3, zero3)
{
    double v[3] = {5.0, -3.0, 7.0};
    MathExtra::zero3(v);
    EXPECT_DOUBLE_EQ(v[0], 0.0);
    EXPECT_DOUBLE_EQ(v[1], 0.0);
    EXPECT_DOUBLE_EQ(v[2], 0.0);
}

TEST(MathExtraVec3, norm3)
{
    double v[3] = {3.0, 4.0, 0.0};
    MathExtra::norm3(v);
    EXPECT_NEAR(v[0], 0.6, 1e-15);
    EXPECT_NEAR(v[1], 0.8, 1e-15);
    EXPECT_NEAR(v[2], 0.0, 1e-15);

    // length should be 1
    double len = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    EXPECT_NEAR(len, 1.0, 1e-15);
}

TEST(MathExtraVec3, norm3_zero_vector)
{
    double v[3] = {0.0, 0.0, 0.0};
    MathExtra::norm3(v);
    // zero vector should remain zero
    EXPECT_DOUBLE_EQ(v[0], 0.0);
    EXPECT_DOUBLE_EQ(v[1], 0.0);
    EXPECT_DOUBLE_EQ(v[2], 0.0);
}

TEST(MathExtraVec3, normalize3)
{
    double v[3] = {0.0, 3.0, 4.0};
    double ans[3];
    MathExtra::normalize3(v, ans);
    EXPECT_NEAR(ans[0], 0.0, 1e-15);
    EXPECT_NEAR(ans[1], 0.6, 1e-15);
    EXPECT_NEAR(ans[2], 0.8, 1e-15);
}

TEST(MathExtraVec3, snormalize3)
{
    double v[3] = {3.0, 4.0, 0.0};
    double ans[3];
    MathExtra::snormalize3(10.0, v, ans);
    EXPECT_NEAR(ans[0], 6.0, 1e-14);
    EXPECT_NEAR(ans[1], 8.0, 1e-14);
    EXPECT_NEAR(ans[2], 0.0, 1e-14);
    EXPECT_NEAR(std::sqrt(ans[0] * ans[0] + ans[1] * ans[1] + ans[2] * ans[2]), 10.0, 1e-14);
}

TEST(MathExtraVec3, negate3)
{
    double v[3] = {1.0, -2.0, 3.0};
    MathExtra::negate3(v);
    EXPECT_DOUBLE_EQ(v[0], -1.0);
    EXPECT_DOUBLE_EQ(v[1], 2.0);
    EXPECT_DOUBLE_EQ(v[2], -3.0);
}

TEST(MathExtraVec3, scale3_inplace)
{
    double v[3] = {1.0, 2.0, 3.0};
    MathExtra::scale3(2.5, v);
    EXPECT_DOUBLE_EQ(v[0], 2.5);
    EXPECT_DOUBLE_EQ(v[1], 5.0);
    EXPECT_DOUBLE_EQ(v[2], 7.5);
}

TEST(MathExtraVec3, scale3_to_ans)
{
    double v[3] = {1.0, 2.0, 3.0};
    double ans[3];
    MathExtra::scale3(-3.0, v, ans);
    EXPECT_DOUBLE_EQ(ans[0], -3.0);
    EXPECT_DOUBLE_EQ(ans[1], -6.0);
    EXPECT_DOUBLE_EQ(ans[2], -9.0);
}

TEST(MathExtraVec3, add3)
{
    double v1[3] = {1.0, 2.0, 3.0};
    double v2[3] = {4.0, 5.0, 6.0};
    double ans[3];
    MathExtra::add3(v1, v2, ans);
    EXPECT_DOUBLE_EQ(ans[0], 5.0);
    EXPECT_DOUBLE_EQ(ans[1], 7.0);
    EXPECT_DOUBLE_EQ(ans[2], 9.0);
}

TEST(MathExtraVec3, scaleadd3_sv1_plus_v2)
{
    double v1[3] = {1.0, 2.0, 3.0};
    double v2[3] = {10.0, 20.0, 30.0};
    double ans[3];
    MathExtra::scaleadd3(2.0, v1, v2, ans);
    EXPECT_DOUBLE_EQ(ans[0], 12.0);
    EXPECT_DOUBLE_EQ(ans[1], 24.0);
    EXPECT_DOUBLE_EQ(ans[2], 36.0);
}

TEST(MathExtraVec3, scaleadd3_s1v1_plus_s2v2)
{
    double v1[3] = {1.0, 2.0, 3.0};
    double v2[3] = {4.0, 5.0, 6.0};
    double ans[3];
    MathExtra::scaleadd3(2.0, v1, 3.0, v2, ans);
    EXPECT_DOUBLE_EQ(ans[0], 14.0);
    EXPECT_DOUBLE_EQ(ans[1], 19.0);
    EXPECT_DOUBLE_EQ(ans[2], 24.0);
}

TEST(MathExtraVec3, sub3)
{
    double v1[3] = {5.0, 7.0, 9.0};
    double v2[3] = {1.0, 2.0, 3.0};
    double ans[3];
    MathExtra::sub3(v1, v2, ans);
    EXPECT_DOUBLE_EQ(ans[0], 4.0);
    EXPECT_DOUBLE_EQ(ans[1], 5.0);
    EXPECT_DOUBLE_EQ(ans[2], 6.0);
}

TEST(MathExtraVec3, len3)
{
    double v[3] = {3.0, 4.0, 0.0};
    EXPECT_DOUBLE_EQ(MathExtra::len3(v), 5.0);

    double v2[3] = {1.0, 2.0, 2.0};
    EXPECT_DOUBLE_EQ(MathExtra::len3(v2), 3.0);
}

TEST(MathExtraVec3, lensq3)
{
    double v[3] = {1.0, 2.0, 3.0};
    EXPECT_DOUBLE_EQ(MathExtra::lensq3(v), 14.0);
}

TEST(MathExtraVec3, distsq3)
{
    double v1[3] = {1.0, 2.0, 3.0};
    double v2[3] = {4.0, 6.0, 3.0};
    EXPECT_DOUBLE_EQ(MathExtra::distsq3(v1, v2), 25.0);
}

TEST(MathExtraVec3, dot3)
{
    double v1[3] = {1.0, 2.0, 3.0};
    double v2[3] = {4.0, 5.0, 6.0};
    EXPECT_DOUBLE_EQ(MathExtra::dot3(v1, v2), 32.0);
}

TEST(MathExtraVec3, dot3_perpendicular)
{
    double v1[3] = {1.0, 0.0, 0.0};
    double v2[3] = {0.0, 1.0, 0.0};
    EXPECT_DOUBLE_EQ(MathExtra::dot3(v1, v2), 0.0);
}

TEST(MathExtraVec3, cross3)
{
    double v1[3] = {1.0, 0.0, 0.0};
    double v2[3] = {0.0, 1.0, 0.0};
    double ans[3];
    MathExtra::cross3(v1, v2, ans);
    EXPECT_DOUBLE_EQ(ans[0], 0.0);
    EXPECT_DOUBLE_EQ(ans[1], 0.0);
    EXPECT_DOUBLE_EQ(ans[2], 1.0);
}

TEST(MathExtraVec3, cross3_general)
{
    double v1[3] = {1.0, 2.0, 3.0};
    double v2[3] = {4.0, 5.0, 6.0};
    double ans[3];
    MathExtra::cross3(v1, v2, ans);
    // (2*6-3*5, 3*4-1*6, 1*5-2*4) = (-3, 6, -3)
    EXPECT_DOUBLE_EQ(ans[0], -3.0);
    EXPECT_DOUBLE_EQ(ans[1], 6.0);
    EXPECT_DOUBLE_EQ(ans[2], -3.0);
}

TEST(MathExtraVec3, cross3_anticommutative)
{
    double v1[3] = {1.0, 2.0, 3.0};
    double v2[3] = {4.0, 5.0, 6.0};
    double c1[3], c2[3];
    MathExtra::cross3(v1, v2, c1);
    MathExtra::cross3(v2, v1, c2);
    EXPECT_DOUBLE_EQ(c1[0], -c2[0]);
    EXPECT_DOUBLE_EQ(c1[1], -c2[1]);
    EXPECT_DOUBLE_EQ(c1[2], -c2[2]);
}

// =========================================================================
// 3x3 matrix operations
// =========================================================================

TEST(MathExtraMat3, zeromat3)
{
    double m[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    MathExtra::zeromat3(m);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            EXPECT_DOUBLE_EQ(m[i][j], 0.0) << "at [" << i << "][" << j << "]";
}

TEST(MathExtraMat3, transpose3)
{
    double m[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    double ans[3][3];
    MathExtra::transpose3(m, ans);
    EXPECT_DOUBLE_EQ(ans[0][0], 1.0);
    EXPECT_DOUBLE_EQ(ans[0][1], 4.0);
    EXPECT_DOUBLE_EQ(ans[0][2], 7.0);
    EXPECT_DOUBLE_EQ(ans[1][0], 2.0);
    EXPECT_DOUBLE_EQ(ans[1][1], 5.0);
    EXPECT_DOUBLE_EQ(ans[1][2], 8.0);
    EXPECT_DOUBLE_EQ(ans[2][0], 3.0);
    EXPECT_DOUBLE_EQ(ans[2][1], 6.0);
    EXPECT_DOUBLE_EQ(ans[2][2], 9.0);
}

TEST(MathExtraMat3, col2mat)
{
    double ex[3] = {1, 4, 7};
    double ey[3] = {2, 5, 8};
    double ez[3] = {3, 6, 9};
    double m[3][3];
    MathExtra::col2mat(ex, ey, ez, m);
    // columns are ex, ey, ez
    EXPECT_DOUBLE_EQ(m[0][0], 1.0);
    EXPECT_DOUBLE_EQ(m[1][0], 4.0);
    EXPECT_DOUBLE_EQ(m[2][0], 7.0);
    EXPECT_DOUBLE_EQ(m[0][1], 2.0);
    EXPECT_DOUBLE_EQ(m[1][1], 5.0);
    EXPECT_DOUBLE_EQ(m[2][1], 8.0);
    EXPECT_DOUBLE_EQ(m[0][2], 3.0);
    EXPECT_DOUBLE_EQ(m[1][2], 6.0);
    EXPECT_DOUBLE_EQ(m[2][2], 9.0);
}

TEST(MathExtraMat3, det3_identity)
{
    double m[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    EXPECT_DOUBLE_EQ(MathExtra::det3(m), 1.0);
}

TEST(MathExtraMat3, det3_known)
{
    double m[3][3] = {{1, 2, 3}, {0, 1, 4}, {5, 6, 0}};
    // det = 1*(0-24) - 2*(0-20) + 3*(0-5) = -24+40-15 = 1
    EXPECT_DOUBLE_EQ(MathExtra::det3(m), 1.0);
}

TEST(MathExtraMat3, det3_singular)
{
    // linearly dependent rows
    double m[3][3] = {{1, 2, 3}, {2, 4, 6}, {1, 1, 1}};
    EXPECT_DOUBLE_EQ(MathExtra::det3(m), 0.0);
}

TEST(MathExtraMat3, plus3)
{
    double m1[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    double m2[3][3] = {{9, 8, 7}, {6, 5, 4}, {3, 2, 1}};
    double ans[3][3];
    MathExtra::plus3(m1, m2, ans);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            EXPECT_DOUBLE_EQ(ans[i][j], 10.0) << "at [" << i << "][" << j << "]";
}

TEST(MathExtraMat3, minus3)
{
    double m1[3][3] = {{10, 20, 30}, {40, 50, 60}, {70, 80, 90}};
    double m2[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    double ans[3][3];
    MathExtra::minus3(m1, m2, ans);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            EXPECT_DOUBLE_EQ(ans[i][j], (i * 3 + j + 1) * 9.0) << "at [" << i << "][" << j << "]";
}

TEST(MathExtraMat3, times3_identity)
{
    double m[3][3]  = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    double id[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double ans[3][3];
    MathExtra::times3(m, id, ans);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            EXPECT_DOUBLE_EQ(ans[i][j], m[i][j]) << "at [" << i << "][" << j << "]";
}

TEST(MathExtraMat3, times3_known)
{
    double A[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    double B[3][3] = {{9, 8, 7}, {6, 5, 4}, {3, 2, 1}};
    double ans[3][3];
    MathExtra::times3(A, B, ans);
    // Row 0: 1*9+2*6+3*3=30, 1*8+2*5+3*2=24, 1*7+2*4+3*1=18
    EXPECT_DOUBLE_EQ(ans[0][0], 30.0);
    EXPECT_DOUBLE_EQ(ans[0][1], 24.0);
    EXPECT_DOUBLE_EQ(ans[0][2], 18.0);
    // Row 1: 4*9+5*6+6*3=84, 4*8+5*5+6*2=69, 4*7+5*4+6*1=54
    EXPECT_DOUBLE_EQ(ans[1][0], 84.0);
    EXPECT_DOUBLE_EQ(ans[1][1], 69.0);
    EXPECT_DOUBLE_EQ(ans[1][2], 54.0);
    // Row 2: 7*9+8*6+9*3=138, 7*8+8*5+9*2=114, 7*7+8*4+9*1=90
    EXPECT_DOUBLE_EQ(ans[2][0], 138.0);
    EXPECT_DOUBLE_EQ(ans[2][1], 114.0);
    EXPECT_DOUBLE_EQ(ans[2][2], 90.0);
}

TEST(MathExtraMat3, transpose_times3)
{
    // transpose_times3(A, B) = A^T * B
    double A[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    double B[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double ans[3][3];
    MathExtra::transpose_times3(A, B, ans);
    // A^T * I = A^T
    EXPECT_DOUBLE_EQ(ans[0][0], 1.0);
    EXPECT_DOUBLE_EQ(ans[0][1], 4.0);
    EXPECT_DOUBLE_EQ(ans[0][2], 7.0);
    EXPECT_DOUBLE_EQ(ans[1][0], 2.0);
    EXPECT_DOUBLE_EQ(ans[1][1], 5.0);
    EXPECT_DOUBLE_EQ(ans[1][2], 8.0);
    EXPECT_DOUBLE_EQ(ans[2][0], 3.0);
    EXPECT_DOUBLE_EQ(ans[2][1], 6.0);
    EXPECT_DOUBLE_EQ(ans[2][2], 9.0);
}

TEST(MathExtraMat3, times3_transpose)
{
    // times3_transpose(A, B) = A * B^T
    double A[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double B[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    double ans[3][3];
    MathExtra::times3_transpose(A, B, ans);
    // I * B^T = B^T
    EXPECT_DOUBLE_EQ(ans[0][0], 1.0);
    EXPECT_DOUBLE_EQ(ans[0][1], 4.0);
    EXPECT_DOUBLE_EQ(ans[0][2], 7.0);
    EXPECT_DOUBLE_EQ(ans[1][0], 2.0);
    EXPECT_DOUBLE_EQ(ans[1][1], 5.0);
    EXPECT_DOUBLE_EQ(ans[1][2], 8.0);
    EXPECT_DOUBLE_EQ(ans[2][0], 3.0);
    EXPECT_DOUBLE_EQ(ans[2][1], 6.0);
    EXPECT_DOUBLE_EQ(ans[2][2], 9.0);
}

TEST(MathExtraMat3, diag_times3)
{
    double d[3]    = {2.0, 3.0, 4.0};
    double m[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    double ans[3][3];
    MathExtra::diag_times3(d, m, ans);
    EXPECT_DOUBLE_EQ(ans[0][0], 2.0);
    EXPECT_DOUBLE_EQ(ans[0][1], 4.0);
    EXPECT_DOUBLE_EQ(ans[0][2], 6.0);
    EXPECT_DOUBLE_EQ(ans[1][0], 12.0);
    EXPECT_DOUBLE_EQ(ans[1][1], 15.0);
    EXPECT_DOUBLE_EQ(ans[1][2], 18.0);
    EXPECT_DOUBLE_EQ(ans[2][0], 28.0);
    EXPECT_DOUBLE_EQ(ans[2][1], 32.0);
    EXPECT_DOUBLE_EQ(ans[2][2], 36.0);
}

TEST(MathExtraMat3, times3_diag)
{
    double m[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    double d[3]    = {2.0, 3.0, 4.0};
    double ans[3][3];
    MathExtra::times3_diag(m, d, ans);
    EXPECT_DOUBLE_EQ(ans[0][0], 2.0);
    EXPECT_DOUBLE_EQ(ans[0][1], 6.0);
    EXPECT_DOUBLE_EQ(ans[0][2], 12.0);
    EXPECT_DOUBLE_EQ(ans[1][0], 8.0);
    EXPECT_DOUBLE_EQ(ans[1][1], 15.0);
    EXPECT_DOUBLE_EQ(ans[1][2], 24.0);
    EXPECT_DOUBLE_EQ(ans[2][0], 14.0);
    EXPECT_DOUBLE_EQ(ans[2][1], 24.0);
    EXPECT_DOUBLE_EQ(ans[2][2], 36.0);
}

TEST(MathExtraMat3, invert3)
{
    // Use a simple invertible matrix: [[2,1,0],[0,3,1],[1,0,2]]
    double m[3][3] = {{2, 1, 0}, {0, 3, 1}, {1, 0, 2}};
    double inv[3][3];
    MathExtra::invert3(m, inv);

    // Verify M * M^-1 = I
    double prod[3][3];
    MathExtra::times3(m, inv, prod);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(prod[i][j], expected, 1e-14) << "at [" << i << "][" << j << "]";
        }
}

TEST(MathExtraMat3, invert3_identity)
{
    double m[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double inv[3][3];
    MathExtra::invert3(m, inv);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(inv[i][j], expected, 1e-14);
        }
}

TEST(MathExtraMat3, matvec_array)
{
    double m[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    double v[3]    = {1.0, 0.0, 0.0};
    double ans[3];
    MathExtra::matvec(m, v, ans);
    EXPECT_DOUBLE_EQ(ans[0], 1.0);
    EXPECT_DOUBLE_EQ(ans[1], 4.0);
    EXPECT_DOUBLE_EQ(ans[2], 7.0);
}

TEST(MathExtraMat3, matvec_columns)
{
    double ex[3] = {1, 4, 7};
    double ey[3] = {2, 5, 8};
    double ez[3] = {3, 6, 9};
    double v[3]  = {1.0, 1.0, 1.0};
    double ans[3];
    MathExtra::matvec(ex, ey, ez, v, ans);
    EXPECT_DOUBLE_EQ(ans[0], 6.0);
    EXPECT_DOUBLE_EQ(ans[1], 15.0);
    EXPECT_DOUBLE_EQ(ans[2], 24.0);
}

TEST(MathExtraMat3, transpose_matvec_array)
{
    double m[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    double v[3]    = {1.0, 1.0, 1.0};
    double ans[3];
    MathExtra::transpose_matvec(m, v, ans);
    // M^T * v: col0 dot v = 1+4+7=12, col1 dot v = 2+5+8=15, col2 dot v = 3+6+9=18
    EXPECT_DOUBLE_EQ(ans[0], 12.0);
    EXPECT_DOUBLE_EQ(ans[1], 15.0);
    EXPECT_DOUBLE_EQ(ans[2], 18.0);
}

TEST(MathExtraMat3, transpose_matvec_columns)
{
    double ex[3] = {1, 2, 3};
    double ey[3] = {4, 5, 6};
    double ez[3] = {7, 8, 9};
    double v[3]  = {1.0, 0.0, 0.0};
    double ans[3];
    MathExtra::transpose_matvec(ex, ey, ez, v, ans);
    // ex . v, ey . v, ez . v
    EXPECT_DOUBLE_EQ(ans[0], 1.0);
    EXPECT_DOUBLE_EQ(ans[1], 4.0);
    EXPECT_DOUBLE_EQ(ans[2], 7.0);
}

TEST(MathExtraMat3, transpose_diag3)
{
    double m[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    double d[3]    = {2.0, 3.0, 4.0};
    double ans[3][3];
    MathExtra::transpose_diag3(m, d, ans);
    // ans = M^T * diag(d)
    EXPECT_DOUBLE_EQ(ans[0][0], 1.0 * 2.0);
    EXPECT_DOUBLE_EQ(ans[0][1], 4.0 * 3.0);
    EXPECT_DOUBLE_EQ(ans[0][2], 7.0 * 4.0);
    EXPECT_DOUBLE_EQ(ans[1][0], 2.0 * 2.0);
    EXPECT_DOUBLE_EQ(ans[1][1], 5.0 * 3.0);
    EXPECT_DOUBLE_EQ(ans[1][2], 8.0 * 4.0);
    EXPECT_DOUBLE_EQ(ans[2][0], 3.0 * 2.0);
    EXPECT_DOUBLE_EQ(ans[2][1], 6.0 * 3.0);
    EXPECT_DOUBLE_EQ(ans[2][2], 9.0 * 4.0);
}

TEST(MathExtraMat3, vecmat)
{
    double v[3]    = {1.0, 2.0, 3.0};
    double m[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    double ans[3];
    MathExtra::vecmat(v, m, ans);
    // v^T * M: 1*1+2*4+3*7=30, 1*2+2*5+3*8=36, 1*3+2*6+3*9=42
    EXPECT_DOUBLE_EQ(ans[0], 30.0);
    EXPECT_DOUBLE_EQ(ans[1], 36.0);
    EXPECT_DOUBLE_EQ(ans[2], 42.0);
}

TEST(MathExtraMat3, scalar_times3)
{
    double m[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    MathExtra::scalar_times3(2.0, m);
    EXPECT_DOUBLE_EQ(m[0][0], 2.0);
    EXPECT_DOUBLE_EQ(m[1][1], 10.0);
    EXPECT_DOUBLE_EQ(m[2][2], 18.0);
    EXPECT_DOUBLE_EQ(m[0][2], 6.0);
}

TEST(MathExtraMat3, outer3)
{
    double v1[3] = {1.0, 2.0, 3.0};
    double v2[3] = {4.0, 5.0, 6.0};
    double ans[3][3];
    MathExtra::outer3(v1, v2, ans);
    EXPECT_DOUBLE_EQ(ans[0][0], 4.0);
    EXPECT_DOUBLE_EQ(ans[0][1], 5.0);
    EXPECT_DOUBLE_EQ(ans[0][2], 6.0);
    EXPECT_DOUBLE_EQ(ans[1][0], 8.0);
    EXPECT_DOUBLE_EQ(ans[1][1], 10.0);
    EXPECT_DOUBLE_EQ(ans[1][2], 12.0);
    EXPECT_DOUBLE_EQ(ans[2][0], 12.0);
    EXPECT_DOUBLE_EQ(ans[2][1], 15.0);
    EXPECT_DOUBLE_EQ(ans[2][2], 18.0);
}

// =========================================================================
// Linear solver (mldivide3)
// =========================================================================

TEST(MathExtraSolver, mldivide3_simple)
{
    // Solve [[2,0,0],[0,3,0],[0,0,4]] * x = [6,9,12] => x = [3,3,3]
    double m[3][3] = {{2, 0, 0}, {0, 3, 0}, {0, 0, 4}};
    double v[3]    = {6.0, 9.0, 12.0};
    double ans[3];
    int ret = MathExtra::mldivide3(m, v, ans);
    EXPECT_EQ(ret, 0);
    EXPECT_NEAR(ans[0], 3.0, 1e-14);
    EXPECT_NEAR(ans[1], 3.0, 1e-14);
    EXPECT_NEAR(ans[2], 3.0, 1e-14);
}

TEST(MathExtraSolver, mldivide3_general)
{
    // Solve [[1,2,3],[0,1,4],[5,6,0]] * x = [1,0,0]
    // det = 1 (from earlier test)
    double m[3][3] = {{1, 2, 3}, {0, 1, 4}, {5, 6, 0}};
    double v[3]    = {1.0, 0.0, 0.0};
    double ans[3];
    int ret = MathExtra::mldivide3(m, v, ans);
    EXPECT_EQ(ret, 0);

    // Verify the solution by multiplying
    double check[3];
    MathExtra::matvec(m, ans, check);
    EXPECT_NEAR(check[0], v[0], 1e-13);
    EXPECT_NEAR(check[1], v[1], 1e-13);
    EXPECT_NEAR(check[2], v[2], 1e-13);
}

TEST(MathExtraSolver, mldivide3_singular)
{
    // Singular matrix should return 1
    double m[3][3] = {{1, 2, 3}, {2, 4, 6}, {1, 1, 1}};
    double v[3]    = {1.0, 2.0, 3.0};
    double ans[3];
    int ret = MathExtra::mldivide3(m, v, ans);
    EXPECT_EQ(ret, 1);
}

// =========================================================================
// Quaternion operations
// =========================================================================

TEST(MathExtraQuat, qnormalize)
{
    double q[4] = {1.0, 2.0, 3.0, 4.0};
    MathExtra::qnormalize(q);
    double norm = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    EXPECT_NEAR(norm, 1.0, 1e-15);
}

TEST(MathExtraQuat, qnormalize_unit)
{
    double q[4] = {1.0, 0.0, 0.0, 0.0};
    MathExtra::qnormalize(q);
    EXPECT_DOUBLE_EQ(q[0], 1.0);
    EXPECT_DOUBLE_EQ(q[1], 0.0);
    EXPECT_DOUBLE_EQ(q[2], 0.0);
    EXPECT_DOUBLE_EQ(q[3], 0.0);
}

TEST(MathExtraQuat, qconjugate)
{
    double q[4] = {0.5, 0.5, 0.5, 0.5};
    double qc[4];
    MathExtra::qconjugate(q, qc);
    EXPECT_DOUBLE_EQ(qc[0], 0.5);
    EXPECT_DOUBLE_EQ(qc[1], -0.5);
    EXPECT_DOUBLE_EQ(qc[2], -0.5);
    EXPECT_DOUBLE_EQ(qc[3], -0.5);
}

TEST(MathExtraQuat, quatquat_identity)
{
    // q * identity = q
    double q[4]  = {0.5, 0.5, 0.5, 0.5};
    double id[4] = {1.0, 0.0, 0.0, 0.0};
    double c[4];
    MathExtra::quatquat(q, id, c);
    for (int i = 0; i < 4; i++)
        EXPECT_NEAR(c[i], q[i], 1e-15) << "component " << i;
}

TEST(MathExtraQuat, quatquat_inverse)
{
    // q * conjugate(q) = identity (for unit quaternion)
    double q[4] = {0.5, 0.5, 0.5, 0.5};
    double qc[4];
    MathExtra::qconjugate(q, qc);
    double c[4];
    MathExtra::quatquat(q, qc, c);
    EXPECT_NEAR(c[0], 1.0, 1e-15);
    EXPECT_NEAR(c[1], 0.0, 1e-15);
    EXPECT_NEAR(c[2], 0.0, 1e-15);
    EXPECT_NEAR(c[3], 0.0, 1e-15);
}

TEST(MathExtraQuat, quatrotvec_identity)
{
    // identity quaternion should not rotate
    double q[4] = {1.0, 0.0, 0.0, 0.0};
    double v[3] = {1.0, 2.0, 3.0};
    double c[3];
    MathExtra::quatrotvec(q, v, c);
    EXPECT_NEAR(c[0], 1.0, 1e-14);
    EXPECT_NEAR(c[1], 2.0, 1e-14);
    EXPECT_NEAR(c[2], 3.0, 1e-14);
}

TEST(MathExtraQuat, quatrotvec_90deg_z)
{
    // 90-degree rotation around z: (1,0,0) -> (0,1,0)
    double angle = MY_PI2;
    double q[4]  = {std::cos(angle / 2.0), 0.0, 0.0, std::sin(angle / 2.0)};
    double v[3]  = {1.0, 0.0, 0.0};
    double c[3];
    MathExtra::quatrotvec(q, v, c);
    EXPECT_NEAR(c[0], 0.0, 1e-14);
    EXPECT_NEAR(c[1], 1.0, 1e-14);
    EXPECT_NEAR(c[2], 0.0, 1e-14);
}

TEST(MathExtraQuat, quatrotvec_preserves_length)
{
    double angle   = 1.23;
    double axis[3] = {1.0, 1.0, 1.0};
    MathExtra::norm3(axis);
    double q[4];
    MathExtra::axisangle_to_quat(axis, angle, q);

    double v[3] = {3.0, 4.0, 5.0};
    double c[3];
    MathExtra::quatrotvec(q, v, c);

    double len_v = MathExtra::len3(v);
    double len_c = MathExtra::len3(c);
    EXPECT_NEAR(len_c, len_v, 1e-13);
}

TEST(MathExtraQuat, axisangle_to_quat)
{
    // 180-degree rotation around z
    double axis[3] = {0.0, 0.0, 1.0};
    double q[4];
    MathExtra::axisangle_to_quat(axis, MY_PI, q);
    EXPECT_NEAR(q[0], std::cos(MY_PI2), 1e-15);
    EXPECT_NEAR(q[1], 0.0, 1e-15);
    EXPECT_NEAR(q[2], 0.0, 1e-15);
    EXPECT_NEAR(q[3], std::sin(MY_PI2), 1e-15);
}

TEST(MathExtraQuat, axisangle_to_quat_zero_angle)
{
    double axis[3] = {1.0, 0.0, 0.0};
    double q[4];
    MathExtra::axisangle_to_quat(axis, 0.0, q);
    EXPECT_NEAR(q[0], 1.0, 1e-15);
    EXPECT_NEAR(q[1], 0.0, 1e-15);
    EXPECT_NEAR(q[2], 0.0, 1e-15);
    EXPECT_NEAR(q[3], 0.0, 1e-15);
}

TEST(MathExtraQuat, vecquat_quatvec_consistency)
{
    // vecquat(a,b) treats a as (0,a); quatvec(b,a) treats a as (0,a)
    // vecquat(a,b) = (0,a)*b and quatvec(b,a) = b*(0,a)
    // These are NOT the same in general, but let's test each works correctly

    double a[3] = {1.0, 0.0, 0.0};
    double b[4] = {1.0, 0.0, 0.0, 0.0}; // identity
    double c1[4], c2[4];

    // vecquat: (0,a)*b where b is identity
    MathExtra::vecquat(a, b, c1);
    // Should give (0, a) since (0,a)*1 = (0,a)
    EXPECT_NEAR(c1[0], 0.0, 1e-15);
    EXPECT_NEAR(c1[1], 1.0, 1e-15);
    EXPECT_NEAR(c1[2], 0.0, 1e-15);
    EXPECT_NEAR(c1[3], 0.0, 1e-15);

    // quatvec: b*(0,a) where b is identity
    MathExtra::quatvec(b, a, c2);
    // Should also give (0, a) since 1*(0,a) = (0,a)
    EXPECT_NEAR(c2[0], 0.0, 1e-15);
    EXPECT_NEAR(c2[1], 1.0, 1e-15);
    EXPECT_NEAR(c2[2], 0.0, 1e-15);
    EXPECT_NEAR(c2[3], 0.0, 1e-15);
}

// =========================================================================
// Quaternion <-> Matrix conversions
// =========================================================================

TEST(MathExtraQuatConvert, quat_to_mat_identity)
{
    double q[4] = {1.0, 0.0, 0.0, 0.0};
    double mat[3][3];
    MathExtra::quat_to_mat(q, mat);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(mat[i][j], expected, 1e-15) << "at [" << i << "][" << j << "]";
        }
}

TEST(MathExtraQuatConvert, quat_to_mat_is_rotation)
{
    // Check that the resulting matrix is a proper rotation (det=1, R^T R = I)
    double angle   = 1.234;
    double axis[3] = {1.0, 2.0, 3.0};
    MathExtra::norm3(axis);
    double q[4];
    MathExtra::axisangle_to_quat(axis, angle, q);

    double mat[3][3];
    MathExtra::quat_to_mat(q, mat);

    EXPECT_NEAR(MathExtra::det3(mat), 1.0, 1e-13);

    double prod[3][3];
    MathExtra::transpose_times3(mat, mat, prod);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(prod[i][j], expected, 1e-13) << "at [" << i << "][" << j << "]";
        }
}

TEST(MathExtraQuatConvert, quat_to_mat_trans_is_transpose)
{
    double angle   = 0.789;
    double axis[3] = {0.0, 1.0, 0.0};
    double q[4];
    MathExtra::axisangle_to_quat(axis, angle, q);

    double mat[3][3], mat_trans[3][3];
    MathExtra::quat_to_mat(q, mat);
    MathExtra::quat_to_mat_trans(q, mat_trans);

    double transposed[3][3];
    MathExtra::transpose3(mat, transposed);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            EXPECT_NEAR(mat_trans[i][j], transposed[i][j], 1e-14)
                << "at [" << i << "][" << j << "]";
}

TEST(MathExtraQuatConvert, quat_to_mat_roundtrip)
{
    // quat -> mat -> quat should give equivalent quaternion
    double angle   = 1.5;
    double axis[3] = {1.0, -1.0, 0.5};
    MathExtra::norm3(axis);
    double q_orig[4];
    MathExtra::axisangle_to_quat(axis, angle, q_orig);

    double mat[3][3];
    MathExtra::quat_to_mat(q_orig, mat);

    double q_back[4];
    MathExtra::mat_to_quat(mat, q_back);

    // Quaternions q and -q represent same rotation
    double sign = (q_orig[0] * q_back[0] >= 0.0) ? 1.0 : -1.0;
    for (int i = 0; i < 4; i++)
        EXPECT_NEAR(q_back[i] * sign, q_orig[i], 1e-13) << "component " << i;
}

TEST(MathExtraQuatConvert, q_to_exyz_identity)
{
    double q[4] = {1.0, 0.0, 0.0, 0.0};
    double ex[3], ey[3], ez[3];
    MathExtra::q_to_exyz(q, ex, ey, ez);
    EXPECT_NEAR(ex[0], 1.0, 1e-15);
    EXPECT_NEAR(ex[1], 0.0, 1e-15);
    EXPECT_NEAR(ex[2], 0.0, 1e-15);
    EXPECT_NEAR(ey[0], 0.0, 1e-15);
    EXPECT_NEAR(ey[1], 1.0, 1e-15);
    EXPECT_NEAR(ey[2], 0.0, 1e-15);
    EXPECT_NEAR(ez[0], 0.0, 1e-15);
    EXPECT_NEAR(ez[1], 0.0, 1e-15);
    EXPECT_NEAR(ez[2], 1.0, 1e-15);
}

TEST(MathExtraQuatConvert, q_to_exyz_and_back)
{
    double angle   = 0.567;
    double axis[3] = {0.0, 0.0, 1.0};
    double q_orig[4];
    MathExtra::axisangle_to_quat(axis, angle, q_orig);

    double ex[3], ey[3], ez[3];
    MathExtra::q_to_exyz(q_orig, ex, ey, ez);

    double q_back[4];
    MathExtra::exyz_to_q(ex, ey, ez, q_back);

    double sign = (q_orig[0] * q_back[0] >= 0.0) ? 1.0 : -1.0;
    for (int i = 0; i < 4; i++)
        EXPECT_NEAR(q_back[i] * sign, q_orig[i], 1e-13) << "component " << i;
}

TEST(MathExtraQuatConvert, quat_to_mat_matches_q_to_exyz)
{
    double angle   = 2.1;
    double axis[3] = {1.0, 0.0, 0.0};
    double q[4];
    MathExtra::axisangle_to_quat(axis, angle, q);

    double mat[3][3];
    MathExtra::quat_to_mat(q, mat);

    double ex[3], ey[3], ez[3];
    MathExtra::q_to_exyz(q, ex, ey, ez);

    // Columns of mat should match ex, ey, ez
    EXPECT_NEAR(mat[0][0], ex[0], 1e-14);
    EXPECT_NEAR(mat[1][0], ex[1], 1e-14);
    EXPECT_NEAR(mat[2][0], ex[2], 1e-14);
    EXPECT_NEAR(mat[0][1], ey[0], 1e-14);
    EXPECT_NEAR(mat[1][1], ey[1], 1e-14);
    EXPECT_NEAR(mat[2][1], ey[2], 1e-14);
    EXPECT_NEAR(mat[0][2], ez[0], 1e-14);
    EXPECT_NEAR(mat[1][2], ez[1], 1e-14);
    EXPECT_NEAR(mat[2][2], ez[2], 1e-14);
}

// =========================================================================
// Rotation matrix builders
// =========================================================================

TEST(MathExtraRotation, BuildRxMatrix)
{
    double R[3][3];
    MathExtra::BuildRxMatrix(R, 0.0);
    // Zero angle => identity
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(R[i][j], expected, 1e-15) << "at [" << i << "][" << j << "]";
        }
}

TEST(MathExtraRotation, BuildRxMatrix_is_rotation)
{
    double R[3][3];
    MathExtra::BuildRxMatrix(R, 0.3);
    EXPECT_NEAR(MathExtra::det3(R), 1.0, 1e-13);

    // R^T * R = I
    double prod[3][3];
    MathExtra::transpose_times3(R, R, prod);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(prod[i][j], expected, 1e-13);
        }
}

TEST(MathExtraRotation, BuildRxMatrix_preserves_x)
{
    double R[3][3];
    MathExtra::BuildRxMatrix(R, 0.5);
    // x-axis should be unchanged by Rx rotation
    double v[3] = {1.0, 0.0, 0.0};
    double ans[3];
    MathExtra::matvec(R, v, ans);
    EXPECT_NEAR(ans[0], 1.0, 1e-14);
    EXPECT_NEAR(ans[1], 0.0, 1e-14);
    EXPECT_NEAR(ans[2], 0.0, 1e-14);
}

TEST(MathExtraRotation, BuildRyMatrix)
{
    double R[3][3];
    MathExtra::BuildRyMatrix(R, 0.0);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(R[i][j], expected, 1e-15);
        }
}

TEST(MathExtraRotation, BuildRyMatrix_preserves_y)
{
    double R[3][3];
    MathExtra::BuildRyMatrix(R, 0.5);
    double v[3] = {0.0, 1.0, 0.0};
    double ans[3];
    MathExtra::matvec(R, v, ans);
    EXPECT_NEAR(ans[0], 0.0, 1e-14);
    EXPECT_NEAR(ans[1], 1.0, 1e-14);
    EXPECT_NEAR(ans[2], 0.0, 1e-14);
}

TEST(MathExtraRotation, BuildRzMatrix)
{
    double R[3][3];
    MathExtra::BuildRzMatrix(R, 0.0);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(R[i][j], expected, 1e-15);
        }
}

TEST(MathExtraRotation, BuildRzMatrix_preserves_z)
{
    double R[3][3];
    MathExtra::BuildRzMatrix(R, 0.5);
    double v[3] = {0.0, 0.0, 1.0};
    double ans[3];
    MathExtra::matvec(R, v, ans);
    EXPECT_NEAR(ans[0], 0.0, 1e-14);
    EXPECT_NEAR(ans[1], 0.0, 1e-14);
    EXPECT_NEAR(ans[2], 1.0, 1e-14);
}

TEST(MathExtraRotation, rotation_generator_x)
{
    // For identity matrix, Gx*I should give the x-rotation generator
    double I[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double ans[3][3];
    MathExtra::rotation_generator_x(I, ans);
    // Generator Lx = [[0,0,0],[0,0,1],[0,-1,0]]
    EXPECT_DOUBLE_EQ(ans[0][0], 0.0);
    EXPECT_DOUBLE_EQ(ans[0][1], 0.0);
    EXPECT_DOUBLE_EQ(ans[0][2], 0.0);
    EXPECT_DOUBLE_EQ(ans[1][0], 0.0);
    EXPECT_DOUBLE_EQ(ans[1][1], 0.0);
    EXPECT_DOUBLE_EQ(ans[1][2], 1.0);
    EXPECT_DOUBLE_EQ(ans[2][0], 0.0);
    EXPECT_DOUBLE_EQ(ans[2][1], -1.0);
    EXPECT_DOUBLE_EQ(ans[2][2], 0.0);
}

TEST(MathExtraRotation, rotation_generator_y)
{
    double I[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double ans[3][3];
    MathExtra::rotation_generator_y(I, ans);
    // Generator Ly = [[0,0,-1],[0,0,0],[1,0,0]]
    EXPECT_DOUBLE_EQ(ans[0][0], 0.0);
    EXPECT_DOUBLE_EQ(ans[0][1], 0.0);
    EXPECT_DOUBLE_EQ(ans[0][2], -1.0);
    EXPECT_DOUBLE_EQ(ans[1][0], 0.0);
    EXPECT_DOUBLE_EQ(ans[1][1], 0.0);
    EXPECT_DOUBLE_EQ(ans[1][2], 0.0);
    EXPECT_DOUBLE_EQ(ans[2][0], 1.0);
    EXPECT_DOUBLE_EQ(ans[2][1], 0.0);
    EXPECT_DOUBLE_EQ(ans[2][2], 0.0);
}

TEST(MathExtraRotation, rotation_generator_z)
{
    double I[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double ans[3][3];
    MathExtra::rotation_generator_z(I, ans);
    // Generator Lz = [[0,1,0],[-1,0,0],[0,0,0]]
    EXPECT_DOUBLE_EQ(ans[0][0], 0.0);
    EXPECT_DOUBLE_EQ(ans[0][1], 1.0);
    EXPECT_DOUBLE_EQ(ans[0][2], 0.0);
    EXPECT_DOUBLE_EQ(ans[1][0], -1.0);
    EXPECT_DOUBLE_EQ(ans[1][1], 0.0);
    EXPECT_DOUBLE_EQ(ans[1][2], 0.0);
    EXPECT_DOUBLE_EQ(ans[2][0], 0.0);
    EXPECT_DOUBLE_EQ(ans[2][1], 0.0);
    EXPECT_DOUBLE_EQ(ans[2][2], 0.0);
}

// =========================================================================
// Volume and other utilities
// =========================================================================

TEST(MathExtraUtil, volume_ellipsoid_sphere)
{
    double shape[3] = {2.0, 2.0, 2.0};
    double vol      = MathExtra::volume_ellipsoid(shape);
    double expected = (4.0 / 3.0) * MY_PI * 8.0;
    EXPECT_NEAR(vol, expected, 1e-10);
}

TEST(MathExtraUtil, volume_ellipsoid_general)
{
    double shape[3] = {1.0, 2.0, 3.0};
    double vol      = MathExtra::volume_ellipsoid(shape);
    double expected = (4.0 / 3.0) * MY_PI * 6.0;
    EXPECT_NEAR(vol, expected, 1e-10);
}

TEST(MathExtraUtil, volume_ellipsoid_sphere_super)
{
    // With blockiness [2,2] and flag_super=0, same as regular ellipsoid
    double shape[3] = {2.0, 2.0, 2.0};
    double block[2] = {2.0, 2.0};
    double vol      = MathExtra::volume_ellipsoid(shape, block, 0);
    double expected = (4.0 / 3.0) * MY_PI * 8.0;
    EXPECT_NEAR(vol, expected, 1e-10);
}

TEST(MathExtraUtil, beta_function)
{
    // Beta(1,1) = 1
    EXPECT_NEAR(MathExtra::beta(1.0, 1.0), 1.0, 1e-14);

    // Beta(x,y) = Gamma(x)*Gamma(y)/Gamma(x+y)
    // Beta(2,3) = 1!*2! / 4! = 2/24 = 1/12
    EXPECT_NEAR(MathExtra::beta(2.0, 3.0), 1.0 / 12.0, 1e-14);

    // Beta(0.5, 0.5) = pi
    EXPECT_NEAR(MathExtra::beta(0.5, 0.5), MY_PI, 1e-12);
}

TEST(MathExtraUtil, multiply_shape_shape)
{
    // Voigt ordering: [0]=xx, [1]=yy, [2]=zz, [3]=yz, [4]=xz, [5]=xy
    // Identity in Voigt: [1,1,1,0,0,0]
    double one[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};
    double two[6] = {2.0, 3.0, 4.0, 0.5, 0.6, 0.7};
    double ans[6];
    MathExtra::multiply_shape_shape(one, two, ans);
    // Identity * X = X
    EXPECT_DOUBLE_EQ(ans[0], 2.0);
    EXPECT_DOUBLE_EQ(ans[1], 3.0);
    EXPECT_DOUBLE_EQ(ans[2], 4.0);
    EXPECT_DOUBLE_EQ(ans[3], 0.5);
    EXPECT_DOUBLE_EQ(ans[4], 0.6);
    EXPECT_DOUBLE_EQ(ans[5], 0.7);
}
