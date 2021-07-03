// clang-format off
/*******************************************************************************
 *  -/_|:|_|_\-
 *
 *  This code is a modification of D.L. Theobald's QCP rotation code.
 *  It has been adapted to calculate the polar decomposition of a 3x3 matrix
 *  Adaption by P.M. Larsen
 *
 *  Original Author(s):          Douglas L. Theobald
 *                                  Department of Biochemistry
 *                                  MS 009
 *                                  Brandeis University
 *                                  415 South St
 *                                  Waltham, MA  02453
 *                                  USA
 *
 *                                  dtheobald@brandeis.edu
 *
 *                                  Pu Liu
 *                                  Johnson & Johnson Pharmaceutical Research and Development, L.L.C.
 *                                  665 Stockton Drive
 *                                  Exton, PA  19341
 *                                  USA
 *
 *                                  pliu24@its.jnj.com
 *
 *
 *        If you use this QCP rotation calculation method in a publication, please
 *        reference:
 *
 *          Douglas L. Theobald (2005)
 *          "Rapid calculation of RMSD using a quaternion-based characteristic
 *          polynomial."
 *          Acta Crystallographica A 61(4):478-480.
 *
 *          Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2009)
 *          "Fast determination of the optimal rotational matrix for macromolecular
 *          superpositions."
 *          Journal of Computational Chemistry 31(7):1561-1563.
 *
 *
 *  Copyright (c) 2009-2013 Pu Liu and Douglas L. Theobald
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without modification, are permitted
 *  provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice, this list
 *        of conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *  * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to
 *        endorse or promote products derived from this software without specific prior written
 *        permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  Source:                 started anew.
 *
 *  Change History:
 *        2009/04/13          Started source
 *        2010/03/28          Modified FastCalcRMSDAndRotation() to handle tiny qsqr
 *                                        If trying all rows of the adjoint still gives too small
 *                                        qsqr, then just return identity matrix. (DLT)
 *        2010/06/30          Fixed prob in assigning A[9] = 0 in InnerProduct()
 *                                        invalid mem access
 *        2011/02/21          Made CenterCoords use weights
 *        2011/05/02          Finally changed CenterCoords declaration in qcprot.h
 *                                        Also changed some functions to static
 *        2011/07/08          put in fabs() to fix taking sqrt of small neg numbers, fp error
 *        2012/07/26          minor changes to comments and main.c, more info (v.1.4)
 *
 *      2016/05/29        QCP method adapted for polar decomposition of a 3x3 matrix,
 *                          for use in Polyhedral Template Matching.
 *
 ******************************************************************************/

#include <cmath>
#include <algorithm>
#include <cstring>
#include "ptm_polar.h"
#include "ptm_quat.h"


namespace ptm {

static void matmul_3x3(double* A, double* x, double* b)
{
        b[0] = A[0] * x[0] + A[1] * x[3] + A[2] * x[6];
        b[3] = A[3] * x[0] + A[4] * x[3] + A[5] * x[6];
        b[6] = A[6] * x[0] + A[7] * x[3] + A[8] * x[6];

        b[1] = A[0] * x[1] + A[1] * x[4] + A[2] * x[7];
        b[4] = A[3] * x[1] + A[4] * x[4] + A[5] * x[7];
        b[7] = A[6] * x[1] + A[7] * x[4] + A[8] * x[7];

        b[2] = A[0] * x[2] + A[1] * x[5] + A[2] * x[8];
        b[5] = A[3] * x[2] + A[4] * x[5] + A[5] * x[8];
        b[8] = A[6] * x[2] + A[7] * x[5] + A[8] * x[8];
}

static double matrix_determinant_3x3(double* A)
{
        return    A[0] * (A[4]*A[8] - A[5]*A[7])
                - A[1] * (A[3]*A[8] - A[5]*A[6])
                + A[2] * (A[3]*A[7] - A[4]*A[6]);
}

static void flip_matrix(double* A)
{
        for (int i=0;i<9;i++)
                A[i] = -A[i];
}

static bool optimal_quaternion(double* A, bool polar, double E0, double* p_nrmsdsq, double* qopt)
{
        const double evecprec = 1e-6;
        const double evalprec = 1e-11;

        double        Sxx = A[0], Sxy = A[1], Sxz = A[2],
                Syx = A[3], Syy = A[4], Syz = A[5],
                Szx = A[6], Szy = A[7], Szz = A[8];

        double        Sxx2 = Sxx * Sxx, Syy2 = Syy * Syy, Szz2 = Szz * Szz,
                Sxy2 = Sxy * Sxy, Syz2 = Syz * Syz, Sxz2 = Sxz * Sxz,
                Syx2 = Syx * Syx, Szy2 = Szy * Szy, Szx2 = Szx * Szx;

        double fnorm_squared = Sxx2 + Syy2 + Szz2 + Sxy2 + Syz2 + Sxz2 + Syx2 + Szy2 + Szx2;

        double SyzSzymSyySzz2 = 2.0 * (Syz * Szy - Syy * Szz);
        double Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;
        double SxzpSzx = Sxz + Szx;
        double SyzpSzy = Syz + Szy;
        double SxypSyx = Sxy + Syx;
        double SyzmSzy = Syz - Szy;
        double SxzmSzx = Sxz - Szx;
        double SxymSyx = Sxy - Syx;
        double SxxpSyy = Sxx + Syy;
        double SxxmSyy = Sxx - Syy;
        double Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

        double C[3];
        C[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
                 + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
                 + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
                 + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
                 + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
                 + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));

        C[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);
        C[2] = -2.0 * fnorm_squared;

        //Newton-Raphson
        double mxEigenV = polar ? sqrt(3 * fnorm_squared) : E0;
        if (mxEigenV > evalprec)
        {
                for (int i=0;i<50;i++)
                {
                        double oldg = mxEigenV;
                        double x2 = mxEigenV*mxEigenV;
                        double b = (x2 + C[2])*mxEigenV;
                        double a = b + C[1];
                        double delta = ((a * mxEigenV + C[0]) / (2 * x2 * mxEigenV + b + a));
                        mxEigenV -= delta;
                        if (fabs(mxEigenV - oldg) < fabs(evalprec * mxEigenV))
                                break;
                }
        }
        else
        {
                mxEigenV = 0.0;
        }

        (*p_nrmsdsq) = std::max(0.0, 2.0 * (E0 - mxEigenV));

        double a11 = SxxpSyy + Szz - mxEigenV;
        double a12 = SyzmSzy;
        double a13 = -SxzmSzx;
        double a14 = SxymSyx;

        double a21 = SyzmSzy;
        double a22 = SxxmSyy - Szz  -mxEigenV;
        double a23 = SxypSyx;
        double a24 = SxzpSzx;

        double a31 = a13;
        double a32 = a23;
        double a33 = Syy - Sxx - Szz - mxEigenV;
        double a34 = SyzpSzy;

        double a41 = a14;
        double a42 = a24;
        double a43 = a34;
        double a44 = Szz - SxxpSyy - mxEigenV;

        double a3344_4334 = a33 * a44 - a43 * a34;
        double a3244_4234 = a32 * a44 - a42 * a34;
        double a3243_4233 = a32 * a43 - a42 * a33;
        double a3143_4133 = a31 * a43 - a41 * a33;
        double a3144_4134 = a31 * a44 - a41 * a34;
        double a3142_4132 = a31 * a42 - a41 * a32;
        double a1324_1423 = a13 * a24 - a14 * a23;
        double a1224_1422 = a12 * a24 - a14 * a22;
        double a1223_1322 = a12 * a23 - a13 * a22;
        double a1124_1421 = a11 * a24 - a14 * a21;
        double a1123_1321 = a11 * a23 - a13 * a21;
        double a1122_1221 = a11 * a22 - a12 * a21;

        double q[4][4];
        q[0][0] =  a12 * a3344_4334 - a13 * a3244_4234 + a14 * a3243_4233;
        q[0][1] = -a11 * a3344_4334 + a13 * a3144_4134 - a14 * a3143_4133;
        q[0][2] =  a11 * a3244_4234 - a12 * a3144_4134 + a14 * a3142_4132;
        q[0][3] = -a11 * a3243_4233 + a12 * a3143_4133 - a13 * a3142_4132;

        q[1][0] =  a22 * a3344_4334 - a23 * a3244_4234 + a24 * a3243_4233;
        q[1][1] = -a21 * a3344_4334 + a23 * a3144_4134 - a24 * a3143_4133;
        q[1][2] =  a21 * a3244_4234 - a22 * a3144_4134 + a24 * a3142_4132;
        q[1][3] = -a21 * a3243_4233 + a22 * a3143_4133 - a23 * a3142_4132;

        q[2][0] =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
        q[2][1] = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
        q[2][2] =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
        q[2][3] = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;

        q[3][0] =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
        q[3][1] = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
        q[3][2] =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
        q[3][3] = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;

        double qsqr[4];
        for (int i=0;i<4;i++)
                qsqr[i] = q[i][0]*q[i][0] + q[i][1]*q[i][1] + q[i][2]*q[i][2] + q[i][3]*q[i][3];

        int bi = 0;
        double max = 0;
        for (int i=0;i<4;i++)
        {
                if (qsqr[i] > max)
                {
                        bi = i;
                        max = qsqr[i];
                }
        }

        bool too_small = false;
        if (qsqr[bi] < evecprec)
        {
                //if qsqr is still too small, return the identity rotation.
                q[bi][0] = 1;
                q[bi][1] = 0;
                q[bi][2] = 0;
                q[bi][3] = 0;
                too_small = true;
        }
        else
        {
                double normq = sqrt(qsqr[bi]);
                q[bi][0] /= normq;
                q[bi][1] /= normq;
                q[bi][2] /= normq;
                q[bi][3] /= normq;
        }

        memcpy(qopt, q[bi], 4 * sizeof(double));
        return !too_small;
}

int polar_decomposition_3x3(double* _A, bool right_sided, double* U, double* P)
{
        double A[9];
        memcpy(A, _A, 9 * sizeof(double));

        double det = matrix_determinant_3x3(A);
        if (det < 0)
                flip_matrix(A);

        double q[4];
        double nrmsdsq = 0;
        optimal_quaternion(A, true, -1, &nrmsdsq, q);
        q[0] = -q[0];
        quaternion_to_rotation_matrix(q, U);

        if (det < 0)
                flip_matrix(U);

        double UT[9] = {U[0], U[3], U[6], U[1], U[4], U[7], U[2], U[5], U[8]};

        if (right_sided)
                matmul_3x3(UT, _A, P);
        else
                matmul_3x3(_A, UT, P);

        return 0;
}

void InnerProduct(double *A, int num, const double (*coords1)[3], double (*coords2)[3], int8_t* permutation)
{
        A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = A[6] = A[7] = A[8] = 0.0;

        for (int i = 0; i < num; ++i)
        {
                double x1 = coords1[i][0];
                double y1 = coords1[i][1];
                double z1 = coords1[i][2];

                double x2 = coords2[permutation[i]][0];
                double y2 = coords2[permutation[i]][1];
                double z2 = coords2[permutation[i]][2];

                A[0] += x1 * x2;
                A[1] += x1 * y2;
                A[2] += x1 * z2;

                A[3] += y1 * x2;
                A[4] += y1 * y2;
                A[5] += y1 * z2;

                A[6] += z1 * x2;
                A[7] += z1 * y2;
                A[8] += z1 * z2;
        }
}

int FastCalcRMSDAndRotation(double *A, double E0, double *p_nrmsdsq, double *q, double* U)
{
        optimal_quaternion(A, false, E0, p_nrmsdsq, q);
        quaternion_to_rotation_matrix(q, U);
        return 0;
}

}

