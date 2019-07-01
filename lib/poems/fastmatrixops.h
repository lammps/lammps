/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: fastmatrixops.h                                         *
 *      AUTHORS: See Author List                                           * 
 *      GRANTS: See Grants List                                            *
 *      COPYRIGHT: (C) 2005 by Authors as listed in Author's List          *
 *      LICENSE: Please see License Agreement                              *
 *      DOWNLOAD: Free at www.rpi.edu/~anderk5                             *
 *      ADMINISTRATOR: Prof. Kurt Anderson                                 *
 *                     Computational Dynamics Lab                          *
 *                     Rensselaer Polytechnic Institute                    *
 *                     110 8th St. Troy NY 12180                           * 
 *      CONTACT:        anderk5@rpi.edu                                    *
 *_________________________________________________________________________*/

#ifndef FASTMATRIXOPS_H
#define FASTMATRIXOPS_H

class ColMatMap;
class ColMatrix;
class Mat3x3;
class Mat4x4;
class Mat6x6;
class Matrix;
class Vect3;
class Vect4;
class Vect6;

void FastCross(Vect3& a, Vect3& b, Vect3& c);
void FastSimpleRotation(Vect3& v, double q, Mat3x3& C);
void FastQuaternions(ColMatrix& q, Mat3x3& C);
void FastQuaternionDerivatives(ColMatrix& q, ColMatrix& omega, ColMatrix& qdot);
void FastLDLT(Matrix& A, Matrix& LD); // C is the LDL^T decomposition of A (SPD)
void FastLDLT(Mat6x6& A, Mat6x6& LD); // C is the LDL^T decomposition of A (SPD)
void FastLDLTSubs(Matrix& LD, Matrix& B, Matrix& C); // Appropriate Forward and Back Substitution
void FastLDLTSubsLH(Matrix& B, Matrix& LD, Matrix& C); // Left handed Forward and Back Substitution
void FastLDLTSubs(Mat6x6& LD, Mat6x6& B, Mat6x6& C); // Appropriate Forward and Back Substitution
void FastLDLTSubs(Mat6x6& LD, Vect6& B, Vect6& C); // Appropriate Forward and Back Substitution
void FastLU(Matrix& A, Matrix& LU, int *indx); // LU is the LU decomposition of A
void FastLU(Mat3x3& A, Mat3x3& LU, int *indx); // LU is the LU decomposition of A
void FastLU(Mat4x4& A, Mat4x4& LU, int *indx); // LU is the LU decomposition of A
void FastLU(Mat6x6& A, Mat6x6& LU, int *indx); // LU is the LU decomposition of A
void FastLUSubs(Matrix& LU, Matrix& B, Matrix& C, int *indx); // Appropriate Forward and Back Substitution
void FastLUSubs(Mat3x3& LU, Matrix& B, Matrix& C, int *indx); // Appropriate Forward and Back Substitution
void FastLUSubs(Mat4x4& LU, Matrix& B, Matrix& C, int *indx); // Appropriate Forward and Back Substitution
void FastLUSubs(Mat6x6& LU, Matrix& B, Matrix& C, int *indx); // Appropriate Forward and Back Substitution
// The following LUSubsLH routine is incomplete at the moment.
void FastLUSubsLH(Matrix& B, Matrix& LU, Matrix& C, int *indx); // Appropriate Forward and Back Subsitution
void FastTripleSum(Vect3& a, Vect3& b, Vect3& c, Vect3& d); // d = a+b+c
void FastTripleSumPPM(Vect3& a, Vect3& b, Vect3& c, Vect3& d); // d = a+b-c

void FastMult(Matrix& A, Matrix& B, Matrix& C);  // C = A*B
void FastTMult(Matrix& A, Matrix& B, Matrix& C);  // C = A^T*B

void FastMult(Mat3x3& A, Vect3& B, Vect3& C);  // C = A*B
void FastMult(Mat3x3& A, ColMatrix& B, Vect3& C); // C = A*B
void FastMult(Mat3x3& A, Vect3& B, ColMatrix& C); // C = A*B
void FastMult(double a, Vect3& B, Vect3& C);  // C = a*B
void FastNegMult(Mat3x3& A, Vect3& B, Vect3& C);  // C = A*B
void FastNegTMult(Mat3x3& A, Vect3& B, Vect3& C);  // C = A*B

void FastMult(Mat4x4& A, Vect4& B, Vect4& C);  // C = A*B
void FastTMult(Mat4x4& A, Vect4& B, Vect4& C);  // C = A^T*B
void FastMult(double a, Vect4& B, Vect4& C);  // C = a*B
void FastNegMult(Mat4x4& A, Vect4& B, Vect4& C);  // C = A*B
void FastNegTMult(Mat4x4& A, Vect4& B, Vect4& C);  // C = A*B

void FastMultT(Matrix& A, Matrix& B, Mat6x6& C);  // C = A*B^T
void FastMult(Mat6x6& A, Matrix& B, Matrix& C);   // C = A*B
void FastTMult(Matrix& A, Vect6& B, ColMatrix& C);// C = A^T*B
void FastMult(Matrix& A, ColMatrix& B, Vect6& C); // C = A*B

void FastMult(Mat3x3& A, Mat3x3& B, Mat3x3& C);  // C = A*B
void FastMultT(Mat3x3& A, Mat3x3& B, Mat3x3& C);  // C = A*B

void FastMult(Mat4x4& A, Mat4x4& B, Mat4x4& C);  // C = A*B
void FastMultT(Mat4x4& A, Mat4x4& B, Mat4x4& C);  // C = A*B

void FastMult(Mat6x6& A, Mat6x6& B, Mat6x6& C);  // C = A*B
void FastMultT(Mat6x6& A, Mat6x6& B, Mat6x6& C);  // C = A*B^T
void FastTMult(Mat6x6& A, Mat6x6& B, Mat6x6& C);  // C = A^T*B

void FastMult(Mat6x6& A, Vect6& B, Vect6& C);  // C = A*B
void FastTMult(Mat6x6& A, Vect6& B, Vect6& C);  // C = A^T*B

void FastAdd(Vect3& A, Vect3& B, Vect3& C);     // C = A+B
void FastAdd(Vect4& A, Vect4& B, Vect3& C);     // C = A+B
void FastAdd(Mat6x6& A, Mat6x6& B, Mat6x6& C);  // C = A+B
void FastAdd(Vect6& A, Vect6& B, Vect6& C);     // C = A+B

void FastSubt(Vect3& A, Vect3& B, Vect3& C);    // C = A-B
void FastSubt(Vect4& A, Vect4& B, Vect4& C);    // C = A-B
void FastSubt(Mat6x6& A, Mat6x6& B, Mat6x6& C); // C = A-B
void FastSubt(Vect6& A, Vect6& B, Vect6& C);    // C = A-B

void FastAssign(ColMatMap& A, ColMatMap& C);  // C = A
void FastAssign(ColMatrix& A, ColMatrix& C);  // C = A
void FastAssign(Vect3& A, Vect3& C);          // C = A
void FastAssign(ColMatrix&A, Vect3& C);
void FastAssign(Vect4& A, Vect4& C);          // C = A
void FastAssignT(Mat3x3& A, Mat3x3& C);       // C = A^T
void FastAssignT(Mat4x4& A, Mat4x4& C);       // C = A^T

#endif
