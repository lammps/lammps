/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: matrixfun.h                                           *
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

#ifndef MATRIXFUN_H
#define MATRIXFUN_H

#include "colmatrix.h"
#include "mat3x3.h"
#include "mat4x4.h"
#include "mat6x6.h"
#include "matrix.h"
#include "rowmatrix.h"
#include "vect3.h"
#include "vect6.h"

class VirtualColMatrix;
class VirtualMatrix;
class VirtualRowMatrix;

// Create a Matrix
VirtualMatrix* NewMatrix(int type);

// Transpose
Matrix T(const VirtualMatrix& A);
Mat3x3 T(const Mat3x3& A);
Mat6x6 T(const Mat6x6& A);
Matrix T(const Vect3& A);
Matrix T(const Vect6& A);
RowMatrix T(const VirtualColMatrix& A);
ColMatrix T(const VirtualRowMatrix& A);

// Symmetric Inverse

Matrix SymInverse(Matrix& A);
Mat6x6 SymInverse(Mat6x6& A);

// Inverse
Matrix Inverse(Matrix& A);
Mat3x3 Inverse(Mat3x3& A);
Mat4x4 Inverse(Mat4x4& A);
Mat6x6 Inverse(Mat6x6& A);

// overloaded addition
Matrix operator+ (const VirtualMatrix &A, const VirtualMatrix &B);	// addition
//Mat3x3 operator+ (const Mat3x3 &A, const Mat3x3 &B);	// addition
//Matrix operator+ (const VirtualMatrix &A, const VirtualMatrix &B);	// addition

// overloaded subtraction
Matrix operator- (const VirtualMatrix &A, const VirtualMatrix &B);	// subtraction

// overloaded matrix multiplication
Matrix operator* (const VirtualMatrix &A, const VirtualMatrix &B);  // multiplication

// overloaded scalar-matrix multiplication
Matrix operator* (const VirtualMatrix &A, double b);			// overloaded *
Matrix operator* (double b, const VirtualMatrix &A);			// overloaded *

// overloaded negative
Matrix operator- (const VirtualMatrix &A);	// negative

Vect3 Cross(Vect3& a, Vect3& b);
Mat3x3 CrossMat(Vect3& a);

Matrix Stack(VirtualMatrix& A, VirtualMatrix& B);
Matrix HStack(VirtualMatrix& A, VirtualMatrix& B);

void Set6DAngularVector(Vect6& v6, Vect3& v3);
void Set6DLinearVector(Vect6& v6, Vect3& v3);

#endif
