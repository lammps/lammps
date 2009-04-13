/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: matrixfun.cpp                                           *
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

#include "matrixfun.h"
#include <math.h>
#include "fastmatrixops.h"
#include <cstdlib>

using namespace std;

//
//  Create a new matrix
//

VirtualMatrix* NewMatrix(int type){
  switch( MatrixType(type) )
    {
      case MATRIX : return new Matrix;
      case COLMATRIX : return new ColMatrix;
      case ROWMATRIX : return new RowMatrix;		
      case MAT3X3 : return new Mat3x3;
		case MAT4X4 : return new Mat4x4;
      case VECT3 : return new Vect3;
		case VECT4 : return new Vect4;
      default  : return 0; // error!
    }
}

//
// Transpose
//

Matrix T(const VirtualMatrix& A){
	int numrows = A.GetNumRows();
	int numcols = A.GetNumCols();
	Matrix C(numcols,numrows);
	for(int i=0;i<numcols;i++)
		for(int j=0;j<numrows;j++)
			C.BasicSet(i,j,A.BasicGet(j,i));
	return C;
}

Mat3x3 T(const Mat3x3& A){
	Mat3x3 C;
	C.elements[0][0] = A.elements[0][0];
	C.elements[1][1] = A.elements[1][1];
	C.elements[2][2] = A.elements[2][2];

	C.elements[0][1] = A.elements[1][0];
	C.elements[0][2] = A.elements[2][0];
	C.elements[1][2] = A.elements[2][1];

	C.elements[1][0] = A.elements[0][1];
	C.elements[2][0] = A.elements[0][2];
	C.elements[2][1] = A.elements[1][2];
	return C;
}

Mat6x6 T(const Mat6x6& A){
	Mat6x6 C;
	int i,j;
	for(i=0;i<6;i++)
		for(j=0;j<6;j++)
			C.elements[i][j] = A.elements[j][i];

	return C;
}

Matrix T(const Vect3& A){
	Matrix C(1,3);
	C.BasicSet(0,0,A.elements[0]);
	C.BasicSet(0,1,A.elements[1]);
	C.BasicSet(0,2,A.elements[2]);

	return C;
}

Matrix T(const Vect6& A){
	Matrix C(1,6);
	C.BasicSet(0,0,A.elements[0]);
	C.BasicSet(0,1,A.elements[1]);
	C.BasicSet(0,2,A.elements[2]);
	C.BasicSet(0,3,A.elements[3]);
	C.BasicSet(0,4,A.elements[4]);
	C.BasicSet(0,5,A.elements[5]);

	return C;
}

RowMatrix T(const VirtualColMatrix &A){
	int numele = A.GetNumRows();
	RowMatrix C(numele);
	for(int i=0;i<numele;i++)
		C.BasicSet(i,A.BasicGet(i));
	return C;
}

ColMatrix T(const VirtualRowMatrix &A){
	int numele = A.GetNumCols();
	ColMatrix C(numele);
	for(int i=0;i<numele;i++)
		C.BasicSet(i,A.BasicGet(i));
	return C;
}

//
// Symmetric Inverse
//

Matrix SymInverse(Matrix &A){
	int r = A.GetNumRows();
	Matrix C(r,r);
	Matrix LD(r,r);
	Matrix I(r,r);
	I.Zeros();
	for(int i=0;i<r;i++) I.BasicSet(i,i,1.0);
	FastLDLT(A,LD);
	FastLDLTSubs(LD,I,C);
	return C;
}

Mat6x6 SymInverse(Mat6x6 &A){
	Mat6x6 C;
	Mat6x6 LD;
	Mat6x6 I;
	I.Zeros();
	for(int i=0;i<6;i++) I.BasicSet(i,i,1.0);
	FastLDLT(A,LD);
	FastLDLTSubs(LD,I,C);
	return C;
}

//
// Inverse
//

Matrix Inverse(Matrix& A){
	int r = A.GetNumRows();
	int indx[10000];
	Matrix LU(r,r);
	Matrix I(r,r);
	Matrix C(r,r);
	I.Zeros();
	for(int i=0;i<r;i++) I.BasicSet(i,i,1.0);
	FastLU(A,LU,indx);
	FastLUSubs(LU,I,C,indx);
	return C;
}

Mat3x3 Inverse(Mat3x3& A){
	int indx[10000];
	Mat3x3 LU;
	Matrix I(3,3);
	Matrix C(3,3);
	I.Zeros();
	for(int i=0;i<3;i++) I.BasicSet(i,i,1.0);
	FastLU(A,LU,indx);
	FastLUSubs(LU,I,C,indx);
	return C;
}

Mat4x4 Inverse(Mat4x4& A){
	int indx[10000];
	Mat4x4 LU;
	Matrix I(4,4);
	Matrix C(4,4);
	I.Zeros();
	for(int i=0;i<4;i++) I.BasicSet(i,i,1.0);
	FastLU(A,LU,indx);
	FastLUSubs(LU,I,C,indx);
	return C;
}

Mat6x6 Inverse(Mat6x6& A){
	int indx[10000];
	Mat6x6 LU;
	Matrix I(6,6);
	Matrix C(6,6);
	I.Zeros();
	for(int i=0;i<6;i++) I.BasicSet(i,i,1.0);
	FastLU(A,LU,indx);
	FastLUSubs(LU,I,C,indx);
	return C;
}

//
// overloaded addition
//

Matrix operator+ (const VirtualMatrix &A, const VirtualMatrix &B){      // addition
	int Arows,Acols,Brows,Bcols;
	Arows = A.GetNumRows();
	Acols = A.GetNumCols();
	Brows = B.GetNumRows();
	Bcols = B.GetNumCols();

	if( !((Arows == Brows) && (Acols == Bcols)) ){
		cerr << "Dimesion mismatch in matrix addition" << endl;
		exit(1);
	}

	Matrix C(Arows,Acols);

	for(int i=0;i<Arows;i++)
		for(int j=0;j<Acols;j++)
			C.BasicSet(i,j,A.BasicGet(i,j) + B.BasicGet(i,j));

	return C;
}

//
// overloaded subtraction
//

Matrix operator- (const VirtualMatrix &A, const VirtualMatrix &B){      // subtraction
	int Arows,Acols,Brows,Bcols;
	Arows = A.GetNumRows();
	Acols = A.GetNumCols();
	Brows = B.GetNumRows();
	Bcols = B.GetNumCols();

	if( !((Arows == Brows) && (Acols == Bcols)) ){
		cerr << "Dimesion mismatch in matrix addition" << endl;
		exit(1);
	}

	Matrix C(Arows,Acols);

	for(int i=0;i<Arows;i++)
		for(int j=0;j<Acols;j++)
			C.BasicSet(i,j,A.BasicGet(i,j) - B.BasicGet(i,j));

	return C;
}

//
// overloaded matrix multiplication
//

Matrix operator* (const VirtualMatrix &A, const VirtualMatrix &B){      // multiplication
	int Arows,Acols,Brows,Bcols;
	Arows = A.GetNumRows();
	Acols = A.GetNumCols();
	Brows = B.GetNumRows();
	Bcols = B.GetNumCols();

	if(Acols != Brows){
		cerr << "Dimesion mismatch in matrix multiplication" << endl;
		exit(1);
	}

	Matrix C(Arows,Bcols);
	C.Zeros();

	for(int i=0;i<Arows;i++)
		for(int j=0;j<Bcols;j++)
			for(int k=0;k<Brows;k++)
				C.BasicIncrement(i,j, A.BasicGet(i,k) * B.BasicGet(k,j) );

	return C;
}

//
// overloaded scalar multiplication
//

Matrix operator* (const VirtualMatrix &A, double b){    // multiplication
	Matrix C = A;
	C *= b;
	return C;
}

Matrix operator* (double b, const VirtualMatrix &A){    // multiplication
	Matrix C = A;  
	C *= b;
	return C;
}

//
// overloaded negative
//

Matrix operator- (const VirtualMatrix& A){  // negative
	int r = A.GetNumRows();
	int c = A.GetNumCols();
	Matrix C(r,c);

	for(int i=0;i<r;i++)
		for(int j=0;j<c;j++)
			C.BasicSet(i,j,-A.BasicGet(i,j));

	return C;
}

//
// Cross product (friend of Vect3)
//

Vect3 Cross(Vect3& a, Vect3& b){
	return CrossMat(a)*b;
}

//
// Cross Matrix (friend of Vect3 & Mat3x3)
//

Mat3x3 CrossMat(Vect3& a){
	Mat3x3 C;
	C.Zeros();
	C.elements[0][1] = -a.elements[2];
	C.elements[1][0] = a.elements[2];
	C.elements[0][2] = a.elements[1];
	C.elements[2][0] = -a.elements[1];
	C.elements[1][2] = -a.elements[0];
	C.elements[2][1] = a.elements[0];

	return C;
}

//
// Stack
//

Matrix Stack(VirtualMatrix& A, VirtualMatrix& B){
	int m,na,nb;
	m = A.GetNumCols();
	if( m != B.GetNumCols()){
		cerr << "Error: cannot stack matrices of differing column dimension" << endl;
		exit(0);
	}
	na = A.GetNumRows();
	nb = B.GetNumRows();

	Matrix C(na+nb,m);

	for(int i=0;i<na;i++){
		for(int j=0;j<m;j++){
			C.BasicSet(i,j,A.BasicGet(i,j));
		}
	}

	for(int i=0;i<nb;i++){
		for(int j=0;j<m;j++){
			C.BasicSet(i+na,j,B.BasicGet(i,j));
		}
	}

	return C;
}

//Hstack

Matrix HStack(VirtualMatrix& A, VirtualMatrix& B){
	int m,na,nb;
	m = A.GetNumRows();
	if( m != B.GetNumRows()){
		cerr << "Error: cannot stack matrices of differing row dimension" << endl;
		exit(0);
	}
	na = A.GetNumCols();
	nb = B.GetNumCols();

	Matrix C(m,na+nb);

	for(int i=0;i<m;i++){
		for(int j=0;j<na;j++){
			C.BasicSet(i,j,A.BasicGet(i,j));
		}
	}

	for(int i=0;i<m;i++){
		for(int j=0;j<nb;j++){
			C.BasicSet(i,j+na,B.BasicGet(i,j));
		}
	}

	return C;
}

//
//

void Set6DAngularVector(Vect6& v6, Vect3& v3){
	v6.elements[0] = v3.elements[0];
	v6.elements[1] = v3.elements[1];
	v6.elements[2] = v3.elements[2];
}

void Set6DLinearVector(Vect6& v6, Vect3& v3){
	v6.elements[3] = v3.elements[0];
	v6.elements[4] = v3.elements[1];
	v6.elements[5] = v3.elements[2];
}

