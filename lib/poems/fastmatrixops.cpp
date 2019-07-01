/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: fastmatrixops.cpp                                       *
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

#include <cmath>

#include "fastmatrixops.h"
#include "colmatmap.h"
#include "colmatrix.h"
#include "mat3x3.h"
#include "mat4x4.h"
#include "mat6x6.h"
#include "matrix.h"
#include "vect3.h"
#include "vect4.h"
#include "vect6.h"

using namespace std;

//
// Cross Product (friend of Vect3)
//

void FastCross(Vect3& a, Vect3& b, Vect3& c){
  c.elements[0] = a.elements[1]*b.elements[2] - a.elements[2]*b.elements[1];
  c.elements[1] = a.elements[2]*b.elements[0] - a.elements[0]*b.elements[2];
  c.elements[2] = a.elements[0]*b.elements[1] - a.elements[1]*b.elements[0];
}

//
//  Simple Rotation (friend of Vect3 and Mat3x3)
//

void FastSimpleRotation(Vect3& v, double q, Mat3x3& C){
  // intermediate quantities
  double cq = cos(q);
  double sq = sin(q);
  double one_m_cq = 1-cq;
  double l12 = v.elements[0]*v.elements[1]*one_m_cq;
  double l13 = v.elements[0]*v.elements[2]*one_m_cq;
  double l23 = v.elements[1]*v.elements[2]*one_m_cq;

  // the transformation
  C.elements[0][0] = v.elements[0]*v.elements[0]*one_m_cq+cq;
  C.elements[0][1] = l12-v.elements[2]*sq;
  C.elements[0][2] = l13+v.elements[1]*sq;
  C.elements[1][0] = l12+v.elements[2]*sq;
  C.elements[1][1] = v.elements[1]*v.elements[1]*one_m_cq+cq;
  C.elements[1][2] = l23-v.elements[0]*sq;
  C.elements[2][0] = l13-v.elements[1]*sq;
  C.elements[2][1] = l23+v.elements[0]*sq;
  C.elements[2][2] = v.elements[2]*v.elements[2]*one_m_cq+cq;
}

//
// Quaternion Functions
//

void FastQuaternions(ColMatrix& q, Mat3x3& C){
  double* e = q.elements;

  // normalize the quaternions
  double length = e[0]*e[0] + e[1]*e[1] + e[2]*e[2] + e[3]*e[3];
  length = sqrt(length);
  e[0] = e[0]/length;
  e[1] = e[1]/length;
  e[2] = e[2]/length;
  e[3] = e[3]/length;

  // compute the transformation
  C.elements[0][0] = e[0]*e[0] + e[1]*e[1] - e[2]*e[2] - e[3]*e[3];
  C.elements[1][1] = e[0]*e[0] - e[1]*e[1] + e[2]*e[2] - e[3]*e[3];
  C.elements[2][2] = e[0]*e[0] - e[1]*e[1] - e[2]*e[2] + e[3]*e[3];

  C.elements[0][1] = 2 * (e[1]*e[2] - e[0]*e[3]);
  C.elements[0][2] = 2 * (e[1]*e[3] + e[0]*e[2]);
  C.elements[1][2] = 2 * (e[2]*e[3] - e[0]*e[1]);

  C.elements[1][0] = 2 * (e[1]*e[2] + e[0]*e[3]);
  C.elements[2][0] = 2 * (e[1]*e[3] - e[0]*e[2]);
  C.elements[2][1] = 2 * (e[2]*e[3] + e[0]*e[1]);
}

void FastQuaternionDerivatives(ColMatrix& q, ColMatrix& omega, ColMatrix& qdot){
  double* w = omega.elements;
  double* e = q.elements;

  qdot.elements[0] = 0.5 * (-w[0]*e[1] - w[1]*e[2] - w[2]*e[3]);
  qdot.elements[1] = 0.5 * ( w[0]*e[0] + w[2]*e[2] - w[1]*e[3]);
  qdot.elements[2] = 0.5 * ( w[1]*e[0] - w[2]*e[1] + w[0]*e[3]);
  qdot.elements[3] = 0.5 * ( w[2]*e[0] + w[1]*e[1] - w[0]*e[2]);
}

void FastInvQuaternions(Mat3x3& C, ColMatrix& q){
}

//
// Inverse
//

// friend of Matrix
//void FastInverse(Matrix& A, Matrix& C){ // C = A^(-1)
//  C.rows[0][0] = 1/A.rows[0][0];
//}

//
// LDL^T Decomposition (from Golub and Van Loan)
//

// friend of Matrix
void FastLDLT(Matrix& A, Matrix& C){ // C is the LD of the LDL^T decomposition of A (SPD)
  double Lv;
  int n = A.numrows;

  for(int j=0;j<n;j++){
    Lv = 0.0;
    for(int i=0;i<j;i++){
      C.rows[i][j] = C.rows[j][i]*C.rows[i][i];
      Lv += C.rows[j][i]*C.rows[i][j];
    }

    C.rows[j][j] = A.rows[j][j] - Lv;
    for(int i=j+1;i<n;i++){
      Lv = 0.0;
      for(int k=0;k<j;k++) Lv += C.rows[i][k]*C.rows[k][j];
      C.rows[i][j] = (A.rows[i][j] - Lv)/C.rows[j][j];
    }
  }
}


// friend of Mat6x6
void FastLDLT(Mat6x6& A, Mat6x6& C){ // C is the LD of the LDL^T decomposition of A (SPD)
  double v[6];
  double Lv;

  for(int j=0;j<6;j++){
    Lv = 0.0;
    for(int i=0;i<j;i++){
      v[i] = C.elements[j][i]*C.elements[i][i];
      Lv += C.elements[j][i]*v[i];
    }

    v[j] = A.elements[j][j] - Lv;
    C.elements[j][j] = v[j];
    for(int i=j+1;i<6;i++){
      Lv = 0.0;
      for(int k=0;k<j;k++) Lv += C.elements[i][k]*v[k];
      C.elements[i][j] = (A.elements[i][j] - Lv)/v[j];
    }
  }
}

// friend of Matrix
void FastLDLTSubs(Matrix& LD, Matrix& B, Matrix& C){
  int n = B.numrows;
  int c = B.numcols;
  double temp;

  for(int k=0;k<c;k++){
    for(int i=0;i<n;i++){
      temp = 0.0;
      for(int j=0;j<i;j++){
        temp += C.rows[j][k] * LD.rows[i][j];
      }
      C.rows[i][k] = B.rows[i][k] - temp;
    }
    for(int i=n-1;i>-1;i--){
      C.rows[i][k] = C.rows[i][k] / LD.rows[i][i];
      temp = 0.0;
      for(int j=n-1;j>i;j--){
        temp += C.rows[j][k] * LD.rows[j][i];
      }
      C.rows[i][k] = C.rows[i][k] - temp;
    }
  }
}

// friend of Matrix
void FastLDLTSubsLH(Matrix& B, Matrix& LD, Matrix& C){
  int n = B.numcols;
  int c = B.numrows;
  double temp;

  for(int k=0;k<c;k++){
    for(int i=0;i<n;i++){
      temp = 0.0;
      for(int j=0;j<i;j++){
        temp += C.rows[k][j] * LD.rows[i][j];
      }
      C.rows[k][i] = B.rows[k][i] - temp;
    }
    for(int i=n-1;i>-1;i--){
      C.rows[k][i] = C.rows[k][i] / LD.rows[i][i];
      temp = 0.0;
      for(int j=n-1;j>i;j--){
        temp += C.rows[k][j] * LD.rows[j][i];
      }
      C.rows[k][i] = C.rows[k][i] - temp;
    }
  }
}

// friend of Mat6x6
void FastLDLTSubs(Mat6x6& LD, Mat6x6& B, Mat6x6& C){
  double temp;

  for(int k=0;k<6;k++){
    for(int i=0;i<6;i++){
      temp = 0.0;
      for(int j=0;j<i;j++){
        temp += C.elements[j][k] * LD.elements[i][j];
      }
      C.elements[i][k] = B.elements[i][k] - temp;
    }
    for(int i=5;i>-1;i--){
      C.elements[i][k] = C.elements[i][k] / LD.elements[i][i];
      temp = 0.0;
      for(int j=5;j>i;j--){
        temp += C.elements[j][k] * LD.elements[j][i];
      }
      C.elements[i][k] = C.elements[i][k] - temp;
    }
  }
}

// friend of Mat6x6 & Vect6
void FastLDLTSubs(Mat6x6& LD, Vect6& B, Vect6& C){
  double temp;

  for(int i=0;i<6;i++){
    temp = 0.0;
    for(int j=0;j<i;j++){
      temp += C.elements[j] * LD.elements[i][j];
    }
    C.elements[i] = B.elements[i] - temp;
  }
  for(int i=5;i>-1;i--){
    C.elements[i] = C.elements[i] / LD.elements[i][i];
    temp = 0.0;
    for(int j=5;j>i;j--){
      temp += C.elements[j] * LD.elements[j][i];
    }
    C.elements[i] = C.elements[i] - temp;
  }
}

// friend of Matrix
void FastLU(Matrix& A, Matrix& LU, int *indx){ // LU is the LU decomposition of A
  int i,imax=0,j,k;
  int n = A.numrows;
  double big, dum, sum, temp;
  double vv[10000];

  LU = A;
  for (i=0;i<n;i++){
    big=0.0;
    for (j=0;j<n;j++){
      temp=fabs(LU.rows[i][j]);
      if (temp > big) big=temp;
    }
    vv[i]=1.0/big;
  }
  for (j=0;j<n;j++){
    for (i=0;i<j;i++){
      sum=LU.rows[i][j];
      for (k=0;k<i;k++) sum -= LU.rows[i][k]*LU.rows[k][j];
      LU.rows[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<n;i++){
      sum=LU.rows[i][j];
      for (k=0;k<j;k++)
        sum -= LU.rows[i][k]*LU.rows[k][j];
      LU.rows[i][j]=sum;
      if ((dum=vv[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k<n;k++) {
        dum=LU.rows[imax][k];
        LU.rows[imax][k]=LU.rows[j][k];
        LU.rows[j][k]=dum;
      }
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    // if (LU.rows[j][j] == 0.0) LU.rows[j][j]=1.0e-20;
    if (j != n-1) {
      dum=1.0/(LU.rows[j][j]);
      for (i=j+1;i<n;i++) LU.rows[i][j] *= dum;
    }
  }
}

// friend of Mat3x3
void FastLU(Mat3x3& A, Mat3x3& LU, int *indx){ // LU is the LU decomposition of A
  int i,imax=0,j,k;
  double big, dum, sum, temp;
  double vv[10000];

  LU = A;
  for (i=0;i<3;i++){
    big=0.0;
    for (j=0;j<3;j++){
      temp=fabs(LU.BasicGet(i,j));
      if (temp > big) big=temp;
    }
    vv[i]=1.0/big;
  }
  for (j=0;j<3;j++){
    for (i=0;i<j;i++){
      sum=LU.BasicGet(i,j);
      for (k=0;k<i;k++) sum -= LU.BasicGet(i,k)*LU.BasicGet(k,j);
      LU.BasicSet(i,j,sum);
    }
    big=0.0;
    for (i=j;i<3;i++){
      sum=LU.BasicGet(i,j);
      for (k=0;k<j;k++)
        sum -= LU.BasicGet(i,k)*LU.BasicGet(k,j);
      LU.BasicSet(i,j,sum);
      if ((dum=vv[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k<3;k++) {
        dum=LU.BasicGet(imax,k);
        LU.BasicSet(imax,k,LU.BasicGet(j,k));
        LU.BasicSet(j,k,dum);
      }
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (j != 3-1) {
      dum=1.0/(LU.BasicGet(j,j));
      for (i=j+1;i<3;i++) LU.BasicSet(i,j,LU.BasicGet(i,j)*dum);
    }
  }
}

// friend of Mat4x4
void FastLU(Mat4x4& A, Mat4x4& LU, int *indx){ // LU is the LU decomposition of A
  int i,imax=0,j,k;
  double big, dum, sum, temp;
  double vv[10000];

  LU = A;
  for (i=0;i<4;i++){
    big=0.0;
    for (j=0;j<4;j++){
      temp=fabs(LU.BasicGet(i,j));
      if (temp > big) big=temp;
      }
    vv[i]=1.0/big;
  }
  for (j=0;j<4;j++){
    for (i=0;i<j;i++){
      sum=LU.BasicGet(i,j);
      for (k=0;k<i;k++) sum -= LU.BasicGet(i,k)*LU.BasicGet(k,j);
      LU.BasicSet(i,j,sum);
    }
    big=0.0;
    for (i=j;i<4;i++){
      sum=LU.BasicGet(i,j);
      for (k=0;k<j;k++)
        sum -= LU.BasicGet(i,k)*LU.BasicGet(k,j);
      LU.BasicSet(i,j,sum);
      if ((dum=vv[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k<4;k++) {
        dum=LU.BasicGet(imax,k);
        LU.BasicSet(imax,k,LU.BasicGet(j,k));
        LU.BasicSet(j,k,dum);
      }
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (j != 4-1) {
      dum=1.0/(LU.BasicGet(j,j));
      for (i=j+1;i<4;i++) LU.BasicSet(i,j,LU.BasicGet(i,j)*dum);
    }
  }
}

// friend of Mat6x6
void FastLU(Mat6x6& A, Mat6x6& LU, int *indx){ // LU is the LU decomposition of A
  int i,imax=0,j,k;
  double big, dum, sum, temp;
  double vv[10000];

  LU = A;
  for (i=0;i<6;i++){
    big=0.0;
    for (j=0;j<6;j++){
      temp=fabs(LU.BasicGet(i,j));
      if (temp > big) big=temp;
      }
    vv[i]=1.0/big;
  }
  for (j=0;j<6;j++){
    for (i=0;i<j;i++){
      sum=LU.BasicGet(i,j);
      for (k=0;k<i;k++) sum -= LU.BasicGet(i,k)*LU.BasicGet(k,j);
      LU.BasicSet(i,j,sum);
    }
    big=0.0;
    for (i=j;i<6;i++){
      sum=LU.BasicGet(i,j);
      for (k=0;k<j;k++)
        sum -= LU.BasicGet(i,k)*LU.BasicGet(k,j);
      LU.BasicSet(i,j,sum);
      if ((dum=vv[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k<6;k++) {
        dum=LU.BasicGet(imax,k);
        LU.BasicSet(imax,k,LU.BasicGet(j,k));
        LU.BasicSet(j,k,dum);
      }
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (j != 6-1) {
      dum=1.0/(LU.BasicGet(j,j));
      for (i=j+1;i<6;i++) LU.BasicSet(i,j,LU.BasicGet(i,j)*dum);
    }
  }
}

// friend of Matrix
void FastLUSubs(Matrix& LU, Matrix& B, Matrix& C, int *indx){ // Appropriate Forward and Back Substitution	
  int i,ip,j,k;
  int n = B.numrows;
  int c = B.numcols;
  double temp;
  
  C = B;
  for (k=0;k<c;k++){
	for (i=0;i<n;i++){
		ip=indx[i];
		temp=C.rows[ip][k];
		C.rows[ip][k]=C.rows[i][k];
		for (j=0;j<i;j++) temp -= LU.rows[i][j]*C.rows[j][k];
		C.rows[i][k]=temp;
	}
	for (i=n-1;i>=0;i--){
		temp=C.rows[i][k];
		for (j=i+1;j<n;j++) temp -= LU.rows[i][j]*C.rows[j][k];
		C.rows[i][k]=temp/LU.rows[i][i];
	}
  }
}

// friend of Matrix and Mat3x3
void FastLUSubs(Mat3x3& LU, Matrix& B, Matrix& C, int *indx){ // Appropriate Forward and Back Substitution	
  int i,ip,j,k;
  int n = B.numrows;
  int c = B.numcols;
  double temp;
  
  C = B;
  for (k=0;k<c;k++){
	for (i=0;i<n;i++){
		ip=indx[i];
		temp=C.rows[ip][k];
		C.rows[ip][k]=C.rows[i][k];
		for (j=0;j<i;j++) temp -= LU.BasicGet(i,j)*C.rows[j][k];
		C.rows[i][k]=temp;
	}
	for (i=n-1;i>=0;i--){
		temp=C.rows[i][k];
		for (j=i+1;j<n;j++) temp -= LU.BasicGet(i,j)*C.rows[j][k];
		C.rows[i][k]=temp/LU.BasicGet(i,i);
	}
  }
}

// friend of Matrix and Mat4x4
void FastLUSubs(Mat4x4& LU, Matrix& B, Matrix& C, int *indx){ // Appropriate Forward and Back Substitution
  int i,ip,j,k;
  int n = B.numrows;
  int c = B.numcols;
  double temp;
  
  C = B;
  for (k=0;k<c;k++){
	for (i=0;i<n;i++){
		ip=indx[i];
		temp=C.rows[ip][k];
		C.rows[ip][k]=C.rows[i][k];
		for (j=0;j<i;j++) temp -= LU.BasicGet(i,j)*C.rows[j][k];
		C.rows[i][k]=temp;
	}
	for (i=n-1;i>=0;i--){
		temp=C.rows[i][k];
		for (j=i+1;j<n;j++) temp -= LU.BasicGet(i,j)*C.rows[j][k];
		C.rows[i][k]=temp/LU.BasicGet(i,i);
	}
  }
}

// friend of Matrix and Mat6x6
void FastLUSubs(Mat6x6& LU, Matrix& B, Matrix& C, int *indx){ // Appropriate Forward and Back Substitution	
  int i,ip,j,k;
  int n = B.numrows;
  int c = B.numcols;
  double temp;
  
  C = B;
  for (k=0;k<c;k++){
	for (i=0;i<n;i++){
		ip=indx[i];
		temp=C.rows[ip][k];
		C.rows[ip][k]=C.rows[i][k];
		for (j=0;j<i;j++) temp -= LU.BasicGet(i,j)*C.rows[j][k];
		C.rows[i][k]=temp;
	}
	for (i=n-1;i>=0;i--){
		temp=C.rows[i][k];
		for (j=i+1;j<n;j++) temp -= LU.BasicGet(i,j)*C.rows[j][k];
		C.rows[i][k]=temp/LU.BasicGet(i,i);
	}
  }
}


// The following LUSubsLH routine is incomplete at the moment.
// friend of Matrix
void FastLUSubsLH(Matrix& LU, Matrix& B, Matrix& C, int *indx){ // Appropriate Forward and Back Substitution	
  int i,ip,j,k;
  int n = B.numcols;
  int c = B.numrows;
  double temp;
  Matrix C_temp;
  
  C_temp = B;
  for (k=0;k<c;k++){
	for (i=0;i<n;i++){
		ip=indx[k];
	    temp=C_temp.rows[ip][i];
		C_temp.rows[ip][i]=C_temp.rows[k][i];
		for (j=0;j<i;j++) temp -= LU.rows[i][j]*C_temp.rows[k][j];
		C_temp.rows[k][i]=temp;
	}
	for (i=n-1;i>=0;i--){
		temp=C_temp.rows[k][i];
		for (j=i+1;j<n;j++) temp -= LU.rows[i][j]*C_temp.rows[k][j];
		C_temp.rows[k][i]=temp/LU.rows[i][i];
	}
  }
  C = C_temp;
}



//
// Triple sum
//

void FastTripleSum(Vect3& a, Vect3& b, Vect3& c, Vect3& d){ // d = a+b+c
  d.elements[0] = a.elements[0]+b.elements[0]+c.elements[0];
  d.elements[1] = a.elements[1]+b.elements[1]+c.elements[1];
  d.elements[2] = a.elements[2]+b.elements[2]+c.elements[2];
}

void FastTripleSumPPM(Vect3& a, Vect3& b, Vect3& c, Vect3& d){ // d = a+b-c
  d.elements[0] = a.elements[0]+b.elements[0]-c.elements[0];
  d.elements[1] = a.elements[1]+b.elements[1]-c.elements[1];
  d.elements[2] = a.elements[2]+b.elements[2]-c.elements[2];
}

//
// Multiplications
//

// friend of matrix
void FastMult(Matrix& A, Matrix& B, Matrix& C){  // C = A*B
  // assumes dimensions are already correct!
  int r = A.numrows;
  int ca = A.numcols;
  int cb = B.numcols;

  int i,j,k;
  for(i=0;i<r;i++)
    for(j=0;j<cb;j++){
      C.rows[i][j] = 0.0;
      for(k=0;k<ca;k++)
        C.rows[i][j] += A.rows[i][k] * B.rows[k][j];
    }
}

// friend of matrix
void FastTMult(Matrix& A, Matrix& B, Matrix& C){  // C = A*B
  // assumes dimensions are already correct!
  int r = A.numcols;
  int ca = A.numrows;
  int cb = B.numcols;

  int i,j,k;
  for(i=0;i<r;i++)
    for(j=0;j<cb;j++){
      C.rows[i][j] = A.rows[0][i] * B.rows[0][j];
      for(k=1;k<ca;k++)
        C.rows[i][j] += A.rows[k][i] * B.rows[k][j];
    }
}

// friend of Mat3x3 & Vect3
void FastMult(Mat3x3& A, Vect3& B, Vect3& C){  // C = A*B
  C.elements[0] = A.elements[0][0]*B.elements[0] + A.elements[0][1]*B.elements[1] + A.elements[0][2]*B.elements[2];
  C.elements[1] = A.elements[1][0]*B.elements[0] + A.elements[1][1]*B.elements[1] + A.elements[1][2]*B.elements[2];
  C.elements[2] = A.elements[2][0]*B.elements[0] + A.elements[2][1]*B.elements[1] + A.elements[2][2]*B.elements[2];
}

// friend of Mat3x3, ColMatrix, & Vect3
void FastMult(Mat3x3& A, ColMatrix& B, Vect3& C){  // C = A*B
  C.elements[0] = A.elements[0][0]*B.elements[0] + A.elements[0][1]*B.elements[1] + A.elements[0][2]*B.elements[2];
  C.elements[1] = A.elements[1][0]*B.elements[0] + A.elements[1][1]*B.elements[1] + A.elements[1][2]*B.elements[2];
  C.elements[2] = A.elements[2][0]*B.elements[0] + A.elements[2][1]*B.elements[1] + A.elements[2][2]*B.elements[2];
}


// friend of Mat3x3, ColMatrix, & Vect3
void FastMult(Mat3x3& A, Vect3& B, ColMatrix& C){  // C = A*B
  C.elements[0] = A.elements[0][0]*B.elements[0] + A.elements[0][1]*B.elements[1] + A.elements[0][2]*B.elements[2];
  C.elements[1] = A.elements[1][0]*B.elements[0] + A.elements[1][1]*B.elements[1] + A.elements[1][2]*B.elements[2];
  C.elements[2] = A.elements[2][0]*B.elements[0] + A.elements[2][1]*B.elements[1] + A.elements[2][2]*B.elements[2];
}

// friend of Mat3x3 & Vect3
void FastTMult(Mat3x3& A, Vect3& B, Vect3& C){  // C = A^T*B
  C.elements[0] = A.elements[0][0]*B.elements[0] + A.elements[1][0]*B.elements[1] + A.elements[2][0]*B.elements[2];
  C.elements[1] = A.elements[0][1]*B.elements[0] + A.elements[1][1]*B.elements[1] + A.elements[2][1]*B.elements[2];
  C.elements[2] = A.elements[0][2]*B.elements[0] + A.elements[1][2]*B.elements[1] + A.elements[2][2]*B.elements[2];
}

// friend of Mat3x3 & Vect3
void FastNegMult(Mat3x3& A, Vect3& B, Vect3& C){  // C = -A*B
  C.elements[0] = -A.elements[0][0]*B.elements[0] - A.elements[0][1]*B.elements[1] - A.elements[0][2]*B.elements[2];
  C.elements[1] = -A.elements[1][0]*B.elements[0] - A.elements[1][1]*B.elements[1] - A.elements[1][2]*B.elements[2];
  C.elements[2] = -A.elements[2][0]*B.elements[0] - A.elements[2][1]*B.elements[1] - A.elements[2][2]*B.elements[2];
}

// friend of Mat3x3 & Vect3
void FastNegTMult(Mat3x3& A, Vect3& B, Vect3& C){  // C = -A^T*B
  C.elements[0] = -A.elements[0][0]*B.elements[0] - A.elements[1][0]*B.elements[1] - A.elements[2][0]*B.elements[2];
  C.elements[1] = -A.elements[0][1]*B.elements[0] - A.elements[1][1]*B.elements[1] - A.elements[2][1]*B.elements[2];
  C.elements[2] = -A.elements[0][2]*B.elements[0] - A.elements[1][2]*B.elements[1] - A.elements[2][2]*B.elements[2];
}

// friend of Vect3
void FastMult(double a, Vect3& B, Vect3& C){  // C = a*B
  C.elements[0] = a*B.elements[0];
  C.elements[1] = a*B.elements[1];
  C.elements[2] = a*B.elements[2];
}

// friend of Mat4x4 & Vect4
void FastMult(Mat4x4& A, Vect4& B, Vect4& C){  // C = A*B
  C.elements[0] = A.elements[0][0]*B.elements[0] + A.elements[0][1]*B.elements[1] + A.elements[0][2]*B.elements[2] + A.elements[0][3]*B.elements[3];
  C.elements[1] = A.elements[1][0]*B.elements[0] + A.elements[1][1]*B.elements[1] + A.elements[1][2]*B.elements[2] + A.elements[1][3]*B.elements[3];
  C.elements[2] = A.elements[2][0]*B.elements[0] + A.elements[2][1]*B.elements[1] + A.elements[2][2]*B.elements[2] + A.elements[2][3]*B.elements[3];
  C.elements[3] = A.elements[3][0]*B.elements[0] + A.elements[3][1]*B.elements[1] + A.elements[3][2]*B.elements[2] + A.elements[3][3]*B.elements[3];
}

// friend of Mat4x4 & Vect4
void FastTMult(Mat4x4& A, Vect4& B, Vect4& C){  // C = A^T*B
  C.elements[0] = A.elements[0][0]*B.elements[0] + A.elements[1][0]*B.elements[1] + A.elements[2][0]*B.elements[2] + A.elements[3][0]*B.elements[3];
  C.elements[1] = A.elements[0][1]*B.elements[0] + A.elements[1][1]*B.elements[1] + A.elements[2][1]*B.elements[2] + A.elements[3][1]*B.elements[3];
  C.elements[2] = A.elements[0][2]*B.elements[0] + A.elements[1][2]*B.elements[1] + A.elements[2][2]*B.elements[2] + A.elements[3][2]*B.elements[3];
  C.elements[3] = A.elements[0][3]*B.elements[0] + A.elements[1][3]*B.elements[1] + A.elements[2][3]*B.elements[2] + A.elements[3][3]*B.elements[3];
}

// friend of Mat4x4 & Vect4
void FastNegMult(Mat4x4& A, Vect4& B, Vect4& C){  // C = -A*B
  C.elements[0] = -A.elements[0][0]*B.elements[0] - A.elements[0][1]*B.elements[1] - A.elements[0][2]*B.elements[2] - A.elements[0][3]*B.elements[3];
  C.elements[1] = -A.elements[1][0]*B.elements[0] - A.elements[1][1]*B.elements[1] - A.elements[1][2]*B.elements[2] - A.elements[1][3]*B.elements[3];
  C.elements[2] = -A.elements[2][0]*B.elements[0] - A.elements[2][1]*B.elements[1] - A.elements[2][2]*B.elements[2] - A.elements[2][3]*B.elements[3];
  C.elements[3] = -A.elements[3][0]*B.elements[0] - A.elements[3][1]*B.elements[1] - A.elements[3][2]*B.elements[2] - A.elements[3][3]*B.elements[3];
}

// friend of Mat4x4 & Vect4
void FastNegTMult(Mat4x4& A, Vect4& B, Vect4& C){  // C = -A^T*B
  C.elements[0] = -A.elements[0][0]*B.elements[0] - A.elements[1][0]*B.elements[1] - A.elements[2][0]*B.elements[2] - A.elements[3][0]*B.elements[3];
  C.elements[1] = -A.elements[0][1]*B.elements[0] - A.elements[1][1]*B.elements[1] - A.elements[2][1]*B.elements[2] - A.elements[3][1]*B.elements[3];
  C.elements[2] = -A.elements[0][2]*B.elements[0] - A.elements[1][2]*B.elements[1] - A.elements[2][2]*B.elements[2] - A.elements[3][2]*B.elements[3];
  C.elements[3] = -A.elements[0][3]*B.elements[0] - A.elements[1][3]*B.elements[1] - A.elements[2][3]*B.elements[2] - A.elements[3][3]*B.elements[3];
}

// friend of Vect4
void FastMult(double a, Vect4& B, Vect4& C){  // C = a*B
  C.elements[0] = a*B.elements[0];
  C.elements[1] = a*B.elements[1];
  C.elements[2] = a*B.elements[2];
  C.elements[3] = a*B.elements[3];
}

// friend of Matrix & Mat6x6
void FastMultT(Matrix& A, Matrix& B, Mat6x6& C){  // C = A*B^T
  int i,j,k,n;
  n = A.numcols;

  for(i=0;i<6;i++)
    for(j=0;j<6;j++){
      C.elements[i][j] = 0.0;
      for(k=0;k<n;k++)
        C.elements[i][j] += A.rows[i][k] * B.rows[j][k];
    }
}

// friend Matrix, Vect6, ColMatrix
void FastMult(Matrix& A, ColMatrix& B, Vect6& C){
  int ca = A.numcols;

  int i,k;
  for(i=0;i<6;i++){
    C.elements[i] = 0.0;
    for(k=0;k<ca;k++)
      C.elements[i] += A.rows[i][k] * B.elements[k];
  }
}

// friend of Matrix & Mat6x6
void FastMult(Mat6x6& A, Matrix& B, Matrix& C){  // C = A*B
  // assumes dimensions are already correct!
  int cb = B.numcols;

  int i,j,k;
  for(i=0;i<6;i++)
    for(j=0;j<cb;j++){
      C.rows[i][j] = 0.0;
      for(k=0;k<6;k++)
        C.rows[i][j] += A.elements[i][k] * B.rows[k][j];
    }
}

// friend Matrix, Vect6, ColMatrix
void FastTMult(Matrix& A, Vect6& B, ColMatrix& C){  // C = A^T*B
  int n = C.numrows;
  int i,k;
  for(i=0;i<n;i++){
    C.elements[i] = 0.0;
    for(k=0;k<6;k++)
      C.elements[i] += A.rows[k][i] * B.elements[k];
  }
}

// friend of Mat3x3
void FastMult(Mat3x3& A, Mat3x3& B, Mat3x3& C){  // C = A*B
  C.elements[0][0] = A.elements[0][0]*B.elements[0][0] + A.elements[0][1]*B.elements[1][0] + A.elements[0][2]*B.elements[2][0];
  C.elements[0][1] = A.elements[0][0]*B.elements[0][1] + A.elements[0][1]*B.elements[1][1] + A.elements[0][2]*B.elements[2][1];
  C.elements[0][2] = A.elements[0][0]*B.elements[0][2] + A.elements[0][1]*B.elements[1][2] + A.elements[0][2]*B.elements[2][2];

  C.elements[1][0] = A.elements[1][0]*B.elements[0][0] + A.elements[1][1]*B.elements[1][0] + A.elements[1][2]*B.elements[2][0];
  C.elements[1][1] = A.elements[1][0]*B.elements[0][1] + A.elements[1][1]*B.elements[1][1] + A.elements[1][2]*B.elements[2][1];
  C.elements[1][2] = A.elements[1][0]*B.elements[0][2] + A.elements[1][1]*B.elements[1][2] + A.elements[1][2]*B.elements[2][2];

  C.elements[2][0] = A.elements[2][0]*B.elements[0][0] + A.elements[2][1]*B.elements[1][0] + A.elements[2][2]*B.elements[2][0];
  C.elements[2][1] = A.elements[2][0]*B.elements[0][1] + A.elements[2][1]*B.elements[1][1] + A.elements[2][2]*B.elements[2][1];
  C.elements[2][2] = A.elements[2][0]*B.elements[0][2] + A.elements[2][1]*B.elements[1][2] + A.elements[2][2]*B.elements[2][2];
}

// friend of Mat3x3
void FastMultT(Mat3x3& A, Mat3x3& B, Mat3x3& C){  // C = A*B^T
  C.elements[0][0] = A.elements[0][0]*B.elements[0][0] + A.elements[0][1]*B.elements[0][1] + A.elements[0][2]*B.elements[0][2];
  C.elements[0][1] = A.elements[0][0]*B.elements[1][0] + A.elements[0][1]*B.elements[1][1] + A.elements[0][2]*B.elements[1][2];
  C.elements[0][2] = A.elements[0][0]*B.elements[2][0] + A.elements[0][1]*B.elements[2][1] + A.elements[0][2]*B.elements[2][2];

  C.elements[1][0] = A.elements[1][0]*B.elements[0][0] + A.elements[1][1]*B.elements[0][1] + A.elements[1][2]*B.elements[0][2];
  C.elements[1][1] = A.elements[1][0]*B.elements[1][0] + A.elements[1][1]*B.elements[1][1] + A.elements[1][2]*B.elements[1][2];
  C.elements[1][2] = A.elements[1][0]*B.elements[2][0] + A.elements[1][1]*B.elements[2][1] + A.elements[1][2]*B.elements[2][2];

  C.elements[2][0] = A.elements[2][0]*B.elements[0][0] + A.elements[2][1]*B.elements[0][1] + A.elements[2][2]*B.elements[0][2];
  C.elements[2][1] = A.elements[2][0]*B.elements[1][0] + A.elements[2][1]*B.elements[1][1] + A.elements[2][2]*B.elements[1][2];
  C.elements[2][2] = A.elements[2][0]*B.elements[2][0] + A.elements[2][1]*B.elements[2][1] + A.elements[2][2]*B.elements[2][2];
}

// friend of Mat4x4
void FastMult(Mat4x4& A, Mat4x4& B, Mat4x4& C){  // C = A*B
  C.elements[0][0] = A.elements[0][0]*B.elements[0][0] + A.elements[0][1]*B.elements[1][0] + A.elements[0][2]*B.elements[2][0] + A.elements[0][3]*B.elements[3][0];
  C.elements[0][1] = A.elements[0][0]*B.elements[0][1] + A.elements[0][1]*B.elements[1][1] + A.elements[0][2]*B.elements[2][1] + A.elements[0][3]*B.elements[3][1];
  C.elements[0][2] = A.elements[0][0]*B.elements[0][2] + A.elements[0][1]*B.elements[1][2] + A.elements[0][2]*B.elements[2][2] + A.elements[0][3]*B.elements[3][2];
  C.elements[0][3] = A.elements[0][0]*B.elements[0][3] + A.elements[0][1]*B.elements[1][3] + A.elements[0][2]*B.elements[2][3] + A.elements[0][3]*B.elements[3][3];

  C.elements[1][0] = A.elements[1][0]*B.elements[0][0] + A.elements[1][1]*B.elements[1][0] + A.elements[1][2]*B.elements[2][0] + A.elements[1][3]*B.elements[3][0];
  C.elements[1][1] = A.elements[1][0]*B.elements[0][1] + A.elements[1][1]*B.elements[1][1] + A.elements[1][2]*B.elements[2][1] + A.elements[1][3]*B.elements[3][1];
  C.elements[1][2] = A.elements[1][0]*B.elements[0][2] + A.elements[1][1]*B.elements[1][2] + A.elements[1][2]*B.elements[2][2] + A.elements[1][3]*B.elements[3][2];
  C.elements[1][3] = A.elements[1][0]*B.elements[0][3] + A.elements[1][1]*B.elements[1][3] + A.elements[1][2]*B.elements[2][3] + A.elements[1][3]*B.elements[3][3];

  C.elements[2][0] = A.elements[2][0]*B.elements[0][0] + A.elements[2][1]*B.elements[1][0] + A.elements[2][2]*B.elements[2][0] + A.elements[2][3]*B.elements[3][0];
  C.elements[2][1] = A.elements[2][0]*B.elements[0][1] + A.elements[2][1]*B.elements[1][1] + A.elements[2][2]*B.elements[2][1] + A.elements[2][3]*B.elements[3][1];
  C.elements[2][2] = A.elements[2][0]*B.elements[0][2] + A.elements[2][1]*B.elements[1][2] + A.elements[2][2]*B.elements[2][2] + A.elements[2][3]*B.elements[3][2];
  C.elements[2][3] = A.elements[2][0]*B.elements[0][3] + A.elements[2][1]*B.elements[1][3] + A.elements[2][2]*B.elements[2][3] + A.elements[2][3]*B.elements[3][3];

  C.elements[3][0] = A.elements[3][0]*B.elements[0][0] + A.elements[3][1]*B.elements[1][0] + A.elements[3][2]*B.elements[2][0] + A.elements[3][3]*B.elements[3][0];
  C.elements[3][1] = A.elements[3][0]*B.elements[0][1] + A.elements[3][1]*B.elements[1][1] + A.elements[3][2]*B.elements[2][1] + A.elements[3][3]*B.elements[3][1];
  C.elements[3][2] = A.elements[3][0]*B.elements[0][2] + A.elements[3][1]*B.elements[1][2] + A.elements[3][2]*B.elements[2][2] + A.elements[3][3]*B.elements[3][2];
  C.elements[3][3] = A.elements[3][0]*B.elements[0][3] + A.elements[3][1]*B.elements[1][3] + A.elements[3][2]*B.elements[2][3] + A.elements[3][3]*B.elements[3][3];
}

// friend of Mat4x4
void FastMultT(Mat4x4& A, Mat4x4& B, Mat4x4& C){  // C = A*B^T
  C.elements[0][0] = A.elements[0][0]*B.elements[0][0] + A.elements[0][1]*B.elements[0][1] + A.elements[0][2]*B.elements[0][2] + A.elements[0][3]*B.elements[0][3];
  C.elements[0][1] = A.elements[0][0]*B.elements[1][0] + A.elements[0][1]*B.elements[1][1] + A.elements[0][2]*B.elements[1][2] + A.elements[0][3]*B.elements[1][3];
  C.elements[0][2] = A.elements[0][0]*B.elements[2][0] + A.elements[0][1]*B.elements[2][1] + A.elements[0][2]*B.elements[2][2] + A.elements[0][3]*B.elements[2][3];
  C.elements[0][3] = A.elements[0][0]*B.elements[3][0] + A.elements[0][1]*B.elements[3][1] + A.elements[0][2]*B.elements[3][2] + A.elements[0][3]*B.elements[3][3];

  C.elements[1][0] = A.elements[1][0]*B.elements[0][0] + A.elements[1][1]*B.elements[0][1] + A.elements[1][2]*B.elements[0][2] + A.elements[1][3]*B.elements[0][3];
  C.elements[1][1] = A.elements[1][0]*B.elements[1][0] + A.elements[1][1]*B.elements[1][1] + A.elements[1][2]*B.elements[1][2] + A.elements[1][3]*B.elements[1][3];
  C.elements[1][2] = A.elements[1][0]*B.elements[2][0] + A.elements[1][1]*B.elements[2][1] + A.elements[1][2]*B.elements[2][2] + A.elements[1][3]*B.elements[2][3];
  C.elements[1][3] = A.elements[1][0]*B.elements[3][0] + A.elements[1][1]*B.elements[3][1] + A.elements[1][2]*B.elements[3][2] + A.elements[1][3]*B.elements[3][3];

  C.elements[2][0] = A.elements[2][0]*B.elements[0][0] + A.elements[2][1]*B.elements[0][1] + A.elements[2][2]*B.elements[0][2] + A.elements[2][3]*B.elements[0][3];
  C.elements[2][1] = A.elements[2][0]*B.elements[1][0] + A.elements[2][1]*B.elements[1][1] + A.elements[2][2]*B.elements[1][2] + A.elements[2][3]*B.elements[1][3];
  C.elements[2][2] = A.elements[2][0]*B.elements[2][0] + A.elements[2][1]*B.elements[2][1] + A.elements[2][2]*B.elements[2][2] + A.elements[2][3]*B.elements[2][3];
  C.elements[2][3] = A.elements[2][0]*B.elements[3][0] + A.elements[2][1]*B.elements[3][1] + A.elements[2][2]*B.elements[3][2] + A.elements[2][3]*B.elements[3][3];

  C.elements[3][0] = A.elements[3][0]*B.elements[0][0] + A.elements[3][1]*B.elements[0][1] + A.elements[3][2]*B.elements[0][2] + A.elements[3][3]*B.elements[0][3];
  C.elements[3][1] = A.elements[3][0]*B.elements[1][0] + A.elements[3][1]*B.elements[1][1] + A.elements[3][2]*B.elements[1][2] + A.elements[3][3]*B.elements[1][3];
  C.elements[3][2] = A.elements[3][0]*B.elements[2][0] + A.elements[3][1]*B.elements[2][1] + A.elements[3][2]*B.elements[2][2] + A.elements[3][3]*B.elements[2][3];
  C.elements[3][3] = A.elements[3][0]*B.elements[3][0] + A.elements[3][1]*B.elements[3][1] + A.elements[3][2]*B.elements[3][2] + A.elements[3][3]*B.elements[3][3];
  C.elements[2][2] = A.elements[2][0]*B.elements[2][0] + A.elements[2][1]*B.elements[2][1] + A.elements[2][2]*B.elements[2][2];
}

// friend of Mat6x6
void FastMult(Mat6x6& A, Mat6x6& B, Mat6x6& C){  // C = A*B
  int i,j,k;
  for(i=0;i<6;i++)
    for(j=0;j<6;j++){
      C.elements[i][j] = 0.0;
      for(k=0;k<6;k++)
        C.elements[i][j] += A.elements[i][k]*B.elements[k][j];
    }
}

// friend of Mat6x6
void FastMultT(Mat6x6& A, Mat6x6& B, Mat6x6& C){  // C = A*B
  int i,j,k;
  for(i=0;i<6;i++)
    for(j=0;j<6;j++){
      C.elements[i][j] = 0.0;
      for(k=0;k<6;k++)
        C.elements[i][j] += A.elements[i][k]*B.elements[j][k];
    }
}

// friend of Mat6x6
void FastTMult(Mat6x6& A, Mat6x6& B, Mat6x6& C){  // C = A^T*B
  int i,j,k;
  for(i=0;i<6;i++)
    for(j=0;j<6;j++){
      C.elements[i][j] = 0.0;
      for(k=0;k<6;k++)
        C.elements[i][j] += A.elements[k][i]*B.elements[k][j];
    }
}

// friend of Mat6x6 & Vect6
void FastMult(Mat6x6& A, Vect6& B, Vect6& C){  // C = A*B
  for(int i=0;i<6;i++)
    C.elements[i] = A.elements[i][0]*B.elements[0] + A.elements[i][1]*B.elements[1] + A.elements[i][2]*B.elements[2] + A.elements[i][3]*B.elements[3] + A.elements[i][4]*B.elements[4] + A.elements[i][5]*B.elements[5];
}

// friend of Mat6x6 & Vect6
void FastTMult(Mat6x6& A, Vect6& B, Vect6& C){  // C = A^T*B
  for(int i=0;i<6;i++)
    C.elements[i] = A.elements[0][i]*B.elements[0] + A.elements[1][i]*B.elements[1] + A.elements[2][i]*B.elements[2] + A.elements[3][i]*B.elements[3] + A.elements[4][i]*B.elements[4] + A.elements[5][i]*B.elements[5];
}

//
// Additions
//

// friend of Vect3
void FastAdd(Vect3& A, Vect3& B, Vect3& C){ // C = A+B
  C.elements[0] = A.elements[0] + B.elements[0];
  C.elements[1] = A.elements[1] + B.elements[1];
  C.elements[2] = A.elements[2] + B.elements[2];
}

// friend of Vect4
void FastAdd(Vect4& A, Vect4& B, Vect4& C){ // C = A+B
  C.elements[0] = A.elements[0] + B.elements[0];
  C.elements[1] = A.elements[1] + B.elements[1];
  C.elements[2] = A.elements[2] + B.elements[2];
  C.elements[3] = A.elements[3] + B.elements[3];
}

void FastAdd(Mat6x6& A, Mat6x6& B, Mat6x6& C){  // C = A+B
  int i,j;
  for(i=0;i<6;i++)
    for(j=0;j<6;j++)
      C.elements[i][j] = A.elements[i][j] + B.elements[i][j];
}

// friend of Vect6
void FastAdd(Vect6& A, Vect6& B, Vect6& C){ // C = A-B
  C.elements[0] = A.elements[0] + B.elements[0];
  C.elements[1] = A.elements[1] + B.elements[1];
  C.elements[2] = A.elements[2] + B.elements[2];
  C.elements[3] = A.elements[3] + B.elements[3];
  C.elements[4] = A.elements[4] + B.elements[4];
  C.elements[5] = A.elements[5] + B.elements[5];
}

//
// Subtractions
//

// friend of Vect3
void FastSubt(Vect3& A, Vect3& B, Vect3& C){ // C = A-B
  C.elements[0] = A.elements[0] - B.elements[0];
  C.elements[1] = A.elements[1] - B.elements[1];
  C.elements[2] = A.elements[2] - B.elements[2];
}

// friend of Vect4
void FastSubt(Vect4& A, Vect4& B, Vect4& C){ // C = A-B
  C.elements[0] = A.elements[0] - B.elements[0];
  C.elements[1] = A.elements[1] - B.elements[1];
  C.elements[2] = A.elements[2] - B.elements[2];
  C.elements[3] = A.elements[3] - B.elements[3];
}

void FastSubt(Mat6x6& A, Mat6x6& B, Mat6x6& C){  // C = A-B
  int i,j;
  for(i=0;i<6;i++)
    for(j=0;j<6;j++)
      C.elements[i][j] = A.elements[i][j] - B.elements[i][j];
}

// friend of Vect6
void FastSubt(Vect6& A, Vect6& B, Vect6& C){ // C = A-B
  C.elements[0] = A.elements[0] - B.elements[0];
  C.elements[1] = A.elements[1] - B.elements[1];
  C.elements[2] = A.elements[2] - B.elements[2];
  C.elements[3] = A.elements[3] - B.elements[3];
  C.elements[4] = A.elements[4] - B.elements[4];
  C.elements[5] = A.elements[5] - B.elements[5];
}

// friend of ColMatMap
void FastAssign(ColMatMap& A, ColMatMap& C){
  for(int i=0;i<C.numrows;i++)
    *(C.elements[i]) = *(A.elements[i]);
}

// friend of ColMatrix
void FastAssign(ColMatrix& A, ColMatrix& C){ //C = A
  for(int i=0;i<C.numrows;i++)
    C.elements[i] = A.elements[i];
}

// friend of Vect3
void FastAssign(Vect3& A, Vect3& C){ //C = A
    C.elements[0] = A.elements[0];
    C.elements[1] = A.elements[1];
    C.elements[2] = A.elements[2];
}

// friend of Vect3 & ColMatrix
void FastAssign(ColMatrix& A, Vect3& C){ //C = A
    C.elements[0] = A.elements[0];
    C.elements[1] = A.elements[1];
    C.elements[2] = A.elements[2];
}

// friend of Vect4
void FastAssign(Vect4& A, Vect4& C){ //C = A
    C.elements[0] = A.elements[0];
    C.elements[1] = A.elements[1];
    C.elements[2] = A.elements[2];
    C.elements[3] = A.elements[3];
}

// friend of Mat3x3
void FastAssignT(Mat3x3& A, Mat3x3& C){
  C.elements[0][0] = A.elements[0][0];
  C.elements[1][1] = A.elements[1][1];
  C.elements[2][2] = A.elements[2][2];

  C.elements[0][1] = A.elements[1][0];
  C.elements[1][0] = A.elements[0][1];

  C.elements[0][2] = A.elements[2][0];
  C.elements[2][0] = A.elements[0][2];

  C.elements[1][2] = A.elements[2][1];
  C.elements[2][1] = A.elements[1][2];
}

// friend of Mat4x4
void FastAssignT(Mat4x4& A, Mat4x4& C){
  C.elements[0][0] = A.elements[0][0];
  C.elements[1][1] = A.elements[1][1];
  C.elements[2][2] = A.elements[2][2];
  C.elements[3][3] = A.elements[3][3];

  C.elements[0][1] = A.elements[1][0];
  C.elements[1][0] = A.elements[0][1];

  C.elements[0][2] = A.elements[2][0];
  C.elements[2][0] = A.elements[0][2];

  C.elements[0][3] = A.elements[3][0];
  C.elements[3][0] = A.elements[0][3];

  C.elements[1][2] = A.elements[2][1];
  C.elements[2][1] = A.elements[1][2];

  C.elements[1][3] = A.elements[3][1];
  C.elements[3][1] = A.elements[1][3];

  C.elements[2][3] = A.elements[3][2];
  C.elements[3][2] = A.elements[2][3];
}
