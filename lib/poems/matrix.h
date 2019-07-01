/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: matrix.h                                                *
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

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>

#include "virtualmatrix.h"

class Mat3x3;
class Mat4x4;
class Mat6x6;
class Vect6;
class ColMatrix;


class Matrix : public VirtualMatrix {
  double **rows;          // row pointer
  double *elements;
public:
  Matrix();
  ~Matrix();
  Matrix(const Matrix& A);  // copy constructor
  Matrix(const VirtualMatrix& A);  // copy constructor
  Matrix(int m, int n);  // size constructor

  double& operator_2int (int i, int j); // array access
  double Get_2int(int i, int j) const;
  void Set_2int(int i, int j, double value);
  double BasicGet_2int(int i, int j) const;
  void BasicSet_2int(int i, int j, double value);
  void BasicIncrement_2int(int i, int j, double value);

  void Const(double value);
  MatrixType GetType() const;
  std::istream& ReadData(std::istream& c);
  std::ostream& WriteData(std::ostream& c) const;

  Matrix& Dim(int m, int n);    // allocate size

  void AssignVM(const VirtualMatrix& A);
  Matrix& operator=(const Matrix& A); // assignment operator
  Matrix& operator=(const VirtualMatrix& A); // overloaded =
  Matrix& operator*=(double b);

  friend void FastLDLT(Matrix& A, Matrix& C);
  friend void FastLDLTSubs(Matrix& LD, Matrix& B, Matrix& C);
  friend void FastLDLTSubsLH(Matrix& B, Matrix& LD, Matrix& C);
  friend void FastLU(Matrix& A, Matrix& LU, int *indx);
  friend void FastLUSubs(Matrix& LU, Matrix& B, Matrix& C, int *indx);
  friend void FastLUSubs(Mat3x3& LU, Matrix& B, Matrix& C, int *indx);
  friend void FastLUSubs(Mat4x4& LU, Matrix& B, Matrix& C, int *indx);
  friend void FastLUSubs(Mat6x6& LU, Matrix& B, Matrix& C, int *indx);
  friend void FastLUSubsLH(Matrix& B, Matrix& LU, Matrix& C, int *indx);
  friend void FastMult(Matrix& A, Matrix& B, Matrix& C);
  friend void FastTMult(Matrix& A, Matrix& B, Matrix& C);
  friend void FastTMult(Matrix& A, Vect6& B, ColMatrix& C);
  friend void FastMult(Mat6x6& A, Matrix& B, Matrix& C);
  friend void FastMult(Matrix& A, ColMatrix& B, Vect6& C);
  friend void FastMultT(Matrix& A, Matrix& B, Mat6x6& C);
  
};

#endif
