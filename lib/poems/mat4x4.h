/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: mat4x4.h                                                *
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

#ifndef MAT4X4_H
#define MAT4X4_H

#include "virtualmatrix.h"
#include "matrix.h"


class Vect4;

class Mat4x4 : public VirtualMatrix  {
  double elements[4][4];
public:
  Mat4x4();
  ~Mat4x4();
  Mat4x4(const Mat4x4& A);  // copy constructor
  Mat4x4(const VirtualMatrix& A);  // copy constructor

  double& operator_2int (int i, int j); // array access
  double Get_2int(int i, int j) const;
  void Set_2int(int i, int j, double value);
  double BasicGet_2int(int i, int j) const;
  void BasicSet_2int(int i, int j, double value);
  void BasicIncrement_2int(int i, int j, double value);

  void Const(double value);
  MatrixType GetType() const;
  std::istream& ReadData(std::istream& c); // input
  std::ostream& WriteData(std::ostream& c) const; // output

  void AssignVM(const VirtualMatrix& A);
  Mat4x4& operator=(const Mat4x4& A); // assignment operator
  Mat4x4& operator=(const VirtualMatrix& A); // overloaded =
  Mat4x4& operator*=(double b);

  Mat4x4& Identity();

  friend Mat4x4 T(const Mat4x4& A);  // a wasteful transpose
  friend Mat4x4 CrossMat(Vect4& a);  // a wasteful cross matrix implementation

  friend void FastLU(Mat4x4& A, Mat4x4& LU, int *indx);
  friend void FastLUSubs(Mat4x4& LU, Matrix& B, Matrix& C, int *indx);
  friend void FastMult(Mat4x4& A, Vect4& B, Vect4& C);
  friend void FastTMult(Mat4x4& A, Vect4& B, Vect4& C);
  friend void FastNegMult(Mat4x4& A, Vect4& B, Vect4& C);
  friend void FastNegTMult(Mat4x4& A, Vect4& B, Vect4& C);
  friend void FastMult(Mat4x4& A, Mat4x4& B, Mat4x4& C);
  friend void FastMultT(Mat4x4& A, Mat4x4& B, Mat4x4& C);
  friend void FastAssignT(Mat4x4& A, Mat4x4& C);
};

#endif
