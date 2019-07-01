/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: mat6x6.h                                                *
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

#ifndef MAT6X6_H
#define MAT6X6_H
#include <iostream>

#include "virtualmatrix.h"

class Matrix;
class Mat3x3;
class Vect6;
class Vect3;

class Mat6x6 : public VirtualMatrix  {
  double elements[6][6];
public:
  Mat6x6();
  ~Mat6x6();
  Mat6x6(const Mat6x6& A);  // copy constructor
  Mat6x6(const VirtualMatrix& A);  // copy constructor

  double& operator_2int(int i, int j); // array access
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
  Mat6x6& operator=(const Mat6x6& A); // assignment operator
  Mat6x6& operator=(const VirtualMatrix& A); // overloaded =
  Mat6x6& operator*=(double b);

  Mat6x6& Identity();

  friend Mat6x6 T(const Mat6x6& A);  // a wasteful transpose

  // fast matrix operations
  friend void FastAdd(Mat6x6& A, Mat6x6& B, Mat6x6& C);
  friend void FastSubt(Mat6x6& A, Mat6x6& B, Mat6x6& C);
  friend void FastMultT(Matrix& A, Matrix& B, Mat6x6& C);
  friend void FastMult(Mat6x6& A, Mat6x6& B, Mat6x6& C);  // C = A*B
  friend void FastMult(Mat6x6& A, Matrix& B, Matrix& C);
  friend void FastMultT(Mat6x6& A, Mat6x6& B, Mat6x6& C);  // C = A*B^T
  friend void FastTMult(Mat6x6& A, Mat6x6& B, Mat6x6& C);  // C = A^T*B
  friend void FastMult(Mat6x6& A, Vect6& B, Vect6& C);
  friend void FastTMult(Mat6x6& A, Vect6& B, Vect6& C);
  friend void FastLDLT(Mat6x6& A, Mat6x6& C);
  friend void FastLDLTSubs(Mat6x6& LD, Mat6x6& B, Mat6x6& C);
  friend void FastLDLTSubs(Mat6x6& LD, Vect6& B, Vect6& C);

  friend void OnPopulateSC(Vect3& gamma, Mat3x3& C, Mat6x6& SC);
  friend void OnPopulateSI(Mat3x3& inertia, double mass, Mat6x6& sI);
};

#endif
