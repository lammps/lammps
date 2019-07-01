/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: mat3x3.h                                                *
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

#ifndef MAT3X3_H
#define MAT3X3_H

#include <iostream>
#include "virtualmatrix.h"

class Vect3;
class Mat6x6;
class ColMatrix;

class Mat3x3 : public VirtualMatrix  {
  double elements[3][3];
public:
  Mat3x3();
  ~Mat3x3();
  Mat3x3(const Mat3x3& A);  // copy constructor
  Mat3x3(const VirtualMatrix& A);  // copy constructor

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
  Mat3x3& operator=(const Mat3x3& A); // assignment operator
  Mat3x3& operator=(const VirtualMatrix& A); // overloaded =
  Mat3x3& operator*=(double b);

  Mat3x3& Identity();

  friend Mat3x3 T(const Mat3x3& A);  // a wasteful transpose
  friend Mat3x3 CrossMat(Vect3& a);  // a wasteful cross matrix implementation

  friend void FastSimpleRotation(Vect3& v, double q, Mat3x3& d);
  friend void FastQuaternions(ColMatrix& q, Mat3x3& C);
  friend void FastInvQuaternions(Mat3x3& C, ColMatrix& q);
  friend void FastLU(Mat3x3& A, Mat3x3& LU, int *indx);
  friend void FastLUSubs(Mat3x3& LU, Mat3x3& B, Mat3x3& C, int *indx);
  friend void FastMult(Mat3x3& A, Vect3& B, Vect3& C);
  friend void FastTMult(Mat3x3& A, Vect3& B, Vect3& C);
  friend void FastNegMult(Mat3x3& A, Vect3& B, Vect3& C);
  friend void FastNegTMult(Mat3x3& A, Vect3& B, Vect3& C);
  friend void FastMult(Mat3x3& A, Mat3x3& B, Mat3x3& C);
  friend void FastMultT(Mat3x3& A, Mat3x3& B, Mat3x3& C);
  friend void FastAssignT(Mat3x3& A, Mat3x3& C);
  friend void FastMult(Mat3x3& A, Vect3& B, ColMatrix& C);  
  
  friend void OnPopulateSC(Vect3& gamma, Mat3x3& C, Mat6x6& SC);
  friend void OnPopulateSI(Mat3x3& inertia, double mass, Mat6x6& sI);

  friend void FastMult(Mat3x3& A, ColMatrix& B, Vect3& C);
  
  friend void EP_Transformation(ColMatrix& q, Mat3x3& C);
  friend void EP_FromTransformation(ColMatrix& q, Mat3x3& C);
  
};

#endif
