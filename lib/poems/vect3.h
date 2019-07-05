/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: vect3.h                                              *
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

#ifndef VECT3_H
#define VECT3_H

#include <iostream>

#include "virtualcolmatrix.h"
#include "virtualmatrix.h"

namespace POEMS {
class Matrix;
class Mat3x3;
class Mat6x6;
class Vect6;
class ColMatrix;

class Vect3 : public VirtualColMatrix  {
  double elements[3];
public:
  Vect3();
  ~Vect3();
  Vect3(const Vect3& A);  // copy constructor
  Vect3(const VirtualMatrix& A);  // copy constructor

  double& operator_1int(int i); // array access
  double Get_1int(int i) const;
  void Set_1int(int i, double value);
  double BasicGet_1int(int i) const;
  void BasicSet_1int(int i, double value);
  void BasicIncrement_1int(int i, double value);

  void Const(double value);
  MatrixType GetType() const;
  std::ostream& WriteData(std::ostream& c) const;
  std::istream& ReadData(std::istream& c);

  void AssignVM(const VirtualMatrix& A);
  Vect3& operator=(const Vect3& A); // assignment operator
  Vect3& operator=(const VirtualMatrix& A); // overloaded =
  Vect3& operator*=(double b);
  Vect3& operator+=(const Vect3& A);
  Vect3& operator-=(const Vect3& A);
  Vect3& operator-();

  friend Matrix T(const Vect3& A);  // a wasteful transpose
  friend Mat3x3 CrossMat(Vect3& a);  // a wasteful cross matrix implementation

  friend void Set6DAngularVector(Vect6& v6, Vect3& v3);
  friend void Set6DLinearVector(Vect6& v6, Vect3& v3);

  // fast matrix functions
  friend void FastAssign(Vect3& a, Vect3& c);
  friend void FastSimpleRotation(Vect3& v, double q, Mat3x3& d);
  friend void FastCross(Vect3& a, Vect3& b, Vect3& c); // cross product axb = c
  friend void FastTripleSum(Vect3& a, Vect3& b, Vect3& c, Vect3& d);
  friend void FastTripleSumPPM(Vect3& a, Vect3& b, Vect3& c, Vect3& d);
  friend void FastMult(Mat3x3& A, Vect3& B, Vect3& C);
  friend void FastTMult(Mat3x3& A, Vect3& B, Vect3& C);
  friend void FastNegMult(Mat3x3& A, Vect3& B, Vect3& C);
  friend void FastNegTMult(Mat3x3& A, Vect3& B, Vect3& C);
  friend void FastMult(double a, Vect3& B, Vect3& C);
  friend void FastAdd(Vect3& A, Vect3& B, Vect3& C);
  friend void FastSubt(Vect3& A, Vect3& B, Vect3& C);
  friend void OnPopulateSVect(Vect3& angular, Vect3& linear, Vect6& sV);
  friend void OnPopulateSC(Vect3& gamma, Mat3x3& C, Mat6x6& SC);

  friend void FastMult(Mat3x3& A, ColMatrix& B, Vect3& C);
  friend void FastAssign(ColMatrix&A, Vect3& C);
  friend void FastMult(Mat3x3& A, Vect3& B, ColMatrix& C);
};
}
#endif
