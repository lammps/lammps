/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: vect4.h                                                 *
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

#ifndef VECT4_H
#define VECT4_H

#include <iostream>

#include "virtualcolmatrix.h"
#include "virtualmatrix.h"

class Matrix;
class Mat4x4;

class Vect4 : public VirtualColMatrix  {
  double elements[4];
public:
  Vect4();
  ~Vect4();
  Vect4(const Vect4& A);  // copy constructor
  Vect4(const VirtualMatrix& A);  // copy constructor

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
  Vect4& operator=(const Vect4& A); // assignment operator
  Vect4& operator=(const VirtualMatrix& A); // overloaded =
  Vect4& operator*=(double b);

  friend Matrix T(const Vect4& A);  // a wasteful transpose
  friend Mat4x4 CrossMat(Vect4& a);  // a wasteful cross matrix implementation

  // fast matrix functions
  friend void FastAssign(Vect4& a, Vect4& c);
  friend void FastSimpleRotation(Vect4& v, double q, Mat4x4& d);
  friend void FastCross(Vect4& a, Vect4& b, Vect4& c); // cross product axb = c
  friend void FastTripleSum(Vect4& a, Vect4& b, Vect4& c, Vect4& d);
  friend void FastTripleSumPPM(Vect4& a, Vect4& b, Vect4& c, Vect4& d);
  friend void FastMult(Mat4x4& A, Vect4& B, Vect4& C);
  friend void FastTMult(Mat4x4& A, Vect4& B, Vect4& C);
  friend void FastNegMult(Mat4x4& A, Vect4& B, Vect4& C);
  friend void FastNegTMult(Mat4x4& A, Vect4& B, Vect4& C);
  friend void FastMult(double a, Vect4& B, Vect4& C);
  friend void FastAdd(Vect4& A, Vect4& B, Vect4& C);
  friend void FastSubt(Vect4& A, Vect4& B, Vect4& C);
};

#endif
