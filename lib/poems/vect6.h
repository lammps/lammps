/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: vect6.h                                                 *
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

#ifndef VECT6_H
#define VECT6_H

#include "virtualcolmatrix.h"

class Matrix;
class Mat6x6;
class ColMatrix;
class Vect3;

class Vect6 : public VirtualColMatrix  {
  double elements[6];
public:
  Vect6();
  ~Vect6();
  Vect6(const Vect6& A);  // copy constructor
  Vect6(const VirtualMatrix& A);  // copy constructor

  double& operator_1int (int i); // array access
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
  Vect6& operator=(const Vect6& A); // assignment operator
  Vect6& operator=(const VirtualMatrix& A); // overloaded =
  Vect6& operator*=(double b);

  friend Matrix T(const Vect6& A);  // a wasteful transpose

  friend void Set6DAngularVector(Vect6& v6, Vect3& v3);
  friend void Set6DLinearVector(Vect6& v6, Vect3& v3);

  // fast matrix operations
  friend void FastAdd(Vect6& A, Vect6& B, Vect6& C);
  friend void FastSubt(Vect6& A, Vect6& B, Vect6& C);
  friend void FastMult(Mat6x6& A, Vect6& B, Vect6& C);
  friend void FastMult(Matrix& A, ColMatrix& B, Vect6& C);
  friend void FastTMult(Mat6x6& A, Vect6& B, Vect6& C);
  friend void FastTMult(Matrix& A, Vect6& B, ColMatrix& C);
  friend void FastLDLTSubs(Mat6x6& LD, Vect6& B, Vect6& C);

  friend void OnPopulateSVect(Vect3& angular, Vect3& linear, Vect6& sV);
};

#endif
