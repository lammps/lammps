/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: colmatrix.h                                             *
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


#ifndef COLMATRIX_H
#define COLMATRIX_H

#include <iostream>
#include "virtualcolmatrix.h"
#include "virtualmatrix.h"

class Matrix;
class Vect6;
class Mat3x3;
class Vect3;

class ColMatrix : public VirtualColMatrix  {
  double* elements;
public:
  ColMatrix();
  ~ColMatrix();
  ColMatrix(const ColMatrix& A);  // copy constructor
  ColMatrix(const VirtualColMatrix& A);  // copy constructor
  ColMatrix(const VirtualMatrix& A);  // copy constructor
  ColMatrix(int m);  // size constructor

  double& operator_1int (int i); // array access
  double Get_1int(int i) const;
  void Set_1int(int i, double value);
  double BasicGet_1int(int i) const;
  void BasicSet_1int(int i, double value);
  void BasicIncrement_1int(int i, double value);

  double* GetElementPointer(int i);
  void Dim(int m);
  void Const(double value);
  MatrixType GetType() const;
  std::istream& ReadData(std::istream& c);
  std::ostream& WriteData(std::ostream& c) const;

  void AssignVM(const VirtualMatrix& A);
  ColMatrix& operator=(const ColMatrix& A); // assignment operator
  ColMatrix& operator=(const VirtualColMatrix& A); // overloaded =
  ColMatrix& operator=(const VirtualMatrix& A); // overloaded =
  ColMatrix& operator*=(double b);

        void Abs();
        void BasicMax(double& value, int& index);
        void BasicMin(double& value, int& index);

  // fast matrix operations
		  friend void FastQuaternions(ColMatrix& q, Mat3x3& C);
		  friend void FastInvQuaternions(Mat3x3& C, ColMatrix& q);
		  friend void FastQuaternionDerivatives(ColMatrix& q, ColMatrix& omega, ColMatrix& qdot);
		  friend void FastTMult(Matrix& A, Vect6& B, ColMatrix& C);
		  friend void FastMult(Matrix& A, ColMatrix& B, Vect6& C);
		  friend void FastAssign(ColMatrix& A, ColMatrix& C);

		  friend void FastMult(Mat3x3& A, ColMatrix& B, Vect3& C);
		  friend void FastMult(Mat3x3& A, Vect3& B, ColMatrix& C);
		  friend void FastAssign(ColMatrix&A, Vect3& C);
  
		  friend void EP_Derivatives(ColMatrix& q, ColMatrix& u, ColMatrix& qdot);
		  friend void EP_Transformation(ColMatrix& q, Mat3x3& C);
		  friend void EP_FromTransformation(ColMatrix& q, Mat3x3& C);
		  friend void EP_Normalize(ColMatrix& q);
		  friend void EPdotdot_udot(ColMatrix& Audot, ColMatrix& Aqdot, ColMatrix& Aq,ColMatrix& Aqddot);
		  friend void qdot_to_u(ColMatrix& q, ColMatrix& u, ColMatrix& qdot);
};

#endif
