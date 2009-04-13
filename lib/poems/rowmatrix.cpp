/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: rowmatrix.cpp                                           *
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

#include "rowmatrix.h"
#include "colmatrix.h"
#include <iostream>
#include <cstdlib>

using namespace std;

RowMatrix::RowMatrix(){
  numcols = 0;
  elements = 0;
}

RowMatrix::~RowMatrix(){
  delete [] elements;
}

RowMatrix::RowMatrix(const RowMatrix& A){  // copy constructor
  numcols = 0;
  elements = 0;
  Dim(A.numcols);
  for(int i=0;i<numcols;i++)
    elements[i] = A.elements[i];
}

RowMatrix::RowMatrix(const VirtualRowMatrix& A){  // copy constructor
  numcols = 0;
  elements = 0;
  Dim(A.GetNumCols());
  for(int i=0;i<numcols;i++)
    elements[i] = A.BasicGet(i);
}

RowMatrix::RowMatrix(const VirtualMatrix& A){  // copy constructor
  if( A.GetNumRows() != 1 ){
    cerr << "error trying to write a 2D matrix to a collumn" << endl;
    exit(1);
  }
  numcols = 0;
  elements = 0;
  Dim(A.GetNumCols());
  for(int i=0;i<numcols;i++)
    elements[i] = A.BasicGet(i,0);
}

RowMatrix::RowMatrix(int n){  // size constructor
  numcols = 0;
  elements = 0;
  Dim(n);
}

void RowMatrix::Dim(int n){
  delete [] elements;
  numcols = n;
  elements = new double [n];
}

void RowMatrix::Const(double value){
  for(int i=0;i<numcols;i++)
    elements[i] = value;
}

double& RowMatrix::operator_1int (int i){ // array access
  if((i>numcols) || (i<1)){
    cerr << "matrix index invalid in operator ()" << endl;
    exit(1);
  }
  return elements[i-1];
}

double RowMatrix::Get_1int(int i) const{
  if((i>numcols) || (i<1)){
    cerr << "matrix index exceeded in Get" << endl;
    exit(1);
  }
  return elements[i-1];
}

void RowMatrix::Set_1int(int i, double value){
  if((i>numcols) || (i<1)){
    cerr << "matrix index exceeded in Set" << endl;
    exit(1);
  }
  elements[i-1] = value;
}

double RowMatrix::BasicGet_1int(int i) const{
  return elements[i];
}

void RowMatrix::BasicSet_1int(int i, double value){
  elements[i] = value;
}

void RowMatrix::BasicIncrement_1int(int i, double value){
  elements[i] += value;
}

MatrixType RowMatrix::GetType() const{
  return ROWMATRIX;
}

istream& RowMatrix::ReadData(istream& c){
  int n;
  c >> n;
  Dim(n);
  for(int i=0;i<n;i++)
    c >> elements[i];

  return c;
}

ostream& RowMatrix::WriteData(ostream& c) const{  //output
  c << numcols << ' ';
  for(int i=0;i<numcols;i++)
    c << elements[i] << ' ';
  return c;
}

void RowMatrix::AssignVM(const VirtualMatrix& A){
  if( A.GetNumRows() != 1 ){
    cerr << "error trying to write a 2D matrix to a collumn" << endl;
    exit(1);
  }
  Dim( A.GetNumCols() );
  for(int i=0;i<numcols;i++)
    elements[i] = A.BasicGet(0,i);
}

RowMatrix& RowMatrix::operator=(const RowMatrix& A){ // assignment operator
  Dim(A.numcols);
  for(int i=0;i<numcols;i++)
    elements[i] = A.elements[i];
  return *this;
}

RowMatrix& RowMatrix::operator=(const VirtualRowMatrix& A){ // overloaded =
  Dim( A.GetNumCols() );
  for(int i=0;i<numcols;i++)
    elements[i] = A.BasicGet(i);
  return *this;
}

RowMatrix& RowMatrix::operator=(const VirtualMatrix& A){ // overloaded =
  if( A.GetNumRows() != 1 ){
    cerr << "error trying to write a 2D matrix to a collumn" << endl;
    exit(1);
  }
  Dim( A.GetNumCols() );
  for(int i=0;i<numcols;i++)
    elements[i] = A.BasicGet(0,i);
  return *this;
}

RowMatrix& RowMatrix::operator*=(double b){
  for(int i=0;i<numcols;i++)
    elements[i] *= b;
  return *this;
}

