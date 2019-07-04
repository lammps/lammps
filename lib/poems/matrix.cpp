/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: matrix.cpp                                              *
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


#include "matrix.h"
#include <iostream>
#include <cstdlib>

using namespace std;
using namespace POEMS;


Matrix::Matrix(){
  numrows = numcols = 0;
  rows = 0;
  elements = 0;
}

Matrix::~Matrix(){
  delete [] rows;
  delete [] elements;
}

Matrix::Matrix(const Matrix& A){
  numrows = numcols = 0;
  rows = 0;
  elements = 0;
  Dim(A.numrows,A.numcols);
  for(int i=0;i<numrows*numcols;i++)
    elements[i] = A.elements[i];
}

Matrix::Matrix(const VirtualMatrix& A){
  numrows = numcols = 0;
  rows = 0;
  elements = 0;
  Dim(A.GetNumRows(),A.GetNumCols());
  for(int i=0;i<numrows;i++)
    for(int j=0;j<numcols;j++)
      rows[i][j] = A.BasicGet(i,j);
}

Matrix::Matrix(int m, int n){
  numrows = numcols = 0;
  rows = 0;
  elements = 0;
  this->Dim(m,n);
}

double& Matrix::operator_2int(int i, int j){
  if((i>numrows) || (j>numcols) || (i*j==0)){
    cerr << "matrix index exceeded in operator ()" << endl;
    exit(1);
  }
  return rows[i-1][j-1];
}

double Matrix::Get_2int(int i, int j) const {
  if((i>numrows) || (j>numcols) || (i*j==0)){
    cerr << "matrix index exceeded in Get" << endl;
    exit(1);
  }
  return rows[i-1][j-1];
}

void Matrix::Set_2int(int i, int j, double value){
  if((i>numrows) || (j>numcols) || (i*j==0)){
    cerr << "matrix index exceeded in Set" << endl;
    exit(1);
  }
  rows[i-1][j-1] = value;
}

double Matrix::BasicGet_2int(int i, int j) const {
  return rows[i][j];
}

void Matrix::BasicSet_2int(int i, int j, double value){
  rows[i][j] = value;
}

void Matrix::BasicIncrement_2int(int i, int j, double value){
  rows[i][j] += value;
}

void Matrix::Const(double value){
  int num = numrows*numcols;
  for(int i=0;i<num;i++) elements[i] = value;
}

MatrixType Matrix::GetType() const{
  return MATRIX;
}

istream& Matrix::ReadData(istream& c){
  int n,m;
  c >> n >> m;
  Dim(n,m);
  for(int i=0;i<numrows;i++)
    for(int j=0;j<numcols;j++){
      c >> rows[i][j];
    }
  return c;
}

ostream& Matrix::WriteData(ostream& c) const{  //output
  c << numrows << ' ' << numcols << ' ';
  for(int i=0;i<numrows;i++)
    for(int j=0;j<numcols;j++){
      c << rows[i][j] << ' ';
    }
  return c;
}

Matrix& Matrix::Dim(int m, int n){  // allocate size
  numrows = m;
  numcols = n;
  delete [] rows;
  delete [] elements;
  elements = new double[n*m];
  rows = new double*[m];
  for(int i=0;i<m;i++)
    rows[i] = &elements[i*numcols];
  return *this;
}

void Matrix::AssignVM(const VirtualMatrix& A){
  Dim( A.GetNumRows(), A.GetNumCols() );
  for(int i=0;i<numrows;i++)
    for(int j=0;j<numcols;j++)
      rows[i][j] = A.BasicGet(i,j);
}

Matrix& Matrix::operator=(const Matrix& A){
  Dim(A.numrows,A.numcols);
  for(int i=0;i<numrows*numcols;i++)
    elements[i] = A.elements[i];
  return *this;
}

Matrix& Matrix::operator=(const VirtualMatrix& A){
  Dim( A.GetNumRows(), A.GetNumCols() );
  for(int i=0;i<numrows;i++)
    for(int j=0;j<numcols;j++)
      rows[i][j] = A.BasicGet(i,j);
  return *this;
}

Matrix& Matrix::operator*=(double b){
  for(int i=0;i<numrows;i++)
    for(int j=0;j<numcols;j++)
      rows[i][j] *= b;
  return *this;
}

