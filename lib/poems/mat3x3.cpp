/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: mat3x3.cpp                                              *
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


#include "mat3x3.h"
#include <cstdlib>

using namespace std;

Mat3x3::Mat3x3(){
  numrows = numcols = 3;
}
Mat3x3::~Mat3x3(){
}

Mat3x3::Mat3x3(const Mat3x3& A){
  numrows = numcols = 3;
  elements[0][0] = A.elements[0][0];elements[0][1] = A.elements[0][1];elements[0][2] = A.elements[0][2];
  elements[1][0] = A.elements[1][0];elements[1][1] = A.elements[1][1];elements[1][2] = A.elements[1][2];
  elements[2][0] = A.elements[2][0];elements[2][1] = A.elements[2][1];elements[2][2] = A.elements[2][2];
}

Mat3x3::Mat3x3(const VirtualMatrix& A){
  numrows = numcols = 3;

  // error check
  if( (A.GetNumRows() != 3) || (A.GetNumCols() != 3) ){
    cerr << "illegal matrix size" << endl;
    exit(0);
  }

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      elements[i][j] = A.BasicGet(i,j);
}

double& Mat3x3::operator_2int(int i, int j){ // array access
  return elements[i-1][j-1];
}

double Mat3x3::Get_2int(int i, int j) const{
  return elements[i-1][j-1];
}

void Mat3x3::Set_2int(int i, int j, double value){
  elements[i-1][j-1] = value;
}

double Mat3x3::BasicGet_2int(int i, int j) const{
  return elements[i][j];
}

void Mat3x3::BasicSet_2int(int i, int j, double value){
  elements[i][j] = value;
}

void Mat3x3::BasicIncrement_2int(int i, int j, double value){
  elements[i][j] += value;
}

void Mat3x3::Const(double value){
  elements[0][0] = elements[0][1] = elements[0][2] = value;
  elements[1][0] = elements[1][1] = elements[1][2] = value;
  elements[2][0] = elements[2][1] = elements[2][2] = value;
}

MatrixType Mat3x3::GetType() const{
  return MAT3X3;
}

istream& Mat3x3::ReadData(istream& c){  //input
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      c >> elements[i][j];
  return c;
}

ostream& Mat3x3::WriteData(ostream& c) const{  //output
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      c << elements[i][j] << ' ';
  return c;
}

void Mat3x3::AssignVM(const VirtualMatrix& A){
  // error check
  if( (A.GetNumRows() != 3) || (A.GetNumCols() != 3) ){
    cerr << "illegal matrix size" << endl;
    exit(0);
  }

  for(int i=0;i<numrows;i++)
    for(int j=0;j<numcols;j++)
      elements[i][j] = A.BasicGet(i,j);
}

Mat3x3& Mat3x3::operator=(const Mat3x3& A){ // assignment operator
  elements[0][0] = A.elements[0][0];elements[0][1] = A.elements[0][1];elements[0][2] = A.elements[0][2];
  elements[1][0] = A.elements[1][0];elements[1][1] = A.elements[1][1];elements[1][2] = A.elements[1][2];
  elements[2][0] = A.elements[2][0];elements[2][1] = A.elements[2][1];elements[2][2] = A.elements[2][2];
  return *this;
}

Mat3x3& Mat3x3::operator=(const VirtualMatrix& A){ // overloaded =
  // error check
  if( (A.GetNumRows() != 3) || (A.GetNumCols() != 3) ){
    cerr << "illegal matrix size" << endl;
    exit(0);
  }

  for(int i=0;i<numrows;i++)
    for(int j=0;j<numcols;j++)
      elements[i][j] = A.BasicGet(i,j);
  return *this;
}

Mat3x3& Mat3x3::operator*=(double b){
  elements[0][0] *= b; elements[0][1] *= b; elements[0][2] *= b;
  elements[1][0] *= b; elements[1][1] *= b; elements[1][2] *= b;
  elements[2][0] *= b; elements[2][1] *= b; elements[2][2] *= b;
  return *this;
}

Mat3x3& Mat3x3::Identity(){
  elements[0][0] = 1.0;
  elements[0][1] = 0.0;
  elements[0][2] = 0.0;

  elements[1][0] = 0.0;
  elements[1][1] = 1.0;
  elements[1][2] = 0.0;

  elements[2][0] = 0.0;
  elements[2][1] = 0.0;
  elements[2][2] = 1.0;

  return *this;
}

