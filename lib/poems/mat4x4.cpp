/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: mat4x4.cpp                                              *
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

#include "mat4x4.h"
#include <cstdlib>

using namespace std;

Mat4x4::Mat4x4(){
  numrows = numcols = 4;
}
Mat4x4::~Mat4x4(){
}

Mat4x4::Mat4x4(const Mat4x4& A){
  numrows = numcols = 4;
  elements[0][0] = A.elements[0][0];elements[0][1] = A.elements[0][1];elements[0][2] = A.elements[0][2]; elements[0][3] = A.elements[0][3];
  elements[1][0] = A.elements[1][0];elements[1][1] = A.elements[1][1];elements[1][2] = A.elements[1][2]; elements[1][3] = A.elements[1][3];
  elements[2][0] = A.elements[2][0];elements[2][1] = A.elements[2][1];elements[2][2] = A.elements[2][2]; elements[2][3] = A.elements[2][3];
  elements[3][0] = A.elements[3][0];elements[3][1] = A.elements[3][1];elements[3][2] = A.elements[3][2]; elements[3][3] = A.elements[3][3];
}

Mat4x4::Mat4x4(const VirtualMatrix& A){
  numrows = numcols = 4;

  // error check
  if( (A.GetNumRows() != 4) || (A.GetNumCols() != 4) ){
    cerr << "illegal matrix size" << endl;
    exit(0);
  }

  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++)
      elements[i][j] = A.BasicGet(i,j);
}

double& Mat4x4::operator_2int(int i, int j){ // array access
  return elements[i-1][j-1];
}

double Mat4x4::Get_2int(int i, int j) const{
  return elements[i-1][j-1];
}

void Mat4x4::Set_2int(int i, int j, double value){
  elements[i-1][j-1] = value;
}

double Mat4x4::BasicGet_2int(int i, int j) const{
  return elements[i][j];
}

void Mat4x4::BasicSet_2int(int i, int j, double value){
  elements[i][j] = value;
}

void Mat4x4::BasicIncrement_2int(int i, int j, double value){
  elements[i][j] += value;
}

void Mat4x4::Const(double value){
  elements[0][0] = elements[0][1] = elements[0][2] = elements[0][3] = value;
  elements[1][0] = elements[1][1] = elements[1][2] = elements[1][3] = value;
  elements[2][0] = elements[2][1] = elements[2][2] = elements[2][3] = value;
  elements[3][0] = elements[3][1] = elements[3][2] = elements[3][3] = value;
}

MatrixType Mat4x4::GetType() const{
  return MAT4X4;
}

istream& Mat4x4::ReadData(istream& c){  //input
  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++)
      c >> elements[i][j];
  return c;
}

ostream& Mat4x4::WriteData(ostream& c) const{  //output
  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++)
      c << elements[i][j] << ' ';
  return c;
}

void Mat4x4::AssignVM(const VirtualMatrix& A){
  // error check
  if( (A.GetNumRows() != 4) || (A.GetNumCols() != 4) ){
    cerr << "illegal matrix size" << endl;
    exit(0);
  }

  for(int i=0;i<numrows;i++)
    for(int j=0;j<numcols;j++)
      elements[i][j] = A.BasicGet(i,j);
}

Mat4x4& Mat4x4::operator=(const Mat4x4& A){ // assignment operator
  elements[0][0] = A.elements[0][0];elements[0][1] = A.elements[0][1];elements[0][2] = A.elements[0][2]; elements[0][3] = A.elements[0][3];
  elements[1][0] = A.elements[1][0];elements[1][1] = A.elements[1][1];elements[1][2] = A.elements[1][2]; elements[1][3] = A.elements[1][3];
  elements[2][0] = A.elements[2][0];elements[2][1] = A.elements[2][1];elements[2][2] = A.elements[2][2]; elements[2][3] = A.elements[2][3];
  elements[3][0] = A.elements[3][0];elements[3][1] = A.elements[3][1];elements[3][2] = A.elements[3][2]; elements[3][3] = A.elements[3][3];
  return *this;
}

Mat4x4& Mat4x4::operator=(const VirtualMatrix& A){ // overloaded =
  // error check
  if( (A.GetNumRows() != 4) || (A.GetNumCols() != 4) ){
    cerr << "illegal matrix size" << endl;
    exit(0);
  }

  for(int i=0;i<numrows;i++)
    for(int j=0;j<numcols;j++)
      elements[i][j] = A.BasicGet(i,j);
  return *this;
}

Mat4x4& Mat4x4::operator*=(double b){
  elements[0][0] *= b; elements[0][1] *= b; elements[0][2] *= b; elements[0][3] *= b;
  elements[1][0] *= b; elements[1][1] *= b; elements[1][2] *= b; elements[1][3] *= b;
  elements[2][0] *= b; elements[2][1] *= b; elements[2][2] *= b; elements[2][3] *= b;
  elements[3][0] *= b; elements[3][1] *= b; elements[3][2] *= b; elements[3][3] *= b;
  return *this;
}

Mat4x4& Mat4x4::Identity(){
  elements[0][0] = 1.0;
  elements[0][1] = 0.0;
  elements[0][2] = 0.0;
  elements[0][3] = 0.0;

  elements[1][0] = 0.0;
  elements[1][1] = 1.0;
  elements[1][2] = 0.0;
  elements[1][3] = 0.0;

  elements[2][0] = 0.0;
  elements[2][1] = 0.0;
  elements[2][2] = 1.0;
  elements[2][3] = 0.0;

  elements[3][0] = 0.0;
  elements[3][1] = 0.0;
  elements[3][2] = 0.0;
  elements[3][3] = 1.0;

  return *this;
}

