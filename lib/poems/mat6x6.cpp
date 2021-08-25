/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: mat6x6.cpp                                              *
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

#include "mat6x6.h"
#include <cstdlib>

using namespace std;

Mat6x6::Mat6x6(){
  numrows = numcols = 6;
}
Mat6x6::~Mat6x6(){
}

Mat6x6::Mat6x6(const Mat6x6& A){
  numrows = numcols = 6;
  int i,j;
  for(i=0;i<6;i++)
    for(j=0;j<6;j++)
      elements[i][j] = A.elements[i][j];
}

Mat6x6::Mat6x6(const VirtualMatrix& A){
  numrows = numcols = 6;

  // error check
  if( (A.GetNumRows() != 6) || (A.GetNumCols() != 6) ){
    cerr << "illegal matrix size" << endl;
    exit(0);
  }

  for(int i=0;i<6;i++)
    for(int j=0;j<6;j++)
      elements[i][j] = A.BasicGet(i,j);
}

double& Mat6x6::operator_2int(int i, int j){ // array access
  return elements[i-1][j-1];
}

double Mat6x6::Get_2int(int i, int j) const{
  return elements[i-1][j-1];
}

void Mat6x6::Set_2int(int i, int j, double value){
  elements[i-1][j-1] = value;
}

double Mat6x6::BasicGet_2int(int i, int j) const{
  return elements[i][j];
}

void Mat6x6::BasicSet_2int(int i, int j, double value){
  elements[i][j] = value;
}

void Mat6x6::BasicIncrement_2int(int i, int j, double value){
  elements[i][j] += value;
}

void Mat6x6::Const(double value){
  int i,j;
  for(i=0;i<6;i++)
    for(j=0;j<6;j++)
      elements[i][j] = value;
}

MatrixType Mat6x6::GetType() const{
  return MAT6X6;
}

istream& Mat6x6::ReadData(istream& c){  //input
  for(int i=0;i<6;i++)
    for(int j=0;j<6;j++)
      c >> elements[i][j];
  return c;
}

ostream& Mat6x6::WriteData(ostream& c) const{  //output
  for(int i=0;i<6;i++)
    for(int j=0;j<6;j++)
      c << elements[i][j] << ' ';
  return c;
}

void Mat6x6::AssignVM(const VirtualMatrix& A){
  // error check
  if( (A.GetNumRows() != 6) || (A.GetNumCols() != 6) ){
    cerr << "illegal matrix size" << endl;
    exit(0);
  }

  for(int i=0;i<6;i++)
    for(int j=0;j<6;j++)
      elements[i][j] = A.BasicGet(i,j);
}

Mat6x6& Mat6x6::operator=(const Mat6x6& A){ // assignment operator
  int i,j;
  for(i=0;i<6;i++)
    for(j=0;j<6;j++)
      elements[i][j] = A.elements[i][j];
  return *this;
}

Mat6x6& Mat6x6::operator=(const VirtualMatrix& A){ // overloaded =
  // error check
  if( (A.GetNumRows() != 6) || (A.GetNumCols() != 6) ){
    cerr << "illegal matrix size" << endl;
    exit(0);
  }

  for(int i=0;i<6;i++)
    for(int j=0;j<6;j++)
      elements[i][j] = A.BasicGet(i,j);
  return *this;
}

Mat6x6& Mat6x6::operator*=(double b){
  int i,j;
  for(i=0;i<6;i++)
    for(j=0;j<6;j++)
      elements[i][j] *= b;
  return *this;
}

Mat6x6& Mat6x6::Identity(){
  int i,j;
  for(i=0;i<5;i++){
    elements[i][i] = 1.0;
    for(j=i+1;j<6;j++){
      elements[i][j] = 0.0;
      elements[j][i] = 0.0;
    }
  }
  elements[i][i] = 1.0;
  return *this;
}

