/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: vect4.cpp                                               *
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

#include "vect4.h"
#include <cstdlib>

using namespace std;

Vect4::Vect4(){
  numrows = 4; numcols = 1;
}
Vect4::~Vect4(){
}

Vect4::Vect4(const Vect4& A){  // copy constructor
  numrows = 4; numcols = 1;

  elements[0] = A.elements[0];
  elements[1] = A.elements[1];
  elements[2] = A.elements[2];
  elements[3] = A.elements[3];
}

Vect4::Vect4(const VirtualMatrix& A){  // copy constructor
  numrows = 4; numcols = 1;

  // error check
  if( (A.GetNumRows() != 4) || (A.GetNumCols() != 1) ){
    cerr << "illegal matrix size" << endl;
    exit(0);
  }

  for(int i=0;i<4;i++)
    elements[i] = A.BasicGet(i,0);
}

double& Vect4::operator_1int (int i){ // array access
  return elements[i-1];
}

double Vect4::Get_1int(int i) const{
  return elements[i-1];
}

void Vect4::Set_1int(int i, double value){
  elements[i-1] = value;
}

double Vect4::BasicGet_1int(int i) const{
  return elements[i];
}

void Vect4::BasicSet_1int(int i, double value){
  elements[i] = value;
}

void Vect4::BasicIncrement_1int(int i, double value){
  elements[i] += value;
}

void Vect4::Const(double value){
  elements[0] = value;
  elements[1] = value;
  elements[2] = value;
  elements[3] = value;
}

MatrixType Vect4::GetType() const{
  return VECT4;
}

istream& Vect4::ReadData(istream& c){  //input
  for(int i=0;i<4;i++)
    c >> elements[i];
  return c;
}

ostream& Vect4::WriteData(ostream& c) const{  //output
  for(int i=0;i<4;i++)
    c << elements[i] << ' ';
  return c;
}

void Vect4::AssignVM(const VirtualMatrix& A){
  // error check
  if( (A.GetNumRows() != 4) || (A.GetNumCols() != 1) ){
    cerr << "illegal matrix size" << endl;
    exit(0);
  }

  for(int i=0;i<numrows;i++)
    elements[i] = A.BasicGet(i,0);
}

Vect4& Vect4::operator=(const Vect4& A){ // assignment operator
  elements[0] = A.elements[0];
  elements[1] = A.elements[1];
  elements[2] = A.elements[2];
  elements[3] = A.elements[3];
  return *this;
}

Vect4& Vect4::operator=(const VirtualMatrix& A){ // overloaded =
  // error check
  if( (A.GetNumRows() != 4) || (A.GetNumCols() != 1) ){
    cerr << "illegal matrix size" << endl;
    exit(0);
  }

  for(int i=0;i<numrows;i++)
    elements[i] = A.BasicGet(i,0);
  return *this;
}

Vect4& Vect4::operator*=(double b){
  elements[0] *= b;
  elements[1] *= b;
  elements[2] *= b;
  elements[3] *= b;
  return *this;
}

