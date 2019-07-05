/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: vect3.cpp                                               *
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

#include "vect3.h"

#include <cstdlib>
#include <iostream>

using namespace std;
using namespace POEMS;


Vect3::Vect3(){
  numrows = 3; numcols = 1;
}
Vect3::~Vect3(){
}

Vect3::Vect3(const Vect3& A){  // copy constructor
  numrows = 3; numcols = 1;

  elements[0] = A.elements[0];
  elements[1] = A.elements[1];
  elements[2] = A.elements[2];
}

Vect3::Vect3(const VirtualMatrix& A){  // copy constructor
  numrows = 3; numcols = 1;

  // error check
  if( (A.GetNumRows() != 3) || (A.GetNumCols() != 1) ){
    cerr << "illegal matrix size" << endl;
    exit(0);
  }

  for(int i=0;i<3;i++)
    elements[i] = A.BasicGet(i,0);
}

double& Vect3::operator_1int (int i){ // array access
	if(i<1 || i>3){
		cerr << "matrix index invalid in operator ()" << endl;
		exit(1);
	}
  return elements[i-1];
}

double Vect3::Get_1int(int i) const{
  return elements[i-1];
}

void Vect3::Set_1int(int i, double value){
  elements[i-1] = value;
}

double Vect3::BasicGet_1int(int i) const{
  return elements[i];
}

void Vect3::BasicSet_1int(int i, double value){
  elements[i] = value;
}

void Vect3::BasicIncrement_1int(int i, double value){
  elements[i] += value;
}

void Vect3::Const(double value){
  elements[0] = value;
  elements[1] = value;
  elements[2] = value;
}

MatrixType Vect3::GetType() const{
  return VECT3;
}

istream& Vect3::ReadData(istream& c){  //input
  for(int i=0;i<3;i++)
    c >> elements[i];
  return c;
}

ostream& Vect3::WriteData(ostream& c) const{  //output
  for(int i=0;i<3;i++)
    c << elements[i] << ' ';
  return c;
}

void Vect3::AssignVM(const VirtualMatrix& A){
  // error check
  if( (A.GetNumRows() != 3) || (A.GetNumCols() != 1) ){
    cerr << "illegal matrix size" << endl;
    exit(0);
  }

  for(int i=0;i<numrows;i++)
    elements[i] = A.BasicGet(i,0);
}

Vect3& Vect3::operator=(const Vect3& A){ // assignment operator
  elements[0] = A.elements[0];
  elements[1] = A.elements[1];
  elements[2] = A.elements[2];
  return *this;
}

Vect3& Vect3::operator=(const VirtualMatrix& A){ // overloaded =
  // error check
  if( (A.GetNumRows() != 3) || (A.GetNumCols() != 1) ){
    cerr << "illegal matrix size" << endl;
    exit(0);
  }

  for(int i=0;i<numrows;i++)
    elements[i] = A.BasicGet(i,0);
  return *this;
}

Vect3& Vect3::operator*=(double b){
  elements[0] *= b;
  elements[1] *= b;
  elements[2] *= b;
  return *this;
}

Vect3& Vect3::operator+=(const Vect3& A){
  elements[0] += A.elements[0];
  elements[1] += A.elements[1];
  elements[2] += A.elements[2];
  return *this;
}

Vect3& Vect3::operator-(){
  elements[0] = -elements[0];
  elements[1] = -elements[1];
  elements[2] = -elements[2];
  return *this;
}

Vect3& Vect3::operator-=(const Vect3& A){
  elements[0] -= A.elements[0];
  elements[1] -= A.elements[1];
  elements[2] -= A.elements[2];
  return *this;
}
