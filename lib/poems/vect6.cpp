/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: vect6.cpp                                               *
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

#include "vect6.h"
#include <cstdlib>

using namespace std;

Vect6::Vect6(){
  numrows = 6; numcols = 1;
}
Vect6::~Vect6(){
}

Vect6::Vect6(const Vect6& A){  // copy constructor
  numrows = 6; numcols = 1;

  elements[0] = A.elements[0];
  elements[1] = A.elements[1];
  elements[2] = A.elements[2];
  elements[3] = A.elements[3];
  elements[4] = A.elements[4];
  elements[5] = A.elements[5];
}

Vect6::Vect6(const VirtualMatrix& A){  // copy constructor
  numrows = 6; numcols = 1;

  // error check
  if( (A.GetNumRows() != 6) || (A.GetNumCols() != 1) ){
    cerr << "illegal matrix size" << endl;
    exit(0);
  }

  for(int i=0;i<6;i++)
    elements[i] = A.BasicGet(i,0);
}

double& Vect6::operator_1int (int i){ // array access
	if(i<1 || i>6){
		cerr << "matrix index invalid in operator ()" << endl;
		exit(1);
	}
	return elements[i-1];
}

double Vect6::Get_1int(int i) const{
  return elements[i-1];
}

void Vect6::Set_1int(int i, double value){
  elements[i-1] = value;
}

double Vect6::BasicGet_1int(int i) const{
  return elements[i];
}

void Vect6::BasicSet_1int(int i, double value){
  elements[i] = value;
}

void Vect6::BasicIncrement_1int(int i, double value){
  elements[i] += value;
}

void Vect6::Const(double value){
  elements[0] = value;
  elements[1] = value;
  elements[2] = value;
  elements[3] = value;
  elements[4] = value;
  elements[5] = value;
}

MatrixType Vect6::GetType() const{
  return VECT6;
}

istream& Vect6::ReadData(istream& c){  //input
  for(int i=0;i<6;i++)
    c >> elements[i];
  return c;
}

ostream& Vect6::WriteData(ostream& c) const{  //output
  for(int i=0;i<6;i++)
    c << elements[i] << ' ';
  return c;
}

void Vect6::AssignVM(const VirtualMatrix& A){
  // error check
  if( (A.GetNumRows() != 6) || (A.GetNumCols() != 1) ){
    cerr << "illegal matrix size" << endl;
    exit(0);
  }

  for(int i=0;i<numrows;i++)
    elements[i] = A.BasicGet(i,0);
}

Vect6& Vect6::operator=(const Vect6& A){ // assignment operator
  elements[0] = A.elements[0];
  elements[1] = A.elements[1];
  elements[2] = A.elements[2];
  elements[3] = A.elements[3];
  elements[4] = A.elements[4];
  elements[5] = A.elements[5];
  return *this;
}

Vect6& Vect6::operator=(const VirtualMatrix& A){ // overloaded =
  // error check
  if( (A.GetNumRows() != 6) || (A.GetNumCols() != 1) ){
    cerr << "illegal matrix size" << endl;
    exit(0);
  }

  for(int i=0;i<numrows;i++)
    elements[i] = A.BasicGet(i,0);
  return *this;
}

Vect6& Vect6::operator*=(double b){
  elements[0] *= b;
  elements[1] *= b;
  elements[2] *= b;
  elements[3] *= b;
  elements[4] *= b;
  elements[5] *= b;
  return *this;
}

