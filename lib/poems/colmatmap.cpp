/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: colmatmap.cpp                                           *
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

#include "colmatmap.h"
#include "colmatrix.h"
#include <iostream>
#include <cstdlib>

using namespace std;

ColMatMap::ColMatMap(){
	numrows = 0;
	elements = 0;
}

ColMatMap::~ColMatMap(){
	delete [] elements;
}

ColMatMap::ColMatMap(const ColMatMap& A){  // copy constructor
	numrows = 0;
	elements = 0;
	Dim(A.numrows);
	for(int i=0;i<numrows;i++)
		elements[i] = A.elements[i];
}

ColMatMap::ColMatMap(ColMatrix& A){  // copy constructor
	numrows = 0;
	elements = 0;
	Dim(A.GetNumRows());
	for(int i=0;i<numrows;i++)
		elements[i] = A.GetElementPointer(i);
}

/*
ColMatrix::ColMatrix(const VirtualMatrix& A){  // copy constructor
  if( A.GetNumCols() != 1 ){
    cerr << "error trying to write a 2D matrix to a collumn" << endl;
    exit(1);
  }
  numrows = 0;
  elements = 0;
  Dim(A.GetNumRows());
  for(int i=0;i<numrows;i++)
    elements[i] = A.BasicGet(i,0);
}
*/

ColMatMap::ColMatMap(int m){  // size constructor
	numrows = 0;
	elements = 0;
	Dim(m);
}

void ColMatMap::Dim(int m){
	delete [] elements;
	numrows = m;
	elements = new double* [m];
}

void ColMatMap::Const(double value){
	for(int i=0;i<numrows;i++)
		*(elements[i]) = value;
}

MatrixType ColMatMap::GetType() const{
	return COLMATMAP;
}

ostream& ColMatMap::WriteData(ostream& c) const{  //output
	c << numrows << ' ';
	for(int i=0;i<numrows;i++)
		c << *(elements[i]) << ' ';
	return c;
}

double& ColMatMap::operator_1int (int i){ // array access
	if((i>numrows) || (i<1)){
		cerr << "matrix index invalid in operator ()" << endl;
		exit(1);
	}
	return *(elements[i-1]);
}

double ColMatMap::Get_1int(int i) const{
	if((i>numrows) || (i<1)){
		cerr << "matrix index exceeded in Get" << endl;
		exit(1);
	}
	return *(elements[i-1]);
}

void ColMatMap::Set_1int(int i, double value){
	if((i>numrows) || (i<1)){
		cerr << "matrix index exceeded in Set" << endl;
		exit(1);
	}
	*(elements[i-1]) = value;
}

double ColMatMap::BasicGet_1int(int i) const{
	return *(elements[i]);
}

void ColMatMap::SetElementPointer(int i, double* p){
	elements[i] = p;
}

double* ColMatMap::GetElementPointer(int i){
	return elements[i];
}

void ColMatMap::BasicSet_1int(int i, double value){
	*(elements[i]) = value;
}

void ColMatMap::BasicIncrement_1int(int i, double value){
	*(elements[i]) += value;
}

void ColMatMap::AssignVM(const VirtualMatrix& A){
	if(A.GetNumRows() != numrows){
		cerr << "dimension mismatch in ColMatMap assignment" << endl;
		exit(0);
	}
	if( A.GetNumCols() != 1 ){
		cerr << "error trying to write a 2D matrix to a collumn" << endl;
		exit(1);
	}
	for(int i=0;i<numrows;i++)
		*(elements[i]) = A.BasicGet(i,0);
}

ColMatMap& ColMatMap::operator=(const ColMatMap& A){ // assignment operator
	if(A.numrows != numrows){
		cerr << "dimension mismatch in ColMatMap assignment" << endl;
		exit(0);
	}
	for(int i=0;i<numrows;i++)
		*(elements[i]) = *(A.elements[i]);
	return *this;
}

ColMatMap& ColMatMap::operator=(const ColMatrix& A){ // assignment operator
	if(A.GetNumRows() != numrows){
		cerr << "dimension mismatch in ColMatMap assignment" << endl;
		exit(0);
	}
	for(int i=0;i<numrows;i++)
		*(elements[i]) = A.BasicGet(i);
	return *this;
}

ColMatMap& ColMatMap::operator=(const VirtualColMatrix& A){ // overloaded =
	if(A.GetNumRows() != numrows){
		cerr << "dimension mismatch in ColMatMap assignment" << endl;
		exit(0);
	}
	for(int i=0;i<numrows;i++)
		*(elements[i]) = A.BasicGet(i);
	return *this;
}

ColMatMap& ColMatMap::operator=(const VirtualMatrix& A){ // overloaded =
	if(A.GetNumRows() != numrows){
		cerr << "dimension mismatch in ColMatMap assignment" << endl;
		exit(0);
	}
	if( A.GetNumCols() != 1 ){
		cerr << "error trying to write a 2D matrix to a collumn" << endl;
		exit(1);
	}
	for(int i=0;i<numrows;i++)
		*(elements[i]) = A.BasicGet(i,0);
	return *this;
}

ColMatMap& ColMatMap::operator*=(double b){
	for(int i=0;i<numrows;i++)
		*(elements[i]) *= b;
	return *this;
}

