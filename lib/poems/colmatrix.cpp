/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: colmatrix.cpp                                           *
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

#include "colmatrix.h"

#include <iostream>
#include <cstdlib>

using namespace std;
using namespace POEMS;


ColMatrix::ColMatrix(){
	numrows = 0;
	elements = 0;
}

ColMatrix::~ColMatrix(){
	delete [] elements;
}

ColMatrix::ColMatrix(const ColMatrix& A){  // copy constructor
	numrows = 0;
	elements = 0;
	Dim(A.numrows);
	for(int i=0;i<numrows;i++)
		elements[i] = A.elements[i];
}

ColMatrix::ColMatrix(const VirtualColMatrix& A){  // copy constructor
	numrows = 0;
	elements = 0;
	Dim(A.GetNumRows());
	for(int i=0;i<numrows;i++)
		elements[i] = A.BasicGet(i);
}

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

ColMatrix::ColMatrix(int m){  // size constructor
	numrows = 0;
	elements = 0;
	Dim(m);
}

void ColMatrix::Dim(int m){
	delete [] elements;
	numrows = m;
	elements = new double [m];
}

void ColMatrix::Const(double value){
	for(int i=0;i<numrows;i++)
		elements[i] = value;
}

MatrixType ColMatrix::GetType() const{
	return COLMATRIX;
}

istream& ColMatrix::ReadData(istream& c){
	int n;
	c >> n;
	Dim(n);
	for(int i=0;i<n;i++)
		c >> elements[i];
	return c;
}

ostream& ColMatrix::WriteData(ostream& c) const{  //output
	c << numrows << ' ';
	for(int i=0;i<numrows;i++)
		c << elements[i] << ' ';
	return c;
}

double& ColMatrix::operator_1int (int i){ // array access
	if((i>numrows) || (i<1)){
		cerr << "matrix index invalid in operator ()" << endl;
		exit(1);
	}
	return elements[i-1];
}

double ColMatrix::Get_1int(int i) const{
	if((i>numrows) || (i<1)){
		cerr << "matrix index exceeded in Get" << endl;
		exit(1);
	}
	return elements[i-1];
}

void ColMatrix::Set_1int(int i, double value){
	if((i>numrows) || (i<1)){
		cerr << "matrix index exceeded in Set" << endl;
		exit(1);
	}
	elements[i-1] = value;
}

double ColMatrix::BasicGet_1int(int i) const{
	return elements[i];
}

double* ColMatrix::GetElementPointer(int i){
	return &(elements[i]);
}

void ColMatrix::BasicSet_1int(int i, double value){
	elements[i] = value;
}

void ColMatrix::BasicIncrement_1int(int i, double value){
	elements[i] += value;
}

void ColMatrix::AssignVM(const VirtualMatrix& A){
	if( A.GetNumCols() != 1 ){
		cerr << "error trying to write a 2D matrix to a collumn" << endl;
		exit(1);
	}
	Dim( A.GetNumRows() );
	for(int i=0;i<numrows;i++)
		elements[i] = A.BasicGet(i,0);
}

ColMatrix& ColMatrix::operator=(const ColMatrix& A){ // assignment operator
	Dim(A.numrows);
	for(int i=0;i<numrows;i++)
		elements[i] = A.elements[i];
	return *this;
}

ColMatrix& ColMatrix::operator=(const VirtualColMatrix& A){ // overloaded =
	Dim( A.GetNumRows() );
	for(int i=0;i<numrows;i++)
		elements[i] = A.BasicGet(i);
	return *this;
}

ColMatrix& ColMatrix::operator=(const VirtualMatrix& A){ // overloaded =
	if( A.GetNumCols() != 1 ){
		cerr << "error trying to write a 2D matrix to a collumn" << endl;
		exit(1);
	}
	Dim( A.GetNumRows() );
	for(int i=0;i<numrows;i++)
		elements[i] = A.BasicGet(i,0);
	return *this;
}

ColMatrix& ColMatrix::operator*=(double b){
	for(int i=0;i<numrows;i++)
		elements[i] *= b;
	return *this;
}

void ColMatrix::Abs(){
	for(int i=0;i<numrows;i++)
		elements[i] = std::abs(elements[i]);
}

void ColMatrix::BasicMax(double& value, int& index){
	value = elements[0];
	index = 0;
	for(int i=1;i<numrows;i++)
		if(elements[i] > value){
			value = elements[i];
			index = i;
		}
}

void ColMatrix::BasicMin(double& value, int& index){
	value = elements[0];
	index = 0;
	for(int i=1;i<numrows;i++)
		if(elements[i] < value){
			value = elements[i];
			index = i;
		}
}

