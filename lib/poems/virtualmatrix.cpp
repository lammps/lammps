/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: virtualmatrix.cpp                                       *
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


#include "virtualmatrix.h"
#include "matrixfun.h"
#include <cstdlib>

using namespace std;

VirtualMatrix::VirtualMatrix(){
  numrows = numcols = 0;
}

VirtualMatrix::~VirtualMatrix(){
}

int VirtualMatrix::GetNumRows() const {
  return numrows;
}

int VirtualMatrix::GetNumCols() const {
  return numcols;
}

double& VirtualMatrix::operator() (int i, int j){ // array access
	return operator_2int(i,j);
}

double VirtualMatrix::Get(int i, int j) const{
	return Get_2int(i,j);
}

void VirtualMatrix::Set(int i, int j, double value){
	Set_2int(i,j,value);
}

double VirtualMatrix::BasicGet(int i, int j) const{
	return BasicGet_2int(i,j);
}

void VirtualMatrix::BasicSet(int i, int j, double value){
	BasicSet_2int(i,j,value);
}

void VirtualMatrix::BasicIncrement(int i, int j, double value){
	BasicIncrement_2int(i,j,value);
}

double& VirtualMatrix::operator() (int i){
	return operator_1int(i);
}

double VirtualMatrix::Get(int i) const{
	return Get_1int(i);
}

void VirtualMatrix::Set(int i, double value){
	Set_1int(i, value);
}

double VirtualMatrix::BasicGet(int i) const{
	return BasicGet_1int(i);
}

void VirtualMatrix::BasicSet(int i, double value){
	BasicSet_1int(i, value);
}

void VirtualMatrix::BasicIncrement(int i, double value){
	BasicIncrement_1int(i, value);
}

double& VirtualMatrix::operator_1int (int i) {
	cerr << "Error: single dimensional access is not defined for matrices of type " << GetType() << endl;
	exit(0);
	return *(new double);
}

double VirtualMatrix::Get_1int(int i) const {
	cerr << "Error: single dimensional access is not defined for matrices of type " << GetType() << endl;
	exit(0);
	return 0.0;
}

void VirtualMatrix::Set_1int(int i, double value){
	cerr << "Error: single dimensional access is not defined for matrices of type " << GetType() << endl;
	exit(0);
}

double VirtualMatrix::BasicGet_1int(int i) const {
	cerr << "Error: single dimensional access is not defined for matrices of type " << GetType() << endl;
	exit(0);
	return 0.0;
}

void VirtualMatrix::BasicSet_1int(int i, double value) {
	cerr << "Error: single dimensional access is not defined for matrices of type " << GetType() << endl;
	exit(0);
}

void VirtualMatrix::BasicIncrement_1int(int i, double value){
	cerr << "Error: single dimensional access is not defined for matrices of type " << GetType() << endl;
	exit(0);
}

void VirtualMatrix::Zeros(){
  Const(0.0);
}

void VirtualMatrix::Ones(){
  Const(1.0);
}

ostream& VirtualMatrix::WriteData(ostream& c) const {
  cerr << "Error: no output definition for matrices of type " << GetType() << endl;
  exit(0);
}

istream& VirtualMatrix::ReadData(istream& c){
	cerr << "Error: no input definition for matrices of type " << GetType() << endl;
	exit(0);
}

//
// operators and functions
//

ostream& operator<< (ostream& c, const VirtualMatrix& A){  //output
  c << A.GetType() << ' ';
  A.WriteData(c);
  c << endl;
  return c;
}

istream& operator>> (istream& c, VirtualMatrix& A){  //input
  VirtualMatrix* vm;
  int matrixtype;
  c >> matrixtype;

  if( MatrixType(matrixtype) == A.GetType() ) A.ReadData(c);
  else{
    // issue a warning?
    cerr << "Warning: During matrix read expected type " << A.GetType() << " and got type " << matrixtype << endl;
    vm = NewMatrix(matrixtype);
    if(!vm){
      cerr << "Error: unable to instantiate matrix of type " << matrixtype << endl;
      exit(0);
    }
    vm->ReadData(c);
    A.AssignVM(*vm);
    delete vm;
  }

  return c;
}

