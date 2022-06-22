/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: virtualrowmatrix.cpp                                    *
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


#include "virtualrowmatrix.h"
#include <iostream>
#include <cstdlib>

using namespace std;

VirtualRowMatrix::VirtualRowMatrix(){
  numrows = 1;
}

double& VirtualRowMatrix::operator_2int (int i, int j){
  if(i!=1){
    cerr << "matrix index invalid in operator ()" << endl;
    exit(1);
  }
  return (*this).operator_1int(j);
}

double VirtualRowMatrix::Get_2int(int i, int j) const{
  if(i!=1){
    cerr << "Subscript out of bounds for row matrix" << endl;
    exit(1);
  }
  return Get_1int(j);
}

void VirtualRowMatrix::Set_2int(int i, int j, double value){
  if(i!=1){
    cerr << "Subscript out of bounds for row matrix" << endl;
    exit(1);
  }
  Set_1int(j,value);
}

double VirtualRowMatrix::BasicGet_2int(int i, int j) const{
  return BasicGet_1int(j);
}

void VirtualRowMatrix::BasicSet_2int(int i, int j, double value){
  BasicSet_1int(j,value);
}

void VirtualRowMatrix::BasicIncrement_2int(int i, int j, double value){
  BasicIncrement_1int(j,value);
}


