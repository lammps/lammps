/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: virtualcolmatrix.cpp                                    *
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

#include "virtualcolmatrix.h"
#include <iostream>
#include <cstdlib>

using namespace std;

VirtualColMatrix::VirtualColMatrix(){
  numcols = 1;
}

double& VirtualColMatrix::operator_2int(int i, int j){
  if(j!=1){
    cerr << "matrix index invalid in operator ()" << endl;
    exit(1);
  }
  return (*this).operator_1int(i);
}

double VirtualColMatrix::Get_2int(int i, int j) const{
  if(j!=1){
    cerr << "Subscript out of bounds for collumn matrix" << endl;
    exit(1);
  }
  return Get_1int(i);
}

void VirtualColMatrix::Set_2int(int i, int j, double value){
  if(j!=1){
    cerr << "Subscript out of bounds for collumn matrix" << endl;
    exit(1);
  }
  Set_1int(i,value);
}

double VirtualColMatrix::BasicGet_2int(int i, int j) const{
  return BasicGet_1int(i);
}

void VirtualColMatrix::BasicSet_2int(int i, int j, double value){
  BasicSet_1int(i,value);
}

void VirtualColMatrix::BasicIncrement_2int(int i, int j, double value){
  BasicIncrement_1int(i,value);
}

