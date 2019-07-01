/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: fixedpoint.cpp                                          *
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


#include "fixedpoint.h"

#include <iostream>
#include <iomanip>

#include "virtualmatrix.h"

using namespace std;

FixedPoint::FixedPoint(){
}
FixedPoint::~FixedPoint(){
}

FixedPoint::FixedPoint(double x, double y, double z){
  position(1) = x;
  position(2) = y;
  position(3) = z;
}

FixedPoint::FixedPoint(Vect3& v){
  position = v;
}

PointType FixedPoint::GetType(){
  return FIXEDPOINT;
}

bool FixedPoint::ReadInPointData(std::istream& in){
  in >> position;
  return true;
}

void FixedPoint::WriteOutPointData(std::ostream& out){  
  out << setprecision(16) << position;
}

Vect3 FixedPoint::GetPoint(){
  return position;
}
