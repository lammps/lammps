/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: point.cpp                                                  *
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


#include "points.h"

Point::Point(){
  position.Zeros();
}
Point::~Point(){
}

bool Point::ReadIn(std::istream& in){
  return ReadInPointData(in);
}

void Point::WriteOut(std::ostream& out){
  out << GetType() << ' ' << GetName() << ' ';
  WriteOutPointData(out);
}

Point* NewPoint(int type){
  switch( PointType(type) )
  {
    case FIXEDPOINT :  // A Fixed Point
      return new FixedPoint();
    default  :  // error
      return 0;
  }
}
