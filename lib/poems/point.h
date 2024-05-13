/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: point.h                                                  *
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

#ifndef POINT_H
#define POINT_H

#include <iostream>
#include "poemsobject.h"
#include "vect3.h"


// emumerated type
enum PointType {
  FIXEDPOINT = 0
};

class Point : public POEMSObject {
public:
  Vect3 position;

  Point();
  bool ReadIn(std::istream& in);
  void WriteOut(std::ostream& out);

  virtual ~Point();
  virtual PointType GetType() = 0;
  virtual Vect3 GetPoint() = 0;
  virtual bool ReadInPointData(std::istream& in) = 0;
  virtual void WriteOutPointData(std::ostream& out) = 0;
};

// global point functions
Point* NewPoint(int type);

#endif
