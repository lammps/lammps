/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: inertialframe.h                                         *
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


#ifndef INERTIALFRAME_H
#define INERTIALFRAME_H

#include "body.h"


class InertialFrame : public Body  {
  Vect3 gravity;
public:
  InertialFrame();
  ~InertialFrame();
  BodyType GetType();
  Vect3 GetGravity();
  void SetGravity(Vect3& g);
  bool ReadInBodyData(std::istream& in);
  void WriteOutBodyData(std::ostream& out);
};

#endif
