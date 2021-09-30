/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: body.h                                                  *
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


#ifndef BODY_H
#define BODY_H

#include "poemslist.h"
#include <iostream>
#include "poemsobject.h"

#include "matrices.h"



// emumerated type
enum BodyType {
  INERTIALFRAME = 0,
  PARTICLE = 1,
  RIGIDBODY = 2
};

class Point;
class Joint;
class CompBody;

class Body : public POEMSObject {
public:
  double mass;
  Mat3x3 inertia;

  Vect3 r;
  Vect3 v;
  Vect3 v_k;
  Vect3 a;
  Vect3 a_t;
  Mat3x3 n_C_k;
  Vect3 omega;
  Vect3 omega_k;
  Vect3 alpha;
  Vect3 alpha_t;
  double KE;


  List<Joint> joints;
  List<Point> points;

  Body();

  bool ReadIn(std::istream& in);
  void WriteOut(std::ostream& out);
  bool ReadInPoints(std::istream& in);
  void WriteOutPoints(std::ostream& out);
  Point* GetPoint(int p);
  void AddJoint(Joint* joint);
  void AddPoint(Point* point);

  virtual bool ReadInBodyData(std::istream& in) = 0;
  virtual void WriteOutBodyData(std::ostream& out) = 0;
  virtual ~Body();
  virtual BodyType GetType() = 0;
};

// global body functions
Body* NewBody(int type);

#endif
