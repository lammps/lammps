/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: revolutejoint.h                                         *
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

#ifndef REVOLUTEJOINT_H
#define REVOLUTEJOINT_H

#include <iostream>
#include "joint.h"
#include "vect3.h"
#include "matrix.h"

namespace POEMS {
class VirtualMatrix;

class RevoluteJoint : public Joint  {
  Vect3 axis_pk; // unit vector in body1 basis
  Vect3 axis_k;  // unit vector in body2 basis
public: 
  RevoluteJoint();
  ~RevoluteJoint();
  JointType GetType();
  void SetAxisK(VirtualMatrix& axis);
  void SetAxisPK(VirtualMatrix& axis);
  bool ReadInJointData(std::istream& in);
  void WriteOutJointData(std::ostream& out);
  Matrix GetForward_sP();
  Matrix GetBackward_sP();
  void UpdateForward_sP( Matrix& sP);
  void UpdateBackward_sP( Matrix& sP);
  void ComputeLocalTransform();
  void ForwardKinematics();
  void BackwardKinematics();
};
}
#endif
