/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: body23joint.h                                                               *
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

#ifndef BODY23JOINT_H
#define BODY23JOINT_H

#include "joint.h"
#include "vect3.h"
#include "mat3x3.h"



class Body23Joint : public Joint  {
  Matrix const_sP;
public:
  Body23Joint();
  ~Body23Joint();
  JointType GetType();
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

#endif
