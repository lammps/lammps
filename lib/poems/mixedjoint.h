/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: mixedjoint.h		                                       *
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

#ifndef MIXEDJOINT_H
#define MIXEDJOINT_H

#include "joint.h"


class MixedJoint : public Joint{
  Matrix const_sP;
  int numrots;
  int numtrans;
  Vect6 dofs;
public: 
  MixedJoint();
  ~MixedJoint();
  
  JointType GetType();
  bool ReadInJointData(std::istream& in);
  void WriteOutJointData(std::ostream& out);
  void ComputeLocalTransform();
  void SetsP(Matrix& sPr, Vect6& temp_dofs, int i, int j);
  Matrix GetForward_sP();
  Matrix GetBackward_sP();
  void UpdateForward_sP( Matrix& sP);
  void UpdateBackward_sP( Matrix& sP);
  void ForwardKinematics();
  void BackwardKinematics();
};

#endif
