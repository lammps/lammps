/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: joint.h                                                 *
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

#ifndef JOINT_H
#define JOINT_H

#include "poemsobject.h"
#include <iostream>
#include "matrices.h"

enum JointType {
  XYZJOINT = 0,
  FREEBODYJOINT = 1,
  REVOLUTEJOINT = 2,
  PRISMATICJOINT = 3,
  SPHERICALJOINT = 4,
  BODY23JOINT = 5,
  MIXEDJOINT = 6		  
};

class Body;
class Point;

class Joint : public POEMSObject {
protected:
  Body* body1;
  Body* body2;
  Point* point1;
  Point* point2;

  ColMatrix qo;     // generalized coordinates (initial value)
  ColMatrix uo;     // generalized speeds (initial value)
  ColMatrix qdoto;     // generalized speeds (initial value)

  ColMatrix q;      // generalized coordinates
  ColMatrix u;      // generalized speeds
  ColMatrix qdot;   // generalized coordinate derivatives
  ColMatrix udot;   // generalized speed derivatives
  ColMatrix qdotdot;
  
  Mat3x3 pk_C_ko;   // transformation relationship for q = 0

  Mat3x3 pk_C_k;  // local transform
  Mat3x3 k_C_pk;

  Vect3 r12;
  Vect3 r21;

public:

  Joint();
  virtual ~Joint();
  virtual JointType GetType() = 0;

  void SetBodies(Body* b1, Body* b2);
  void SetPoints(Point* p1, Point* p2);

  int GetBodyID1();
  int GetBodyID2();
  int GetPointID1();
  int GetPointID2();

  ColMatrix* GetQ();
  ColMatrix* GetU();
  ColMatrix* GetQdot();
  ColMatrix* GetUdot();
  ColMatrix* GetQdotdot();
  /*ColMatrix* GetAcc();
  ColMatrix* GetAng();*/

  void DimQandU(int i);
  void DimQandU(int i, int j);

  Body* GetBody1();
  Body* GetBody2();
  Body* OtherBody(Body* body);

  Vect3* GetR12();
  Vect3* GetR21();
  Mat3x3* Get_pkCk();
  Mat3x3* Get_kCpk();

  //void SetInitialState(VirtualMatrix& q, VirtualMatrix& u);  
  void SetInitialState(ColMatrix& q, ColMatrix& u);  
  void SetZeroOrientation(VirtualMatrix& C);
  void ResetQdot();
  void ResetQ();
  bool ReadIn(std::istream& in);  
  void WriteOut(std::ostream& out);
  
  virtual void WriteOutJointData(std::ostream& out) = 0;
  virtual bool ReadInJointData(std::istream& in) = 0;
  virtual Matrix GetForward_sP();
  virtual Matrix GetBackward_sP();
  virtual void UpdateForward_sP(Matrix& sP);
  virtual void UpdateBackward_sP(Matrix& sP);
  virtual void ComputeForwardTransforms();
  virtual void ComputeBackwardTransforms();
  virtual void ComputeLocalTransform()=0;
  virtual void ComputeForwardGlobalTransform();
  virtual void ComputeBackwardGlobalTransform();
  virtual void ForwardKinematics()=0;
  virtual void BackwardKinematics()=0;
};

// global joint functions
Joint* NewJoint(int type);

#endif
