/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: revolutejoint.cpp                                        *
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

#include "revolutejoint.h"
#include "point.h"
#include "matrixfun.h"
#include "body.h"
#include "fastmatrixops.h"

RevoluteJoint::RevoluteJoint(){
  DimQandU(1);
  Vect3 axis;
  axis.Zeros();
  axis(3) = 1;
  SetAxisPK(axis);
}

RevoluteJoint::~RevoluteJoint(){
}

JointType RevoluteJoint::GetType(){
  return REVOLUTEJOINT;
}

void RevoluteJoint::SetAxisK(VirtualMatrix& axis){
  axis_k = axis;
  axis_pk = pk_C_ko*axis_k;
}

void RevoluteJoint::SetAxisPK(VirtualMatrix& axis){
  axis_pk = axis;
  axis_k = T(pk_C_ko)*axis_pk;
}

bool RevoluteJoint::ReadInJointData(std::istream& in){
  Vect3 axis;
  in >> axis;
  SetAxisPK(axis);
  return true;
}

void RevoluteJoint::WriteOutJointData(std::ostream& out){
  out << axis_pk;
}

Matrix RevoluteJoint::GetForward_sP(){
  Vect3 v_kk;

  // v_kk = axis x r
  FastCross(point2->position,axis_k,v_kk);

  // sP = [axis;v_kk]
  return Stack(axis_k,v_kk);
}

void RevoluteJoint::UpdateForward_sP( Matrix& sP){
  // sP is constant, do nothing.
}

Matrix RevoluteJoint::GetBackward_sP(){
  Vect3 v_kk;

  // v_kk = axis x r
  FastCross(point1->position,axis_pk,v_kk);

  // sP = [axis;v_kk]
  return -Stack(axis_pk,v_kk);;
}

void RevoluteJoint::UpdateBackward_sP( Matrix& sP){
  // sP is constant, do nothing.
}

void RevoluteJoint::ComputeLocalTransform(){
  Mat3x3 ko_C_k;
  FastSimpleRotation(axis_k,q.BasicGet(0),ko_C_k);
  // pk_C_k = pk_C_ko * ko_C_k;
  FastMult(pk_C_ko,ko_C_k,pk_C_k);
}

void RevoluteJoint::ForwardKinematics(){
  Vect3 result1,result2,result3,result4,result5;
  Vect3 pk_w_k;

  // orientations
  ComputeForwardTransforms();

  // compute position vector r12
  //r12 = point1->position - pk_C_k * point2->position;
  FastMult(pk_C_k,point2->position,result1);
  FastSubt(point1->position,result1,r12);// Jacks comment: needs flipping!!!

  // compute position vector r21
  FastNegMult(k_C_pk,r12,r21);

  // compute global location
  // body2->r = body1->r + body1->n_C_k * r12;
  FastMult(body1->n_C_k,r12,result1);
  FastAdd(body1->r,result1,body2->r);

  // compute qdot (for revolute joint qdot = u)
  // qdot = u
  FastAssign(u,qdot);

  // angular velocities
  //body2->omega = body1->omega + body1->n_C_k * axis_pk * u;
  //pk_w_k = axis_k * u;
  //body2->omega_k = T(pk_C_k) * body1->omega_k + pk_w_k;
  double u_double = u.BasicGet(0);
  FastMult(u_double,axis_pk,result1);
  FastMult(body1->n_C_k,result1,result2);
  FastAdd(body1->omega,result2,body2->omega);
  FastMult(u_double,axis_k,pk_w_k);
  FastTMult(pk_C_k,body1->omega_k,result1);
  FastAdd(result1,pk_w_k,body2->omega_k);

  // compute velocities
  FastCross(body1->omega_k,r12,result1);
  FastCross(point2->position,pk_w_k,result2);
  FastAdd(body1->v_k,result1,result3);
  FastTMult(pk_C_k,result3,result4);
  FastAdd(result2,result4,body2->v_k);
  FastMult(body2->n_C_k,body2->v_k,body2->v);

  // compute state explicit angular acceleration
  FastCross(body2->omega_k,pk_w_k,result1);
  FastTMult(pk_C_k,body1->alpha_t,result2);
  FastAdd(result1,result2,body2->alpha_t);

  // compute state explicit acceleration
  FastCross(body1->alpha_t,point1->position,result1);
  FastCross(body1->omega_k,point1->position,result2);
  FastCross(body1->omega_k,result2,result3);
  FastTripleSum(body1->a_t,result1,result3,result4);
  FastTMult(pk_C_k,result4,result5);

  FastCross(point2->position,body2->alpha_t,result1);
  FastCross(point2->position,body2->omega_k,result2);
  FastCross(body2->omega_k,result2,result3);

  FastTripleSum(result5,result1,result3,body2->a_t);
}

void RevoluteJoint::BackwardKinematics(){
  Vect3 result1,result2,result3,result4,result5;
  Vect3 k_w_pk;

  // orientations
  ComputeBackwardTransforms();

  // compute position vector r21
  //r21 = point2->position - k_C_pk * point1->position;
  FastMult(k_C_pk,point1->position,result1);
  FastSubt(point2->position,result1,r21);

  // compute position vector r21
  FastNegMult(pk_C_k,r21,r12);

  // compute global location
  // body1->r = body2->r + body2->n_C_k * r21;
  FastMult(body2->n_C_k,r21,result1);
  FastAdd(body2->r,result1,body1->r);

  // compute qdot (for revolute joint qdot = u)
  // qdot = u
  FastAssign(u,qdot);

  // angular velocities
  //body1->omega = body2->omega - body2->n_C_k * axis_k * u;
  //k_w_pk = - axis_pk * u;
  //body1->omega_k = pk_C_k * body2->omega_k + k_w_pk;
  double u_double = u.BasicGet(0);
  FastMult(-u_double,axis_k,result1);
  FastMult(body2->n_C_k,result1,result2);
  FastAdd(body2->omega,result2,body1->omega);
  FastMult(-u_double,axis_pk,k_w_pk);
  FastMult(pk_C_k,body2->omega_k,result1);
  FastAdd(result1,k_w_pk,body1->omega_k);

  // compute velocities
  FastCross(body2->omega_k,r21,result1);
  FastCross(point1->position,k_w_pk,result2);
  FastAdd(body2->v_k,result1,result3);
  FastMult(pk_C_k,result3,result4);
  FastAdd(result2,result4,body1->v_k);
  FastMult(body1->n_C_k,body1->v_k,body1->v);

  // compute state explicit angular acceleration
  FastCross(body1->omega_k,k_w_pk,result1);
  FastMult(pk_C_k,body2->alpha_t,result2);
  FastAdd(result1,result2,body1->alpha_t);

  // compute state explicit acceleration
  FastCross(body2->alpha_t,point2->position,result1);
  FastCross(body2->omega_k,point2->position,result2);
  FastCross(body2->omega_k,result2,result3);
  FastTripleSum(body2->a_t,result1,result3,result4);
  FastMult(pk_C_k,result4,result5);

  FastCross(point1->position,body1->alpha_t,result1);
  FastCross(point1->position,body1->omega_k,result2);
  FastCross(body1->omega_k,result2,result3);

  FastTripleSum(result5,result1,result3,body1->a_t);
}
