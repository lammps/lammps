/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: prismaticjoint.cpp                                      *
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

#include "prismaticjoint.h"
#include "point.h"
#include "matrixfun.h"
#include "body.h"
#include "fastmatrixops.h"

PrismaticJoint::PrismaticJoint(){
  q.Dim(1);
  qdot.Dim(1);
  u.Dim(1);
  udot.Dim(1);
}
PrismaticJoint::~PrismaticJoint(){
}

JointType PrismaticJoint::GetType(){
  return PRISMATICJOINT;
}

bool PrismaticJoint::ReadInJointData(std::istream& in){
  in >> axis_pk;
  axis_k = T(pk_C_ko)*axis_pk;

  // init the constant transforms
  pk_C_k = pk_C_ko;
  k_C_pk = T(pk_C_k);

  return true;
}

void PrismaticJoint::WriteOutJointData(std::ostream& out){
  out << axis_pk;
}

Matrix PrismaticJoint::GetForward_sP(){
  Vect3 zero;
  zero.Zeros();

  // sP = [zero;axis]
  return Stack(zero,axis_k);
}

void PrismaticJoint::UpdateForward_sP( Matrix& sP){
  // sP is constant, do nothing.
}

Matrix PrismaticJoint::GetBackward_sP(){
  Vect3 zero;
  zero.Zeros();

  // sP = [zero;axis]
  return -Stack(zero,axis_pk);
}

void PrismaticJoint::UpdateBackward_sP( Matrix& sP){
  // sP is constant, do nothing.
}

void PrismaticJoint::ComputeForwardTransforms(){
  ComputeForwardGlobalTransform();
}

void PrismaticJoint::ComputeBackwardTransforms(){
  ComputeBackwardGlobalTransform();
}

void PrismaticJoint::ComputeLocalTransform(){
  // the transform is constant, do nothing
}

void PrismaticJoint::ForwardKinematics(){
  Vect3 result1,result2,result3;
  Vect3 d_pk;

  // orientations
  ComputeForwardTransforms();

  // compute position vector r12
  //r12 = point1->position + axis_pk * q - pk_C_k * point2->position;
  FastMult(pk_C_k,point2->position,result1);
  FastMult(q.BasicGet(0),axis_pk,d_pk);
  FastTripleSumPPM(point1->position,d_pk,result1,r12);

  // compute position vector r21
  FastNegMult(k_C_pk,r12,r21);

  // compute global location
  // body2->r = body1->r + body1->n_C_k * r12;
  FastMult(body1->n_C_k,r12,result1);
  FastAdd(body1->r,result1,body2->r);

  // compute qdot (for Prismatic joint qdot = u)
  // qdot = u
  FastAssign(u,qdot);

  // angular velocities
  //body2->omega = body1->omega;
  //body2->omega_k = T(pk_C_k) * body1->omega_k;
  FastAssign(body1->omega,body2->omega);
  FastMult(k_C_pk,body1->omega_k,body2->omega_k);

  // compute velocities
  Vect3 pk_v_k;
  Vect3 wxgamma;
  FastMult(u.BasicGet(0),axis_k,pk_v_k);
  FastMult(k_C_pk,body1->v_k,result1);
  FastCross(body2->omega_k,r12,wxgamma);
  FastTripleSum(result1,pk_v_k,wxgamma,body2->v_k);
  FastMult(body2->n_C_k,body2->v_k,body2->v);

  // compute state explicit angular acceleration
  FastMult(k_C_pk,body1->alpha_t,body2->alpha_t);

  // compute state explicit acceleration
  FastCross(r21,body1->alpha_t,result1);
  FastAdd(body1->a_t,result1,result2);
  FastMult(k_C_pk,result2,result1);

  FastCross(body2->omega_k,pk_v_k,result2);
  FastMult(2.0,result2,result3);
  FastCross(body2->omega_k,wxgamma,result2);
  FastTripleSum(result1,result2,result3,body2->a_t);
}

void PrismaticJoint::BackwardKinematics(){
  Vect3 result1,result2,result3;
  Vect3 d_k;

  // orientations
  ComputeBackwardTransforms();

  // compute position vector r21
  //r21 = point2->position + axis_k * q - k_C_pk * point1->position;
  FastMult(k_C_pk,point1->position,result1);
  FastMult(-q.BasicGet(0),axis_k,d_k);
  FastTripleSumPPM(point2->position,d_k,result1,r21);

  // compute position vector r12
  FastNegMult(pk_C_k,r21,r12);

  // compute global location
  // body1->r = body2->r + body2->n_C_k * r21;
  FastMult(body2->n_C_k,r21,result1);
  FastAdd(body2->r,result1,body1->r);

  // compute qdot (for Prismatic joint qdot = u)
  // qdot = u
  FastAssign(u,qdot);

  // angular velocities
  //body1->omega = body2->omega;
  //body1->omega_k = pk_C_k * body2->omega_k;
  FastAssign(body2->omega,body1->omega);
  FastMult(pk_C_k,body2->omega_k,body1->omega_k);

  // compute velocities
  Vect3 k_v_pk;
  Vect3 wxgamma;
  FastMult(-u.BasicGet(0),axis_pk,k_v_pk);
  FastMult(pk_C_k,body2->v_k,result1);
  FastCross(body1->omega_k,r21,wxgamma);
  FastTripleSum(result1,k_v_pk,wxgamma,body1->v_k);
  FastMult(body1->n_C_k,body1->v_k,body1->v);

  // compute state explicit angular acceleration
  FastMult(pk_C_k,body2->alpha_t,body1->alpha_t);

  // compute state explicit acceleration
  FastCross(r12,body2->alpha_t,result1);
  FastAdd(body2->a_t,result1,result2);
  FastMult(pk_C_k,result2,result1);

  FastCross(body1->omega_k,k_v_pk,result2);
  FastMult(2.0,result2,result3);
  FastCross(body1->omega_k,wxgamma,result2);
  FastTripleSum(result1,result2,result3,body1->a_t);
}
