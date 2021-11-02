/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: sphericaljoint.cpp                                      *
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


#include "sphericaljoint.h"
#include "point.h"
#include "matrixfun.h"
#include "body.h"
#include "fastmatrixops.h"
#include "norm.h"
#include "eulerparameters.h"
#include "matrices.h"
#include <iomanip>


SphericalJoint::SphericalJoint(){
  DimQandU(4,3);
}
SphericalJoint::~SphericalJoint() = default;

JointType SphericalJoint::GetType(){
  return SPHERICALJOINT;
}

bool SphericalJoint::ReadInJointData(std::istream& in){
  return true;
}

void SphericalJoint::WriteOutJointData(std::ostream& out){
}

Matrix SphericalJoint::GetForward_sP(){
  Mat3x3 sPa,sPl;
  Matrix sP(6,3);
  sPa.Identity();
  sPl.Zeros();
  Vect3 temp = -(point2->position);

  sPl(1,2) = temp(3);
  sPl(1,3) = -temp(2);

  sPl(2,1) = -temp(3);
  sPl(2,3) = temp(1);

  sPl(3,1) = temp(2);
  sPl(3,2) = -temp(1);

  sP=Stack(sPa,sPl);
  return sP;
}

void SphericalJoint::UpdateForward_sP( Matrix& sP){
    // sP is constant, do nothing.
  // linear part is not constant
}

Matrix SphericalJoint::GetBackward_sP(){
  cout<<" -----------"<<endl;
  cout<<"Am I coming here "<<endl;
  cout<<" -----------"<<endl;
  Mat3x3 sPa,sPl;
  Matrix sP;
  sPa.Identity();
  sPl.Zeros();
  sPl(3,2)=(point2->position(1));
  sPl(2,3)=-(point2->position(1));
  sP=Stack(sPa,sPl);
  return sP;
}

void SphericalJoint::UpdateBackward_sP( Matrix& sP){
  // sP is constant, do nothing.
}

void SphericalJoint::ComputeLocalTransform(){
  Mat3x3 ko_C_k;
  // Obtain the transformation matrix from euler parameters
  EP_Transformation(q, ko_C_k);
  FastMult(pk_C_ko,ko_C_k,pk_C_k);
  }


void SphericalJoint::ForwardKinematics(){
  Vect3 result1,result2,result3,result4,result5;
  Vect3 pk_w_k;

  //cout<<"Check in spherical "<<q<<" "<<qdot<<endl;
  EP_Normalize(q);


  // orientations
  ComputeForwardTransforms();


  //----------------------------------//
  // COMPUTE POSITION VECTOR R12 aka GAMMA

  FastNegMult(pk_C_k,point2->position,result1); // parents basis
  FastAdd(result1,point1->position,r12);

  // compute position vector r21
  FastNegMult(k_C_pk,r12,r21);



  //----------------------------------//
  // COMPUTE GLOBAL LOCATION
  FastMult(body1->n_C_k,(body1->GetPoint(2))->position,result1);
  FastAdd(result1,body1->r,result1);
  FastNegMult(body2->n_C_k,(body2->GetPoint(1))->position,result2);
  FastAdd(result1,result2,body2->r);

  qdot_to_u(q, u, qdot);


  //-----------------------------------
  // angular velocities

  FastAssign(u,pk_w_k);
  FastTMult(pk_C_k,body1->omega_k,result1);
  FastAdd(result1,pk_w_k,body2->omega_k);
  FastMult(body2->n_C_k,body2->omega_k,body2->omega);



  //-----------------------------------

  // compute velocities
  FastCross(body1->omega_k,(body1->GetPoint(2))->position,result1);
  FastAdd(body1->v_k,result1,result2);
  FastTMult(pk_C_k,result2,result1); // In body basis

  FastCross((body2->GetPoint(1))->position,body2->omega_k,result2);
  FastAdd(result1,result2,body2->v_k);    // In body basis
  FastMult(body2->n_C_k,body2->v_k,body2->v);


  //------------------------------------------
  //Compute the KE
  Matrix tempke;
  tempke = T(body2->v)*(body2->v);
  double ke = 0.0;
  ke = body2->mass*tempke(1,1);
  FastMult(body2->inertia,body2->omega_k,result1);
  tempke= T(body2->omega_k)*result1;
  ke = 0.5*ke + 0.5*tempke(1,1);
  body2->KE=ke;

  //-----------------------------------
  // compute state explicit angular acceleration  // Has to be in body basis
  FastTMult(pk_C_k,body1->alpha_t,result2);
  FastCross(body2->omega_k,pk_w_k,result1);
  FastAdd(result1,result2,body2->alpha_t);

  //-----------------------------------
  // compute state explicit acceleration
  // NEED TO DO THIS COMPLETELY IN BODY BASIS

  FastCross(body1->omega_k,(body1->GetPoint(2))->position,result1);
  FastCross(body1->omega_k,result1,result2);
  FastTMult(pk_C_k,result2,result1);

  //FastCross(body2->omega_k,-(body2->GetPoint(1))->position,result3);
  FastCross((body2->GetPoint(1))->position,body2->omega_k,result3);
  FastCross(body2->omega_k,result3,result2);
  FastAdd(result1,result2,result3); //wxwxr in body basis

  FastCross(body1->alpha_t,(body1->GetPoint(2))->position,result4);
  FastTMult(pk_C_k,result4,result5);
  FastAssign(result5,result4);

  FastCross((body2->GetPoint(1))->position,body2->alpha_t,result2);
  FastAdd(result2,result4,result1); //alphaxr in body basis

  FastTMult(pk_C_k,body1->a_t,result2);
  FastTripleSum(result3,result1,result2,body2->a_t);     // in body basis


  //-----------------------------------
}

// NOTE: NOT USING BACKWARDKINEMATICS AT PRESENT
void SphericalJoint::BackwardKinematics(){
  cout<<"what about here "<<endl;
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
  ColMatrix us(3);
  /*us(1)=0;
  us(2)=u(1);
  us(3)=u(2);*/
  EP_Derivatives(q,u,qdot);

  // angular velocities

  FastMult(body2->n_C_k,u,result2);
  FastAdd(body2->omega,result2,body1->omega);
  FastAssign(u,k_w_pk);
  FastMult(pk_C_k,body2->omega_k,result1);
  FastSubt(result1,k_w_pk,body1->omega_k);
  cout<<"The program was here"<<endl;

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
