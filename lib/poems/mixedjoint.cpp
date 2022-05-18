/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: mixedjoint.cpp                                          *
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

#include "mixedjoint.h"
#include "point.h"
#include "matrixfun.h"
#include "body.h"
#include "fastmatrixops.h"
#include "norm.h"
#include "eulerparameters.h"
#include "matrices.h"



MixedJoint::MixedJoint() = default;

MixedJoint::~MixedJoint() = default;

JointType MixedJoint::GetType(){
  return MIXEDJOINT;
}

bool MixedJoint::ReadInJointData(std::istream& in){
  return true;
}

void MixedJoint::WriteOutJointData(std::ostream& out){
}

void MixedJoint::SetsP(Matrix& sPr, Vect6& temp_dofs, int i, int j){
  const_sP = sPr;
  dofs = temp_dofs;
  numrots = i;
  numtrans = j;
  if (numrots < 2)
    DimQandU(numrots+numtrans,numrots+numtrans);
  else
    DimQandU((4+numtrans),(numrots+numtrans));
  cout<<"Check "<<4+numtrans<<" "<<numrots+numtrans<<" "<<i<<" "<<j<<endl;
}

void MixedJoint::ComputeLocalTransform(){
  Mat3x3 ko_C_k;
  EP_Transformation(q, ko_C_k);
  FastMult(pk_C_ko,ko_C_k,pk_C_k);
}

Matrix MixedJoint::GetForward_sP(){
  Mat6x6 temp_sP;
  Matrix sP;

  temp_sP.Zeros();
  Mat3x3 temp0=T(pk_C_k);
  for(int i=1;i<4;i++){
    temp_sP(i,i)=1.0;
    for(int j=1;j<4;j++){
      temp_sP(3+i,3+j)=temp0(i,j);
    }
  }
  sP = temp_sP*const_sP;
  return sP;
}

Matrix MixedJoint::GetBackward_sP(){
  Mat6x6 sP;
  sP.Identity();
  sP =-1.0*sP;
  cout<<"Did I come here in "<<endl;
  return sP;
}


void MixedJoint::UpdateForward_sP( Matrix& sP){
  // do nothing
}

void MixedJoint::UpdateBackward_sP( Matrix& sP){
  // do nothing
}

void MixedJoint::ForwardKinematics(){
  if(numrots > 1)
    EP_Normalize(q);

  // COMMENT STEP1: CALCULATE ORIENTATIONS
  ComputeForwardTransforms();


  //COMMENT STEP2: CALCULATE POSITION VECTORS
  Vect3 result1, result2, result3, result4;
  result1.Zeros();
  for (int k=0; k<3; k++){
    if( dofs(3+k) != 0.0 ){
      if (numrots > 1)
        result1.BasicSet(k,q.BasicGet(4 + k));
      else
        result1.BasicSet(k,q.BasicGet(numrots + k));
    }
  }



  FastAssign(result1,r12);     // r12 in parents basis  i.e. Newtonian
  FastNegMult(k_C_pk,r12,r21);      // r21 in body basis

  FastAssign(r12,body2->r);  // This is right

  //COMMENT STEP3: CALCULATE QDOT
  int pp = 0;
  if (numrots > 1){
    ColMatrix temp_u(3+numtrans);
    qdot_to_u(q,temp_u,qdot);
    for (int k=1;k<=6;k++){
      if(dofs(k) != 0.0){
        u.BasicSet(pp,temp_u.BasicGet(k-1));
        pp = pp+1;
      }
    }
  }
  else u = qdot;



  Vect3 WN; WN.Zeros();
  int p = 0;
  for (int k=0;k<3;k++){
    if(dofs(k+1) != 0.0){
      WN.BasicSet(k,u.BasicGet(p));
      p=p+1;
    }
  }// WN is in body basis


  Vect3 VN; VN.Zeros();
  for (int k=0;k<3;k++){
    if( dofs(3+k+1) != 0.0 ) {
      VN.BasicSet(k,u.BasicGet(p));
      p=p+1;
    }
  }// VN is the vector of translational velocities in Newtonian basis

  FastAssign(WN,body2->omega_k);

  // cout<<"Angular Velocity "<<WN<<endl;
  Vect3 pk_w_k;
  FastMult(body2->n_C_k,WN,pk_w_k);
  FastAssign(pk_w_k,body2->omega);



  //COMMENT STEP5: CALCULATE VELOCITES
  FastAssign(VN,body2->v);
  FastTMult(body2->n_C_k,body2->v,body2->v_k);


  //CALCULATE KE

  Matrix tempke;
  tempke = T(body2->v)*(body2->v);
  double ke = 0.0;
  ke = body2->mass*tempke(1,1);
  FastMult(body2->inertia,body2->omega_k,result1);
  tempke= T(body2->omega_k)*result1;
  ke = 0.5*ke + 0.5*tempke(1,1);
  body2->KE=ke;


  //COMMENT STEP6: CALCULATE STATE EXPLICIT ANGULAR ACCELERATIONS
  body2->alpha_t.Zeros();


  //COMMENT STEP7: CALCULATE STATE EXPLICIT ACCELERATIONS
  body2->a_t.Zeros();

  }

void MixedJoint::BackwardKinematics(){
cout<<"Did I come here "<<endl;

}
