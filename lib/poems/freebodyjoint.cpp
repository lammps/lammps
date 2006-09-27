/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: freebodyjoint.cpp                                       *
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
 
#include "freebodyjoint.h"
#include "point.h"
#include "matrixfun.h"
#include "body.h"
#include "fastmatrixops.h"
#include "norm.h"
#include "eulerparameters.h"
#include "matrices.h"
#include <iomanip>
   

FreeBodyJoint::FreeBodyJoint(){
  DimQandU(7,6);  
}

FreeBodyJoint::~FreeBodyJoint(){
}

JointType FreeBodyJoint::GetType(){
  return FREEBODYJOINT;
}

bool FreeBodyJoint::ReadInJointData(std::istream& in){
  return true;
}

void FreeBodyJoint::WriteOutJointData(std::ostream& out){
}

void FreeBodyJoint::ComputeLocalTransform(){
  Mat3x3 ko_C_k;  
  EP_Transformation(q, ko_C_k);   
  FastMult(pk_C_ko,ko_C_k,pk_C_k);  
}

Matrix FreeBodyJoint::GetForward_sP(){
  Mat6x6 sP;
  //sP.Identity();      
  
  sP.Zeros();  
  Mat3x3 temp0=T(pk_C_k);
  for(int i=1;i<4;i++){
	  sP(i,i)=1.0;
	  for(int j=1;j<4;j++){
		  sP(3+i,3+j)=temp0(i,j);
		  }
	  } 	  	  
  return sP;    
}

Matrix FreeBodyJoint::GetBackward_sP(){
  Mat6x6 sP;
  sP.Identity();        
  sP =-1.0*sP;
  cout<<"Did I come here in "<<endl;
  return sP;    
}


void FreeBodyJoint::UpdateForward_sP( Matrix& sP){
  // do nothing
}
  
void FreeBodyJoint::UpdateBackward_sP( Matrix& sP){
  // do nothing
}  
  
void FreeBodyJoint::ForwardKinematics(){
 //cout<<"Check in freebody "<<q<<" "<<qdot<<endl;
  EP_Normalize(q);

  // COMMENT STEP1: CALCULATE ORIENTATIONS
  ComputeForwardTransforms();  
  
    
  //COMMENT STEP2: CALCULATE POSITION VECTORS
  Vect3 result1, result2, result3, result4; 
  
  result1.BasicSet(0,q.BasicGet(4));
  result1.BasicSet(1,q.BasicGet(5));
  result1.BasicSet(2,q.BasicGet(6));
  
  FastAssign(result1,r12);    
  FastNegMult(k_C_pk,r12,r21);
  
  FastAssign(r12,body2->r);  
    
  //COMMENT STEP3: CALCULATE QDOT  
  qdot_to_u(q, u, qdot);
  
  
  Vect3 WN;
  WN.BasicSet(0,u.BasicGet(0));
  WN.BasicSet(1,u.BasicGet(1));
  WN.BasicSet(2,u.BasicGet(2));
    
  Vect3 VN; 
  VN.BasicSet(0,u.BasicGet(3));
  VN.BasicSet(1,u.BasicGet(4));
  VN.BasicSet(2,u.BasicGet(5));  
    
  FastAssign(WN,body2->omega_k);  
 
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

void FreeBodyJoint::BackwardKinematics(){
cout<<"Did I come here "<<endl;
         
}
