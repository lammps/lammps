/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: onbody.cpp                                              *
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

#include "onbody.h"

#include <cstdlib>
#include <iostream>

#include "body.h"
#include "inertialframe.h"
#include "joint.h"
#include "onfunctions.h"
#include "virtualmatrix.h"
#include "matrixfun.h"
#include "eulerparameters.h"
#include "colmatrix.h"

using namespace std;
using namespace POEMS;



OnBody::OnBody(){
  system_body = 0;
  system_joint = 0;
  parent = 0;

  // these terms have zeros which are NEVER overwritten
  sI.Zeros();
  sSC.Zeros();
}

OnBody::~OnBody(){
  children.DeleteValues();
}

int OnBody::RecursiveSetup (InertialFrame* basebody){
  int ID = 0;
  system_body = basebody;

  // record that the traversal algorithm has been here
  if( system_body->GetID() ) return 0;
  ID++;
  system_body->SetID(ID);  

  // setup inertial frame
  SetupInertialFrame();

  Joint* joint;
  OnBody* child;
  int tempID;

  // loop through children calling the recursive setup function
  ListElement<Joint>* ele = system_body->joints.GetHeadElement();
  while(ele){
    joint = ele->value;
    child = new OnBody;

    tempID = child->RecursiveSetup(ID,this,joint);
    if( tempID ){
      children.Append(child);
      ID = tempID;      
    }
    else delete child;

    ele = ele->next;
  }

  return ID;
}

int OnBody::RecursiveSetup(int ID, OnBody* parentbody, Joint* sys_joint){
  // initialize the topology and system analogs
  parent = parentbody;  
  system_joint = sys_joint;
  system_body = system_joint->OtherBody(parentbody->system_body);
  
    
  // record that the traversal algorithm has been here
  if( system_body->GetID() ) return 0;
  ID++;  

  // perform the local setup operations
  Setup();

  Joint* joint;
  OnBody* child;
  int tempID;

  // loop through children calling the recursive setup function
  ListElement<Joint>* ele = system_body->joints.GetHeadElement();
  while(ele){
    joint = ele->value;
    if(joint != sys_joint){
      child = new OnBody;

      tempID = child->RecursiveSetup(ID,this,joint);
      if( tempID ){
        children.Append(child);
        ID = tempID;
      }
      else delete child;
    }

    ele = ele->next;
  }

  return ID;  
}

void OnBody::SetupInertialFrame(){
  // error check
  if(system_body->GetType() != INERTIALFRAME){
    cerr << "ERROR: attempting to setup inertial frame for non-inertial body" << endl;
    exit(1);
  }
  
  // setup gravity
  Vect3 neg_gravity = -((InertialFrame*) system_body)->GetGravity();  
  sAhat.Zeros();
  Set6DLinearVector(sAhat,neg_gravity);
  
}


void OnBody::Setup(){  
  
  // get the direction of the joint
  if( system_joint->GetBody2() == system_body ) joint_dir = FORWARD;
  else joint_dir = BACKWARD;

  // set the function pointers for the joint direction
  if( joint_dir == FORWARD ){    
    kinfun = &Joint::ForwardKinematics;
    updatesP = &Joint::UpdateForward_sP;
    gamma = system_joint->GetR12(); 
    pk_C_k = system_joint->Get_pkCk();
  } else {
    kinfun = &Joint::BackwardKinematics;
    updatesP = &Joint::UpdateBackward_sP;
    gamma = system_joint->GetR21();
    pk_C_k = system_joint->Get_kCpk();
  }

  // initialize variables and dimensions

  // sI
  OnPopulateSI(system_body->inertia, system_body->mass, sI);
	
  // sP
  if( joint_dir == FORWARD )
    sP = system_joint->GetForward_sP();
  else
    sP = system_joint->GetBackward_sP();

  // dimension these quantities
  sM = T(sP)*sP;
  sMinv = sM;
  sPsMinv = sP;
  sIhatsP = sP;
  

  // get the state and state derivative column matrix pointers
  q = system_joint->GetQ();
  u = system_joint->GetU();
  qdot = system_joint->GetQdot();  
  udot = system_joint->GetUdot();
  qdotdot = system_joint->GetQdotdot();
  }

void OnBody::RecursiveKinematics(){
  OnBody* child;
  // Perform local kinematics recursively outward
  ListElement<OnBody>* ele = children.GetHeadElement();
  while(ele){
    child = ele->value;    
    child->LocalKinematics();
    child->RecursiveKinematics();
    Mat3x3 result=*child->pk_C_k;           
    ele = ele->next;
  }
  
}

void OnBody::RecursiveTriangularization(){
  OnBody* child;

  // Perform local triangularization recursively inward
  ListElement<OnBody>* ele = children.GetHeadElement();
  while(ele){
    child = ele->value;
    child->RecursiveTriangularization();
    //child->LocalTriangularization();
    ele = ele->next;
  }

}

void OnBody::RecursiveForwardSubstitution(){
  OnBody* child;
  // Perform local forward substitution recursively outward
  ListElement<OnBody>* ele = children.GetHeadElement();
  while(ele){
    child = ele->value;
   // child->LocalForwardSubstitution();
    child->RecursiveForwardSubstitution();
    ele = ele->next;
  }
}

void OnBody::LocalKinematics(){
  // do the directional kinematics
  (system_joint->*kinfun)();
  (system_joint->*updatesP)( sP );
  OnPopulateSC( *gamma, *pk_C_k, sSC );      
}

Mat3x3 OnBody::GetN_C_K(){
Mat3x3 nck=system_body->n_C_k;
return nck;
}

 
int OnBody::GetBodyID(){
int ID=system_body->GetID();
return ID;
}

Vect3 OnBody::LocalCart(){  
  (system_joint->*kinfun)();  
  Vect3 RR=system_body->r;  
  return RR;
}

  
  
void OnBody::LocalTriangularization(Vect3& Torque, Vect3& Force){	  

	Vect3 Iw,wIw,Ialpha,wIwIalpha,ma,Bforce,Bforce_ma,Btorque,Btorque_wIwIalpha;
  Iw.Zeros();wIw.Zeros();wIwIalpha.Zeros();ma.Zeros();
	Bforce.Zeros();Bforce_ma.Zeros();Btorque.Zeros();Btorque_wIwIalpha.Zeros();
     
	FastMult(system_body->inertia,system_body->omega_k,Iw);  
	FastCross(Iw,system_body->omega_k,wIw);  
  
	FastMult(system_body->inertia, system_body->alpha_t, Ialpha);    
	FastSubt(wIw,Ialpha,wIwIalpha); 
	FastMult(-system_body->mass,system_body->a_t,ma);    
	
	
	Mat3x3 k_C_n=T(system_body->n_C_k);      
	Bforce=k_C_n*Force;
	Btorque=k_C_n*Torque;  
	
	FastAdd(Bforce,ma,Bforce_ma);       
	FastAdd(Btorque,wIwIalpha,Btorque_wIwIalpha);	  
	OnPopulateSVect(Btorque_wIwIalpha,Bforce_ma,sF);
  
  
	sIhat = sI;    
	sFhat = sF;     
	Mat6x6 sPsMinvsPT;
	Mat6x6 sIhatsPsMsPT;
	Mat6x6 sUsIhatsPsMsPT;
	Mat6x6 sTsIhat;
	Mat6x6 sTsIhatsSCT;
	Vect6 sTsFhat;
	Mat6x6 sU;
	sU.Identity();  
  
	OnBody* child;
	ListElement<OnBody>* ele = children.GetHeadElement();
  
	while(ele){
		child = ele->value;
    
		FastMultT(child->sIhatsP,child->sPsMinv,sIhatsPsMsPT); 
		FastSubt(sU,sIhatsPsMsPT,sUsIhatsPsMsPT);              
		FastMult(child->sSC,sUsIhatsPsMsPT,child->sT);         
        
		FastMult(child->sT,child->sIhat,sTsIhat);
		FastMultT(sTsIhat,child->sSC,sTsIhatsSCT); 
		FastAdd(sIhat,sTsIhatsSCT,sIhat);
    
		FastMult(child->sT,child->sFhat,sTsFhat);
		FastAdd(sFhat,sTsFhat,sFhat);            
		ele = ele->next;
	}  
  
	FastMult(sIhat,sP,sIhatsP);   
	FastTMult(sP,sIhatsP,sM);       
	sMinv=SymInverse(sM);    
	FastMult(sP,sMinv,sPsMinv);	
}

void OnBody::LocalForwardSubstitution(){  
	Vect6 sSCTsAhat;
	Vect6 sIhatsSCTsAhat;
	Vect6 sFsIhatsSCTsAhat;
	Vect6 sPudot;
	int type = system_joint->GetType();
	if(type == 1){
		FastTMult(sSC,parent->sAhat,sSCTsAhat);
		FastMult(sIhat,sSCTsAhat,sIhatsSCTsAhat);       
		FastSubt(sFhat,sIhatsSCTsAhat,sFsIhatsSCTsAhat);
		FastTMult(sPsMinv,sFsIhatsSCTsAhat,*udot);		
  
		ColMatrix result1=*udot;      
		ColMatrix result2=*qdot;
		ColMatrix result3=*q;
		int num=result1.GetNumRows();
		ColMatrix result4(num+1);
		result4.Zeros();		
		EPdotdot_udot(result1, result2, result3, result4);
		FastAssign(result4, *qdotdot);    
		FastMult(sP,*udot,sPudot);      
		FastAdd(sSCTsAhat,sPudot,sAhat);		
	}
	else if (type == 4){		
		FastTMult(sSC,parent->sAhat,sSCTsAhat);
		FastMult(sIhat,sSCTsAhat,sIhatsSCTsAhat);       
		FastSubt(sFhat,sIhatsSCTsAhat,sFsIhatsSCTsAhat);
		FastTMult(sPsMinv,sFsIhatsSCTsAhat,*udot);
		
		ColMatrix result1=*udot;      		
		ColMatrix result2=*qdot;
		ColMatrix result3=*q;
		int num=result1.GetNumRows();
		ColMatrix result4(num+1);
		result4.Zeros();
		
		EPdotdot_udot(result1, result2, result3, result4);
		FastAssign(result4, *qdotdot);    
		FastMult(sP,*udot,sPudot);      
		FastAdd(sSCTsAhat,sPudot,sAhat);
	}	
	else if (type == 5){		
		FastTMult(sSC,parent->sAhat,sSCTsAhat);
		FastMult(sIhat,sSCTsAhat,sIhatsSCTsAhat);       
		FastSubt(sFhat,sIhatsSCTsAhat,sFsIhatsSCTsAhat);
		FastTMult(sPsMinv,sFsIhatsSCTsAhat,*udot);				
		ColMatrix temp_u = *udot;
		
		ColMatrix result1(3);		
		result1(1)=0.0;
		result1(2)=temp_u(1);      		
		result1(3)=temp_u(2); 		
		ColMatrix result2=*qdot;
		ColMatrix result3=*q;
		int num=result1.GetNumRows();
		ColMatrix result4(num+1);
		result4.Zeros();
		
		EPdotdot_udot(result1, result2, result3, result4);		
		FastAssign(result4, *qdotdot);    
		FastMult(sP,*udot,sPudot);      
		FastAdd(sSCTsAhat,sPudot,sAhat);
	}
	else if (type == 6){		
		FastTMult(sSC,parent->sAhat,sSCTsAhat);
		FastMult(sIhat,sSCTsAhat,sIhatsSCTsAhat);       
		FastSubt(sFhat,sIhatsSCTsAhat,sFsIhatsSCTsAhat);
		FastTMult(sPsMinv,sFsIhatsSCTsAhat,*udot);				
		ColMatrix temp_u = *udot;
		int tt = temp_u.GetNumRows() + 1;
		ColMatrix result1(tt);				
		result1(1)=0.0;
		for (int k =2; k<=tt; k++){
			result1(k)=temp_u(k-1);      				
		}			
		ColMatrix result2=*qdot;
		ColMatrix result3=*q;
		int num=result1.GetNumRows();
		ColMatrix result4(num+1);
		result4.Zeros();
		
		EPdotdot_udot(result1, result2, result3, result4);		
		FastAssign(result4, *qdotdot);    
		FastMult(sP,*udot,sPudot);      
		FastAdd(sSCTsAhat,sPudot,sAhat);
	}		
	else{
		cout<<"Joint type not recognized in onbody.cpp LocalForwardSubsitution() "<<type<<endl;
		exit(-1);
	}
	CalculateAcceleration();
	
}


void OnBody::CalculateAcceleration(){
}
