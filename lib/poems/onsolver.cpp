/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: onsolver.cpp                                            *
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


#include "onsolver.h"
#include "system.h"
#include "onbody.h"
#include "body.h"
#include "matrixfun.h"
#include <fstream>


using namespace std;

OnSolver::OnSolver(){
  numbodies = 0;
  bodyarray = 0;
  q=0;u=0; qdot=0; udot=0; qdotdot=0;  
  type = ONSOLVER;
}

OnSolver::~OnSolver(){
  DeleteModel();
}

void OnSolver::DeleteModel(){
  delete [] bodyarray;

  delete [] q;
  delete [] u;

  delete [] qdot;
  delete [] udot; 
  
  delete [] qdotdot;
  
  numbodies = 0;
}

void OnSolver::CreateModel(){
  // delete old model
  DeleteModel();

  // clear system body IDs (primer for traversal algorithm)
  system->ClearBodyIDs();
  

  // error check for inertial frame
  Body* sysbasebody = system->bodies.GetHeadElement()->value;
  if( sysbasebody->GetType() != INERTIALFRAME ){
    cerr << "ERROR: inertial frame not at head of bodies list" << endl;
    exit(1);
  }

  // setup the O(n) spanning tree
  numbodies = inertialframe.RecursiveSetup( (InertialFrame*) sysbasebody );
  if(!numbodies){
    cerr << "ERROR: unable to create O(n) model" << endl;
    exit(1);
  }
  
  
  
  bodyarray = new OnBody* [numbodies];
  
  CreateTopologyArray(0,&inertialframe);	  
  
  CreateStateMatrixMaps();  
}

int OnSolver::CreateTopologyArray(int num, OnBody* body){
  int i = num;
  bodyarray[i] = body;    
  i++;

  
  OnBody* child;
  ListElement<OnBody>* ele = body->children.GetHeadElement();    
  
  while(ele){
    child = ele->value;    
    i = CreateTopologyArray(i,child);
    ele = ele->next;        
  }  
  return i;  
}

void OnSolver::CreateStateMatrixMaps(){
  
	
  int numstates=0;
  for(int i=1;i<numbodies;i++)    
    numstates += bodyarray[i]->q->GetNumRows();       
  
  state.Dim(numstates);
  statedot.Dim(numstates);
  statedoubledot.Dim(numstates); 

  
  int count=0;
  
  for(int i=1;i<numbodies;i++){
	  for(int j=0;j<bodyarray[i]->q->GetNumRows();j++){
		  state.SetElementPointer(count,bodyarray[i]->q->GetElementPointer(j));      
		  statedot.SetElementPointer(count,bodyarray[i]->qdot->GetElementPointer(j));        	
		  statedoubledot.SetElementPointer(count,bodyarray[i]->qdotdot->GetElementPointer(j));       
		  count++;      
	  }        
  }
}


void OnSolver::Solve(double time, Matrix& FF){
	system->SetTime(time);
	for(int i=1;i<numbodies;i++)
		bodyarray[i]->LocalKinematics();
	
	Vect3 Torque; Torque.Zeros();
	Vect3 Force; Force.Zeros();
	
	for(int i=numbodies-1;i>0;i--){
		Torque(1)=FF(1,i);     
		Torque(2)=FF(2,i);
		Torque(3)=FF(3,i);		
		Force(1)=FF(4,i);     
		Force(2)=FF(5,i);
		Force(3)=FF(6,i);   
		bodyarray[i]->LocalTriangularization(Torque,Force);
	}
   	
	for(int i=1;i<numbodies;i++){
		bodyarray[i]->LocalForwardSubstitution();
	}  
}
