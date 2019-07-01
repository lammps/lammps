/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: system.cpp                                              *
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


#include "system.h"

#include <cstddef>
#include <iostream>

#include "body.h"
#include "joint.h"
#include "colmatrix.h"
#include "eulerparameters.h"
#include "fixedpoint.h"
#include "freebodyjoint.h"
#include "inertialframe.h"
#include "mat3x3.h"
#include "matrix.h"
#include "matrixfun.h"
#include "rigidbody.h"
#include "sphericaljoint.h"
#include "vect3.h"
#include "virtualmatrix.h"

class Point;


System::System(){
	mappings = NULL;
}

System::~System(){
  Delete();
}

void System::Delete(){
  delete [] mappings;
  bodies.DeleteValues();
  joints.DeleteValues();
}

int System::GetNumBodies(){
  return bodies.GetNumElements();
}

int * System::GetMappings()
{
	return mappings;
}

void System::AddBody(Body* body){
  bodies.Append(body);
}

void System::AddJoint(Joint* joint){
  joints.Append(joint);
}

void System::SetTime(double t){
  time = t;
}

double System::GetTime(){
  return time;
}

void System::ComputeForces(){	
	// NOT DOING ANYTHING AT THIS TIME
  }  
  
bool System::ReadIn(istream& in){
  int numbodies;
  Body* body;
  int bodytype;
  char bodyname[256];
  int index;  

  // get number of bodies
  in >> numbodies;

  // bodies loop
  for(int i=0;i<numbodies;i++){

    // error check
    in >> index;
    if(index != i){
      cerr << "Error reading bodies" << endl;
      return false;
    }

    in >> bodytype >> bodyname;
    body = NewBody(bodytype);

    // type check
    if(!body){
      cerr << "Unrecognized body type '" << bodytype << "'" << endl;
      return false;
    }

    // add the body
    AddBody(body);

    // set generic body info
    body->ChangeName(bodyname);

    // read in the rest of its data
    if(!body->ReadIn(in)) return false;
  }

  // create a temporary array for fast indexed access
  Body** bodyarray = bodies.CreateArray();

  int numjoints;
  int jointtype;
  char jointname[256];
  Joint* joint;
  int body1, body2;
  int point1, point2;

  // get number of joints
  in >> numjoints;

  // joints loop
  for(int i=0;i<numjoints;i++){

    // error check
    in >> index;
    if(index != i){
      cerr << "Error reading joints" << endl;
      return false;
    }

    in >> jointtype >> jointname;
    joint = NewJoint(jointtype);

    // joint type check
    if(!joint){
      cerr << "Unrecognized joint type '" << jointtype << "'" << endl;
      return false;
    }

    // add the joint
    AddJoint(joint);

    // set the generic joint info
    joint->ChangeName(jointname);

    in >> body1 >> body2;
    if( !(body1<numbodies) || !(body2<numbodies) ){
      cerr << "Body index out of range" << endl;
      delete [] bodyarray;
      return false;
    }
    
    joint->SetBodies(bodyarray[body1], bodyarray[body2]);

    bodyarray[body1]->AddJoint(joint);
    bodyarray[body2]->AddJoint(joint);
    
    in >> point1 >> point2;

    joint->SetPoints(bodyarray[body1]->GetPoint(point1),bodyarray[body2]->GetPoint(point2));

    // read in the rest of its data
    if(!joint->ReadIn(in)){
      delete [] bodyarray;
      return false;
    }
  }

  // delete the temporary array
  delete [] bodyarray;

  return true;
}

void System::WriteOut(ostream& out){
  // number of bodies
  out << bodies.GetNumElements() << endl;

  // bodies loop
  int i = 0;
  Body* body;
  ListElement<Body>* b_ele = bodies.GetHeadElement();
  while(b_ele !=0){
    out << i << ' ';

    body = b_ele->value;

    // set the body ID for later identification
    body->SetID(i);
    
    // write out the data
    body->WriteOut(out);

    i++; b_ele = b_ele->next;
  }

  // number of joints
  out << joints.GetNumElements() << endl;  

  // joints loop
  i = 0;
  Joint* joint;
  ListElement<Joint>* j_ele = joints.GetHeadElement();
  while(j_ele !=0){
    out << i << ' ';
    joint = j_ele->value;

    // set the joint ID for later identification
    joint->SetID(i);

    // write out the data
    joint->WriteOut(out);
    
    i++; j_ele = j_ele->next;
  }
}

void System::ClearBodyIDs(){
  ListElement<Body>* current = bodies.GetHeadElement();

  while(current){
    current->value->SetID(0);
    current = current->next;
  }
}

void System::ClearJointIDs(){
  ListElement<Joint>* current = joints.GetHeadElement();

  while(current){
    current->value->SetID(0);
    current = current->next;
  }
}


void System::Create_DegenerateSystem(int& nfree, int*freelist, double *&masstotal, double **&inertia, double **&xcm, double **&vcm, double **&omega, double **&ex_space, double **&ey_space, double **&ez_space){

 //-------------------------------------------------------------------------//    
  // Declaring Temporary Entities 
  //-------------------------------------------------------------------------//   
     Body* body = NULL;
     Body* prev;
     Body* Inertial;
     Point* origin;
     Joint* joint;
     Point* point_CM; 
     Point* point_p;  
     Point* point_k;  
     Point* point_ch = NULL;
     Vect3 r1,r2,r3,v1,v2,v3; 
     Mat3x3 IM, N, PKCK,PKCN ;
     ColMatrix qo, uo, q, qdot,w;
     
	 mappings = new int[nfree];
	 for(int i = 0; i < nfree; i++)
	 {
		 mappings[i] = freelist[i];		 
	 }
     qo.Dim(4);
     uo.Dim(3);
     q.Dim(4);
     qdot.Dim(4);     
     PKCN.Identity();
     PKCK.Identity();
     w.Dim(3);
	
//-------------------------------------------------------------------------//    
  // Setting up Inertial Frame, gravity and Origin  
  //-------------------------------------------------------------------------//
     Inertial= new InertialFrame;  
     AddBody(Inertial);
     
     Vect3 temp1;  
     temp1.Zeros();  
     ((InertialFrame*) Inertial)->SetGravity(temp1);    
     origin= new FixedPoint(temp1);
     Inertial->AddPoint(origin);    
//-------------------------------------------------------------------------//
	double ** xh1 = new double*[nfree];
	double ** xh2 = new double*[nfree];
	
	for (int i=0; i<nfree; i++){		
 		xh1[i] = new double[3];
		xh2[i] = new double[3];
	}	
	for (int i=0; i<nfree; i++){		
		for (int j=0; j<3; j++){
			xh1[i][j] = xcm[mappings[i]-1][j];
		}
	}
	
//-------------------------------------------------------------------------//    
// Begin looping over each body for recursive kinematics
//-------------------------------------------------------------------------//
    for(int i=0;i<nfree;i++){
	    prev=Inertial;
	    point_p=origin;               
      
     body = new RigidBody;   
     body->mass=masstotal[mappings[i]-1];
     IM(1,1)=inertia[mappings[i]-1][0];
     IM(2,2)=inertia[mappings[i]-1][1];
     IM(3,3)=inertia[mappings[i]-1][2];
     IM(1,2)=0.0;
     IM(1,3)=0.0;
     IM(2,3)=0.0;
     IM(2,1)=IM(1,2);
     IM(3,1)=IM(1,3);
     IM(3,2)=IM(2,3);
     body->inertia = IM;   

//-------------------------------------------------------//                  
    
    
	    for (int k=0;k<3;k++){ 
          r1(k+1)=xh1[i][k]-xcm[mappings[i]-1][k];     
			 r3(k+1) = xcm[mappings[i]-1][k];
			 r3(k+1)=xh2[i][k]-xcm[mappings[i]-1][k]; 	  
	    }
	  
     r2.Zeros(); 
     
     for (int k=1;k<=3;k++){     
          N(k,1)=ex_space[mappings[i]-1][k-1];              
          N(k,2)=ey_space[mappings[i]-1][k-1];              
          N(k,3)=ez_space[mappings[i]-1][k-1];              
         }                  
     
     PKCK=T(N);
     PKCN=T(N);
    
     q.Zeros();     
     EP_FromTransformation(q,N);
                    
     r1=PKCN*r1;
     r3=PKCN*r3;
     
     for (int k=1;k<=3;k++){
     w(k)=omega[mappings[i]-1][k-1];                                   
     }     
     
     Vect3 cart_r, cart_v;
     for (int k=1;k<=3;k++){
	     cart_r(k)=xcm[mappings[i]-1][k-1];	     
	     cart_v(k)=vcm[mappings[i]-1][k-1];
	     }
	     	     
	     w=PKCN*w;     
	     EP_Derivatives(q,w,qdot);
	     
     
//-------------------------------------------------------------------------//    
// Create bodies and joints with associated properties for POEMS
//-------------------------------------------------------------------------//

     point_CM = new FixedPoint(r2);
     point_k = new FixedPoint(r1);
     point_ch = new FixedPoint(r3);
     body->AddPoint(point_CM);
     body->AddPoint(point_k);
     body->AddPoint(point_ch);
     AddBody(body);    
   
     Mat3x3 One;
     One.Identity();	
	  ColMatrix qq=Stack(q,cart_r);
          ColMatrix vv=Stack(qdot,cart_v);          
          joint=new FreeBodyJoint;
          AddJoint(joint);
          joint->SetBodies(prev,body);
          body->AddJoint(joint);
          prev->AddJoint(joint);
          joint->SetPoints(point_p,point_k);
          joint->SetZeroOrientation(One);
          joint->DimQandU(7,6);
          joint->SetInitialState(qq,vv);
          joint->ForwardKinematics();                             
  }
  for(int i = 0; i < nfree; i++) {
	  delete [] xh1[i];
	  delete [] xh2[i];
  }
  delete [] xh1;
  delete [] xh2;  
}


void System::Create_System_LAMMPS(int numbodies, double *mass,double **inertia, double ** xcm, double ** xjoint,double **vcm,double **omega,double **ex_space, double **ey_space, double **ez_space, int b, int * mapping, int count){

	//-------------------------------------------------------------------------//    
  // Declaring Temporary Entities 
	//-------------------------------------------------------------------------//   
  	
	Body* body = NULL;
	Body* prev;
	Body* Inertial;
	Point* origin;
	Joint* joint;
	Point* point_CM;   
	Point* point_p;    
	Point* point_k;    
	Point* point_ch = NULL;  
	Vect3 r1,r2,r3,v1,v2,v3; 
	Mat3x3 IM, N, PKCK,PKCN ;
	ColMatrix qo, uo, q, qdot,w;
	Vect3 cart_r, cart_v;
	mappings = new int[b];
	for(int i = 0; i < b; i++){
		mappings[i] = mapping[i];		 
	}	 
	 
     
	qo.Dim(4);
	uo.Dim(3);
	q.Dim(4);
	qdot.Dim(4);     
	PKCN.Identity();
	PKCK.Identity();
	w.Dim(3);
	
	//-------------------------------------------------------------------------//    
  // Setting up Inertial Frame, gravity and Origin  
	//-------------------------------------------------------------------------//
	Inertial= new InertialFrame;  
	AddBody(Inertial);
     
	Vect3 temp1;  
	temp1.Zeros();  
	((InertialFrame*) Inertial)->SetGravity(temp1);    
	origin= new FixedPoint(temp1);
	Inertial->AddPoint(origin);    
	//-------------------------------------------------------------------------//

	double ** xh1;
	double ** xh2;
		  
	xh1 = new double*[b];		  
	xh2 = new double*[b];
		  
	
	for (int i=0; i<b; i++){			
		xh1[i] = new double[3];
		xh2[i] = new double[3];			
	}	
		  
		  
		  
	for (int j=0; j<3; j++){
		xh1[0][j] = xcm[mapping[0]-1][j];
		xh2[b-1][j] = xcm[mapping[b-1]-1][j];
	}
	
	for (int i=0; i<b-1; i++){
		for (int j=0; j<3; j++){
			xh1[i+1][j] = xjoint[mapping[i]-count-1][j];
		}
	}		  
	
	for (int i=0; i<b-1; i++){
		for (int j=0; j<3; j++){
			xh2[i][j] = xjoint[mapping[i]-count-1][j];
		}
	}

	
	//-------------------------------------------------------------------------//    
// Begin looping over each body for recursive kinematics
	//-------------------------------------------------------------------------//
	for(int i=0;i<b;i++){	    
		if (i == 0){
			prev=Inertial;
			point_p=origin;
		}    
		else{
			prev = body;            
			point_p = point_ch;     
		}   
    
  	  
		body = new RigidBody;  
		body->mass=mass[mapping[i]-1];
		IM(1,1)=inertia[mapping[i]-1][0];
		IM(2,2)=inertia[mapping[i]-1][1];
		IM(3,3)=inertia[mapping[i]-1][2];
		IM(1,2)=0.0;
		IM(1,3)=0.0;
		IM(2,3)=0.0;
		IM(2,1)=IM(1,2);
		IM(3,1)=IM(1,3);
		IM(3,2)=IM(2,3);
		body->inertia = IM;   
     
		//-------------------------------------------------------//
             
		for (int k=0;k<3;k++){
			r1(k+1)=xh1[i][k]-xcm[mapping[i]-1][k]; 
			r3(k+1)=xh2[i][k]-xcm[mapping[i]-1][k]; 
		}
		r2.Zeros();   
     
		for (int k=1;k<=3;k++){     
			N(k,1)=ex_space[mapping[i]-1][k-1];              
			N(k,2)=ey_space[mapping[i]-1][k-1];              
			N(k,3)=ez_space[mapping[i]-1][k-1];              
		}
	 
         
		if (i==0){
			PKCK=T(N);     
			PKCN=T(N);     
	     
			q.Zeros();	     
			EP_FromTransformation(q,N);
	     
			r1=PKCN*r1;
			r3=PKCN*r3;
	     
			for (int k=1;k<=3;k++){     
				w(k)=omega[mappings[i]-1][k-1];
			}	     
	     
			for (int k=1;k<=3;k++){
				cart_r(k)=xcm[mappings[i]-1][k-1];		     
				cart_v(k)=vcm[mappings[i]-1][k-1];
			}
			w=PKCN*w;     
			EP_Derivatives(q,w,qdot);    
		 
		}                 
		else{	     
			PKCK=PKCN*N;   
			PKCN=T(N);     
    
			q.Zeros();
			EP_FromTransformation(q,PKCK);
               
			r1=PKCN*r1;
			r3=PKCN*r3;
     
			for (int k=1;k<=3;k++){     
				w(k)=omega[mapping[i]-1][k-1]-omega[mapping[i-1]-1][k-1];                         
			}
    	
			w=PKCN*w; 
			EP_Derivatives(q, w, qdot);        
		}  
        
	
		//-------------------------------------------------------------------------//    
// Create bodies and joints with associated properties for POEMS
		//-------------------------------------------------------------------------//

		point_CM = new FixedPoint(r2);
		point_k = new FixedPoint(r1);
		point_ch = new FixedPoint(r3);
		body->AddPoint(point_CM);
		body->AddPoint(point_k);
		body->AddPoint(point_ch);
		AddBody(body);    
   
		Mat3x3 One;
		One.Identity();	
		if (i==0){	  
			ColMatrix qq=Stack(q,cart_r);
			ColMatrix vv=Stack(qdot,cart_v);          
			joint=new FreeBodyJoint;
			AddJoint(joint);
			joint->SetBodies(prev,body);
			body->AddJoint(joint);
			prev->AddJoint(joint);
			joint->SetPoints(point_p,point_k);
			joint->SetZeroOrientation(One);
			joint->DimQandU(7,6);
			joint->SetInitialState(qq,vv);  	  
			joint->ForwardKinematics();  
		}    
		else{	     
			joint= new SphericalJoint;
			AddJoint(joint);
			joint->SetBodies(prev,body);
			body->AddJoint(joint);
			prev->AddJoint(joint);
			joint->SetPoints(point_p,point_k);          
			joint->SetZeroOrientation(One);
			joint->DimQandU(4,3);
			joint->SetInitialState(q,qdot);
			joint->ForwardKinematics();  
		}
	}
	for(int i = 0; i < b; i++)
	{
		delete [] xh1[i];
		delete [] xh2[i];
	}
	delete [] xh1;
	delete [] xh2;
  
}
