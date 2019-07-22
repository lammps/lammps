/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: eulerparameters.cpp                                     *
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

#include "eulerparameters.h"
#include <cmath>

using namespace std;

void EP_Derivatives(ColMatrix& q, ColMatrix& u, ColMatrix& qdot){  
  EP_Normalize(q);
  int num=u.GetNumRows();
  if (3<num){
	  for (int i=4; i<=num; i++){
		  qdot.elements[i]=u.elements[i-1];
	  }
  }
  
  qdot.elements[0] = 0.5 *(q.elements[3]*u.elements[0] - q.elements[2]*u.elements[1] + q.elements[1]*u.elements[2]);  
  qdot.elements[1] = 0.5 *(q.elements[2]*u.elements[0] + q.elements[3]*u.elements[1] - q.elements[0]*u.elements[2]);
  qdot.elements[2] = 0.5 *(-q.elements[1]*u.elements[0] + q.elements[0]*u.elements[1] + q.elements[3]*u.elements[2]);
  qdot.elements[3] = -0.5 *(q.elements[0]*u.elements[0] + q.elements[1]*u.elements[1] + q.elements[2]*u.elements[2]);
    
  }

void EP_Transformation(ColMatrix& q, Mat3x3& C){
  EP_Normalize(q);
  
  double q11 = q.elements[0]*q.elements[0];
  double q22 = q.elements[1]*q.elements[1];
  double q33 = q.elements[2]*q.elements[2];
  double q44 = q.elements[3]*q.elements[3];

  double q12 = q.elements[0]*q.elements[1];
  double q13 = q.elements[0]*q.elements[2];
  double q14 = q.elements[0]*q.elements[3];
  double q23 = q.elements[1]*q.elements[2];
  double q24 = q.elements[1]*q.elements[3];
  double q34 = q.elements[2]*q.elements[3];

  C.elements[0][0] = q11 - q22 - q33 + q44;
  C.elements[1][1] = -q11 + q22 - q33 + q44;
  C.elements[2][2] = -q11 - q22 + q33 + q44;

  C.elements[0][1] = 2*(q12 - q34);
  C.elements[1][0] = 2*(q12 + q34);

  C.elements[0][2] = 2*(q13 + q24);
  C.elements[2][0] = 2*(q13 - q24);

  C.elements[1][2] = 2*(q23 - q14);
  C.elements[2][1] = 2*(q23 + q14);
}

void EP_FromTransformation(ColMatrix& q, Mat3x3& C){
  double b[4];

  // condition indicators
  b[0] = C.elements[0][0] - C.elements[1][1] - C.elements[2][2] + 1;
  b[1] = -C.elements[0][0] + C.elements[1][1] - C.elements[2][2] + 1;
  b[2] = -C.elements[0][0] - C.elements[1][1] + C.elements[2][2] + 1;
  b[3] = C.elements[0][0] + C.elements[1][1] + C.elements[2][2] + 1;

  int max = 0;
  for(int i=1;i<4;i++){
    if( b[i] > b[max] ) max = i;
  }

  if( max == 3 ){
    q.elements[3] = 0.5 * sqrt( b[3] );
    q.elements[0] = ( C.elements[2][1] - C.elements[1][2] ) / ( 4 * q.elements[3] );
    q.elements[1] = ( C.elements[0][2] - C.elements[2][0] ) / ( 4 * q.elements[3] );
    q.elements[2] = ( C.elements[1][0] - C.elements[0][1] ) / ( 4 * q.elements[3] );
    return;
  }

  if( max == 0 ){
    q.elements[0] = 0.5 * sqrt( b[0] );
    q.elements[1] = ( C.elements[0][1] + C.elements[1][0] ) / ( 4 * q.elements[0] );
    q.elements[2] = ( C.elements[0][2] + C.elements[2][0] ) / ( 4 * q.elements[0] );
    q.elements[3] = ( C.elements[2][1] - C.elements[1][2] ) / ( 4 * q.elements[0] );
    return;
  }

  if( max == 1 ){
    q.elements[1] = 0.5 * sqrt( b[1] );
    q.elements[0] = ( C.elements[0][1] + C.elements[1][0] ) / ( 4 * q.elements[1] );
    q.elements[2] = ( C.elements[1][2] + C.elements[2][1] ) / ( 4 * q.elements[1] );
    q.elements[3] = ( C.elements[0][2] - C.elements[2][0] ) / ( 4 * q.elements[1] );
    return;
  }

  if( max == 2 ){
    q.elements[2] = 0.5 * sqrt( b[2] );
    q.elements[0] = ( C.elements[0][2] + C.elements[2][0] ) / ( 4 * q.elements[2] );
    q.elements[1] = ( C.elements[1][2] + C.elements[2][1] ) / ( 4 * q.elements[2] );
    q.elements[3] = ( C.elements[1][0] - C.elements[0][1] ) / ( 4 * q.elements[2] );
    return;
  }
  EP_Normalize(q);
}

void EP_Normalize(ColMatrix& q){
  double one = 1.0/sqrt(q.elements[0]*q.elements[0] + q.elements[1]*q.elements[1] + q.elements[2]*q.elements[2] + q.elements[3]*q.elements[3]);  
  q.elements[0] = one*q.elements[0];
  q.elements[1] = one*q.elements[1];
  q.elements[2] = one*q.elements[2];
  q.elements[3] = one*q.elements[3];
}




void EPdotdot_udot(ColMatrix& Audot, ColMatrix& Aqdot, ColMatrix& Aq, ColMatrix& Aqddot){  
  int num=Audot.GetNumRows();
  if (3<num){
	  for (int i=4; i<=num; i++){
		  Aqddot.elements[i]=Audot.elements[i-1];
	  }
  }
    
  double AA;
  AA=Aqdot.elements[0]*Aqdot.elements[0]+Aqdot.elements[1]*Aqdot.elements[1]+Aqdot.elements[2]*Aqdot.elements[2]+Aqdot.elements[3]*Aqdot.elements[3];
  
  Aqddot.elements[0] = 0.5 *(Aq.elements[3]*Audot.elements[0] - Aq.elements[2]*Audot.elements[1] + Aq.elements[1]*Audot.elements[2]-2*Aq.elements[0]*AA);
  
  Aqddot.elements[1] = 0.5 *(Aq.elements[2]*Audot.elements[0] + Aq.elements[3]*Audot.elements[1] - Aq.elements[0]*Audot.elements[2]-2*Aq.elements[1]*AA);
  
  Aqddot.elements[2] = 0.5 *(-Aq.elements[1]*Audot.elements[0] + Aq.elements[0]*Audot.elements[1] + Aq.elements[3]*Audot.elements[2]-2*Aq.elements[2]*AA);
  
  Aqddot.elements[3] = -0.5 *(Aq.elements[0]*Audot.elements[0] + Aq.elements[1]*Audot.elements[1] + Aq.elements[2]*Audot.elements[2]+2*Aq.elements[3]*AA);
  
}


void qdot_to_u(ColMatrix& q, ColMatrix& u, ColMatrix& qdot){
	EP_Normalize(q);
	int num=qdot.GetNumRows();
      	if (4<num){
      		for (int i=5; i<=num; i++){
      			u.elements[i-2]=qdot.elements[i-1];
      		}
      	}		
      	u.elements[0]=2*(q.elements[3]*qdot.elements[0]+q.elements[2]*qdot.elements[1]-q.elements[1]*qdot.elements[2]-q.elements[0]*qdot.elements[3]);
      	u.elements[1]=2*(-q.elements[2]*qdot.elements[0]+q.elements[3]*qdot.elements[1]+q.elements[0]*qdot.elements[2]-q.elements[1]*qdot.elements[3]);
      	u.elements[2]=2*(q.elements[1]*qdot.elements[0]-q.elements[0]*qdot.elements[1]+q.elements[3]*qdot.elements[2]-q.elements[2]*qdot.elements[3]);
  
}


