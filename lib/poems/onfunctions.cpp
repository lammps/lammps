/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: onfunction.cpp                                          *
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
 

#include "onfunctions.h"
#include "matrixfun.h"
#include <iostream>
using namespace std;

// friend of Vect3 & Vect6
void OnPopulateSVect(Vect3& angular, Vect3& linear, Vect6& sV){
  sV.elements[0] = angular.elements[0];
  sV.elements[1] = angular.elements[1];
  sV.elements[2] = angular.elements[2];
  sV.elements[3] = linear.elements[0];
  sV.elements[4] = linear.elements[1];
  sV.elements[5] = linear.elements[2];
}

// friend of Vect3, Mat3x3, & Mat6x6
void OnPopulateSC(Vect3& gamma, Mat3x3& C, Mat6x6& SC){
  // the block diagonals
  
  // the gamma cross with transform
  Mat3x3 temp; Mat3x3 temp2;
  SC.Zeros();
  temp.Zeros();  
  temp2.Zeros();
  //FastTMult(C,gamma,temp);     
  temp(1,2)= -gamma(3); temp(1,3)= gamma(2); temp(2,1)= gamma(3);
  temp(2,3)= -gamma(1); temp(3,1)= -gamma(2); temp(3,2)= gamma(1);  
  FastMult(temp,C,temp2);
  
  SC(1,4)=temp2(1,1);   SC(2,4)=temp2(2,1);   SC(3,4)=temp2(3,1);
  SC(1,5)=temp2(1,2);   SC(2,5)=temp2(2,2);   SC(3,5)=temp2(3,2);
  SC(1,6)=temp2(1,3);   SC(2,6)=temp2(2,3);   SC(3,6)=temp2(3,3);
  
  SC(1,1)=C(1,1);   SC(2,1)=C(2,1);   SC(3,1)=C(3,1);
  SC(1,2)=C(1,2);   SC(2,2)=C(2,2);   SC(3,2)=C(3,2);
  SC(1,3)=C(1,3);   SC(2,3)=C(2,3);   SC(3,3)=C(3,3);

  SC(4,4)=C(1,1);   SC(5,4)=C(2,1);   SC(6,4)=C(3,1);
  SC(4,5)=C(1,2);   SC(5,5)=C(2,2);   SC(6,5)=C(3,2);
  SC(4,6)=C(1,3);   SC(5,6)=C(2,3);   SC(6,6)=C(3,3);  
  
  }

// friend of Mat3x3 & Mat6x6
void OnPopulateSI(Mat3x3& inertia, double mass, Mat6x6& sI){
  
	sI(4,4)=mass; sI(5,5)=mass; sI(6,6)=mass;
	sI(1,1)=inertia(1,1); sI(1,2)=inertia(1,2); sI(1,3)=inertia(1,3);
	sI(2,1)=inertia(2,1); sI(2,2)=inertia(2,2); sI(2,3)=inertia(2,3);
	sI(3,1)=inertia(3,1); sI(3,2)=inertia(3,2); sI(3,3)=inertia(3,3);	
}
