/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: system.h                                                *
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

 
#ifndef SYSTEM_H
#define SYSTEM_H


#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>   
#include <iomanip>  

#include "poemslist.h"
#include "matrices.h"
#include "workspace.h"
#include "matrixfun.h"
#include "onsolver.h"
#include "system.h"
#include "inertialframe.h"
#include "rigidbody.h"
#include "revolutejoint.h"
#include "fixedpoint.h"
#include "freebodyjoint.h"
#include "sphericaljoint.h"
#include "body23joint.h"
#include "mixedjoint.h"
#include "eulerparameters.h"
#include "matrices.h"
#include "norm.h"


     class Body;
     class Joint;

     class System{
	 private:
		int * mappings;
		
     public:
       double time;
       List<Body> bodies;
       List<Joint> joints;
  
       System();
       ~System();
       void Delete();
       
       int GetNumBodies();
       
	   int * GetMappings();

       void AddBody(Body* body);
       
       void AddJoint(Joint* joint);
       
       void SetTime(double t);
       
       double GetTime();
       
       void ComputeForces();
       
       bool ReadIn(std::istream& in);
       
       void WriteOut(std::ostream& out);
       
       void ClearBodyIDs();
       
       void ClearJointIDs();              

       void Create_System_LAMMPS(int numbodies, double *mass,double **inertia, double ** xcm, double ** xjoint,double **vh1,double **omega,double **ex_space, double **ey_space, double **ez_space, int b, int * mapping, int count);
       
       void Create_DegenerateSystem(int& nfree, int*freelist, double *&masstotal, double **&inertia, double **&xcm, double **&vcm, double **&omega, double **&ex_space, double **&ey_space, double **&ez_space);

};

#endif
