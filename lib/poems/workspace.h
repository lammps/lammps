/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: workspace.h                                             *
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


#ifndef WORKSPACE_H
#define WORKSPACE_H

#include "matrices.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>   
#include <iomanip>  
#include <vector>


class System;
class Solver;

struct SysData{
	System * system;
	int solver;
	int integrator;
};

class Workspace {
	SysData * system; // the multibody systems data
	int currentIndex;
	int maxAlloc;
	
public:
     Workspace();
     ~Workspace();
     
     double Thalf;
     double Tfull;
     double ConFac;
     double KE_val;
	  int FirstTime;
     
     bool LoadFile(char* filename);
     
     bool SaveFile(char* filename, int index = -1);

     System* GetSystem(int index = -1);
     
     void AddSolver(Solver* s, int index = -1);        
    

     void LobattoOne(double **&xcm, double **&vcm,double **&omega,double **&torque, double **&fcm, double **&ex_space, double **&ey_space, double **&ez_space);	
     
     void LobattoTwo(double **&vcm,double **&omega,double **&torque, double **&fcm);	 
     
 
     bool MakeSystem(int& nbody, double *&masstotal, double **&inertia, double **&xcm, double **&vcm, double **&omega, double **&ex_space, double **&ey_space, double **&ez_space, int &njoint, int **&jointbody, double **&xjoint, int& nfree, int*freelist, double dthalf, double dtv, double tempcon, double KE);
     																																							
	  
	  bool SaveSystem(int& nbody, double *&masstotal, double **&inertia, double **&xcm, double **&xjoint, double **&vcm, double **&omega, double **&ex_space, double **&ey_space, double **&ez_space, double **&acm, double **&alpha, double **&torque, double **&fcm, int **&jointbody, int &njoint);
	  																		
	 bool MakeDegenerateSystem(int& nfree, int*freelist, double *&masstotal, double **&inertia, double **&xcm, double **&vcm, double **&omega, double **&ex_space, double **&ey_space, double **&ez_space);
     int getNumberOfSystems();
     
     void SetLammpsValues(double dtv, double dthalf, double tempcon);
     void SetKE(int temp, double SysKE);
	  
	  void RKStep(double **&xcm, double **&vcm,double **&omega,double **&torque, double **&fcm, double **&ex_space, double **&ey_space, double **&ez_space);	
	  
	  void WriteFile(char* filename);

private:
	void allocateNewSystem(); //helper function to handle vector resizing and such for the array of system pointers
};

#endif
