/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: onbody.h                                              *
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

#ifndef ONBODY_H
#define ONBODY_H

#include "poemslist.h"
#include "matrix.h"
#include "vect6.h"
#include "mat6x6.h"

// emumerated type
enum Direction {
  BACKWARD = 0,
  FORWARD= 1
};

class Body;
class InertialFrame;
class Joint;
class OnSolver;

class OnBody {
       Body* system_body;
       Joint* system_joint;
       OnBody* parent;
       List<OnBody> children;

       Direction joint_dir;
       void (Joint::*kinfun)();          // kinematics function
       void (Joint::*updatesP)(Matrix&); // sP update function
       Vect3* gamma;                     // pointer to gamma vector
       Mat3x3* pk_C_k;                   // pointer to transformation


       Mat6x6 sI;      // spatial inertias
       Mat6x6 sIhat;   // recursive spatial inertias
       Mat6x6 sSC;     // spatial shift
       Mat6x6 sT;      // spatial triangularization

       Vect6 sF;       // spatial forces
       Vect6 sFhat;    // recursive spatial forces
       Vect6 sAhat;    // recursive spatial acceleration

       Matrix sP;      // spatial partial velocities
       Matrix sM;      // triangularized mass matrix diagonal elements
       Matrix sMinv;   // inverse of sM
       Matrix sPsMinv;
       Matrix sIhatsP;

       // states and state derivatives
       ColMatrix* q;
       ColMatrix* u;
       ColMatrix* qdot;
       ColMatrix* udot;
       ColMatrix* qdotdot;

       ColMatrix* r;
       ColMatrix* acc;
       ColMatrix* ang;

  // friend classes
  friend class OnSolver;


public:
       OnBody();
       ~OnBody();
       int RecursiveSetup(InertialFrame* basebody);
       int RecursiveSetup(int ID, OnBody* parentbody, Joint* sys_joint);
       void RecursiveKinematics();
       void RecursiveTriangularization();
       void RecursiveForwardSubstitution();
       Mat3x3 GetN_C_K();
       Vect3 LocalCart();
       int GetBodyID();
       void CalculateAcceleration();
       void Setup();
       void SetupInertialFrame();
       void LocalKinematics();
                 void LocalTriangularization(Vect3& Torque, Vect3& Force);
       void LocalForwardSubstitution();
};

#endif
