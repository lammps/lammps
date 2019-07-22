/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: solver.h                                                *
 *      AUTHORS: See Author List                                           * 
 *      GRANTS: See Grants List                                            *
 *      COPYRIGHT: (C) 2005 by Authors as listed in Author's List          *
 *      LICENSE: Please see License Agreement                              *
 *      DOWNLOAD: Free at www.rpi.edu/~anderk5                             *
 *      ADMINISTRATOR: Prof. Kurt Anderson                                 *
 *                     Computational Dynamics Lab                          *
 *                     Rensselaer Polytechnic Institute                    *
 *                     110 8th St. Troy NY 12180                           * 
 *      CONTACT:       anderk5@rpi.edu                                     *
 *_________________________________________________________________________*/

#ifndef SOLVER_H
#define SOLVER_H
#include <fstream>
#include "colmatmap.h"
#include "matrices.h"
#include "defines.h"

class System;

class Solver{
protected:
  System* system;

  
  double time;
  ColMatMap state;
  ColMatMap statedot;
  ColMatMap statedoubledot;  

  SolverType type;

virtual void ComputeForces();

public:
  Solver();
  virtual ~Solver();
  void SetSystem(System* s);
  SolverType GetSolverType();
  static Solver * GetSolver(SolverType solverToMake);

  virtual void DeleteModel() = 0;
  virtual void CreateModel() = 0;
  virtual void Solve(double time, Matrix& FF) = 0;  
  
    
    
  ColMatMap* GetState();
  ColMatMap* GetStateDerivative();
  ColMatMap* GetStateDerivativeDerivative();
};

#endif
