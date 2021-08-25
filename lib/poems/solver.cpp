/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: solver.cpp                                              *
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

#include "solver.h"
#include "system.h"
#include "matrices.h"

Solver::Solver(){

}

Solver::~Solver(){
}

void Solver::SetSystem(System* s){
  system = s;
  CreateModel();
}

void Solver::ComputeForces(){
system->ComputeForces();
}

SolverType Solver::GetSolverType()
{
  return type;
}

Solver * Solver::GetSolver(SolverType solverToMake) //returning a pointer to a new Solver object of the appropriate type
{
  switch((int)solverToMake)
  {
    case ONSOLVER: return new OnSolver();
    default: return nullptr;
  }
}

ColMatMap* Solver::GetState(){
  return &state;
}

ColMatMap* Solver::GetStateDerivative(){
  return &statedot;
}

ColMatMap* Solver::GetStateDerivativeDerivative(){
  return &statedoubledot;
}
