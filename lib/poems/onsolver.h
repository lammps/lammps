/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: onsolver.h                                              *
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


#ifndef ONSOLVER_H
#define ONSOLVER_H

#include "solver.h"
#include "onbody.h"

namespace POEMS {
class ColMatrix;
class Matrix;

class OnSolver : public Solver {
  OnBody inertialframe;
  int numbodies;
  OnBody** bodyarray;
  ColMatrix** q;
  ColMatrix** qdot;
                 ColMatrix** qdotdot;
  ColMatrix** u;
  ColMatrix** udot;



  void DeleteModel();
  int CreateTopologyArray(int i, OnBody* body);
  void CreateStateMatrixMaps();
            void GetType();
public:
  OnSolver();
  ~OnSolver();
  void CreateModel();
  void Solve(double time, Matrix& FF);
};
}
#endif
