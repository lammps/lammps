/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: nrom.cpp                                               *
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

#include <cmath>
#include "norm.h"

double Magnitude(ColMatrix& A){
  double G;
  G = 0;
  for (int i=1;i<=A.GetNumRows();i++) G += A.Get(i)*A.Get(i);
  G = sqrt(G);
  return G;
}

double Magnitude(RowMatrix& A){
  double G;
  G = 0;
  for (int i=1;i<=A.GetNumCols();i++) G += A.Get(i)*A.Get(i);
  G = sqrt(G);
  return G;
}

double Magnitude(Vect3& A){
  double G;
  G = 0;
  for (int i=1;i<=3;i++) G += A.Get(i)*A.Get(i);
  G = sqrt(G);
  return G;
}

double Magnitude(Vect4& A){
  double G;
  G = 0;
  for (int i=1;i<=4;i++) G += A.Get(i)*A.Get(i);
  G = sqrt(G);
  return G;
}

double Magnitude(Vect6& A){
  double G;
  G = 0;
  for (int i=1;i<=6;i++) G += A.Get(i)*A.Get(i);
  G = sqrt(G);
  return G;
}
