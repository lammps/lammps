/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: inertialframe.cpp                                       *
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

#include "inertialframe.h"
#include "fixedpoint.h"

using namespace std;

InertialFrame::InertialFrame(){
  gravity.Zeros();
  n_C_k.Identity();

  r.Zeros();
  v.Zeros();
  v_k.Zeros();
  a.Zeros();
  omega.Zeros();
  omega_k.Zeros();
  alpha.Zeros();
  alpha_t.Zeros();
}
InertialFrame::~InertialFrame(){
}

BodyType InertialFrame::GetType(){
  return INERTIALFRAME;
}

Vect3 InertialFrame::GetGravity(){
  return gravity;
}

void InertialFrame::SetGravity(Vect3& g){
  gravity = g;
}

bool InertialFrame::ReadInBodyData(istream& in){
  in >> gravity;
  return true;
}

void InertialFrame::WriteOutBodyData(ostream& out){
  out << gravity;
}
