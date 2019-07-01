/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: eulerparameters.h                                       *
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

#ifndef EULERPARAMETERS_H
#define EULERPARAMETERS_H

class ColMatrix;
class Mat3x3;

void EP_Derivatives(ColMatrix& q, ColMatrix& u, ColMatrix& qdot);

void EP_Transformation(ColMatrix& q, Mat3x3& C);

void EP_FromTransformation(ColMatrix& q, Mat3x3& C);

void EP_Normalize(ColMatrix& q);

void EPdotdot_udot(ColMatrix& Audot, ColMatrix& Aqdot, ColMatrix& Aq,ColMatrix& Aqddot);

void qdot_to_u(ColMatrix& q, ColMatrix& u, ColMatrix& qdot);

#endif

