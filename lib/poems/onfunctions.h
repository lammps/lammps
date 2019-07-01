/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: onfunction.h                                            *
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

#ifndef ONFUNCTIONS_H
#define ONFUNCTIONS_H

class Mat3x3;
class Mat6x6;
class Vect3;
class Vect6;

void OnPopulateSVect(Vect3& angular, Vect3& linear, Vect6& sV);
void OnPopulateSC(Vect3& gamma, Mat3x3& C, Mat6x6& SC);
void OnPopulateSI(Mat3x3& inertia, double mass, Mat6x6& sI);


void Create_Map(int MM);
int ICELL(int IX,int IY,int IZ, int MM);

#endif
