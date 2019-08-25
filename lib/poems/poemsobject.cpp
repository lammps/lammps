/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: poemsobject.cpp                                         *
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
 

#include "poemsobject.h"
#include <cstring>

POEMSObject::POEMSObject(){
  name = 0;
  ChangeName((const char*)"unnamed");
  ID = -1;
}

POEMSObject::~POEMSObject(){
  delete [] name;
}

void POEMSObject::ChangeName(const char* newname){
  delete [] name;
  name = new char[strlen(newname)+1];
  strcpy(name,newname);
}

char* POEMSObject::GetName(){
  return name;
}

int POEMSObject::GetID(){
  return ID;
}

void POEMSObject::SetID(int id){
  ID = id;
}
