/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: body.cpp                                                *
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
 
#include "bodies.h"
#include "point.h"

using namespace std;

Body::Body()
{
  inertia.Zeros();
  mass = 0;
  a_t.Zeros();
}

Body::~Body()
{
  points.DeleteValues();
}

bool Body::ReadIn(istream& in){
  return (ReadInBodyData(in) && ReadInPoints(in));
}

void Body::WriteOut(ostream& out){
  // write out header <type> <name>
  out << GetType() << ' ' << GetName() << endl;

  // write out body specific data
  WriteOutBodyData(out);

  // write out points
  WriteOutPoints(out);

}

bool Body::ReadInPoints(istream& in){
  // get numfixed points
  int numpoints;
  int index;
  Point* point;
  int pointtype;
  char pointname[256];

  in >> numpoints;

  for(int i=points.GetNumElements();i<numpoints;i++){
    // error check
    in >> index;
    if(index != i){
      cerr << "Invalid file format" << endl;
      return false;
    }

    in >> pointtype >> pointname;
    point = NewPoint(pointtype);

    if(!point){
      cerr << "Unrecognized point type '" << pointtype << endl;
      return false;
    }

    // add the point
    AddPoint(point);

    // set generic point info
    point->ChangeName(pointname);

    // read in the rest of its data
    if(!point->ReadIn(in)) return false;
  }
  return true;
}

void Body::WriteOutPoints(std::ostream& out){
  int numpoints = points.GetNumElements();

  out << numpoints << endl;

  // list element pointer
  ListElement<Point>* ele = points.GetHeadElement();

  // set the ID of each point
  for(int i=0;i<numpoints;i++){
    ele->value->SetID(i);
    out << i << ' ';
    ele->value->WriteOut(out);

    // increment list element pointer
    ele = ele->next;
  }
  out << endl;
}

Point* Body::GetPoint(int p) {
  return points(p);
}

void Body::AddJoint(Joint* joint){
  joints.Append(joint);
}

void Body::AddPoint(Point* point){
  points.Append(point);
}

//
// global body functions
//

Body* NewBody(int type){
  switch( BodyType(type) )
    {
      case INERTIALFRAME :  // The inertial reference frame
        return new InertialFrame;
      case RIGIDBODY :  // A Rigid Body
        return new RigidBody;
      case PARTICLE :  // A Particle
        return new Particle;
      default  :  // error
        return 0;
    }
}
