/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: joint.cpp                                              *
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


#include "joints.h"
#include "body.h"
#include "point.h"
#include <string>
#include "matrixfun.h"
#include "fastmatrixops.h"
#include <iomanip>

using namespace std;

Joint::Joint(){
  body1 = body2 = nullptr;
  point1 = point2 = nullptr;
  pk_C_ko.Identity();
  pk_C_k.Identity();
}

Joint::~Joint() = default;

void Joint::SetBodies(Body* b1, Body* b2){
  body1 = b1;
  body2 = b2;
}

void Joint::SetPoints(Point* p1, Point* p2){
  point1 = p1;
  point2 = p2;
}

int Joint::GetBodyID1(){
  return body1->GetID();
}

int Joint::GetBodyID2(){
  return body2->GetID();
}

int Joint::GetPointID1(){
  return point1->GetID();
}

int Joint::GetPointID2(){
  return point2->GetID();
}

bool Joint::ReadIn(std::istream& in){
  in >>setprecision(20)>> qo >> setprecision(20)>>qdoto >> setprecision(20)>>pk_C_ko;
  q = qo;
  qdot = qdoto;
  EP_Normalize(q);

  return ReadInJointData(in);
}

void Joint::ResetQdot(){
  //EP_Derivatives(q,u,qdot);
  qdot_to_u(q,u,qdot);
}

void Joint::ResetQ(){
  EP_Normalize(q);
}

void Joint::SetInitialState(ColMatrix& a, ColMatrix& adot){
  if( (qo.GetNumRows() != a.GetNumRows()) || (qdoto.GetNumRows() != adot.GetNumRows()) ){
    cout<<qo.GetNumRows()<<"  "<<a.GetNumRows()<<" "<<qdoto.GetNumRows()<<" "<<adot.GetNumRows()<<endl;
    cerr << "ERROR::Illegal matrix size for initial condition" << endl;
    exit(1);
  }
  qo = a;
  qdoto = adot;
  EP_Normalize(qo);
  q=qo;           //////////////////////////Check this ...added May 14 2005
  qdot=qdoto;     //////////////////////////Check this ...added May 14 2005
}

void Joint::SetZeroOrientation(VirtualMatrix& C){
  pk_C_ko = C;
}


void Joint::WriteOut(std::ostream& out){
  out << GetType() << ' ' << GetName() << ' ';
  out << GetBodyID1() << ' ' << GetBodyID2() << ' ';
  out << GetPointID1() << ' ' << GetPointID2() << endl;
  out <<setprecision(16)<< qo <<setprecision(16)<< qdoto <<setprecision(16)<< pk_C_ko;
  WriteOutJointData(out);
  out << endl;
}

ColMatrix* Joint::GetQ(){
  return &q;
}

ColMatrix* Joint::GetU(){
  return &u;
}

ColMatrix* Joint::GetQdot(){
  return &qdot;
}

ColMatrix* Joint::GetUdot(){
  return &udot;
}

ColMatrix* Joint::GetQdotdot(){
  return &qdotdot;
}


void Joint::DimQandU(int i){
  DimQandU(i,i);
}

void Joint::DimQandU(int i, int j){
  qo.Dim(i);
  q.Dim(i);
  qdot.Dim(i);
  qdoto.Dim(i);
  uo.Dim(j);
  u.Dim(j);
  udot.Dim(j);
  qdotdot.Dim(i);

  // zero them
  qo.Zeros();
  qdoto.Zeros();
  q.Zeros();
  qdot.Zeros();
  uo.Zeros();
  u.Zeros();
  udot.Zeros();
  qdotdot.Zeros();

}

Body* Joint::GetBody1(){
  return body1;
}

Body* Joint::GetBody2(){
  return body2;
}

Body* Joint::OtherBody(Body* body){
  if(body1 == body) return body2;
  if(body2 == body) return body1;
  return nullptr;
}

Vect3* Joint::GetR12(){
  return &r12;
}

Vect3* Joint::GetR21(){
  return &r21;
}

Mat3x3* Joint::Get_pkCk(){
  return &pk_C_k;
}


Mat3x3* Joint::Get_kCpk(){
  return &k_C_pk;
}

Matrix Joint::GetForward_sP(){
  cerr << "ERROR: Forward Spatial Partial Velocity is not supported for joint type " << GetType() << endl;
  exit(0);
}

Matrix Joint::GetBackward_sP(){
  cerr << "ERROR: Backward Spatial Partial Velocity is not supported for joint type " << GetType() << endl;
  exit(0);
}

void Joint::UpdateForward_sP(Matrix& sP){
  cerr << "WARNING: Using default Update sP procedure" << endl;
  sP = GetForward_sP();
}

void Joint::UpdateBackward_sP(Matrix& sP){
  cerr << "WARNING: Using default Update sP procedure" << endl;
  sP = GetBackward_sP();
}

void Joint::ComputeForwardTransforms(){
  ComputeLocalTransform();
//  k_C_pk = T(pk_C_k);
  FastAssignT(pk_C_k, k_C_pk);
  ComputeForwardGlobalTransform();
}

void Joint::ComputeBackwardTransforms(){
  ComputeLocalTransform();
  // k_C_pk = pk_C_k^T
  FastAssignT(pk_C_k, k_C_pk);
  ComputeBackwardGlobalTransform();
}

void Joint::ComputeForwardGlobalTransform(){
  // body2->n_C_k = body1->n_C_k * pk_C_k;
  FastMult(body1->n_C_k,pk_C_k,body2->n_C_k);
  }

void Joint::ComputeBackwardGlobalTransform(){
  // body1->n_C_k = body2->n_C_k * T(pk_C_k);
  FastMultT(body2->n_C_k,pk_C_k,body1->n_C_k);
}


//
// global joint functions
//

Joint* NewJoint(int type){
  switch( JointType(type) )
  {
    case FREEBODYJOINT : return new FreeBodyJoint;
    case REVOLUTEJOINT : return new RevoluteJoint;
    case PRISMATICJOINT : return new PrismaticJoint;
    case SPHERICALJOINT : return new SphericalJoint;
   case BODY23JOINT : return new Body23Joint;
   case MIXEDJOINT : return new MixedJoint;
    default : return nullptr; // error
  }
}
