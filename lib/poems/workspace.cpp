/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: workspace.cpp                                           *
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


#include "workspace.h"
#include "system.h"
#include "solver.h"
#include "SystemProcessor.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>


using namespace std;

void Workspace::allocateNewSystem() {
  currentIndex++;
  if(currentIndex < maxAlloc)
  {
    system[currentIndex].system = new System;
  }
  else
  {
    maxAlloc = (maxAlloc + 1) * 2;

    SysData * tempPtrSys = new SysData[maxAlloc];
    for(int i = 0; i < currentIndex; i++)
    {
      tempPtrSys[i] = system[i];
    }
    delete [] system;
    system = tempPtrSys;
    system[currentIndex].system = new System;
  }
}

Workspace::Workspace(){
  currentIndex = -1;
  maxAlloc = 0;
  system = nullptr;
}

Workspace::~Workspace(){
  for(int i = 0; i <= currentIndex; i++){
    delete system[i].system;
  }
  delete [] system;
}


bool Workspace::LoadFile(char* filename){
  bool ans;
  ifstream file;
  file.open(filename, ifstream::in);
  if( !file.is_open() ){
    cerr << "File '" << filename << "' not found." << endl;
    return false;
  }
  allocateNewSystem();
  ans = system[currentIndex].system->ReadIn(file);
  file.close();

  return ans;
}

void Workspace::SetLammpsValues(double dtv, double dthalf, double tempcon){
  Thalf = dthalf;
  Tfull = dtv;
  ConFac = tempcon;
}


bool Workspace::MakeSystem(int& nbody, double *&masstotal, double **&inertia, double **&xcm, double **&vcm, double **&omega, double **&ex_space, double **&ey_space, double **&ez_space, int &njoint, int **&jointbody, double **&xjoint, int& nfree, int*freelist, double dthalf, double dtv, double tempcon, double KE){

  SetLammpsValues(dtv, dthalf, tempcon);

if(njoint){
  SystemProcessor sys;
  sys.processArray(jointbody, njoint);
  List<POEMSChain> * results = sys.getSystemData();

  int numsyschains = results->GetNumElements();
  int headvalue = 0;
  List<POEMSChain> * newresults = results;
  ListElement<POEMSChain> * tempNode = results->GetHeadElement();
  int stop = 1;
  int counter = 1;
  for(int n = 1; n<=numsyschains; n++){
    while(stop){
      if ( (*(tempNode->value->listOfNodes.GetHeadElement()->value) == (headvalue+1) ) || (*(tempNode->value->listOfNodes.GetTailElement()->value) == (headvalue+1) ) ) {
      newresults->Append(tempNode->value);
      headvalue = headvalue + tempNode->value->listOfNodes.GetNumElements();
      tempNode = results->GetHeadElement();
      stop = 0;
      counter ++;
    }
    else{
      tempNode = tempNode->next;
    }
    }
    stop=1;
  }

  ListElement<POEMSChain> * TNode = newresults->GetHeadElement();
  ListElement<POEMSChain> * TTNode = newresults->GetHeadElement();
  for(int kk=1; kk<=numsyschains; kk++){
    if(kk!=numsyschains)
      TTNode = TNode->next;
    newresults->Remove(TNode);
    if(kk!=numsyschains)
      TNode = TTNode;
  }
  ListElement<POEMSChain> * NodeValue = newresults->GetHeadElement();
  int count = 0;
  int * array;
  int ** arrayFromChain;
  int numElementsInSystem;
  int ttk = 0;


  while(NodeValue != nullptr) {
    array = new int[NodeValue->value->listOfNodes.GetNumElements()];
    arrayFromChain = NodeValue->value->listOfNodes.CreateArray();
    numElementsInSystem = NodeValue->value->listOfNodes.GetNumElements();
    for(counter = 0; counter < numElementsInSystem; counter++){
      array[counter] = *arrayFromChain[counter];
    }

    SetKE(1,KE);
    allocateNewSystem();
    system[currentIndex].system->Create_System_LAMMPS(nbody,masstotal,inertia,xcm,xjoint,vcm,omega,ex_space,ey_space,ez_space, numElementsInSystem, array, count);

    system[currentIndex].solver = ONSOLVER;
    ttk = ttk + 1;
    count = ttk;
    delete [] array;
    delete [] arrayFromChain;
    NodeValue= NodeValue->next;
  }
}
if(nfree){
  MakeDegenerateSystem(nfree,freelist,masstotal,inertia,xcm,vcm,omega,ex_space,ey_space,ez_space);
}
  return true;
}


bool Workspace::MakeDegenerateSystem(int& nfree, int*freelist, double *&masstotal, double **&inertia, double **&xcm, double **&vcm, double **&omega, double **&ex_space, double **&ey_space, double **&ez_space){
  allocateNewSystem();
  system[currentIndex].system->Create_DegenerateSystem(nfree,freelist,masstotal,inertia,xcm,vcm,omega,ex_space,ey_space,ez_space);
  system[currentIndex].solver = ONSOLVER;
  return true;
}

bool Workspace::SaveFile(char* filename, int index){
  if(index < 0){
    index = currentIndex;
  }
  ofstream file;

  file.open(filename, ofstream::out);

  if( !file.is_open() ){
    cerr << "File '" << filename << "' could not be opened." << endl;
    return false;
  }
  if(index >= 0 && index <= currentIndex){
    system[index].system->WriteOut(file);
  }
  else {
    cerr << "Error, requested system index " << index << ", minimum index 0 and maximum index " << currentIndex << endl;
  }
  file.close();
  return true;
}


System* Workspace::GetSystem(int index){
  if(index <= currentIndex){
    if(index >= 0){
      return system[index].system;
    }
    else{
      return system[currentIndex].system;
    }
  }
  else{
    return nullptr;
  }
}

void Workspace::AddSolver(Solver* s, int index){
  if(currentIndex >= index){
    if(index >= 0){
      system[index].solver = (int)s->GetSolverType();
    }
    else{
      system[currentIndex].solver = (int)s->GetSolverType();
    }
  }
  else{
    cout << "Error adding solver to index " << index << endl;
  }
}

int Workspace::getNumberOfSystems(){
  return currentIndex + 1;
}



void Workspace::SetKE(int temp, double SysKE){
KE_val = SysKE;
FirstTime =temp;
}


void Workspace::LobattoOne(double **&xcm, double **&vcm,double **&omega,double **&torque, double **&fcm, double **&ex_space, double **&ey_space, double **&ez_space){

  int numsys = currentIndex;
  int numbodies;
  double time = 0.0;
  int * mappings;
  double SysKE=0.0;

  for (int i = 0; i <= numsys; i++){
    mappings = system[i].system->GetMappings();
    numbodies = system[i].system->GetNumBodies() - 1;
    Matrix FF(6,numbodies);
    FF.Zeros();
    for (int j=1; j<=numbodies; j++){
      FF(1,j)  = torque[mappings[j - 1]-1][0]*ConFac;
      FF(2,j)  = torque[mappings[j - 1]-1][1]*ConFac;
      FF(3,j)  = torque[mappings[j - 1]-1][2]*ConFac;

      FF(4,j)  = fcm[mappings[j - 1]-1][0]*ConFac;
      FF(5,j)  = fcm[mappings[j - 1]-1][1]*ConFac;
      FF(6,j)  = fcm[mappings[j - 1]-1][2]*ConFac;
    }

    //------------------------------------//
    // Get a solver and solve that system.
    Solver * theSolver = Solver::GetSolver((SolverType)system[i].solver);
    theSolver->SetSystem(system[i].system);
    theSolver->Solve(time, FF);


    theSolver->Solve(time, FF);
    ColMatrix tempx = *((*theSolver).GetState());
    ColMatrix tempv = *((*theSolver).GetStateDerivative());
    ColMatrix tempa = *((*theSolver).GetStateDerivativeDerivative());


    for(int numint =0 ; numint<3; numint++){
      theSolver->Solve(time, FF);
      tempa = *((*theSolver).GetStateDerivativeDerivative());
      *((*theSolver).GetStateDerivative())= tempv + Thalf*tempa;
    }

    ColMatrix TempV= *((*theSolver).GetStateDerivative());
    *((*theSolver).GetState())= tempx + Tfull*TempV;

    int numjoints = system[i].system->joints.GetNumElements();
    for(int k = 0; k < numjoints; k++){
      system[i].system->joints(k)->ForwardKinematics();
    }

    for(int k = 0; k < numbodies; k++){
      Vect3 temp1 =system[i].system->bodies(k+1)->r;
      Vect3 temp2 =system[i].system->bodies(k+1)->v;
      Vect3 temp3 =system[i].system->bodies(k+1)->omega;
      Mat3x3 temp4 =system[i].system->bodies(k+1)->n_C_k;
      for(int m = 0; m < 3; m++){
        xcm[mappings[k]-1][m]   = temp1(m+1);
        vcm[mappings[k]-1][m]   = temp2(m+1);
        omega[mappings[k]-1][m] = temp3(m+1);
        ex_space[mappings[k]-1][m] = temp4(m+1,1);
        ey_space[mappings[k]-1][m] = temp4(m+1,2);
        ez_space[mappings[k]-1][m] = temp4(m+1,3);
      }
      SysKE = SysKE + system[i].system->bodies(k+1)->KE;
    }
    delete theSolver;
  }
}

void Workspace::LobattoTwo(double **&vcm,double **&omega,double **&torque, double **&fcm){
  int numsys = currentIndex;
  int numbodies;
  double time = 0.0;
  int * mappings;
  double SysKE =0.0;
  for (int i = 0; i <= numsys; i++){
    mappings = system[i].system->GetMappings();
    numbodies = system[i].system->GetNumBodies() - 1;
    Matrix FF(6,numbodies);

    for (int j=1; j<=numbodies; j++){
      FF(1,j)  = torque[mappings[j - 1]-1][0]*ConFac;
      FF(2,j)  = torque[mappings[j - 1]-1][1]*ConFac;
      FF(3,j)  = torque[mappings[j - 1]-1][2]*ConFac;
      FF(4,j)  = fcm[mappings[j - 1]-1][0]*ConFac;
      FF(5,j)  = fcm[mappings[j - 1]-1][1]*ConFac;
      FF(6,j)  = fcm[mappings[j - 1]-1][2]*ConFac;
    }

    //------------------------------------//
    // Get a solver and solve that system.
    Solver * theSolver = Solver::GetSolver((SolverType)system[i].solver);
    theSolver->SetSystem(system[i].system);
    theSolver->Solve(time, FF);
    ColMatrix tempv = *((*theSolver).GetStateDerivative());
    ColMatrix tempa = *((*theSolver).GetStateDerivativeDerivative());
    *((*theSolver).GetStateDerivative()) = tempv + Thalf*tempa;


    int numjoints = system[i].system->joints.GetNumElements();
    for(int k = 0; k < numjoints; k++){
      system[i].system->joints(k)->ForwardKinematics();
    }

    for(int k = 0; k < numbodies; k++){
      Vect3 temp1 =system[i].system->bodies(k+1)->r;
      Vect3 temp2 =system[i].system->bodies(k+1)->v;
      Vect3 temp3 =system[i].system->bodies(k+1)->omega;
      SysKE = SysKE + system[i].system->bodies(k+1)->KE;
      for(int m = 0; m < 3; m++){
        vcm[mappings[k]-1][m]   = temp2(m+1);
        omega[mappings[k]-1][m] = temp3(m+1);
      }
    }
    delete theSolver;
  }
}



  void Workspace::RKStep(double **&xcm, double **&vcm,double **&omega,double **&torque, double **&fcm, double **&ex_space, double **&ey_space, double **&ez_space){

    double a[6];
    double b[6][6];
    double c[6];
    //double e[6];

    a[1] = 0.2;
    a[2] = 0.3;
    a[3] = 0.6;
    a[4] = 1.0;
    a[5] = 0.875;

    b[1][0] = 0.2;
    b[2][0] = 0.075;           b[2][1] = 0.225;
    b[3][0] = 0.3;             b[3][1] = -0.9;        b[3][2] = 1.2;
    b[4][0] = -11.0/54.0;      b[4][1] = 2.5;         b[4][2] = -70.0/27.0;    b[4][3] = 35.0/27.0;
    b[5][0] = 1631.0/55296.0; b[5][1] = 175.0/512.0; b[5][2] = 575.0/13824.0; b[5][3] = 44275.0/110592.0; b[5][4] = 253.0/4096.0;

    c[0] = 37.0/378.0;
    c[1] = 0.0;
    c[2] = 250.0/621.0;
    c[3] = 125.0/594.0;
    c[4] = 0.0;
    c[5] = 512.0/1771.0;

    int numsys = currentIndex;
    int numbodies;
    double time = 0.0;
    int * mappings;
    double SysKE =0.0;


    for (int i = 0; i <= numsys; i++){
      mappings = system[i].system->GetMappings();

      numbodies = system[i].system->GetNumBodies() - 1;
      Matrix FF(6,numbodies);
      for (int j=1; j<=numbodies; j++){
        FF(1,j)  = ConFac*torque[mappings[j - 1]-1][0];
        FF(2,j)  = ConFac*torque[mappings[j - 1]-1][1];
        FF(3,j)  = ConFac*torque[mappings[j - 1]-1][2];

        FF(4,j)  = ConFac*fcm[mappings[j - 1]-1][0];
        FF(5,j)  = ConFac*fcm[mappings[j - 1]-1][1];
        FF(6,j)  = ConFac*fcm[mappings[j - 1]-1][2];
      }


      //------------------------------------//
    // Get a solver and solve that system.
      Solver * theSolver = Solver::GetSolver((SolverType)system[i].solver);
      theSolver->SetSystem(system[i].system);
      theSolver->Solve(time, FF);

      ColMatrix initial_x;
      ColMatrix initial_xdot;
      ColMatMap* x;
      ColMatMap* xdot;
      ColMatMap* xdotdot;

      x = theSolver->GetState();
      xdot = theSolver->GetStateDerivative();
      xdotdot=theSolver->GetStateDerivativeDerivative();

      initial_x = *x;
      initial_xdot = *xdot;
      ColMatrix f[6];
      ColMatrix ff[6];

      ff[0] = initial_xdot;
      f[0] = *xdotdot;

      for(int ii=1;ii<6;ii++){
        time = a[ii] * Tfull;
        (*x) = initial_x;
        (*xdot) = initial_xdot;
        for(int j=0;j<ii;j++){
          (*x) = (*x) + (b[ii][j]*Tfull)*ff[j];
          (*xdot) = (*xdot) + (b[ii][j]*Tfull)*f[j];
        }
        theSolver->Solve(time,FF);
        f[ii] = (*xdotdot);
        ff[ii] = (*xdot);
      }


      (*x) = initial_x + (c[0]*Tfull)*ff[0] + (c[2]*Tfull)*ff[2] + (c[3]*Tfull)*ff[3] + (c[5]*Tfull)*ff[5];

      (*xdot) = initial_xdot + (c[0]*Tfull)*f[0] + (c[2]*Tfull)*f[2] + (c[3]*Tfull)*f[3] + (c[5]*Tfull)*f[5];


      int numjoints = system[i].system->joints.GetNumElements();
      for(int k = 0; k < numjoints; k++){
        system[i].system->joints(k)->ForwardKinematics();
      }


      for(int k = 0; k < numbodies; k++){
        Vect3 temp1 =system[i].system->bodies(k+1)->r;
        Vect3 temp2 =system[i].system->bodies(k+1)->v;
        Vect3 temp3 =system[i].system->bodies(k+1)->omega;
        Mat3x3 temp4 =system[i].system->bodies(k+1)->n_C_k;
        SysKE = SysKE + system[i].system->bodies(k+1)->KE;
        for(int m = 0; m < 3; m++){
          xcm[mappings[k]-1][m]   = temp1(m+1);
          vcm[mappings[k]-1][m]   = temp2(m+1);
          omega[mappings[k]-1][m] = temp3(m+1);

          ex_space[mappings[k]-1][m] = temp4(m+1,1);
          ey_space[mappings[k]-1][m] = temp4(m+1,2);
          ez_space[mappings[k]-1][m] = temp4(m+1,3);
        }
      }
      delete theSolver;
    }
  }


void Workspace::WriteFile(char * filename){
  int numsys = currentIndex;
  int numbodies;


  for (int i = 0; i <= numsys; i++){
    numbodies = system[i].system->GetNumBodies() - 1;
    ofstream outfile;
    outfile.open(filename,ofstream::out | ios::app);
    outfile << numbodies<<endl;
    outfile << "Atoms "<<endl;
    for(int k = 0; k < numbodies; k++){
      Vect3 temp1 =system[i].system->bodies(k+1)->r;
      outfile<<1<<"\t"<<temp1(1)<<"\t"<<temp1(2)<<"\t"<<temp1(3)<<endl;
    }
    outfile.close();
  }
  }
