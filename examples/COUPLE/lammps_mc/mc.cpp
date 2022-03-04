/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
------------------------------------------------------------------------- */

// MC code used with LAMMPS in client/server mode
// MC is the client, LAMMPS is the server

// Syntax: mc infile mode modearg
//         mode = file, zmq
//         modearg = filename for file, localhost:5555 for zmq

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "mc.h"
#include "random_park.h"

#include "cslib.h"
using namespace CSLIB_NS;

void error(const char *);
CSlib *cs_create(char *, char *);

#define MAXLINE 256

/* ---------------------------------------------------------------------- */

// main program

int main(int narg, char **arg)
{
  if (narg != 4) {
    error("Syntax: mc infile mode modearg");
    exit(1);
  }

  // initialize CSlib

  CSlib *cs = cs_create(arg[2],arg[3]);

  // create MC class and perform run

  MC *mc = new MC(arg[1],cs);
  mc->run();

  // final MC stats

  int naccept = mc->naccept;
  int nattempt = mc->nattempt;

  printf("------ MC stats ------\n");
  printf("MC attempts = %d\n",nattempt);
  printf("MC accepts = %d\n",naccept);
  printf("Acceptance ratio = %g\n",1.0*naccept/nattempt);

  // clean up

  delete cs;
  delete mc;
}

/* ---------------------------------------------------------------------- */

void error(const char *str)
{
  printf("ERROR: %s\n",str);
  exit(1);
}

/* ---------------------------------------------------------------------- */

CSlib *cs_create(char *mode, char *arg)
{
  CSlib *cs = new CSlib(0,mode,arg,NULL);

  // initial handshake to agree on protocol

  cs->send(0,1);
  cs->pack_string(1,(char *) "mc");

  int msgID,nfield;
  int *fieldID,*fieldtype,*fieldlen;
  msgID = cs->recv(nfield,fieldID,fieldtype,fieldlen);

  return cs;
}

// ----------------------------------------------------------------------
// MC class
// ----------------------------------------------------------------------

MC::MC(char *mcfile, void *cs_caller)
//MC::MC(char *mcfile, CSlib *cs_caller)
{
  cs_void = cs_caller;

  // setup MC params

  options(mcfile);

  // random # generator

  random = new RanPark(seed);
}

/* ---------------------------------------------------------------------- */

MC::~MC()
{
  free(x);
  delete random;
}

/* ---------------------------------------------------------------------- */

void MC::run()
{
  int iatom,accept,msgID,nfield;
  double pe_initial,pe_final,edelta;
  double dx,dy,dz;
  double xold[3],xnew[3];
  int *fieldID,*fieldtype,*fieldlen;

  enum{NATOMS=1,EINIT,DISPLACE,ACCEPT,RUN};

  CSlib *cs = (CSlib *) cs_void;

  // one-time request for atom count from MD
  // allocate 1d coord buffer

  cs->send(NATOMS,0);

  msgID = cs->recv(nfield,fieldID,fieldtype,fieldlen);
  natoms = cs->unpack_int(1);

  x = (double *) malloc(3*natoms*sizeof(double));

  // loop over MC moves

  naccept = nattempt = 0;

  for (int iloop = 0; iloop < nloop; iloop++) {

    // request current energy from MD
    // recv energy, coords from MD

    cs->send(EINIT,0);

    msgID = cs->recv(nfield,fieldID,fieldtype,fieldlen);
    pe_initial = cs->unpack_double(1);
    double *x = (double *) cs->unpack(2);

    // perform simple MC event
    // displace a single atom by random amount

    iatom = (int) natoms*random->uniform();
    xold[0] = x[3*iatom+0];
    xold[1] = x[3*iatom+1];
    xold[2] = x[3*iatom+2];

    dx = 2.0*delta*random->uniform() - delta;
    dy = 2.0*delta*random->uniform() - delta;
    dz = 2.0*delta*random->uniform() - delta;

    xnew[0] = xold[0] + dx;
    xnew[1] = xold[1] + dx;
    xnew[2] = xold[2] + dx;

    // send atom ID and its new coords to MD
    // recv new energy

    cs->send(DISPLACE,2);
    cs->pack_int(1,iatom+1);
    cs->pack(2,4,3,xnew);

    msgID = cs->recv(nfield,fieldID,fieldtype,fieldlen);
    pe_final = cs->unpack_double(1);

    // decide whether to accept/reject MC event

    if (pe_final <= pe_initial) accept = 1;
    else if (temperature == 0.0) accept = 0;
    else if (random->uniform() > 
             exp(natoms*(pe_initial-pe_final)/temperature)) accept = 0;
    else accept = 1;
  
    nattempt++;
    if (accept) naccept++;

    // send accept (1) or reject (0) flag to MD

    cs->send(ACCEPT,1);
    cs->pack_int(1,accept);

    msgID = cs->recv(nfield,fieldID,fieldtype,fieldlen);
    
    // send dynamics timesteps

    cs->send(RUN,1);
    cs->pack_int(1,ndynamics);

    msgID = cs->recv(nfield,fieldID,fieldtype,fieldlen);
  }

  // send exit message to MD

  cs->send(-1,0);
  msgID = cs->recv(nfield,fieldID,fieldtype,fieldlen);
}

/* ---------------------------------------------------------------------- */

void MC::options(char *filename)
{
  // default params

  nsteps = 0;
  ndynamics = 100;
  delta = 0.1;
  temperature = 1.0;
  seed = 12345;

  // read and parse file

  FILE *fp = fopen(filename,"r");
  if (fp == NULL) error("Could not open MC file");

  char line[MAXLINE];
  char *keyword,*value;
  char *eof = fgets(line,MAXLINE,fp);

  while (eof) {
    if (line[0] == '#') {                     // comment line
      eof = fgets(line,MAXLINE,fp);
      continue;
    }

    value = strtok(line," \t\n\r\f");
    if (value == NULL) {                      // blank line
      eof = fgets(line,MAXLINE,fp);
      continue;
    }

    keyword = strtok(NULL," \t\n\r\f");
    if (keyword == NULL) error("Missing keyword in MC file");

    if (strcmp(keyword,"nsteps") == 0) nsteps = atoi(value);
    else if (strcmp(keyword,"ndynamics") == 0) ndynamics = atoi(value);
    else if (strcmp(keyword,"delta") == 0) delta = atof(value);
    else if (strcmp(keyword,"temperature") == 0) temperature = atof(value);
    else if (strcmp(keyword,"seed") == 0) seed = atoi(value);
    else error("Unknown param in MC file");

    eof = fgets(line,MAXLINE,fp);
  }

  // derived params

  nloop = nsteps/ndynamics;
}
