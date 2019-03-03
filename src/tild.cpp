<<<<<<< HEAD
/* ------------------------------------------------------------------
    TILD - Theoretically Informed Langevan Dynamics
    This replicates the TILD coe done by the Riggleman group, 
    previously known as Dynamical Mean Field Theory. 
    
    Copyright (2019) Christian Tabedzki and Zachariah Vicars.
    tabedzki@seas.upenn.edu zvicars@seas.upenn.edu
-------------------------------------------------------------------- */

#include <mpi.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "tild.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include <iostream>

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.00001

/* ------------------------------------------------------------ */

TILD::TILD(LAMMPS *lmp) : KSpace(lmp)  {
    if (screen) fprintf(screen,"TILD construction...\n");
    if (logfile) fprintf(logfile,"TILD construction...\n");
  };

TILD::~TILD(){
  return;
}

void TILD::settings(int narg, char **arg)
{
  std::cout<<"help me "<< std::endl;
}

void TILD::init()
{
  if (comm->me == 0) {
    if (screen) fprintf(screen,"TILD initialization ...\n");
    if (logfile) fprintf(logfile,"TILD initialization ...\n");
  }
}


void TILD::setup(){
  return;
}

void TILD::compute(int i1, int i2){
  return;
}

void TILD::compute_group_group(int, int, int){
  return;
}

double TILD::memory_usage(){
  return 0;
}

void TILD::eik_dot_r(){
  return;
}

void TILD::allocate(){
  return ;
}
=======
/* ------------------------------------------------------------------
    TILD - Theoretically Informed Langevan Dynamics
    This replicates the TILD coe done by the Riggleman group, 
    previously known as Dynamical Mean Field Theory. 
    
    Copyright (2019) Christian Tabedzki and Zachariah Vicars.
    tabedzki@seas.upenn.edu zvicars@seas.upenn.edu
-------------------------------------------------------------------- */

#include <mpi.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "tild.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include <iostream>

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.00001

/* ------------------------------------------------------------ */

TILD::TILD(LAMMPS *lmp) : KSpace(lmp)  {
    if (screen) fprintf(screen,"TILD construction...\n");
    if (logfile) fprintf(logfile,"TILD construction...\n");
  };

TILD::~TILD(){
  return;
}

void TILD::settings(int narg, char **arg)
{
  std::cout<<"help me "<< std::endl;
}

void TILD::init()
{
  if (comm->me == 0) {
    if (screen) fprintf(screen,"TILD initialization ...\n");
    if (logfile) fprintf(logfile,"TILD initialization ...\n");
  }
}


void TILD::setup(){
  return;
}

void TILD::compute(int i1, int i2){
  return;
}

void TILD::compute_group_group(int, int, int){
  return;
}

double TILD::memory_usage(){
  return 0;
}

void TILD::eik_dot_r(){
  return;
}

void TILD::allocate(){
  return ;
}
>>>>>>> 4670323e5724f284949fe7a1644765c3cdf10afc
