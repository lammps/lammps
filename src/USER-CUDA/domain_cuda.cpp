/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author (triclinic) : Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "stdio.h"
#include "math.h"
#include "domain_cuda.h"
#include "style_region.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_deform.h"
#include "region.h"
#include "lattice.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include "cuda.h"
#include "domain_cu.h"

using namespace LAMMPS_NS;

#define BIG   1.0e20
#define SMALL 1.0e-4
#define DELTA 1

enum{NO_REMAP,X_REMAP,V_REMAP};                   // same as fix_deform.cpp

/* ----------------------------------------------------------------------
   default is periodic 
------------------------------------------------------------------------- */

DomainCuda::DomainCuda(LAMMPS *lmp) : Domain(lmp)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");
}

/* ---------------------------------------------------------------------- */

void DomainCuda::init()
{
  cuda->accelerator(0,NULL);
  Domain::init();

  if(not cuda->finished_run)
  {
    cuda->setDomainParams();
    Cuda_Domain_Init(&cuda->shared_data);
  }
}

/* ----------------------------------------------------------------------
   set global box params
   assumes boxlo/hi and triclinic tilts are already set
------------------------------------------------------------------------- */

void DomainCuda::set_global_box()
{
  Domain::set_global_box();

  if(not cuda->finished_run)
  {
    cuda->setDomainParams();
  }
}

/* ----------------------------------------------------------------------
   set lamda box params, only need be done one time
   assumes global box is defined and proc assignment has been made by comm
   for uppermost proc, insure subhi = 1.0 (in case round-off occurs)
------------------------------------------------------------------------- */

void DomainCuda::set_lamda_box()
{
  Domain::set_lamda_box();

  if(not cuda->finished_run)
  {
    cuda->setDomainParams();
  }
}

/* ----------------------------------------------------------------------
   set local subbox params
   assumes global box is defined and proc assignment has been made
   for uppermost proc, insure subhi = boxhi (in case round-off occurs)
------------------------------------------------------------------------- */

void DomainCuda::set_local_box()
{
  Domain::set_local_box();

  if(not cuda->finished_run)
  {
   // cuda->setDomainParams();
    //Cuda_Domain_Init(&cuda->shared_data);
  }
}

/* ----------------------------------------------------------------------
   reset global & local boxes due to global box boundary changes
   if shrink-wrapped, determine atom extent and reset boxlo/hi
   if shrink-wrapped and triclinic, perform shrink-wrap in box coords
------------------------------------------------------------------------- */

void DomainCuda::reset_box()
{
  if (nonperiodic == 2) {

    // convert back to box coords for shrink-wrap operation

    if (triclinic) lamda2x(atom->nlocal);

    // compute extent of atoms on this proc

    double extent[3][2],all[3][2];

    extent[2][0] = extent[1][0] = extent[0][0] = BIG;
    extent[2][1] = extent[1][1] = extent[0][1] = -BIG;

    double **x = atom->x;
    int nlocal = atom->nlocal;

    if (cuda->finished_setup&&(!cuda->oncpu))
      {
        extent[0][0]=cuda->extent[0];
        extent[0][1]=cuda->extent[1];
        extent[1][0]=cuda->extent[2];
        extent[1][1]=cuda->extent[3];
        extent[2][0]=cuda->extent[4];
        extent[2][1]=cuda->extent[5];
      }
    else
      for (int i = 0; i < nlocal; i++) {
        extent[0][0] = MIN(extent[0][0],x[i][0]);
        extent[0][1] = MAX(extent[0][1],x[i][0]);
        extent[1][0] = MIN(extent[1][0],x[i][1]);
        extent[1][1] = MAX(extent[1][1],x[i][1]);
        extent[2][0] = MIN(extent[2][0],x[i][2]);
        extent[2][1] = MAX(extent[2][1],x[i][2]);
      }

    // compute extent across all procs
    // flip sign of MIN to do it in one Allreduce MAX

    extent[0][0] = -extent[0][0];
    extent[1][0] = -extent[1][0];
    extent[2][0] = -extent[2][0];

    MPI_Allreduce(extent,all,6,MPI_DOUBLE,MPI_MAX,world);

    // in shrink-wrapped dims, set box by atom extent
    // if minimum set, enforce min box size settings

    if (xperiodic == 0) {
      if (boundary[0][0] == 2) boxlo[0] = -all[0][0] - SMALL;
      else if (boundary[0][0] == 3) boxlo[0] = MIN(-all[0][0]-SMALL,minxlo);
      if (boundary[0][1] == 2) boxhi[0] = all[0][1] + SMALL;
      else if (boundary[0][1] == 3) boxhi[0] = MAX(all[0][1]+SMALL,minxhi);
      if (boxlo[0] > boxhi[0]) error->all(FLERR,"Illegal simulation box");
    }
    if (yperiodic == 0) {
      if (boundary[1][0] == 2) boxlo[1] = -all[1][0] - SMALL;
      else if (boundary[1][0] == 3) boxlo[1] = MIN(-all[1][0]-SMALL,minylo);
      if (boundary[1][1] == 2) boxhi[1] = all[1][1] + SMALL;
      else if (boundary[1][1] == 3) boxhi[1] = MAX(all[1][1]+SMALL,minyhi);
      if (boxlo[1] > boxhi[1]) error->all(FLERR,"Illegal simulation box");
    }
    if (zperiodic == 0) {
      if (boundary[2][0] == 2) boxlo[2] = -all[2][0] - SMALL;
      else if (boundary[2][0] == 3) boxlo[2] = MIN(-all[2][0]-SMALL,minzlo);
      if (boundary[2][1] == 2) boxhi[2] = all[2][1] + SMALL;
      else if (boundary[2][1] == 3) boxhi[2] = MAX(all[2][1]+SMALL,minzhi);
      if (boxlo[2] > boxhi[2]) error->all(FLERR,"Illegal simulation box");
    }
  }

  set_global_box();
  set_local_box();

  if(not cuda->finished_run)
  {
    cuda->setDomainParams();
    Cuda_Domain_Init(&cuda->shared_data);
  }

  // if shrink-wrapped, convert to lamda coords for new box
  // must re-invoke pbc() b/c x2lamda result can be outside 0,1 due to roundoff

  if (nonperiodic == 2 && triclinic) {
    x2lamda(atom->nlocal);
    pbc();
  }
}

/* ----------------------------------------------------------------------
   enforce PBC and modify box image flags for each atom
   called every reneighboring and by other commands that change atoms
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   if fix deform, remap velocity of fix group atoms by box edge velocities
   for triclinic, atoms must be in lamda coords (0-1) before pbc is called
   image = 10 bits for each dimension
   increment/decrement in wrap-around fashion
------------------------------------------------------------------------- */

void DomainCuda::pbc()
{
  if(cuda->finished_setup&&(!cuda->oncpu))
  {
  	cuda->setDomainParams();
    Cuda_Domain_PBC(&cuda->shared_data, deform_vremap, deform_groupbit,cuda->extent);
    return;
  }
  
  Domain::pbc();
}


/* ----------------------------------------------------------------------
   convert triclinic 0-1 lamda coords to box coords for all N atoms
   x = H lamda + x0;
------------------------------------------------------------------------- */

void DomainCuda::lamda2x(int n)
{
  if(cuda->finished_setup&&(!cuda->oncpu))
  {
    Cuda_Domain_lamda2x(&cuda->shared_data,n);
    return;
  }

  Domain::lamda2x(n);
}

/* ----------------------------------------------------------------------
   convert box coords to triclinic 0-1 lamda coords for all N atoms
   lamda = H^-1 (x - x0)
------------------------------------------------------------------------- */

void DomainCuda::x2lamda(int n)
{
  if(cuda->finished_setup&&(!cuda->oncpu))
  {
    Cuda_Domain_x2lamda(&cuda->shared_data,n);
    return;
  }
 
  Domain::x2lamda(n);
}
