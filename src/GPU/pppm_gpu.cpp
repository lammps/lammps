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
   Contributing authors: Roy Pollock (LLNL), Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "pppm_gpu.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "domain.h"
#include "fft3d_wrap.h"
#include "remap_wrap.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXORDER 7
#define OFFSET 16384
#define SMALL 0.00001
#define LARGE 10000.0
#define EPS_HOC 1.0e-7

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PPPM_GPU::PPPM_GPU(LAMMPS *lmp, int narg, char **arg) : PPPM(lmp, narg, arg)
{
}

/* ----------------------------------------------------------------------
   free all memory 
------------------------------------------------------------------------- */

PPPM_GPU::~PPPM_GPU()
{
}

/* ----------------------------------------------------------------------
   called once before run 
------------------------------------------------------------------------- */

void PPPM_GPU::init()
{
  PPPM::init();
}

/* ----------------------------------------------------------------------
   adjust PPPM_GPU coeffs, called initially and whenever volume has changed 
------------------------------------------------------------------------- */

void PPPM_GPU::setup()
{
  PPPM::setup();
}

/* ----------------------------------------------------------------------
   compute the PPPM_GPU long-range force, energy, virial 
------------------------------------------------------------------------- */

void PPPM_GPU::compute(int eflag, int vflag)
{
  int i;

  // convert atoms from box to lamda coords
  
  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  // extend size of per-atom arrays if necessary

  if (atom->nlocal > nmax) {
    memory->destroy_2d_int_array(part2grid);
    nmax = atom->nmax;
    part2grid = memory->create_2d_int_array(nmax,3,"pppm:part2grid");
  }

  energy = 0.0;
  if (vflag) for (i = 0; i < 6; i++) virial[i] = 0.0;

  // find grid points for all my particles
  // map my particle charge onto my local 3d density grid

  particle_map();
  make_rho();

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  brick2fft();

  // compute potential gradient on my FFT grid and
  //   portion of e_long on this proc's FFT grid
  // return gradients (electric fields) in 3d brick decomposition
  
  poisson(eflag,vflag);

  // all procs communicate E-field values to fill ghost cells
  //   surrounding their 3d bricks

  fillbrick();

  // calculate the force on my particles

  fieldforce();

  // sum energy across procs and add in volume-dependent term

  if (eflag) {
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;
   
    energy *= 0.5*volume;
    energy -= g_ewald*qsqsum/1.772453851 +
      0.5*PI*qsum*qsum / (g_ewald*g_ewald*volume);
    energy *= qqrd2e*scale;
  }

  // sum virial across procs

  if (vflag) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = 0.5*qqrd2e*scale*volume*virial_all[i];
  }

  // 2d slab correction

  if (slabflag) slabcorr(eflag);

  // convert atoms back from lamda to box coords
  
  if (triclinic) domain->lamda2x(atom->nlocal);
}

