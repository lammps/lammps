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
   Contributing authors: Mike Brown (SNL), wmbrown@sandia.gov
                         Peng Wang (Nvidia), penwang@nvidia.com
                         Paul Crozier (SNL), pscrozi@sandia.gov
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "pair_lj_cut_gpu.h"
#include "math_extra.h"
#include "atom.h"
#include "domain.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "integrate.h"
#include "memory.h"
#include "error.h"
#include "universe.h"

#include <string>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#ifndef WINDLL

// External functions from cuda library for force decomposition

bool lj_gpu_init(int &ij_size, const int ntypes, double **cutsq,double **sigma, 
                  double **epsilon, double **host_lj1, double **host_lj2, 
                  double **host_lj3, double **host_lj4, double **offset, 
                  double *special_lj, double *boxlo, double *boxhi, double cell_len, double skin,
                  const int max_nbors, const int gpu_id);
void lj_gpu_clear();
double lj_gpu_cell(double **force, double *virial, double **host_x, int *host_type, const int inum, const int nall, 
		   const int ago, const bool eflag, const bool vflag, 
		   const double *boxlo, const double *boxhi);
void lj_gpu_name(const int gpu_id, const int max_nbors, char * name);
void lj_gpu_time();
double lj_gpu_bytes();

#else
#include <windows.h>

typedef bool (*_lj_gpu_init)(int &ij_size, const int ntypes, double **cutsq,double **sigma, 
                  double **epsilon, double **host_lj1, double **host_lj2, 
                  double **host_lj3, double **host_lj4, double **offset, 
                  double *special_lj, double *boxlo, double *boxhi, double cell_len, double skin,
                  const int max_nbors, const int gpu_id);
typedef void (*_lj_gpu_clear)();
typedef double (*_lj_gpu_cell)(double **force, double *virial, double **host_x, int *host_type, const int inum, const int nall, 
		   const int ago, const bool eflag, const bool vflag, 
		   const double *boxlo, const double *boxhi);
typedef void (*_lj_gpu_name)(const int gpu_id, const int max_nbors, char * name);
typedef void (*_lj_gpu_time)();
typedef double (*_lj_gpu_bytes)();

_lj_gpu_init lj_gpu_init;
_lj_gpu_clear lj_gpu_clear;
_lj_gpu_cell lj_gpu_cell;
_lj_gpu_name lj_gpu_name;
_lj_gpu_time lj_gpu_time;
_lj_gpu_bytes lj_gpu_bytes;

#endif

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCutGPU::PairLJCutGPU(LAMMPS *lmp) : PairLJCut(lmp), multi_gpu_mode(0)
{
  ij_new=NULL;
  respa_enable = 0;

#ifdef WINDLL
  HINSTANCE hinstLib = LoadLibrary(TEXT("gpu.dll"));
  if (hinstLib == NULL) {
    printf("\nUnable to load gpu.dll\n");
    exit(1);
  }

  lj_gpu_init=(_lj_gpu_init)GetProcAddress(hinstLib,"lj_gpu_init");
  lj_gpu_clear=(_lj_gpu_clear)GetProcAddress(hinstLib,"lj_gpu_clear");
  lj_gpu_cell=(_lj_gpu_cell)GetProcAddress(hinstLib,"lj_gpu_cell");
  lj_gpu_name=(_lj_gpu_name)GetProcAddress(hinstLib,"lj_gpu_name");
  lj_gpu_time=(_lj_gpu_time)GetProcAddress(hinstLib,"lj_gpu_time");
  lj_gpu_bytes=(_lj_gpu_bytes)GetProcAddress(hinstLib,"lj_gpu_bytes");
#endif

}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairLJCutGPU::~PairLJCutGPU()
{
  printf("\n\n-------------------------------------");
  printf("--------------------------------\n");
  printf("      GPU Time Stamps: ");
  printf("\n-------------------------------------");
  printf("--------------------------------\n");
  lj_gpu_time();
  printf("-------------------------------------");
  printf("--------------------------------\n\n");
  lj_gpu_clear();
  if (ij_new!=NULL)
    delete [] ij_new;
}

/* ---------------------------------------------------------------------- */

void PairLJCutGPU::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  // compute forces on GPU
  eng_vdwl = lj_gpu_cell(atom->f, virial, atom->x, atom->type, atom->nlocal, atom->nlocal + atom->nghost, 
			 neighbor->ago, eflag, vflag, domain->boxlo, domain->boxhi);

  if (vflag_fdotr) virial_compute();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJCutGPU::settings(int narg, char **arg)
{
  // strip off GPU keyword/value and send remaining args to parent

  if (narg < 2) error->all("Illegal pair_style command");

  // set multi_gpu_mode to one/node for multiple gpus on 1 node
  // -- param is starting gpu id
  // set multi_gpu_mode to one/gpu to select the same gpu id on every node
  // -- param is id of gpu
  // set multi_gpu_mode to multi/gpu to get ma
  // -- param is number of gpus per node

  if (strcmp(arg[0],"one/node") == 0)
    multi_gpu_mode = ONE_NODE;
  else if (strcmp(arg[0],"one/gpu") == 0)
    multi_gpu_mode = ONE_GPU;
  else if (strcmp(arg[0],"multi/gpu") == 0)
    multi_gpu_mode = MULTI_GPU;
  else error->all("Illegal pair_style command");

  multi_gpu_param = atoi(arg[1]);

  if (multi_gpu_mode == MULTI_GPU && multi_gpu_param < 1)
    error->all("Illegal pair_style command");

  PairLJCut::settings(narg-2,&arg[2]);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutGPU::init_style()
{
  if (force->pair_match("gpu",0) == NULL)
    error->all("Cannot use pair hybrid with multiple GPU pair styles");

  // set the GPU ID

  int my_gpu=comm->me+multi_gpu_param;
  int ngpus=universe->nprocs;
  if (multi_gpu_mode==ONE_GPU) {
    my_gpu=multi_gpu_param;
    ngpus=1;
  } else if (multi_gpu_mode==MULTI_GPU) {
    ngpus=multi_gpu_param;
    my_gpu=comm->me%ngpus;
    if (ngpus>universe->nprocs)
      ngpus=universe->nprocs;
  }

  cut_respa=NULL;

  // Repeat cutsq calculation because done after call to init_style
  double cut;
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      cut = init_one(i,j);
      cutsq[i][j] = cutsq[j][i] = cut*cut;
    }

  // use the max cutoff length as the cell length
  double maxcut = -1.0;
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      if (cutsq[i][j] > maxcut) maxcut = cutsq[i][j];

  // for this problem, adding skin results in better perf
  // this may be a parameter in the future
  double cell_len = sqrt(maxcut) + neighbor->skin;

  if (!lj_gpu_init(ij_size, atom->ntypes+1, cutsq, sigma, epsilon, lj1, lj2, lj3, 
                 lj4, offset, force->special_lj, domain->boxlo, domain->boxhi, 
                 cell_len, neighbor->skin, neighbor->oneatom, my_gpu))
    error->one("At least one process could not allocate a CUDA-enabled gpu");

  if (ij_new!=NULL)
    delete [] ij_new;
  ij_new=new int[ij_size];
  
  if (force->newton_pair) 
    error->all("Cannot use newton pair with GPU lj/cut pair style");

  if (comm->me == 0 && screen) {
    printf("\n-------------------------------------");
    printf("-------------------------------------\n");
    printf("- Using GPGPU acceleration for LJ-Cut:\n");
    printf("-------------------------------------");
    printf("-------------------------------------\n");

    for (int i=0; i<ngpus; i++) {
      int gpui=my_gpu;
      if (multi_gpu_mode==ONE_NODE)
        gpui=i+multi_gpu_param;
      else if (multi_gpu_mode==MULTI_GPU)
        gpui=i;
      char gpu_string[500];
      lj_gpu_name(gpui,neighbor->oneatom,gpu_string);
      printf("GPU %d: %s\n",gpui,gpu_string);   
    }
    printf("-------------------------------------");
    printf("-------------------------------------\n\n");
  }
}

/* ---------------------------------------------------------------------- */

double PairLJCutGPU::memory_usage()
{
  double bytes=Pair::memory_usage()+ij_size*sizeof(int);
  return bytes+lj_gpu_bytes();
}
