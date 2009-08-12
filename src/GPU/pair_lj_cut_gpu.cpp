/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   cetain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "pair_lj_cut_gpu.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "integrate.h"
#include "memory.h"
#include "error.h"
#include "neigh_request.h"
#include "universe.h"

#include <string>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

// External functions from cuda library for force decomposition

int * lj_gpu_init(int &ij_size, const int ntypes, double **cutsq, 
                  double **sigma, double **epsilon, double **host_lj1, 
                  double **host_lj2, double **host_lj3, double **host_lj4, 
                  double **offset, double *special_lj, const int max_nbors, 
                  const int gpu_id);
void lj_gpu_clear();
bool lj_gpu_reset_nbors(const int nall, const int inum, int *ilist, 
                        const int *numj);
void lj_gpu_nbors(const int num_ij);
void lj_gpu_atom(double **host_x, const int *host_type, const bool rebuild);
void lj_gpu(const bool eflag, const bool vflag, const bool rebuild);
double lj_gpu_forces(double **f, const int *ilist, const bool eflag, 
                     const bool vflag, const bool eflag_atom, 
                     const bool vflag_atom, double *eatom, double **vatom,
                     double *virial);
std::string lj_gpu_name(const int gpu_id, const int max_nbors);
void lj_gpu_time();
int lj_gpu_num_devices();
double lj_gpu_bytes();

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCutGPU::PairLJCutGPU(LAMMPS *lmp) : PairLJCut(lmp), multi_gpu_mode(0)
{
  ij_new=NULL;
  respa_enable = 0;
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairLJCutGPU::~PairLJCutGPU()
{
  if (comm->me == 0 && screen) {
    printf("\n\n-------------------------------------");
    printf("--------------------------------\n");
    printf("      GPU Time Stamps: ");
    printf("\n-------------------------------------");
    printf("--------------------------------\n");
    lj_gpu_time();
    printf("Procs: %d\n",comm->nprocs);
    printf("-------------------------------------");
    printf("--------------------------------\n\n");
  }
  lj_gpu_clear();
  if (ij_new!=NULL) {
    ij_new=NULL;
    delete [] ij_new;
  }
}

/* ---------------------------------------------------------------------- */

void PairLJCutGPU::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;
  if (vflag_atom) 
    error->all("Per-atom virial not available with GPU Gay-Berne.");

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int inum = list->inum;
  int *ilist = list->ilist;
  
  bool rebuild=false;
  if (neighbor->ncalls > last_neighbor) {
    last_neighbor=neighbor->ncalls;
    rebuild=true;
  }
  
  // copy nbors to GPU
  if (rebuild)
    if (!lj_gpu_reset_nbors(nall, inum, ilist, list->numneigh))
      error->one("Total # of atoms exceed maximum allowed per GPGPU.\n");
  
  // copy atom data to GPU
  lj_gpu_atom(atom->x,atom->type,rebuild);

  int i,j,ii,jj,jnum;
  double factor_lj;
  int *jlist;

  if (rebuild==true) {
    int num_ij = 0;

    // loop over neighbors of my atoms
    int *ijp=ij_new;
    for (ii = 0; ii<inum; ii++) {
      i = ilist[ii];
      jlist = list->firstneigh[i];
      jnum = list->numneigh[i];
      
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];

        *ijp=j;
        ijp++;
        num_ij++;
          
        if (num_ij==ij_size) {
          memcpy(ij,ij_new,num_ij*sizeof(int));
          lj_gpu_nbors(num_ij);
          ijp=ij_new;
          num_ij=0;
        }
      }
    }
    if (num_ij>0) {
      memcpy(ij,ij_new,num_ij*sizeof(int));
      lj_gpu_nbors(num_ij);
    }
  }
  
  lj_gpu(eflag,vflag,rebuild);
  eng_vdwl=lj_gpu_forces(atom->f,ilist,eflag,vflag, eflag_atom, vflag_atom, 
                         eatom, vatom, virial);

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
  
  int irequest = neighbor->request(this);
  cut_respa=NULL;

  // Repeat cutsq calculation because done after call to init_style
  double cut;
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      cut = init_one(i,j);
      cutsq[i][j] = cutsq[j][i] = cut*cut;
    }

  ij=lj_gpu_init(ij_size, atom->ntypes+1, cutsq, sigma, epsilon, lj1, lj2, lj3, 
                 lj4, offset, force->special_lj, neighbor->oneatom, my_gpu);
  if (ij==0)
    error->one("AT LEAST ONE PROCESS COULD NOT ALLOCATE A CUDA-ENABLED GPU.");
    
  if (ij_new!=NULL)
    delete [] ij_new;
  ij_new=new int[ij_size];
  
  last_neighbor = -1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  if (force->newton_pair) 
    error->all("Cannot use newton with GPU LJCut pair style.");

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
      std::string gpu_string=lj_gpu_name(gpui,neighbor->oneatom);
      printf("GPU %d: %s\n",gpui,gpu_string.c_str());  
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
