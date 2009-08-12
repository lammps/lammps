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
#include "pair_gayberne_gpu.h"
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

#ifdef GB_GPU_OMP
#include "omp.h"
#endif

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

// External functions from cuda library for atom decomposition

int * gb_gpu_init(int &ij_size, const int ntypes, const double gamma,
                  const double upsilon, const double mu, double **shape,
                  double **well, double **cutsq, double **sigma, 
                  double **epsilon, double *host_lshape, int **form,
                  double **host_lj1, double **host_lj2, double **host_lj3, 
                  double **host_lj4, double **offset, double *special_lj, 
                  const int max_nbors, const int thread, const int gpu_id);
void gb_gpu_clear(const int thread);
int * gb_gpu_reset_nbors(const int nall, const int nlocal, const int inum, 
                         int *ilist, const int *numj, const int *type,
                         const int thread, bool &success);
void gb_gpu_nbors(const int num_ij, const bool eflag, const int thread);
void gb_gpu_atom(double **host_x, double **host_quat, const int *host_type, 
                 const bool rebuild, const int thread);
void gb_gpu_gayberne(const bool eflag, const bool vflag, const bool rebuild, 
                     const int thread);
double gb_gpu_forces(double **f, double **tor, const int *ilist,
                     const bool eflag, const bool vflag, const bool eflag_atom,
                     const bool vflag_atom, double *eatom, double **vatom,
                     double *virial, const int thread);
std::string gb_gpu_name(const int i, const int max_nbors);
void gb_gpu_time(const int thread);
int gb_gpu_num_devices();
double gb_gpu_bytes();

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGayBerneGPU::PairGayBerneGPU(LAMMPS *lmp) : PairGayBerne(lmp), my_thread(0),
                                                omp_chunk(0), nthreads(1),
                                                multi_gpu_mode(ONE_NODE),
                                                multi_gpu_param(0)
{
  ij_new[0]=NULL;
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairGayBerneGPU::~PairGayBerneGPU()
{
  if (comm->me == 0 && screen) {
    printf("\n\n-------------------------------------");
    printf("--------------------------------\n");
    printf("      GPU Time Stamps: ");
    printf("\n-------------------------------------");
    printf("--------------------------------\n");
    gb_gpu_time(my_thread);
    printf("Procs: %d\n",comm->nprocs);
    printf("-------------------------------------");
    printf("--------------------------------\n\n");
  }
  #pragma omp parallel
  {
    #ifdef GB_GPU_OMP
    int my_thread=omp_get_thread_num();
    #endif
    gb_gpu_clear(my_thread);
    if (ij_new[my_thread]!=NULL) {
      ij_new[my_thread]=NULL;
      delete [] ij_new[my_thread];
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairGayBerneGPU::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;
  if (vflag_atom) 
    error->all("Per-atom virial not available with GPU Gay-Berne.");

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int inum = list->inum;
  
  bool rebuild=false;
  if (neighbor->ncalls > last_neighbor) {
    last_neighbor=neighbor->ncalls;
    rebuild=true;
  }
  
  #pragma omp parallel
  {
    
  bool success=true;
  #ifdef GB_GPU_OMP
  int my_thread=omp_get_thread_num();
  if (rebuild) {
    omp_chunk=static_cast<int>(ceil(static_cast<double>(inum)/nthreads));
    if (my_thread==nthreads-1)
      thread_inum[my_thread]=inum-(nthreads-1)*omp_chunk;
    else
      thread_inum[my_thread]=omp_chunk;
    olist[my_thread]=gb_gpu_reset_nbors(nall, atom->nlocal, 
                                        thread_inum[my_thread],  
                                        list->ilist+omp_chunk*my_thread, 
                                        list->numneigh, atom->type, my_thread,
                                        success);
  }
  #else
  if (rebuild)
    olist[my_thread]=gb_gpu_reset_nbors(nall, atom->nlocal, inum, list->ilist, 
                                        list->numneigh, atom->type, my_thread,
                                        success);
  #endif
  if (!success)
    error->one("Total # of atoms exceeds maximum allowed per GPGPU.\n");
  
  // copy atom data to GPU
  gb_gpu_atom(atom->x,atom->quat,atom->type,rebuild,my_thread);

  int i,j,ii,jj,jnum;
  double factor_lj;
  int *jlist;

  if (rebuild==true) {
    int num_ij = 0;

    // loop over neighbors of my atoms
    int *ijp=ij_new[my_thread];
    #ifdef GB_GPU_OMP
    int mgo=my_thread*omp_chunk;
    int mgot=mgo+thread_inum[my_thread];
    #else
    int mgo=0, mgot=inum;
    #endif
    for (ii = mgo; ii<mgot; ii++) {
      i = olist[my_thread][ii];
      jlist = list->firstneigh[i];
      jnum = list->numneigh[i];
      
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];

        *ijp=j;
        ijp++;
        num_ij++;
          
        if (num_ij==ij_size) {
          memcpy(ij[my_thread],ij_new[my_thread],num_ij*sizeof(int));
          gb_gpu_nbors(num_ij,eflag,my_thread);
          ijp=ij_new[my_thread];
          num_ij=0;
        }
      }
    }
    if (num_ij>0) {
      memcpy(ij[my_thread],ij_new[my_thread],num_ij*sizeof(int));
      gb_gpu_nbors(num_ij,eflag,my_thread);
    }
  }
  
  gb_gpu_gayberne(eflag,vflag,rebuild,my_thread);
  double lvirial[6];
  for (int i=0; i<6; i++) lvirial[i]=0.0;
  double my_eng=gb_gpu_forces(atom->f,atom->torque,olist[my_thread],eflag,vflag,
                              eflag_atom, vflag_atom, eatom, vatom, lvirial,
                              my_thread);
  #pragma omp critical
  {
    eng_vdwl+=my_eng;
    virial[0]+=lvirial[0];
    virial[1]+=lvirial[1];
    virial[2]+=lvirial[2];
    virial[3]+=lvirial[3];
    virial[4]+=lvirial[4];
    virial[5]+=lvirial[5];
  }
  
  } //End parallel
  
  if (vflag_fdotr) virial_compute();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGayBerneGPU::settings(int narg, char **arg)
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

  PairGayBerne::settings(narg-2,&arg[2]);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGayBerneGPU::init_style()
{
  if (force->pair_match("gpu",0) == NULL)
    error->all("Cannot use pair hybrid with multiple GPU pair styles");
  if (!atom->quat_flag || !atom->torque_flag || !atom->avec->shape_type)
    error->all("Pair gayberne requires atom attributes quat, torque, shape");
  if (atom->radius_flag)
    error->all("Pair gayberne cannot be used with atom attribute diameter");

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

  // per-type shape precalculations

  for (int i = 1; i <= atom->ntypes; i++) {
    if (setwell[i]) {
      double *one = atom->shape[i];
      shape[i][0] = one[0]*one[0];
      shape[i][1] = one[1]*one[1];
      shape[i][2] = one[2]*one[2];
      lshape[i] = (one[0]*one[1]+one[2]*one[2])*sqrt(one[0]*one[1]);
    }
  }

  // Repeat cutsq calculation because done after call to init_style
  double cut;
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      cut = init_one(i,j);
      cutsq[i][j] = cutsq[j][i] = cut*cut;
    }

  // If compiled with OpenMP and only 1 proc, try to use multiple GPUs w/threads
  #ifdef GB_GPU_OMP
  if (multi_gpu_mode!=ONE_GPU)
    nthreads=ngpus=gb_gpu_num_devices();
  else
    nthreads=ngpus=1;
  if (nthreads>MAX_GPU_THREADS)
    nthreads=MAX_GPU_THREADS;
  omp_set_num_threads(nthreads);
  #endif
    
  #pragma omp parallel firstprivate(my_gpu)
  {
    #ifdef GB_GPU_OMP
    int my_thread = omp_get_thread_num();
    if (multi_gpu_mode!=ONE_GPU)
      my_gpu=my_thread;
    if (multi_gpu_mode==ONE_NODE)
      my_gpu+=multi_gpu_param;
    #endif
    
    ij[my_thread]=gb_gpu_init(ij_size, atom->ntypes+1, gamma, upsilon, mu, 
                              shape, well, cutsq, sigma, epsilon, lshape, form,
                              lj1, lj2, lj3, lj4, offset, force->special_lj, 
                              neighbor->oneatom, my_thread, my_gpu);
    if (ij[my_thread]==0)
      error->one("AT LEAST ONE PROCESS COULD NOT ALLOCATE A CUDA-ENABLED GPU.");
    
    if (ij_new[my_thread]!=NULL)
      delete [] ij_new[my_thread];
    ij_new[my_thread]=new int[ij_size];
  }
  
  last_neighbor = -1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  if (force->newton_pair) 
    error->all("Cannot use newton with GPU GayBerne pair style.");

  if (comm->me == 0 && screen) {
    printf("\n-------------------------------------");
    printf("-------------------------------------\n");
    printf("- Using GPGPU acceleration for Gay-Berne:\n");
    printf("-------------------------------------");
    printf("-------------------------------------\n");

    for (int i=0; i<ngpus; i++) {
      int gpui=my_gpu;
      if (multi_gpu_mode==ONE_NODE)
        gpui=i+multi_gpu_param;
      else if (multi_gpu_mode==MULTI_GPU)
        gpui=i;
      std::string gpu_string=gb_gpu_name(gpui,neighbor->oneatom);
      printf("GPU %d: %s\n",gpui,gpu_string.c_str());  
    }
    printf("-------------------------------------");
    printf("-------------------------------------\n\n");
  }
}

/* ---------------------------------------------------------------------- */

double PairGayBerneGPU::memory_usage()
{
  double bytes=Pair::memory_usage()+nthreads*ij_size*sizeof(int);
  return bytes+gb_gpu_bytes();
}
