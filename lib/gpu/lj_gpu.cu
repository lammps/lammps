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
#include <iostream>
#include <cassert>
#include <string.h>
#include "cudatimer.h"
#include "lj_tex.h"
#include "neigh.h"
#include "cell.h"
#include "lj_gpu_kernel.h"

#ifdef WINDLL
#define EXTERN extern "C" __declspec(dllexport) 
#else
#define EXTERN 
#endif

static float h_boxlo[3], h_boxhi[3];
static float cell_size;
static float *energy = NULL, *d_energy = NULL;
static float3 *d_force = NULL, *f_temp = NULL, *v_temp = NULL, *d_virial = NULL;
static float4 *d_pos = NULL, *temp_pos = NULL;
static int *d_type = NULL;
static int ncellx, ncelly, ncellz;

static neigh_list_gpu d_neigh_list;
static cell_list_gpu d_cell_list;

#define TIMING(x) 

// ---------------------------------------------------------------------------
// Return string with GPU info
// ---------------------------------------------------------------------------
EXTERN void lj_gpu_name(const int id, const int max_nbors, char * name) 
{
  struct cudaDeviceProp prop;
  CUDA_SAFE_CALL( cudaGetDeviceProperties(&prop, id) );
#ifdef _WIN32
  strcpy_s(name, strlen(prop.name)+1, prop.name);
#else
  strncpy(name, prop.name, strlen(prop.name)+1);
#endif
}

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
EXTERN bool lj_gpu_init(int &ij_size, const int ntypes, 
			double **cutsq,double **sigma, 
			 double **epsilon, double **host_lj1, double **host_lj2, 
			 double **host_lj3, double **host_lj4, double **offset, 
			 double *special_lj, double *boxlo, double *boxhi, 
			 double cellsize, double skin,
			 const int max_nbors, const int gpu_id) 
{
  int num_devices;

  /* get device count */
  CUDA_SAFE_CALL( cudaGetDeviceCount(&num_devices) );
  if (num_devices == 0) {
    printf("NO CUDA-capable GPU detected.\n");
    exit(1);
  }

  if (gpu_id > num_devices) {
    printf("gpu_id %d is larger than the number of GPUs %d\n", 
	   gpu_id, num_devices);
    exit(1);
  }

  /* set CUDA device to the specified GPU */
  cudaThreadExit();
  CUDA_SAFE_CALL( cudaSetDevice(gpu_id) );
  
  ij_size=0;

  cell_size = cellsize;
  ncellx = ceil(((boxhi[0] - boxlo[0]) + 2.0*cell_size) / cell_size);
  ncelly = ceil(((boxhi[1] - boxlo[1]) + 2.0*cell_size) / cell_size);
  ncellz = ceil(((boxhi[2] - boxlo[2]) + 2.0*cell_size) / cell_size);
   
  for (int i = 0; i < 3; i++) {
    h_boxhi[i] = boxhi[i];
    h_boxlo[i] = boxlo[i];
  }

  init_force_const(ntypes, cutsq, host_lj1, host_lj2, host_lj3, host_lj4, offset);

  init_cell_list_const(cellsize, skin, boxlo, boxhi);

  return true;
}

// ---------------------------------------------------------------------------
// Clear memory on host and device
// ---------------------------------------------------------------------------
EXTERN void lj_gpu_clear() {

  free(energy);
  free(v_temp);
  CUDA_SAFE_CALL( cudaFreeHost(f_temp) );
  if (d_force) CUDA_SAFE_CALL( cudaFree(d_force) );
  if (d_energy) CUDA_SAFE_CALL( cudaFree(d_energy) );
  if (d_virial) CUDA_SAFE_CALL( cudaFree(d_virial) );
  if (d_pos) CUDA_SAFE_CALL( cudaFree(d_pos) );
  if (d_type) CUDA_SAFE_CALL( cudaFree(d_type) );
  if (temp_pos) CUDA_SAFE_CALL( cudaFreeHost(temp_pos) );
  clear_neigh_list_gpu(d_neigh_list);
  clear_cell_list_gpu(d_cell_list);

  if (useCache) {
    unbind_pos();
    unbind_type();
  }


  //LJMF.clear();
}


template <class numtyp, class acctyp>
double _lj_gpu_neigh(double **force, double *virial,
		     double **host_x, int *host_type, const int inum, 
		     const int nall, const int ago, const bool eflag, const bool vflag, 
		     const double *boxlo, const double *boxhi)
{

  double evdwl=0.0;

  static int first_call = 1;
  
  TIMING( static CUDATimer cuTimer );  
  TIMING( static CTimer cTimer );
  TIMING( static CTimer cTimer2 );
  
  double *atom_pos = host_x[0];

  static int szTailList = inum*32;
  
  TIMING( cTimer.Start() );
  TIMING( cTimer2.Start() );
   
  /* MPI communication just happened, reallocate space using new inum & nall
     FIXME: this is costly: ~ total kernel time! Use a DIY GPU memory allocator.*/

  if (first_call || ago == 0) {

    if (!first_call) {
      if (useCache) {
	unbind_pos();
	unbind_type();
      }
      
      CUDA_SAFE_CALL( cudaFree(d_force) );
      CUDA_SAFE_CALL( cudaFree(d_energy) );
      CUDA_SAFE_CALL( cudaFree(d_virial) );
      CUDA_SAFE_CALL( cudaFree(d_pos) );
      CUDA_SAFE_CALL( cudaFree(d_type) );

      clear_neigh_list_gpu(d_neigh_list);

      CUDA_SAFE_CALL( cudaFreeHost(f_temp) );
      CUDA_SAFE_CALL( cudaFreeHost(temp_pos) );

      free(energy);
      free(v_temp);
    }

    CUDA_SAFE_CALL( cudaMalloc((void**)&d_force,     inum*sizeof(float3)) );
    CUDA_SAFE_CALL( cudaMalloc((void**)&d_energy,    inum*sizeof(float)) );
    CUDA_SAFE_CALL( cudaMalloc((void**)&d_virial,    inum*3*sizeof(float3)) );
    CUDA_SAFE_CALL( cudaMalloc((void**)&d_pos, nall*sizeof(float4)) );
    CUDA_SAFE_CALL( cudaMalloc((void**)&d_type, nall*sizeof(int)) );
    
    init_neigh_list_gpu(d_neigh_list, inum, NEIGH_BIN_SIZE, szTailList);

    CUDA_SAFE_CALL( cudaMallocHost((void**)&temp_pos, nall*sizeof(float4)) );
    CUDA_SAFE_CALL( cudaMallocHost((void**)&f_temp,   inum*sizeof(float3)) );

    energy    = (float*) malloc(inum*sizeof(float));
    v_temp    = (float3*)malloc(inum*2*sizeof(float3));

    if (useCache) {
      bind_pos(d_pos, nall);
      bind_type(d_type, nall);
    }

    first_call = 0;
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
    CUDA_SAFE_CALL( cudaGetLastError() );
    CUDA_SAFE_CALL( cudaMemcpy(d_type, host_type, nall*sizeof(int), 
			       cudaMemcpyHostToDevice) );

  }

  TIMING( static double mallocTime = 0. );
  TIMING( mallocTime += cTimer2.GetET() );
  TIMING( printf("malloc time = %f ms\n", mallocTime*1e3) );

  TIMING( cTimer2.Start() );
  for (int i = 0; i < 3*nall; i+=3) { 
    temp_pos[i/3] = make_float4(atom_pos[i], atom_pos[i+1], atom_pos[i+2], 0.f);
  }

  TIMING( static double copyTime = 0. );
  TIMING( copyTime += cTimer2.GetET() );
  TIMING( printf("position copy time = %f ms\n", copyTime*1e3) );

  
  TIMING( cTimer2.Start() );
  CUDA_SAFE_CALL( cudaMemcpy(d_pos, temp_pos, nall*sizeof(float4), cudaMemcpyHostToDevice) );

  TIMING( static double h2dTime = 0. );
  TIMING( h2dTime += cTimer2.GetET() );
  TIMING( printf("h2d copy time = %f ms\n", h2dTime*1e3) );

  TIMING( cTimer2.Start() );
  if (ago == 0) {
    build_neigh_list_gpu(d_pos,
			 d_neigh_list,
			 h_boxlo, h_boxhi, cell_size,
			 inum, nall);
  }
  TIMING( static double neighTime = 0. );
  TIMING( neighTime += cTimer2.GetET() );
  TIMING( printf("Neigh List time = %f ms\n", neighTime*1e3) );

  TIMING( cTimer2.Start() );
  calc_lj_neigh_gpu(d_force, d_energy, d_virial,
		    d_pos, d_type,
		    d_neigh_list,
		    inum, nall,
		    eflag, vflag);
  TIMING( static double forceTime = 0. );
  TIMING( forceTime += cTimer2.GetET() );
  TIMING( printf("Force time = %f ms\n", forceTime*1e3) );
  TIMING( printf("GPU kernel time = %f ms\n", (forceTime + neighTime)*1e3) );


  TIMING( cTimer2.Start() );
  CUDA_SAFE_CALL( cudaMemcpy(f_temp, d_force, inum*sizeof(float3), cudaMemcpyDeviceToHost) );
  TIMING( static double d2hTime = 0. );
  TIMING( d2hTime += cTimer2.GetET() );
  TIMING( printf("d2h copy time = %f ms\n", d2hTime*1e3) );
  TIMING( printf("GPU-CPU data transfer time = %f ms\n", (h2dTime+d2hTime)*1e3) );

  TIMING( cTimer2.Start() );

  for (int i = 0; i < inum; i++) {
    force[i][0] += f_temp[i].x;
    force[i][1] += f_temp[i].y;
    force[i][2] += f_temp[i].z;
  }

  if (eflag) {
    CUDA_SAFE_CALL( cudaMemcpy(energy, d_energy, 
			       inum*sizeof(float), cudaMemcpyDeviceToHost) );
    for (int i = 0; i < inum; i++) {
      evdwl += energy[i];
    }
    evdwl *= 0.5f;
  }
  
  if (vflag) {
    CUDA_SAFE_CALL( cudaMemcpy(v_temp, d_virial, inum*2*sizeof(float3), 
			       cudaMemcpyDeviceToHost) ); 
    for (int i = 0; i < inum; i++) {
      virial[0] += v_temp[2*i].x;
      virial[1] += v_temp[2*i].y;
      virial[2] += v_temp[2*i].z;
      virial[3] += v_temp[2*i+1].x;
      virial[4] += v_temp[2*i+1].y;
      virial[5] += v_temp[2*i+1].z;
    }
    for (int i = 0; i < 6; i++) 
      virial[i] *= 0.5f;
  }


  TIMING( static double postTime = 0. );
  TIMING( postTime += cTimer2.GetET() );
  TIMING( printf("postprocess Time = %f ms\n", postTime*1e3) );
  TIMING( printf("Data process time = %f ms\n", (postTime+copyTime)*1e3) );

  TIMING( static double totalTime = 0. );
  TIMING( totalTime += cTimer.GetET() );
  TIMING( printf("lj_gpu time = %f ms\n", totalTime*1e3) );

  return evdwl;
 
}

EXTERN double lj_gpu_neigh(double **force, double *virial, 
			  double **host_x, int *host_type, 
			  const int inum, const int nall, 
			  const int ago, const bool eflag, const bool vflag, 
			  const double *boxlo, const double *boxhi) 
{
  return _lj_gpu_neigh<float,float>(force, virial, 
				    host_x, host_type, inum, nall, 
				    ago, eflag, vflag, boxlo, boxhi);
}


template <class numtyp, class acctyp>
double _lj_gpu_cell(double **force, double *virial,
		    double **host_x, int *host_type, const int inum, 
		    const int nall, const int ago, 
		    const bool eflag, const bool vflag, 
		    const double *boxlo, const double *boxhi)
{
  
  double evdwl=0.0;

  static int ncell = ncellx*ncelly*ncellz;

  static int first_call = 1;

  // allocate memory on CPU and GPU
  if (first_call || ago == 0) {
    if (!first_call) {
     if (useCache) {
	unbind_pos();
	unbind_type();
      }

      free(energy);
      free(v_temp);
      
      CUDA_SAFE_CALL( cudaFree(d_force) );
      CUDA_SAFE_CALL( cudaFree(d_energy) );
      CUDA_SAFE_CALL( cudaFree(d_virial) );

      CUDA_SAFE_CALL( cudaFree(d_pos) );
      CUDA_SAFE_CALL( cudaFree(d_type) );
      CUDA_SAFE_CALL( cudaFreeHost(f_temp) );
      CUDA_SAFE_CALL( cudaFreeHost(temp_pos) );

      clear_cell_list_gpu(d_cell_list);
    }

    energy    = (float*) malloc(inum*sizeof(float));
    v_temp    = (float3*)malloc(inum*2*sizeof(float3));


    cudaMalloc((void**)&d_force,     inum*sizeof(float3));
    cudaMalloc((void**)&d_energy,    inum*sizeof(float));
    cudaMalloc((void**)&d_virial,    inum*3*sizeof(float3));

    CUDA_SAFE_CALL( cudaMalloc((void**)&d_pos, nall*sizeof(float4)) );
    CUDA_SAFE_CALL( cudaMalloc((void**)&d_type, nall*sizeof(int)) );

    CUDA_SAFE_CALL( cudaMallocHost((void**)&f_temp,   inum*sizeof(float3)) );
    CUDA_SAFE_CALL( cudaMallocHost((void**)&temp_pos, nall*sizeof(float4)) );

    init_cell_list_gpu(d_cell_list, nall, ncell);

    CUDA_SAFE_CALL( cudaMemcpy(d_type, host_type, nall*sizeof(int), 
			       cudaMemcpyHostToDevice) );

    if (useCache) {
      bind_pos(d_pos, nall);
      bind_type(d_type, nall);
    }

    first_call = 0;
  }

  /* build cell-list on GPU */
  double *atom_pos = host_x[0];
  for (int i = 0; i < 3*nall; i+=3) { 
    temp_pos[i/3] = make_float4(atom_pos[i], atom_pos[i+1], atom_pos[i+2], 0.f);
  }
  CUDA_SAFE_CALL( cudaMemcpy(d_pos, temp_pos, nall*sizeof(float4), 
			     cudaMemcpyHostToDevice) );
  if (ago == 0) {
    build_cell_list_gpu(d_pos, d_cell_list, h_boxlo, h_boxhi, 
			cell_size, inum, nall);
  }

  calc_lj_cell_gpu(d_force, d_energy, d_virial,
		   d_pos, d_type, d_cell_list,
		   inum, nall, ncellx, 
		   ncelly, ncellz, cell_size,
		   eflag, vflag);

  CUDA_SAFE_CALL( cudaMemcpy(f_temp, d_force, inum*sizeof(float3), 
			     cudaMemcpyDeviceToHost) );

  for (int i = 0; i < inum; i++) {
    force[i][0] += f_temp[i].x;
    force[i][1] += f_temp[i].y;
    force[i][2] += f_temp[i].z;
  }
  
  if (eflag) {
    CUDA_SAFE_CALL( cudaMemcpy(energy, d_energy, 
			       inum*sizeof(float), cudaMemcpyDeviceToHost) );
    for (int i = 0; i < inum; i++) {
      evdwl += energy[i];
    }
    evdwl *= 0.5f;
  }
  
  if (vflag) {
    CUDA_SAFE_CALL( cudaMemcpy(v_temp, d_virial, inum*2*sizeof(float3), 
			       cudaMemcpyDeviceToHost) ); 
    for (int i = 0; i < inum; i++) {
      virial[0] += v_temp[2*i].x;
      virial[1] += v_temp[2*i].y;
      virial[2] += v_temp[2*i].z;
      virial[3] += v_temp[2*i+1].x;
      virial[4] += v_temp[2*i+1].y;
      virial[5] += v_temp[2*i+1].z;
    }
    for (int i = 0; i < 6; i++) 
      virial[i] *= 0.5f;
  }

  return evdwl; 
}

EXTERN double lj_gpu_cell(double **force, double *virial, 
			  double **host_x, int *host_type, 
			  const int inum, const int nall, 
			  const int ago, const bool eflag, const bool vflag, 
			  const double *boxlo, const double *boxhi) 
{
  return _lj_gpu_cell<float,float>(force, virial, 
				   host_x, host_type, inum, nall, 
				   ago, eflag, vflag, boxlo, boxhi);
}

EXTERN void lj_gpu_time() {
  /*  cout.precision(4);
  cout << "Atom copy:     " << LJMF.time_atom.total_seconds() << " s.\n";
  cout << "Neighbor copy: " << LJMF.time_nbor.total_seconds() << " s.\n";
  cout << "LJ calc:       " << LJMF.time_pair.total_seconds() << " s.\n";*/
  //cout << "Answer copy:   " << LJMF.time_answer.total_seconds() << " s.\n";
}

EXTERN int lj_gpu_num_devices() {
  int num_devices;
  CUDA_SAFE_CALL( cudaGetDeviceCount(&num_devices) );
  return num_devices;
}

EXTERN double lj_gpu_bytes() {
  return 0.0;
}
