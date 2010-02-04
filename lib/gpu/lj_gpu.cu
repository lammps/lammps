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
#include "nvc_macros.h"
#include "nvc_timer.h"
#include "nvc_device.h"
#include "pair_gpu_texture.h"
#include "pair_gpu_cell.h"
#include "lj_gpu_memory.cu"
#include "lj_gpu_kernel.h"

#ifdef WINDLL
#include <windows.h>
BOOL APIENTRY DllMain(HANDLE hModule, DWORD dwReason, LPVOID lpReserved)
{
    return TRUE;
}
#endif

#ifdef WINDLL
#define EXTERN extern "C" __declspec(dllexport) 
#else
#define EXTERN 
#endif
using namespace std;

static LJ_GPU_Memory<PRECISION,ACC_PRECISION> LJMF;
#define LJMT LJ_GPU_Memory<numtyp,acctyp>

static float kernelTime = 0.0;
static int ncell1D;
static float *energy, *d_energy;
static float3 *d_force, *f_temp, *v_temp, *d_virial;
static cell_list cell_list_gpu;

// ---------------------------------------------------------------------------
// Convert something to a string
// ---------------------------------------------------------------------------
#include <sstream>

template <class t>
inline string lj_gpu_toa(const t& in) {
  ostringstream o;
  o.precision(2);
  o << in;
  return o.str();
}

// ---------------------------------------------------------------------------
// Return string with GPU info
// ---------------------------------------------------------------------------
EXTERN void lj_gpu_name(const int id, const int max_nbors, char * name) {
  string sname=LJMF.gpu.name(id)+", "+
              lj_gpu_toa(LJMF.gpu.cores(id))+" cores, "+
              lj_gpu_toa(LJMF.gpu.gigabytes(id))+" GB, "+
              lj_gpu_toa(LJMF.gpu.clock_rate(id))+" GHZ";
  strcpy(name,sname.c_str());
}

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
EXTERN bool lj_gpu_init(int &ij_size, const int ntypes, double **cutsq,double **sigma, 
			 double **epsilon, double **host_lj1, double **host_lj2, 
			 double **host_lj3, double **host_lj4, double **offset, 
			 double *special_lj, double *boxlo, double *boxhi, 
			 double cell_size, double skin,
			 const int max_nbors, const int gpu_id) {
  LJMF.gpu.init();
  if (LJMF.gpu.num_devices()==0)
    return false;                   

  ij_size=IJ_SIZE;

  bool ret = LJMF.init(ij_size, ntypes, cutsq, sigma, epsilon, host_lj1, host_lj2, 
		       host_lj3, host_lj4, offset, special_lj, max_nbors, gpu_id);

  ncell1D = ceil(((boxhi[0] - boxlo[0]) + 2.0*cell_size) / cell_size);

  init_cell_list_const(cell_size, skin, boxlo, boxhi);

  return ret;
}

// ---------------------------------------------------------------------------
// Clear memory on host and device
// ---------------------------------------------------------------------------
EXTERN void lj_gpu_clear() {
  free(energy);
  free(v_temp);
  cudaFreeHost(f_temp);
  cudaFree(d_force);
  cudaFree(d_energy);
  cudaFree(d_virial);
  clear_cell_list(cell_list_gpu);

  LJMF.clear();
}


// ---------------------------------------------------------------------------
// Calculate energies and forces for all ij interactions
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void _lj_gpu(LJMT &ljm, const bool eflag, const bool vflag, const bool rebuild){
  // Compute the block size and grid size to keep all cores busy
  const int BX=BLOCK_1D;

  int GX=static_cast<int>(ceil(static_cast<double>(ljm.atom.inum())/BX));

  ljm.time_pair.start();

  if (ljm.shared_types)
    kernel_lj_fast<numtyp,acctyp><<<GX,BX,0,ljm.pair_stream>>>
           (ljm.special_lj.begin(), ljm.nbor.dev_nbor.begin(), 
            ljm.nbor.ij.begin(), ljm.nbor.dev_nbor.row_size(), 
            ljm.atom.ans.begin(), ljm.atom.ans.row_size(), eflag,
            vflag, ljm.atom.inum(), ljm.atom.nall());
  else
    kernel_lj<numtyp,acctyp><<<GX,BX,0,ljm.pair_stream>>>
           (ljm.special_lj.begin(), ljm.nbor.dev_nbor.begin(), 
            ljm.nbor.ij.begin(), ljm.nbor.dev_nbor.row_size(), 
            ljm.atom.ans.begin(), ljm.atom.ans.row_size(), eflag, 
            vflag, ljm.atom.inum(), ljm.atom.nall());
	    ljm.time_pair.stop();
}

EXTERN void lj_gpu(const bool eflag, const bool vflag, const bool rebuild) {
  _lj_gpu<PRECISION,ACC_PRECISION>(LJMF, eflag,vflag,rebuild);
}

template <class numtyp, class acctyp>
double _lj_gpu_cell(LJMT &ljm, double **force, double *virial,
		    double **host_x, int *host_type, const int inum, 
		    const int nall, const int ago, const bool eflag, const bool vflag, 
		    const double *boxlo, const double *boxhi)
{
  ljm.atom.nall(nall);
  ljm.atom.inum(inum);

  ljm.nbor.time_nbor.start();
  ljm.nbor.time_nbor.stop();

  double evdwl=0.0;

  static int buffer = CELL_SIZE;
  static int ncell = (int)pow((float)ncell1D,3);

  static int first_call = 1;

  // allocate memory on CPU and GPU
  if (first_call) {
    energy    = (float*) malloc(inum*sizeof(float));
    v_temp    = (float3*)malloc(inum*2*sizeof(float3));
    cudaMallocHost((void**)&f_temp,   inum*sizeof(float3));

    cudaMalloc((void**)&d_force,     inum*sizeof(float3));
    cudaMalloc((void**)&d_energy,    inum*sizeof(float));
    cudaMalloc((void**)&d_virial,    inum*3*sizeof(float3));

    init_cell_list(cell_list_gpu, nall, ncell, buffer);

    first_call = 0;
  }

  if (!first_call && ago == 0) {
    free(energy);
    free(v_temp);
    cudaFreeHost(f_temp);
    cudaFree(d_force);
    cudaFree(d_energy);
    cudaFree(d_virial);

    energy    = (float*) malloc(inum*sizeof(float));
    v_temp    = (float3*)malloc(inum*2*sizeof(float3));
    cudaMallocHost((void**)&f_temp,   inum*sizeof(float3));

    cudaMalloc((void**)&d_force,     inum*sizeof(float3));
    cudaMalloc((void**)&d_energy,    inum*sizeof(float));
    cudaMalloc((void**)&d_virial,    inum*3*sizeof(float3));

    clear_cell_list(cell_list_gpu);
    init_cell_list(cell_list_gpu, nall, ncell, buffer);
  }

  // build cell-list on GPU
  ljm.atom.time_atom.start();
  build_cell_list(host_x[0], host_type, cell_list_gpu, 
		  ncell, ncell1D, buffer, inum, nall, ago);
  ljm.atom.time_atom.stop();

  ljm.time_pair.start();

#ifdef TIMING
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
#endif

  // call the cell-list force kernel
  const int BX=BLOCK_1D;
  dim3 GX(ncell1D, ncell1D*ncell1D);
  if (eflag == 0 && vflag == 0) {
    kernel_lj_cell<false,false><<<GX, BX, 0>>>
      (d_force, d_energy, d_virial, 
       cell_list_gpu.pos, 
       cell_list_gpu.idx, 
       cell_list_gpu.type, 
       cell_list_gpu.natom,
       inum, nall, ncell);
  } else {
    kernel_lj_cell<true,true><<<GX, BX, 3*sizeof(float)*MAX_SHARED_TYPES*MAX_SHARED_TYPES>>>
      (d_force, d_energy, d_virial, 
       cell_list_gpu.pos, 
       cell_list_gpu.idx, 
       cell_list_gpu.type, 
       cell_list_gpu.natom,
       inum, nall, ncell);
  }

#ifdef TIMING
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  float kTime;
  cudaEventElapsedTime(&kTime, start, stop);
  kernelTime += kTime;
  printf("kernelTime = %f, eflag=%d, vflag=%d\n", kTime, eflag, vflag);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
#endif

  // copy results from GPU to CPU
  cudaMemcpy(f_temp, d_force, inum*sizeof(float3), cudaMemcpyDeviceToHost);
  if (eflag) {
    cudaMemcpy(energy, d_energy, inum*sizeof(float), cudaMemcpyDeviceToHost);
    for (int i = 0; i < inum; i++) {
      evdwl += energy[i];
    }
    evdwl *= 0.5f;
  }
  if (vflag) {
    cudaMemcpy(v_temp, d_virial, inum*2*sizeof(float3), cudaMemcpyDeviceToHost);
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

  for (int i = 0; i < inum; i++) {
    force[i][0] += f_temp[i].x;
    force[i][1] += f_temp[i].y;
    force[i][2] += f_temp[i].z;
  }

  ljm.time_pair.stop();

  ljm.atom.time_atom.add_to_total();
  ljm.nbor.time_nbor.add_to_total();
  ljm.time_pair.add_to_total();


  return evdwl;
 
}

EXTERN double lj_gpu_cell(double **force, double *virial, double **host_x, int *host_type, const int inum, const int nall, 
		   const int ago, const bool eflag, const bool vflag, 
		   const double *boxlo, const double *boxhi) 
{
  return _lj_gpu_cell<PRECISION,ACC_PRECISION>(LJMF, force, virial, host_x, host_type, inum, nall, 
					       ago, eflag, vflag, boxlo, boxhi);
}

template <class numtyp, class acctyp>
double _lj_gpu_n2(LJMT &ljm, double **force, double *virial,
		  double **host_x, int *host_type, const int inum, const int nall, const bool eflag, const bool vflag,
		  const double *boxlo, const double *boxhi)
{
  ljm.atom.nall(nall);
  ljm.atom.inum(inum);


  ljm.nbor.time_nbor.start();
  ljm.nbor.time_nbor.stop();

  
  double evdwl=0.0;

#ifdef NOUSE
  static int first_call = 1;

  if (first_call) {
    energy    = (float*)    malloc(inum*sizeof(float));
    v_temp    = (float3*)   malloc(inum*2*sizeof(float3));
    cudaMallocHost((void**)&f_temp,   inum*sizeof(float3));
    cudaMallocHost((void**)&pos_temp, nall*sizeof(float3));
    cudaMalloc((void**)&d_force,     inum*sizeof(float3));
    cudaMalloc((void**)&d_energy,    inum*sizeof(float));
    cudaMalloc((void**)&d_virial,    inum*3*sizeof(float3));
    cudaMalloc((void**)&d_pos,       nall*sizeof(float3));
    cudaMalloc((void**)&d_type,      nall*sizeof(int));
    first_call = 0;
  }


  ljm.atom.time_atom.start();
  double *atom_pos = host_x[0];
  for (int i = 0; i < 3*nall; i+=3) { 
    pos_temp[i/3] = make_float3(atom_pos[i], atom_pos[i+1], atom_pos[i+2]);
  }
  cudaMemcpy(d_pos, pos_temp, nall*sizeof(float3), cudaMemcpyHostToDevice);
  cudaMemcpy(d_type, host_type, nall*sizeof(int),  cudaMemcpyHostToDevice);

  ljm.atom.time_atom.stop();

  ljm.time_pair.start();
  
  // Compute the block size and grid size to keep all cores busy
  const int BX=BLOCK_1D;
  dim3 GX(static_cast<int>(ceil(static_cast<double>(inum)/BX)));

#ifdef TIMING
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
#endif

  // N^2 force kernel
  kernel_lj_n2<numtyp, acctyp><<<GX, BX>>>(d_force, d_energy, d_virial, 
					   d_pos, d_type, eflag, vflag, inum, nall);

#ifdef TIMING
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  float kTime;
  cudaEventElapsedTime(&kTime, start, stop);
  kernelTime += kTime;
  printf("kernelTime = %f, eflag=%d, vflag=%d\n", kTime, eflag, vflag);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
#endif

  // copy results from GPU to CPU
  cudaMemcpy(f_temp, d_force, inum*sizeof(float3), cudaMemcpyDeviceToHost);
  if (eflag) {
    cudaMemcpy(energy, d_energy, inum*sizeof(float), cudaMemcpyDeviceToHost);
    for (int i = 0; i < inum; i++) {
      evdwl += energy[i];
    }
    evdwl *= 0.5f;
  }
  if (vflag) {
    cudaMemcpy(v_temp, d_virial, inum*2*sizeof(float3), cudaMemcpyDeviceToHost);
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

  for (int i = 0; i < inum; i++) {
    force[i][0] += f_temp[i].x;
    force[i][1] += f_temp[i].y;
    force[i][2] += f_temp[i].z;
  }
#endif
  ljm.time_pair.stop();

  ljm.atom.time_atom.add_to_total();
  ljm.nbor.time_nbor.add_to_total();
  ljm.time_pair.add_to_total();

  return evdwl;
}

EXTERN double lj_gpu_n2(double **force, double *virial, double **host_x, int *host_type, const int inum, const int nall, 
		 const bool eflag, const bool vflag,
		 const double *boxlo, const double *boxhi) 
{
  return _lj_gpu_n2<PRECISION,ACC_PRECISION>(LJMF, force, virial, host_x, host_type, inum, nall, 
					  eflag, vflag, boxlo, boxhi);
}

EXTERN void lj_gpu_time() {
  cout.precision(4);
  cout << "Atom copy:     " << LJMF.atom.time_atom.total_seconds() << " s.\n";
  cout << "Neighbor copy: " << LJMF.nbor.time_nbor.total_seconds() << " s.\n";
  cout << "LJ calc:       " << LJMF.time_pair.total_seconds() << " s.\n";
  cout << "Answer copy:   " << LJMF.atom.time_answer.total_seconds() << " s.\n";
}

EXTERN int lj_gpu_num_devices() {
  return LJMF.gpu.num_devices();
}

EXTERN double lj_gpu_bytes() {
  return LJMF.host_memory_usage();
}
