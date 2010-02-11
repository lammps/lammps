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

#ifndef PAIR_GPU_CELL_H
#define PAIR_GPU_CELL_H

#ifdef WINDLL
#include <windows.h>
#endif

#ifdef WINDLL
#define EXTERN extern "C" __declspec(dllexport) 
#else
#define EXTERN 
#endif
using namespace std;

static float kernelTime = 0.0;
static int ncellx, ncelly, ncellz;
static float *energy, *d_energy;
static float3 *d_force, *f_temp, *v_temp, *d_virial;


typedef struct {
  float3 *pos;
  unsigned int *idx;
  int *type;
  int *natom;
} cell_list;

static cell_list cell_list_gpu;

__global__ void kernel_set_cell_list(unsigned int *cell_idx);
__global__ void kernel_build_cell_list(float3 *cell_list, 
				       unsigned int *cell_idx, 
				       int *cell_type, 
				       int *cell_atom,
				       float3 *pos, 
				       int *type, 
				       const int inum, 
				       const int nall);
__global__ void kernel_test_rebuild(float3 *cell_list, int *cell_atom, int *rebuild);
__global__ void kernel_copy_list(float3 *cell_list, 
				 unsigned int *cell_idx, 
				 int *cell_atom, 
				 float3 *pos);
__global__ void kernel_test_overflow(int *cell_atom, int *overflow, const int ncell);
void sortBlocks(unsigned int *keys, float3 *values1, int *values2, const int size);

void init_cell_list_const(double cell_size, double skin,
			 double *boxlo, double *boxhi);
void init_cell_list(cell_list &cell_list_gpu, 
		   const int nall,
		   const int ncell, 
		   const int buffer);

void build_cell_list(double *atom_pos, int *atom_type, 
		    cell_list &cell_list_gpu, 
		    const int ncell, const int ncellx, const int ncelly, const int ncellz, 
		    const int buffer, const int inum, const int nall, const int ago);

void clear_cell_list(cell_list &cell_list_gpu);

#endif
