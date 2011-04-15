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
   Contributing authors: Peng Wang (Nvidia), penwang@nvidia.com
                         Mike Brown (ORNL), brownw@ornl.gov
------------------------------------------------------------------------- */

#ifdef NV_KERNEL

#include "nv_kernel_def.h"
texture<float4> neigh_tex;

#ifdef _DOUBLE_DOUBLE
__inline double4 fetch_pos(const int i, const double4 *pos)
{
  return pos[i];
}
#else
__inline float4 fetch_pos(const int& i, const float4 *pos)
{
  return tex1Dfetch(neigh_tex, i);
}
#endif

#else

#define fetch_pos(i,y) x_[i]
#define BLOCK_NBOR_BUILD 64

#endif

#ifdef _DOUBLE_DOUBLE
#define numtyp double
#define numtyp4 double4
#endif

#ifdef _SINGLE_DOUBLE
#define numtyp float
#define numtyp4 float4
#endif

#ifndef numtyp
#define numtyp float
#define numtyp4 float4
#endif

#define BLOCK_CELL_2D 8

#define SBBITS 30

#define SBBITS 30

__kernel void transpose(int *out, int *in, int columns_in, int rows_in)
{
	__local float block[BLOCK_CELL_2D][BLOCK_CELL_2D+1];
	
	unsigned ti=THREAD_ID_X;
	unsigned tj=THREAD_ID_Y;
	unsigned bi=BLOCK_ID_X;
	unsigned bj=BLOCK_ID_Y;
	
	unsigned i=bi*BLOCK_CELL_2D+ti;
	unsigned j=bj*BLOCK_CELL_2D+tj;
	if ((i<columns_in) && (j<rows_in))
		block[tj][ti]=in[j*columns_in+i];

	__syncthreads();

	i=bj*BLOCK_CELL_2D+ti;
	j=bi*BLOCK_CELL_2D+tj;
	if ((i<rows_in) && (j<columns_in))
		out[j*rows_in+i] = block[ti][tj];
}

__kernel void calc_cell_id(numtyp4 *pos, unsigned *cell_id, int *particle_id,
                           numtyp boxlo0, 
                           numtyp boxlo1, numtyp boxlo2, numtyp boxhi0, 
                           numtyp boxhi1, numtyp boxhi2, numtyp cell_size, 
                           int ncellx, int ncelly, int nall) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;

  if (i < nall) {
    numtyp4 p = fetch_pos(i,pos); //pos[i];

    p.x -= boxlo0;
    p.y -= boxlo1;
    p.z -= boxlo2;
    
    p.x = fmaxf(p.x, -cell_size);
    p.x = fminf(p.x, boxhi0-boxlo0+cell_size);
    p.y = fmaxf(p.y, -cell_size);
    p.y = fminf(p.y, boxhi1-boxlo1+cell_size);
    p.z = fmaxf(p.z, -cell_size);
    p.z = fminf(p.z, boxhi2-boxlo2+cell_size);
    
    unsigned int id = (unsigned int)(p.x/cell_size + 1.0) 
      + (unsigned int)(p.y/cell_size + 1.0) * ncellx
      + (unsigned int)(p.z/cell_size + 1.0) * ncellx * ncelly;
    
    cell_id[i] = id;
    particle_id[i] = i;
  }
}

__kernel void kernel_calc_cell_counts(unsigned *cell_id,
                                      int *cell_counts, int nall, int ncell) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < nall) {
    int id = cell_id[idx];

    // handle boundary cases
    if (idx == 0) {
      for (int i = 0; i < id + 1; i++) 
        cell_counts[i] = 0;
    }
    if (idx == nall - 1) {
      for (int i = id+1; i <= ncell; i++) 
        cell_counts[i] = nall;
    }

    if (idx > 0 && idx < nall) {
      int id_l = cell_id[idx-1];
      if (id != id_l) {
        for (int i = id_l+1; i <= id; i++) 
          cell_counts[i] = idx;
      }
    }
  }
}

__kernel void calc_neigh_list_cell(numtyp4 *pos,
				     int *cell_particle_id, 
				     int *cell_counts,
				     int *nbor_list,
				     int *host_nbor_list,
				     int *host_numj, 
				     int neigh_bin_size, 
				     numtyp cell_size,
				     int ncellx, int ncelly, int ncellz,
				     int inum, int nt, int nall)
{
  int tid = threadIdx.x;
  int ix = blockIdx.x;
  int iy = blockIdx.y % ncelly;
  int iz = blockIdx.y / ncelly;
	  
  int icell = ix + iy*ncellx + iz*ncellx*ncelly;

  __shared__ int cell_list_sh[BLOCK_NBOR_BUILD];
  __shared__ numtyp4 pos_sh[BLOCK_NBOR_BUILD];

  int icell_begin = cell_counts[icell];
  int icell_end = cell_counts[icell+1];

  int nborz0 = max(iz-1,0), nborz1 = min(iz+1, ncellz-1),
      nbory0 = max(iy-1,0), nbory1 = min(iy+1, ncelly-1),
      nborx0 = max(ix-1,0), nborx1 = min(ix+1, ncellx-1);

  numtyp4 diff;
  numtyp r2;
  for (int ii = 0; ii < ceil((numtyp)(icell_end - icell_begin)/blockDim.x); ii++) {
    int i = icell_begin + tid + ii*blockDim.x;
    int pid_i = nall, pid_j, stride;
    numtyp4 atom_i, atom_j;
    int cnt = 0;    
    int *neigh_counts, *neigh_list;
    
    if (i < icell_end)
      pid_i = cell_particle_id[i];

    if (pid_i < nt) {
      atom_i = fetch_pos(pid_i,pos); //pos[pid_i];
    }
    if (pid_i < inum) {
      stride=inum;
      neigh_counts=nbor_list+stride+pid_i;
      neigh_list=neigh_counts+stride;
      nbor_list[pid_i]=pid_i;
    } else {
      stride=1;
    	neigh_counts=host_numj+pid_i-inum;
      neigh_list=host_nbor_list+(pid_i-inum)*neigh_bin_size;
    }
    
    // loop through neighbors

    for (int nborz = nborz0; nborz <= nborz1; nborz++) {
      for (int nbory = nbory0; nbory <= nbory1; nbory++) {
        for (int nborx = nborx0; nborx <= nborx1; nborx++) {
	
          int jcell = nborx + nbory*ncellx + nborz*ncellx*ncelly;
		
          int jcell_begin = cell_counts[jcell];
          int jcell_end = cell_counts[jcell+1];
          int num_atom_cell = jcell_end - jcell_begin;
	  
          // load jcell to shared memory
          int num_iter = (int)ceil((numtyp)num_atom_cell/BLOCK_NBOR_BUILD);

          for (int k = 0; k < num_iter; k++) {
            int end_idx = min(BLOCK_NBOR_BUILD, num_atom_cell-k*BLOCK_NBOR_BUILD);
	    
            if (tid < end_idx) {
              pid_j =  cell_particle_id[tid+k*BLOCK_NBOR_BUILD+jcell_begin];
              cell_list_sh[tid] = pid_j;
              atom_j = fetch_pos(pid_j,pos); //[pid_j];
              pos_sh[tid].x = atom_j.x;
              pos_sh[tid].y = atom_j.y;
              pos_sh[tid].z = atom_j.z;
            }
            __syncthreads();
	    
            if (pid_i < nt) {
	    
              for (int j = 0; j < end_idx; j++) {
                int pid_j = cell_list_sh[j]; // gather from shared memory
                diff.x = atom_i.x - pos_sh[j].x;
                diff.y = atom_i.y - pos_sh[j].y;
                diff.z = atom_i.z - pos_sh[j].z;
		
                r2 = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
                if (r2 < cell_size*cell_size && r2 > 1e-5) {
                  if (cnt < neigh_bin_size) {
                    *neigh_list = pid_j;
                    neigh_list+=stride;
                  }
                  cnt++;
                }		
              }
            }
	          __syncthreads();
	        } // for (k)
        }
      }
    }
    if (pid_i < nt)
      *neigh_counts = cnt;
  } // for (i)
}

__kernel void kernel_special(__global int *dev_nbor, 
                             __global int *host_nbor_list, 
                             __global int *host_numj, __global int *tag,
                             __global int *nspecial, __global int *special,
                             int inum, int nt, int nall, int max_nbors) {
  // ii indexes the two interacting particles in gi
  int ii=GLOBAL_ID_X;

  if (ii<nt) {
    int stride;
    __global int *list, *list_end;
    
    int n1=nspecial[ii*3];
    int n2=nspecial[ii*3+1];
    int n3=nspecial[ii*3+2];

    int numj;
    if (ii < inum) {
      stride=inum;
      list=dev_nbor+stride+ii;
      numj=*list;
      list+=stride;
    } else {
      stride=1;
      list=host_nbor_list+(ii-inum)*max_nbors;
      numj=host_numj[ii-inum];
    }
    list_end=list+numj*stride;
  
    for ( ; list<list_end; list+=stride) {
      int nbor=*list;
      int jtag=tag[nbor];

      int offset=ii;
      for (int i=0; i<n3; i++) {
        if (special[offset]==jtag) {
          int which = 1;
          if (i>=n1)
            which++;
          if (i>=n2)
            which++;
          nbor=nbor ^ (which << SBBITS);
          *list=nbor;
        }
        offset+=nt;
      }
    }
  } // if ii
}
