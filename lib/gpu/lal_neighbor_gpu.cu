// **************************************************************************
//                               neighbor_gpu.cu
//                             -------------------
//                              Peng Wang (Nvidia)
//                           W. Michael Brown (ORNL)
//
//  Device code for handling GPU generated neighbor lists
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : 
//    email                : penwang@nvidia.com, brownw@ornl.gov
// ***************************************************************************/

#ifdef NV_KERNEL
#include "lal_preprocessor.h"
texture<float4> neigh_tex;
#ifndef _DOUBLE_DOUBLE
ucl_inline float4 fetch_pos(const int& i, const float4 *pos) 
  { return tex1Dfetch(neigh_tex, i); }
#endif

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

#endif



__kernel void transpose(__global int *out, __global int *in, int columns_in, 
                        int rows_in)
{
	__local int block[BLOCK_CELL_2D][BLOCK_CELL_2D+1];
	
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

__kernel void calc_neigh_list_cell(__global numtyp4 *x_, 
                                   __global int *cell_particle_id, 
                                   __global int *cell_counts, 
                                   __global int *nbor_list,
                                   __global int *host_nbor_list, 
                                   __global int *host_numj, 
                                   int neigh_bin_size, numtyp cell_size,
                                   int ncellx, int ncelly, int ncellz,
                                   int inum, int nt, int nall, int t_per_atom)
{
  int tid = THREAD_ID_X;
  int ix = BLOCK_ID_X;
  int iy = BLOCK_ID_Y % ncelly;
  int iz = BLOCK_ID_Y / ncelly;
  int bsx = BLOCK_SIZE_X;
	  
  int icell = ix + iy*ncellx + iz*ncellx*ncelly;

  __local int cell_list_sh[BLOCK_NBOR_BUILD];
  __local numtyp4 pos_sh[BLOCK_NBOR_BUILD];

  int icell_begin = cell_counts[icell];
  int icell_end = cell_counts[icell+1];

  int nborz0 = max(iz-1,0), nborz1 = min(iz+1, ncellz-1),
      nbory0 = max(iy-1,0), nbory1 = min(iy+1, ncelly-1),
      nborx0 = max(ix-1,0), nborx1 = min(ix+1, ncellx-1);

  numtyp4 diff;
  numtyp r2;
  int cap=ucl_ceil((numtyp)(icell_end - icell_begin)/bsx);
  for (int ii = 0; ii < cap; ii++) {
    int i = icell_begin + tid + ii*bsx;
    int pid_i = nall, pid_j, stride;
    numtyp4 atom_i, atom_j;
    int cnt = 0;    
    __global int *neigh_counts, *neigh_list;
    
    if (i < icell_end)
      pid_i = cell_particle_id[i];

    if (pid_i < nt) {
      atom_i = fetch_pos(pid_i,x_); //pos[pid_i];
    }
    if (pid_i < inum) {
      stride=inum;
      neigh_counts=nbor_list+stride+pid_i;
      neigh_list=neigh_counts+stride+pid_i*(t_per_atom-1);
      stride=stride*t_per_atom-t_per_atom;
      nbor_list[pid_i]=pid_i;
    } else {
      stride=0;
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
          int num_iter = ucl_ceil((numtyp)num_atom_cell/bsx);

          for (int k = 0; k < num_iter; k++) {
            int end_idx = min(bsx, num_atom_cell-k*bsx);
	    
            if (tid < end_idx) {
              pid_j =  cell_particle_id[tid+k*bsx+jcell_begin];
              cell_list_sh[tid] = pid_j;
              atom_j = fetch_pos(pid_j,x_); //[pid_j];
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
                  cnt++;
                  if (cnt <= neigh_bin_size) {
                    *neigh_list = pid_j;
                    neigh_list++;
                    if ((cnt & (t_per_atom-1))==0)
                      neigh_list=neigh_list+stride;
                  }
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
                             int inum, int nt, int max_nbors, int t_per_atom) {
  int tid=THREAD_ID_X;
  int ii=fast_mul((int)BLOCK_ID_X,(int)(BLOCK_SIZE_X)/t_per_atom);
  ii+=tid/t_per_atom;
  int offset=tid & (t_per_atom-1);

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
      list+=stride+fast_mul(ii,t_per_atom-1);
      stride=fast_mul(inum,t_per_atom);
      int njt=numj/t_per_atom;
      list_end=list+fast_mul(njt,stride)+(numj & (t_per_atom-1));
      list+=offset;
    } else {
      stride=1;
      list=host_nbor_list+(ii-inum)*max_nbors;
      numj=host_numj[ii-inum];
      list_end=list+fast_mul(numj,stride);
    }
  
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

