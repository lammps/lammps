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

#if defined(__APPLE__)
#if _GLIBCXX_ATOMIC_BUILTINS == 1
#undef _GLIBCXX_ATOMIC_BUILTINS
#endif // _GLIBCXX_ATOMIC_BUILTINS
#endif // __APPLE__

#include <assert.h>
#include "lj_gpu_memory.h"
#include "pair_gpu_cell.h"

static __constant__ float d_boxlo[3];
static __constant__ float d_boxhi[3];
static __constant__ float d_cell_size[1];
static __constant__ float d_skin[1];

void init_cell_list_const(double cell_size, double skin,
			  double *boxlo, double *boxhi)
{
  float cell_size1 = cell_size;
  float skin1 = skin;
  float boxlo1[3], boxhi1[3];
  for (int i = 0; i < 3; i++) {
    boxlo1[i] = boxlo[i];
    boxhi1[i] = boxhi[i];
  }

  cudaMemcpyToSymbol(d_cell_size, &cell_size1, sizeof(float),   
		     0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_boxlo,     boxlo1,      3*sizeof(float), 
		     0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_boxhi,     boxhi1,      3*sizeof(float), 
		     0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_skin,      &skin1,       sizeof(float),   
		     0, cudaMemcpyHostToDevice); 
}

__global__ void kernel_set_cell_list(unsigned int *cell_idx)
{
  unsigned int gid = threadIdx.x + blockIdx.x*blockDim.x;
  cell_idx[gid] = BIG_NUMBER;
}

// build the cell list
__global__ void kernel_build_cell_list(float3 *cell_list, 
				       unsigned int *cell_idx, 
				       int *cell_type, 
				       int *cell_atom,
				       float3 *pos, 
				       int *type, 
				       const int inum, 
				       const int nall,
				       const int cell_size)
{
  unsigned int gid = threadIdx.x + blockIdx.x*blockDim.x;
  float cSize = d_cell_size[0];
  int ncellx = ceil(((d_boxhi[0] - d_boxlo[0]) + 2.0f*cSize) / cSize);
  int ncelly = ceil(((d_boxhi[1] - d_boxlo[1]) + 2.0f*cSize) / cSize);
  int ncellz = ceil(((d_boxhi[2] - d_boxlo[2]) + 2.0f*cSize) / cSize);

  if (gid < nall) {
    float3 p = pos[gid];
    p.x = fmaxf(p.x, d_boxlo[0]-cSize);
    p.x = fminf(p.x, d_boxhi[0]+cSize);
    p.y = fmaxf(p.y, d_boxlo[1]-cSize);
    p.y = fminf(p.y, d_boxhi[1]+cSize);
    p.z = fmaxf(p.z, d_boxlo[2]-cSize);
    p.z = fminf(p.z, d_boxhi[2]+cSize);

    int cell_id = (int)(p.x/cSize + 1.0) + (int)(p.y/cSize + 1.0) * ncellx
    		    + (int)(p.z/cSize + 1.0) * ncellx * ncelly;

    int atom_pos = atomicAdd(&cell_atom[cell_id], 1);
    int pid = cell_id*cell_size + atom_pos;

    cell_list[pid] = pos[gid];
    cell_type[pid] = type[gid];
    cell_idx [pid] = gid;
    
  }
}

__global__ void kernel_test_rebuild(float3 *cell_list, int *cell_atom, int *rebuild)
{

  float cSize = d_cell_size[0];
  int ncellx = ceil(((d_boxhi[0] - d_boxlo[0]) + 2.0f*cSize) / cSize);
  int ncelly = ceil(((d_boxhi[1] - d_boxlo[1]) + 2.0f*cSize) / cSize);
  int ncellz = ceil(((d_boxhi[2] - d_boxlo[2]) + 2.0f*cSize) / cSize);

  // calculate 3D block idx from 2d block
  int bx = blockIdx.x;
  int by = blockIdx.y % ncelly;
  int bz = blockIdx.y / ncelly;

  int tid = threadIdx.x;

  // compute cell idx from 3D block idx
  int cid = bx + INT_MUL(by, ncellx) + INT_MUL(bz, INT_MUL(ncellx,ncelly));
  int pbase = INT_MUL(cid,blockDim.x); // atom position id in cell list

  float skin = d_skin[0];
  float lowx = d_boxlo[0] + (bx-1)*cSize - 0.5*skin;
  float hix  = lowx + cSize + skin;
  float lowy = d_boxlo[1] + (by-1)*cSize - 0.5*skin;
  float hiy  = lowy + cSize + skin;
  float lowz = d_boxlo[2] + (bz-1)*cSize - 0.5*skin;
  float hiz  = lowz + cSize + skin;

  for (int i = tid; i < cell_atom[cid]; i += blockDim.x) {
    int pid = pbase + i;
    float3 p = cell_list[pid];
    p.x = fmaxf(p.x, d_boxlo[0]-cSize);
    p.x = fminf(p.x, d_boxhi[0]+cSize);
    p.y = fmaxf(p.y, d_boxlo[1]-cSize);
    p.y = fminf(p.y, d_boxhi[1]+cSize);
    p.z = fmaxf(p.z, d_boxlo[2]-cSize);
    p.z = fminf(p.z, d_boxhi[2]+cSize);

    if (p.x < lowx || p.x > hix || p.y < lowy || p.y > hiy || p.z < lowz || p.z > hiz) {
      *rebuild = 1;   
    }
  }

}


__global__ void kernel_test_overflow(int *cell_atom, int *overflow, const int ncell)
{
  unsigned int gid = threadIdx.x + blockIdx.x*blockDim.x;

  if (gid < ncell) {
    if (cell_atom[gid] > blockDim.x) 
      *overflow = 1;
  }
}

__global__ void kernel_copy_list(float3 *cell_list, unsigned int *cell_idx, int *cell_atom, float3 *pos)
{
  float cSize = d_cell_size[0];
  int ncellx = ceil(((d_boxhi[0] - d_boxlo[0]) + 2.0f*cSize) / cSize);
  int ncelly = ceil(((d_boxhi[1] - d_boxlo[1]) + 2.0f*cSize) / cSize);
  int ncellz = ceil(((d_boxhi[2] - d_boxlo[2]) + 2.0f*cSize) / cSize);

  // calculate 3D block idx from 2d block
  int bx = blockIdx.x;
  int by = blockIdx.y % ncelly;
  int bz = blockIdx.y / ncelly;

  int tid = threadIdx.x;

  // compute cell idx from 3D block idx
  int cid = bx + INT_MUL(by, ncellx) + INT_MUL(bz, INT_MUL(ncellx,ncelly));
  int pbase = INT_MUL(cid,blockDim.x); // atom position id in cell list

  for (int i = tid; i < cell_atom[cid]; i += blockDim.x) {
    int pid = pbase + i;
    cell_list[pid] = pos[cell_idx[pid]];
  }

}


__global__ void radixSortBlocks(unsigned int *keys, float3 *values1, int *values2, unsigned int nbits, unsigned int startbit); 



#ifdef __DEVICE_EMULATION__
#define __SYNC __syncthreads();
#else
#define __SYNC 
#endif


#define WARP_SIZE 32

template<class T, int maxlevel>
__device__ T scanwarp(T val, T* sData)
{
    // The following is the same as 2 * RadixSort::WARP_SIZE * warpId + threadInWarp = 
    // 64*(threadIdx.x >> 5) + (threadIdx.x & (RadixSort::WARP_SIZE - 1))
    int idx = 2 * threadIdx.x - (threadIdx.x & (WARP_SIZE - 1));
    sData[idx] = 0;
    idx += WARP_SIZE;
    sData[idx] = val;          __SYNC

#ifdef __DEVICE_EMULATION__
	T t = sData[idx -  1]; __SYNC 
        sData[idx] += t;       __SYNC
        t = sData[idx -  2];   __SYNC 
        sData[idx] += t;       __SYNC
        t = sData[idx -  4];   __SYNC 
        sData[idx] += t;       __SYNC
        t = sData[idx -  8];   __SYNC 
        sData[idx] += t;       __SYNC
        t = sData[idx - 16];   __SYNC 
        sData[idx] += t;       __SYNC
#else
        if (0 <= maxlevel) { sData[idx] += sData[idx - 1]; } __SYNC
        if (1 <= maxlevel) { sData[idx] += sData[idx - 2]; } __SYNC
        if (2 <= maxlevel) { sData[idx] += sData[idx - 4]; } __SYNC
        if (3 <= maxlevel) { sData[idx] += sData[idx - 8]; } __SYNC
        if (4 <= maxlevel) { sData[idx] += sData[idx -16]; } __SYNC
#endif

        return sData[idx] - val;  // convert inclusive -> exclusive
}

__device__ unsigned int scan(unsigned int idata)
{    
    extern  __shared__  unsigned int ptr[];
    
    unsigned int idx = threadIdx.x;
    
    unsigned int val = idata;
    
    val = scanwarp<unsigned int, 4>(val, ptr);
    __syncthreads();

    if ((idx & (WARP_SIZE - 1)) == WARP_SIZE - 1)
    {
        ptr[idx >> 5] = val + idata;
    }
    __syncthreads();

#ifndef __DEVICE_EMULATION__
    if (idx < WARP_SIZE)
#endif
    {
        ptr[idx] = scanwarp<unsigned int, 2>(ptr[idx], ptr);
    }
    __syncthreads();

    val += ptr[idx >> 5];

    return val;
}


__device__ unsigned int rank(unsigned int preds)
{
    unsigned int address = scan(preds);  

    __shared__ unsigned int numtrue;
    if (threadIdx.x == blockDim.x - 1)
    {
        numtrue = address + preds;
    }
    __syncthreads();

    unsigned int rank;
    unsigned int idx = threadIdx.x;
    rank = (preds) ? address : numtrue + idx - address;

    return rank;
}

template<int blockSize>
__device__ void radixSortBlock(unsigned int *key, float3 *value1, int *value2, unsigned int nbits, unsigned int startbit)
{
  extern __shared__ unsigned int sMem1[];
  __shared__ float sMem2[blockSize];
  __shared__ int sMem3[blockSize];

  int tid = threadIdx.x;

  for(unsigned int shift = startbit; shift < (startbit + nbits); ++shift) {
    unsigned int lsb;
    lsb = !(((*key) >> shift) & 0x1);

    unsigned int r;
		
    r = rank(lsb);

    // This arithmetic strides the ranks across 4 CTA_SIZE regions
    sMem1[r] = *key;
    __syncthreads();

    // The above allows us to read without 4-way bank conflicts:
    *key = sMem1[tid];    
    __syncthreads();

    sMem2[r] = (*value1).x;
    __syncthreads();
    (*value1).x = sMem2[tid];
    __syncthreads();

    sMem2[r] = (*value1).y;
    __syncthreads();
    (*value1).y = sMem2[tid];
    __syncthreads();

    sMem2[r] = (*value1).z;
    __syncthreads();
    (*value1).z = sMem2[tid];
    __syncthreads();

    sMem3[r] = *value2;
    __syncthreads();
    *value2 = sMem3[tid];
    __syncthreads();

  }

}

__global__ void radixSortBlocks(unsigned int *keys, 
				float3 *values1, 
				int *values2, 
				unsigned int nbits, 
				unsigned int startbit)
{

  extern __shared__ unsigned int sMem[];

  int gid = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int key;
  float3 value1;
  int value2;
  key = keys[gid];
  value1 = values1[gid];
  value2 = values2[gid];
  __syncthreads();

  if (blockDim.x == 64) 
    radixSortBlock<64>(&key, &value1, &value2, nbits, startbit);
  else if (blockDim.x == 128) 
    radixSortBlock<128>(&key, &value1, &value2, nbits, startbit);
  else if (blockDim.x == 256)
    radixSortBlock<256>(&key, &value1, &value2, nbits, startbit);

  keys[gid] = key;
  values1[gid] = value1;
  values2[gid] = value2;
}

void sortBlocks(unsigned int *keys, float3 *values1, int *values2, const int size, int cell_size)
{
  int i = 0;
  const unsigned int bitSize = sizeof(unsigned int)*8;
  const unsigned int bitStep = 4;
  const int gSize = size/cell_size;
  while (bitSize > i*bitStep) {
    radixSortBlocks<<<gSize, cell_size, 2*cell_size*sizeof(unsigned int)>>>(keys, values1, values2, bitStep, i*bitStep);
    i++;
  }
}

static float3 *d_pos, *pos_temp;
static int *d_type;
static int *d_overflow, *d_rebuild;

void init_cell_list(cell_list &cell_list_gpu, 
		   const int nall,
		   const int ncell, 
		   const int buffer)
{
  cudaMalloc((void**)&(cell_list_gpu.pos), ncell*buffer*sizeof(float3));
  cudaMalloc((void**)&(cell_list_gpu.idx),  ncell*buffer*sizeof(unsigned int));
  cudaMalloc((void**)&(cell_list_gpu.type), ncell*buffer*sizeof(int));
  cudaMalloc((void**)&(cell_list_gpu.natom), ncell*sizeof(int));

  cudaMallocHost((void**)&pos_temp, nall*sizeof(float3));
  cudaMalloc((void**)&d_pos,       nall*sizeof(float3));
  cudaMalloc((void**)&d_type,      nall*sizeof(int));
  cudaMalloc((void**)&d_overflow, sizeof(int));
  cudaMalloc((void**)&d_rebuild, sizeof(int));

  cudaMemset(cell_list_gpu.natom, 0, ncell*sizeof(int));
  cudaMemset(cell_list_gpu.pos, 0, ncell*buffer*sizeof(float3));
}

void clear_cell_list(cell_list &cell_list_gpu)
{
  cudaFree(cell_list_gpu.pos);
  cudaFree(cell_list_gpu.idx);
  cudaFree(cell_list_gpu.natom);
  cudaFree(cell_list_gpu.type);

  cudaFreeHost(pos_temp);
  cudaFree(d_pos);
  cudaFree(d_type);
  cudaFree(d_overflow);
  cudaFree(d_rebuild);
}


void build_cell_list(double *atom_pos, int *atom_type, 
		     cell_list &cell_list_gpu, 
		     const int ncell, const int ncellx, const int ncelly, const int ncellz, 
		     const int buffer, const int inum, const int nall, const int ago)
{

  cudaError_t err;				     

  cudaMemset(d_overflow, 0, sizeof(int));
  cudaMemset(d_rebuild, 0, sizeof(int));

  // copy position and type to GPU
  for (int i = 0; i < 3*nall; i+=3) { 
    pos_temp[i/3] = make_float3(atom_pos[i], atom_pos[i+1], atom_pos[i+2]);
  }
  cudaMemcpy(d_pos, pos_temp, nall*sizeof(float3), cudaMemcpyHostToDevice);
  cudaMemcpy(d_type, atom_type, nall*sizeof(int),  cudaMemcpyHostToDevice);

  static int first_build = 1;
  int rebuild = 0;

  // copy the last built cell-list and test whether it needs to be rebuilt
  if (!first_build) {
    
    dim3 grid(ncellx, ncelly*ncellz);
    kernel_copy_list<<<grid, buffer>>>(cell_list_gpu.pos, 
				 cell_list_gpu.idx, 
				 cell_list_gpu.natom, d_pos);
    cudaMemset(d_rebuild, 0, sizeof(int));
    kernel_test_rebuild<<<grid, buffer>>>(cell_list_gpu.pos, 
					 cell_list_gpu.natom,
					 d_rebuild);
    cudaMemcpy(&rebuild, d_rebuild, sizeof(int), cudaMemcpyDeviceToHost);
    
    err = cudaGetLastError();
    assert(err == cudaSuccess);
  }

  if (ago == 0) rebuild = 1;
  
  // build cell-list for the first time
  if (first_build || rebuild) {
    first_build = 0;
    // cout << "Building cell list..." << endl;
    cudaMemset(cell_list_gpu.natom, 0, ncell*sizeof(int));
    // initialize d_cell_idx for the sorting routine
    kernel_set_cell_list<<<ncell, buffer>>>(cell_list_gpu.idx);
    
    // build cell list
    dim3 blockDim(128);
    dim3 gridDim(static_cast<int>(ceil(static_cast<double>(nall)/blockDim.x)));
    kernel_build_cell_list<<<gridDim, blockDim>>>(cell_list_gpu.pos, 
						  cell_list_gpu.idx, 
						  cell_list_gpu.type, 
						  cell_list_gpu.natom, 
						  d_pos, d_type, inum, nall, buffer);
    err = cudaGetLastError();
    assert(err == cudaSuccess);
    // check cell list overflow
    int overflow = 0;
    int gDimCell = static_cast<int>(ceil(static_cast<double>(ncell)/buffer));
    kernel_test_overflow<<<gDimCell, buffer>>>(cell_list_gpu.natom, 
					       d_overflow, ncell);
    cudaMemcpy(&overflow, d_overflow, sizeof(int), cudaMemcpyDeviceToHost);
     
    if (overflow > 0) {
      printf("\n BLOCK_1D too small for cell list, please increase it!");
      printf("\n BLOCK_1D = %d",BLOCK_1D);
      printf("\n ncell = %d",ncell);
      printf("\n gDimCell = %d",gDimCell);
      printf("\n overflow = %d \n",overflow);
      exit(0);
    }
    
    // sort atoms in every cell by atom index to avoid floating point associativity problem.
    sortBlocks(cell_list_gpu.idx, cell_list_gpu.pos, 
	       cell_list_gpu.type, ncell*buffer, buffer);

    cudaThreadSynchronize();
    err = cudaGetLastError();
    assert(err == cudaSuccess);
  }

}
