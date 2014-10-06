/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

extern __shared__ ENERGY_CFLOAT sharedmem[];

static inline __device__ void PairVirialCompute_A_Kernel(int eflag, int vflag, int coulflag = 0)
{
  __syncthreads();
  ENERGY_CFLOAT* shared = sharedmem;

  if(eflag) {
    reduceBlock(shared);
    shared += blockDim.x;

    if(coulflag) {
      reduceBlock(shared);
      shared += blockDim.x;
    }
  }

  if(vflag) {
    reduceBlock(shared + 0 * blockDim.x);
    reduceBlock(shared + 1 * blockDim.x);
    reduceBlock(shared + 2 * blockDim.x);
    reduceBlock(shared + 3 * blockDim.x);
    reduceBlock(shared + 4 * blockDim.x);
    reduceBlock(shared + 5 * blockDim.x);
  }

  if(threadIdx.x == 0) {
    shared = sharedmem;
    ENERGY_CFLOAT* buffer = (ENERGY_CFLOAT*) _buffer;

    if(eflag) {
      buffer[blockIdx.x * gridDim.y + blockIdx.y] = ENERGY_F(0.5) * shared[0];
      shared += blockDim.x;
      buffer += gridDim.x * gridDim.y;

      if(coulflag) {
        buffer[blockIdx.x * gridDim.y + blockIdx.y] = ENERGY_F(0.5) * shared[0];
        shared += blockDim.x;
        buffer += gridDim.x * gridDim.y;
      }
    }

    if(vflag) {
      buffer[blockIdx.x * gridDim.y + blockIdx.y + 0 * gridDim.x * gridDim.y] = ENERGY_F(0.5) * shared[0 * blockDim.x];
      buffer[blockIdx.x * gridDim.y + blockIdx.y + 1 * gridDim.x * gridDim.y] = ENERGY_F(0.5) * shared[1 * blockDim.x];
      buffer[blockIdx.x * gridDim.y + blockIdx.y + 2 * gridDim.x * gridDim.y] = ENERGY_F(0.5) * shared[2 * blockDim.x];
      buffer[blockIdx.x * gridDim.y + blockIdx.y + 3 * gridDim.x * gridDim.y] = ENERGY_F(0.5) * shared[3 * blockDim.x];
      buffer[blockIdx.x * gridDim.y + blockIdx.y + 4 * gridDim.x * gridDim.y] = ENERGY_F(0.5) * shared[4 * blockDim.x];
      buffer[blockIdx.x * gridDim.y + blockIdx.y + 5 * gridDim.x * gridDim.y] = ENERGY_F(0.5) * shared[5 * blockDim.x];
    }
  }

  __syncthreads();
}

__global__ void MY_AP(PairVirialCompute_reduce)(int n)
{
  sharedmem[threadIdx.x] = ENERGY_F(0.0);
  ENERGY_CFLOAT sum = ENERGY_F(0.0);
  ENERGY_CFLOAT* buf = (ENERGY_CFLOAT*) _buffer;
  buf = &buf[blockIdx.x * n];
  //if(blockIdx.x==2) buf=&buf[n];

  for(int i = 0; i < n; i += blockDim.x) {
    sharedmem[threadIdx.x] = (i + threadIdx.x < n) ? buf[i + threadIdx.x] : ENERGY_F(0.0);
    __syncthreads();
    reduceBlock(sharedmem);

    if(threadIdx.x == 0) sum += sharedmem[0];
  }

  if(threadIdx.x == 0) {
    if(gridDim.x == 1) { //evdwl
      _eng_vdwl[0] += sum;
    }

    if(gridDim.x == 2) { //evdwl + ecoul only
      if(blockIdx.x == 0)
        _eng_vdwl[0] += sum;
      else
        _eng_coul[0] += sum;
    }

    if(gridDim.x == 6) { //virial
      _virial[blockIdx.x] += sum;
    }

    if(gridDim.x == 7) { //evdwl+virial
      if(blockIdx.x == 0)
        _eng_vdwl[0] += sum;
      else _virial[blockIdx.x - 1] += sum;
    }

    if(gridDim.x == 8) { //evdwl+ecoul+virial
      if(blockIdx.x == 0)
        _eng_vdwl[0] += sum;
      else if(blockIdx.x == 1)
        _eng_coul[0] += sum;
      else
        _virial[blockIdx.x - 2] += sum;
    }
  }
}
