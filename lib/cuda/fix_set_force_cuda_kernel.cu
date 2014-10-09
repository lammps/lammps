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

extern __shared__ F_CFLOAT sharedmem[];


__global__ void Cuda_FixSetForceCuda_PostForce_Kernel(int groupbit, F_CFLOAT xvalue, F_CFLOAT yvalue, F_CFLOAT zvalue, int flagx, int flagy, int flagz)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  sharedmem[threadIdx.x] = 0;
  sharedmem[threadIdx.x + blockDim.x] = 0;
  sharedmem[threadIdx.x + 2 * blockDim.x] = 0;

  if(i < _nlocal)
    if(_mask[i] & groupbit) {
      sharedmem[threadIdx.x] = _f[i];
      sharedmem[threadIdx.x + blockDim.x] = _f[i + 1 * _nmax];
      sharedmem[threadIdx.x + 2 * blockDim.x] = _f[i + 2 * _nmax];

      if(flagx) _f[i] = xvalue;

      if(flagy) _f[i + 1 * _nmax] = yvalue;

      if(flagz) _f[i + 2 * _nmax] = zvalue;
    }


  reduceBlock(sharedmem);
  reduceBlock(&sharedmem[blockDim.x]);
  reduceBlock(&sharedmem[2 * blockDim.x]);
  F_CFLOAT* buffer = (F_CFLOAT*)_buffer;

  if(threadIdx.x == 0) {
    buffer[blockIdx.x * gridDim.y + blockIdx.y] = sharedmem[0];
    buffer[blockIdx.x * gridDim.y + blockIdx.y + gridDim.x * gridDim.y] = sharedmem[blockDim.x];
    buffer[blockIdx.x * gridDim.y + blockIdx.y + 2 * gridDim.x * gridDim.y] = sharedmem[2 * blockDim.x];
  }
}


__global__ void Cuda_FixSetForceCuda_Reduce_FOriginal(int n, F_CFLOAT* foriginal)
{
  int i = 0;
  sharedmem[threadIdx.x] = 0;
  F_CFLOAT myforig = 0.0;
  F_CFLOAT* buf = (F_CFLOAT*)_buffer;
  buf = &buf[blockIdx.x * n];

  while(i < n) {
    sharedmem[threadIdx.x] = 0;

    if(i + threadIdx.x < n)
      sharedmem[threadIdx.x] = buf[i + threadIdx.x];

    __syncthreads();
    reduceBlock(sharedmem);
    i += blockDim.x;

    if(threadIdx.x == 0)
      myforig += sharedmem[0];
  }

  if(threadIdx.x == 0)
    foriginal[blockIdx.x] = myforig;
}

