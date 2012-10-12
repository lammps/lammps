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

extern __shared__ ENERGY_FLOAT sharedmem[];


__global__ void Cuda_ComputeTempPartialCuda_Scalar_Kernel(int groupbit, int xflag, int yflag, int zflag)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  sharedmem[threadIdx.x] = 0;

  if(i < _nlocal) {
    if(_rmass_flag) {
      if(_mask[i] & groupbit)
        sharedmem[threadIdx.x] = (_v[i] * _v[i] * xflag + _v[i + _nmax] * _v[i + _nmax] * yflag + _v[i + 2 * _nmax] * _v[i + 2 * _nmax] * zflag) * _rmass[i];
    } else {
      if(_mask[i] & groupbit)
        sharedmem[threadIdx.x] = (_v[i] * _v[i] * xflag + _v[i + _nmax] * _v[i + _nmax] * yflag + _v[i + 2 * _nmax] * _v[i + 2 * _nmax] * zflag) * (_mass[_type[i]]);
    }
  }

  reduceBlock(sharedmem);
  ENERGY_FLOAT* buffer = (ENERGY_FLOAT*) _buffer;

  if(threadIdx.x == 0) {
    buffer[blockIdx.x * gridDim.y + blockIdx.y] = sharedmem[0];
  }
}

__global__ void Cuda_ComputeTempPartialCuda_Vector_Kernel(int groupbit, int xflag, int yflag, int zflag)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  sharedmem[threadIdx.x] = 0;
  sharedmem[threadIdx.x + blockDim.x] = 0;
  sharedmem[threadIdx.x + 2 * blockDim.x] = 0;
  sharedmem[threadIdx.x + 3 * blockDim.x] = 0;
  sharedmem[threadIdx.x + 4 * blockDim.x] = 0;
  sharedmem[threadIdx.x + 5 * blockDim.x] = 0;

  if(i < _nlocal)
    if(_mask[i] & groupbit) {
      V_FLOAT massone;

      if(_rmass_flag) massone = _rmass[i];
      else massone = _mass[_type[i]];

      sharedmem[threadIdx.x] = massone * _v[i] * _v[i] * xflag;
      sharedmem[threadIdx.x + blockDim.x] = massone * _v[i + _nmax] * _v[i + _nmax] * yflag;
      sharedmem[threadIdx.x + 2 * blockDim.x] = massone * _v[i + 2 * _nmax] * _v[i + 2 * _nmax] * zflag;
      sharedmem[threadIdx.x + 3 * blockDim.x] = massone * _v[i] * _v[i + _nmax] * xflag * yflag;
      sharedmem[threadIdx.x + 4 * blockDim.x] = massone * _v[i] * _v[i + 2 * _nmax] * xflag * zflag;
      sharedmem[threadIdx.x + 5 * blockDim.x] = massone * _v[i + _nmax] * _v[i + 2 * _nmax] * yflag * zflag;
    }

  reduceBlock(sharedmem);
  reduceBlock(&sharedmem[blockDim.x]);
  reduceBlock(&sharedmem[2 * blockDim.x]);
  reduceBlock(&sharedmem[3 * blockDim.x]);
  reduceBlock(&sharedmem[4 * blockDim.x]);
  reduceBlock(&sharedmem[5 * blockDim.x]);
  ENERGY_FLOAT* buffer = (ENERGY_FLOAT*) _buffer;

  if(threadIdx.x == 0) {
    buffer[blockIdx.x * gridDim.y + blockIdx.y] = sharedmem[0];
    buffer[blockIdx.x * gridDim.y + blockIdx.y + gridDim.x * gridDim.y] = sharedmem[blockDim.x];
    buffer[blockIdx.x * gridDim.y + blockIdx.y + 2 * gridDim.x * gridDim.y] = sharedmem[2 * blockDim.x];
    buffer[blockIdx.x * gridDim.y + blockIdx.y + 3 * gridDim.x * gridDim.y] = sharedmem[3 * blockDim.x];
    buffer[blockIdx.x * gridDim.y + blockIdx.y + 4 * gridDim.x * gridDim.y] = sharedmem[4 * blockDim.x];
    buffer[blockIdx.x * gridDim.y + blockIdx.y + 5 * gridDim.x * gridDim.y] = sharedmem[5 * blockDim.x];
  }
}


__global__ void Cuda_ComputeTempPartialCuda_Reduce_Kernel(int n, ENERGY_FLOAT* t)
{
  int i = 0;
  sharedmem[threadIdx.x] = 0;
  ENERGY_FLOAT myforig = 0.0;
  ENERGY_FLOAT* buf = (ENERGY_FLOAT*) _buffer;
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
    t[blockIdx.x] = myforig;
}

__global__ void Cuda_ComputeTempPartialCuda_RemoveBiasAll_Kernel(int groupbit, int xflag, int yflag, int zflag, V_FLOAT* vbiasall)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < _nlocal)
    if(_mask[i] & groupbit) {
      if(!xflag) {
        vbiasall[i] = _v[i];
        _v[i] = V_F(0.0);
      }

      if(!yflag) {
        vbiasall[i + _nmax] = _v[i + _nmax];
        _v[i + _nmax] = V_F(0.0);
      }

      if(!zflag) {
        vbiasall[i + 2 * _nmax] = _v[i + 2 * _nmax];
        _v[i + 2 * _nmax] = V_F(0.0);
      }
    }
}

__global__ void Cuda_ComputeTempPartialCuda_RestoreBiasAll_Kernel(int groupbit, int xflag, int yflag, int zflag, V_FLOAT* vbiasall)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < _nlocal)
    if(_mask[i] & groupbit) {
      if(!xflag) {
        _v[i] += vbiasall[i];
      }

      if(!yflag) {
        _v[i + _nmax] += vbiasall[i + _nmax];
      }

      if(!zflag) {
        _v[i + 2 * _nmax] += vbiasall[i + 2 * _nmax];
      }
    }
}
