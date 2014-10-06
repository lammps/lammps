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

__global__ void Cuda_CommCuda_PackComm_Kernel(int* sendlist, int n, int maxlistlength, int iswap, X_CFLOAT dx, X_CFLOAT dy, X_CFLOAT dz, void* buffer)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  int* list = sendlist + iswap * maxlistlength;

  if(i < n) {
    int j = list[i];

    if(j > _nmax) _flag[0] = 1;

    ((X_CFLOAT*) buffer)[i] = _x[j] + dx;
    ((X_CFLOAT*) buffer)[i + 1 * n] = _x[j + _nmax] + dy;
    ((X_CFLOAT*) buffer)[i + 2 * n] = _x[j + 2 * _nmax] + dz;
  }
}

__global__ void Cuda_CommCuda_PackCommVel_Kernel(int* sendlist, int n, int maxlistlength, int iswap, X_CFLOAT dx, X_CFLOAT dy, X_CFLOAT dz, void* buffer)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  int* list = sendlist + iswap * maxlistlength;

  if(i < n) {
    int j = list[i];

    if(j > _nmax) _flag[0] = 1;

    ((X_CFLOAT*) buffer)[i] = _x[j] + dx;
    ((X_CFLOAT*) buffer)[i + 1 * n] = _x[j + _nmax] + dy;
    ((X_CFLOAT*) buffer)[i + 2 * n] = _x[j + 2 * _nmax] + dz;
    ((X_CFLOAT*) buffer)[i + 3 * n] = _v[j];
    ((X_CFLOAT*) buffer)[i + 4 * n] = _v[j + _nmax];
    ((X_CFLOAT*) buffer)[i + 5 * n] = _v[j + 2 * _nmax];
  }
}

__global__ void Cuda_CommCuda_PackComm_Self_Kernel(int* sendlist, int n, int maxlistlength, int iswap, X_CFLOAT dx, X_CFLOAT dy, X_CFLOAT dz, int first)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  int* list = sendlist + iswap * maxlistlength;

  if(i < n) {
    int j = i;
    j = list[i];

    _x[i + first] = _x[j] + dx;
    _x[i + first + _nmax] = _x[j + _nmax] + dy;
    _x[i + first + 2 * _nmax] = _x[j + 2 * _nmax] + dz;
  }
}

__global__ void Cuda_CommCuda_PackCommVel_Self_Kernel(int* sendlist, int n, int maxlistlength, int iswap, X_CFLOAT dx, X_CFLOAT dy, X_CFLOAT dz, int first)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  int* list = sendlist + iswap * maxlistlength;

  if(i < n) {
    int j = i;
    j = list[i];

    _x[i + first] = _x[j] + dx;
    _x[i + first + _nmax] = _x[j + _nmax] + dy;
    _x[i + first + 2 * _nmax] = _x[j + 2 * _nmax] + dz;
    _v[i + first] = _v[j];
    _v[i + first + _nmax] = _v[j + _nmax];
    _v[i + first + 2 * _nmax] = _v[j + 2 * _nmax];
  }
}

__global__ void Cuda_CommCuda_UnpackComm_Kernel(int n, int first, void* buffer)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < n) {
    _x[i + first] = ((X_CFLOAT*) buffer)[i];
    _x[i + first + _nmax] = ((X_CFLOAT*) buffer)[i + 1 * n];
    _x[i + first + 2 * _nmax] = ((X_CFLOAT*) buffer)[i + 2 * n];
  }
}


__global__ void Cuda_CommCuda_UnpackCommVel_Kernel(int n, int first, void* buffer)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < n) {
    _x[i + first] = ((X_CFLOAT*) buffer)[i];
    _x[i + first + _nmax] = ((X_CFLOAT*) buffer)[i + 1 * n];
    _x[i + first + 2 * _nmax] = ((X_CFLOAT*) buffer)[i + 2 * n];
    _v[i + first] = ((X_CFLOAT*) buffer)[i + 3 * n];
    _v[i + first + _nmax] = ((X_CFLOAT*) buffer)[i + 4 * n];
    _v[i + first + 2 * _nmax] = ((X_CFLOAT*) buffer)[i + 5 * n];
  }
}

__global__ void Cuda_CommCuda_PackReverse_Kernel(int n, int first)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < n) {
    ((F_CFLOAT*) _buffer)[i] = _f[i + first];
    ((F_CFLOAT*) _buffer)[i + n] = _f[i + first + _nmax];
    ((F_CFLOAT*) _buffer)[i + 2 * n] = _f[i + first + 2 * _nmax];
  }

}

__global__ void Cuda_CommCuda_UnpackReverse_Kernel(int* sendlist, int n, int maxlistlength, int iswap)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  int* list = sendlist + iswap * maxlistlength;

  if(i < n) {
    int j = list[i];
    _f[j] += ((F_CFLOAT*)_buffer)[i];
    _f[j + _nmax] += ((F_CFLOAT*) _buffer)[i + n];
    _f[j + 2 * _nmax] += ((F_CFLOAT*) _buffer)[i + 2 * n];
  }

}

__global__ void Cuda_CommCuda_UnpackReverse_Self_Kernel(int* sendlist, int n, int maxlistlength, int iswap, int first)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  int* list = sendlist + iswap * maxlistlength;

  if(i < n) {
    int j = list[i];

    _f[j] += _f[i + first];
    _f[j + _nmax] += _f[i + first + _nmax];
    _f[j + 2 * _nmax] += _f[i + first + 2 * _nmax];
  }

}

extern __shared__ int shared[];

__global__ void Cuda_CommCuda_BuildSendlist_Single(int bordergroup, int ineed, int atom_nfirst,
    int nfirst, int nlast, int dim, int iswap, X_CFLOAT* slablo, X_CFLOAT* slabhi, int* sendlist, int maxlistlength)
{
  int* list = sendlist + iswap * maxlistlength;
  X_CFLOAT lo = slablo[iswap];
  X_CFLOAT hi = slabhi[iswap];
  bool add = false;

  if(!bordergroup || ineed >= 2) {
    int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x + nfirst;

    if(i < nlast)
      if(_x[i + dim * _nmax] >= lo && _x[i + dim * _nmax] <= hi) {
        add = true;
      }

    shared[threadIdx.x] = add ? 1 : 0;

    __syncthreads();

    int nsend = 0;

    if(threadIdx.x == 0) {
      for(int k = 0; k < blockDim.x; k++) {
        if(shared[k]) {
          nsend++;
          shared[k] = nsend;
        }
      }

      shared[blockDim.x] = atomicAdd((int*) _buffer, nsend);
    }

    __syncthreads();

    nsend = shared[blockDim.x] + shared[threadIdx.x] - 1;

    if(add && nsend < maxlistlength)
      list[nsend] = i;


  } else {

    int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

    if(i < atom_nfirst)
      if(_x[i + dim * _nmax] >= lo && _x[i + dim * _nmax] <= hi) {
        add = true;
      }

    shared[threadIdx.x] = add ? 1 : 0;

    __syncthreads();

    int nsend = 0;

    if(threadIdx.x == 0) {
      for(int k = 0; k < blockDim.x; k++) {
        if(shared[k]) {
          nsend++;
          shared[k] = nsend;
        }
      }

      shared[blockDim.x] = atomicAdd((int*) _buffer, nsend);
    }

    __syncthreads();

    nsend = shared[blockDim.x] + shared[threadIdx.x] - 1;

    if(add && nsend < maxlistlength)
      list[nsend] = i;

    __syncthreads();

    add = false;
    i += _nlocal;

    if(i < nlast)
      if(_x[i + dim * _nmax] >= lo && _x[i + dim * _nmax] <= hi) {
        add = true;
      }

    shared[threadIdx.x] = add ? 1 : 0;

    __syncthreads();

    nsend = 0;

    if(threadIdx.x == 0) {
      for(int k = 0; k < blockDim.x; k++) {
        if(shared[k]) {
          nsend++;
          shared[k] = nsend;
        }
      }

      shared[blockDim.x] = atomicAdd((int*) _buffer, nsend);
    }

    __syncthreads();

    nsend = shared[blockDim.x] + shared[threadIdx.x] - 1;

    if(add && nsend < maxlistlength)
      list[nsend] = i;

  }
}


__global__ void Cuda_CommCuda_BuildSendlist_Multi(int bordergroup, int ineed, int atom_nfirst
    , int nfirst, int nlast, int dim, int iswap, X_CFLOAT* multilo, X_CFLOAT* multihi, int* sendlist, int maxlistlength)
{
  int* list = sendlist + iswap * maxlistlength;
  X_CFLOAT* mlo = &multilo[iswap * _cuda_ntypes];
  X_CFLOAT* mhi = &multihi[iswap * _cuda_ntypes];
  int itype = 0;
  bool add = false;

  if(!bordergroup || ineed >= 2) {
    int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x + nfirst;

    if(i < nlast) {
      itype = _type[i];

      if(_x[i + dim * _nmax] >= mlo[itype] && _x[i + dim * _nmax] <= mhi[itype]) {
        add = true;
      }
    }

    shared[threadIdx.x] = add ? 1 : 0;

    __syncthreads();

    int nsend = 0;

    if(threadIdx.x == 0) {
      for(int k = 0; k < blockDim.x; k++) {
        if(shared[k]) {
          nsend++;
          shared[k] = nsend;
        }
      }

      shared[blockDim.x] = atomicAdd((int*) _buffer, nsend);
    }

    __syncthreads();

    nsend = shared[blockDim.x] + shared[threadIdx.x] - 1;

    if(add && nsend < maxlistlength)
      list[nsend] = i;


  } else {

    int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

    if(i < atom_nfirst) {
      itype = _type[i];

      if(_x[i + dim * _nmax] >= mlo[itype] && _x[i + dim * _nmax] <= mhi[itype]) {
        add = true;
      }
    }

    shared[threadIdx.x] = add ? 1 : 0;

    __syncthreads();

    int nsend = 0;

    if(threadIdx.x == 0) {
      for(int k = 0; k < blockDim.x; k++) {
        if(shared[k]) {
          nsend++;
          shared[k] = nsend;
        }
      }

      shared[blockDim.x] = atomicAdd((int*) _buffer, nsend);
    }

    __syncthreads();

    nsend = shared[blockDim.x] + shared[threadIdx.x] - 1;

    if(add && nsend < maxlistlength)
      list[nsend] = i;

    __syncthreads();

    add = false;
    i += _nlocal;

    if(i < nlast) {
      itype = _type[i];

      if(_x[i + dim * _nmax] >= mlo[itype] && _x[i + dim * _nmax] <= mhi[itype]) {
        add = true;
      }
    }

    shared[threadIdx.x] = add ? 1 : 0;

    __syncthreads();

    nsend = 0;

    if(threadIdx.x == 0) {
      for(int k = 0; k < blockDim.x; k++) {
        if(shared[k]) {
          nsend++;
          shared[k] = nsend;
        }
      }

      shared[blockDim.x] = atomicAdd((int*) _buffer, nsend);
    }

    __syncthreads();

    nsend = shared[blockDim.x] + shared[threadIdx.x] - 1;

    if(add && nsend < maxlistlength)
      list[nsend] = i;

  }
}
