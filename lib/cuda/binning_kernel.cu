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

// load some variables from shared cuda data into device's constant memory:
__device__ __constant__ X_FLOAT rez_bin_size[3];
__device__ __constant__ unsigned* bin_error_count;

__device__ __constant__ int cuda_dummy_type;
__device__ __constant__ unsigned binned_size_all;
__device__ __constant__ X_FLOAT outside[3];

__global__ void PreBinning_Kernel()
{
  const unsigned bin = gridDim.y * blockIdx.x + blockIdx.y;

  if(bin < gridDim.x * gridDim.y) { // TODO: suspected always to be true
    _binned_type[blockDim.x * bin + threadIdx.x] = cuda_dummy_type;

    const int i = 3 * blockDim.x * bin + threadIdx.x;
    X_FLOAT* binned_x = _binned_x + i;
    *binned_x = _subhi[0] + outside[0] * (1 + i);
    binned_x += blockDim.x;
    *binned_x = _subhi[1] + outside[1] * (1 + i);
    binned_x += blockDim.x;
    *binned_x = _subhi[2] + outside[2] * (1 + i);
    _binned_tag[i] = -1;
  }
}

__global__ void Binning_Kernel(X_FLOAT* x, X_FLOAT* binned_x, int q_flag, int offset, int rmass_flag)
{
  const unsigned i = blockDim.x * blockIdx.x + threadIdx.x + offset;

  int binatoms = _natoms;

  if(offset == 0) binatoms = _nlocal ;

  if(i < binatoms) {
    // copy atom position from global device memory to local register
    // in this 3 steps to get as much coalesced access as possible
    X_FLOAT my_xX, my_xY, my_xZ;
    x += i;
    my_xX = *x;
    x += _nmax;
    my_xY = *x;
    x += _nmax;
    my_xZ = *x;
    //my_xX=x[i];
    //my_xY=x[i+_nmax];
    //my_xZ=x[i+2*_nmax];


    // calculate flat bin index
    int bx = __float2int_rd(rez_bin_size[0] * (my_xX - _sublo[0])) + 2;
    int by = __float2int_rd(rez_bin_size[1] * (my_xY - _sublo[1])) + 2;
    int bz = __float2int_rd(rez_bin_size[2] * (my_xZ - _sublo[2])) + 2;

    bx -= bx * negativCUDA(1.0f * bx);
    bx -= (bx - _bin_dim.x + 1) * negativCUDA(1.0f * _bin_dim.x - 1.0f - 1.0f * bx);
    by -= by * negativCUDA(1.0f * by);
    by -= (by - _bin_dim.y + 1) * negativCUDA(1.0f * _bin_dim.y - 1.0f - 1.0f * by);
    bz -= bz * negativCUDA(1.0f * bz);
    bz -= (bz - _bin_dim.z + 1) * negativCUDA(1.0f * _bin_dim.z - 1.0f - 1.0f * bz);


    const unsigned j = _bin_dim.z * (_bin_dim.y * bx + by) + bz;

    // add new atom to bin, get bin-array position
    const unsigned k = atomicAdd(& _bin_count_all[j], 1);

    if(offset == 0) atomicAdd(& _bin_count_local[j], 1);

    if(k < _bin_nmax) {
      // copy register values back to global device memory
      unsigned pos = 3 * _bin_nmax * j + k;
      _binpos[i] = pos;
      binned_x += pos;
      *binned_x = my_xX;
      binned_x += _bin_nmax;
      *binned_x = my_xY;
      binned_x += _bin_nmax;
      *binned_x = my_xZ;

      // also copy velocity and force accordingly

      binned_x  = _binned_v + pos;
      x  = _v + i;
      *binned_x = *x;
      binned_x += _bin_nmax;
      x += _nmax;
      *binned_x = *x;
      binned_x += _bin_nmax;
      x += _nmax;
      *binned_x = *x;

      binned_x  = _binned_f + pos;
      x  = _f + i;
      *binned_x = *x;
      binned_x += _bin_nmax;
      x += _nmax;
      *binned_x = *x;
      binned_x += _bin_nmax;
      x += _nmax;
      *binned_x = *x;

      pos = _bin_nmax * j + k;
      _binned_type [pos] = _type[i];
      _binned_tag  [pos] = _tag[i];

      if(rmass_flag)
        _binned_rmass[pos] = _rmass[i];

      if(q_flag)
        _binned_q    [pos] = _q[i];
    } else {
      // normally, this should not happen:
      int errorn = atomicAdd(bin_error_count, 1);
      MYEMUDBG(printf("# CUDA: Binning_Kernel: WARNING: atom %i ignored, no place left in bin %u\n", i, j);)
    }
  }
}

__global__ void ReverseBinning_Kernel(X_FLOAT* x, X_FLOAT* binned_x, int q_flag)
{
  const unsigned i = blockDim.x * blockIdx.x + threadIdx.x;

  if(i < _nlocal) {
    unsigned bin_pos3 = _binpos[i];
    unsigned bin_pos = bin_pos3 / (3 * _bin_nmax);
    bin_pos *= _bin_nmax;
    bin_pos += bin_pos3 - bin_pos * 3;

    binned_x  = _binned_x + bin_pos3;
    x  = x + i;
    *x = *binned_x;
    binned_x += _bin_nmax;
    x += _nmax;
    *x = *binned_x;
    binned_x += _bin_nmax;
    x += _nmax;
    *x = *binned_x;

    binned_x  = _binned_v + bin_pos3;
    x  = _v + i;
    *x = *binned_x;
    binned_x += _bin_nmax;
    x += _nmax;
    *x = *binned_x;
    binned_x += _bin_nmax;
    x += _nmax;
    *x = *binned_x;

    binned_x  = _binned_f + bin_pos3;
    x  = _f + i;
    *x = *binned_x;
    binned_x += _bin_nmax;
    x += _nmax;
    *x = *binned_x;
    binned_x += _bin_nmax;
    x += _nmax;
    *x = *binned_x;


    _type[i] = _binned_type[bin_pos];
    _tag[i] = _binned_tag[bin_pos];

    if(q_flag) _q[i] = _binned_q[bin_pos];
  }
}
