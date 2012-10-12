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

extern __shared__ X_FLOAT sharedmem[];

#define BIG 1e10
__global__ void Domain_PBC_Kernel(int deform_remap, int deform_groupbit, int box_change)
{
  int idim, otherdims;
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  X_FLOAT lo[3];
  X_FLOAT hi[3];
  X_FLOAT* period;

  if(_triclinic == 0) {
    lo[0] = _boxlo[0];
    lo[1] = _boxlo[1];
    lo[2] = _boxlo[2];

    hi[0] = _boxhi[0];
    hi[1] = _boxhi[1];
    hi[2] = _boxhi[2];
    period = _prd;
  } else {
    lo[0] = _boxlo_lamda[0];
    lo[1] = _boxlo_lamda[1];
    lo[2] = _boxlo_lamda[2];

    hi[0] = _boxhi_lamda[0];
    hi[1] = _boxhi_lamda[1];
    hi[2] = _boxhi_lamda[2];
    period = _prd_lamda;
  }


  X_FLOAT tmpx = X_F(0.5) * (hi[0] + lo[0]);
  X_FLOAT tmpy = X_F(0.5) * (hi[1] + lo[1]);
  X_FLOAT tmpz = X_F(0.5) * (hi[2] + lo[2]);

  X_FLOAT* buf = (X_FLOAT*) _buffer;
  buf += blockIdx.x * gridDim.y + blockIdx.y;
  buf[0] = tmpx;
  buf += gridDim.x * gridDim.y;
  buf[0] = tmpx;
  buf += gridDim.x * gridDim.y;
  buf[0] = tmpy;
  buf += gridDim.x * gridDim.y;
  buf[0] = tmpy;
  buf += gridDim.x * gridDim.y;
  buf[0] = tmpz;
  buf += gridDim.x * gridDim.y;
  buf[0] = tmpz;

  if(i < _nlocal) {

    if(_periodicity[0]) {
      if(_x[i] < lo[0]) {
        _x[i] += period[0];

        if(deform_remap && _mask[i] & deform_groupbit) _v[i] += _h_rate[0];

        idim = _image[i] & 1023;
        otherdims = _image[i] ^ idim;
        idim--;
        idim &= 1023;
        _image[i] = otherdims | idim;
      }

      if(_x[i] >= hi[0]) {
        _x[i] -= period[0];
        _x[i] = MAX(_x[i], lo[0]);

        if(deform_remap && _mask[i] & deform_groupbit) _v[i] -= _h_rate[0];

        idim = _image[i] & 1023;
        otherdims = _image[i] ^ idim;
        idim++;
        idim &= 1023;
        _image[i] = otherdims | idim;
      }
    }

    if(_periodicity[1]) {
      if(_x[i + _nmax] < lo[1]) {
        _x[i + _nmax] += period[1];

        if(deform_remap && _mask[i] & deform_groupbit) {
          _v[i] += _h_rate[5];
          _v[i + _nmax] += _h_rate[1];
        }

        idim = (_image[i] >> 10) & 1023;
        otherdims = _image[i] ^ (idim << 10);
        idim--;
        idim &= 1023;
        _image[i] = otherdims | (idim << 10);
      }

      if(_x[i + _nmax] >= hi[1]) {
        _x[i + _nmax] -= period[1];
        _x[i + _nmax] = MAX(_x[i + _nmax], lo[1]);

        if(deform_remap && _mask[i] & deform_groupbit) {
          _v[i] -= _h_rate[5];
          _v[i + _nmax] -= _h_rate[1];
        }

        idim = (_image[i] >> 10) & 1023;
        otherdims = _image[i] ^ (idim << 10);
        idim++;
        idim &= 1023;
        _image[i] = otherdims | (idim << 10);
      }
    }

    if(_periodicity[2]) {
      if(_x[i + 2 * _nmax] < lo[2]) {
        _x[i + 2 * _nmax] += period[2];

        if(deform_remap && _mask[i] & deform_groupbit) {
          _v[i] += _h_rate[4];
          _v[i + _nmax] += _h_rate[3];
          _v[i + 2 * _nmax] += _h_rate[2];
        }

        idim = _image[i] >> 20;
        otherdims = _image[i] ^ (idim << 20);
        idim--;
        idim &= 1023;
        _image[i] = otherdims | (idim << 20);
      }

      if(_x[i + 2 * _nmax] >= hi[2]) {
        _x[i + 2 * _nmax] -= period[2];
        _x[i + 2 * _nmax] = MAX(_x[i + 2 * _nmax], lo[2]);

        if(deform_remap && _mask[i] & deform_groupbit) {
          _v[i] -= _h_rate[4];
          _v[i + _nmax] -= _h_rate[3];
          _v[i + 2 * _nmax] -= _h_rate[2];
        }

        idim = _image[i] >> 20;
        otherdims = _image[i] ^ (idim << 20);
        idim++;
        idim &= 1023;
        _image[i] = otherdims | (idim << 20);
      }
    }

    if(box_change) {
      tmpx = _x[i];
      tmpy = _x[i + _nmax];
      tmpz = _x[i + 2 * _nmax];


    }
  }

  __syncthreads();

  if(box_change) {
    X_FLOAT minx = BIG;
    X_FLOAT maxx = -BIG;
    X_FLOAT miny = BIG;
    X_FLOAT maxy = -BIG;
    X_FLOAT minz = BIG;
    X_FLOAT maxz = -BIG;

    if(not _periodicity[0]) {
      sharedmem[threadIdx.x] = tmpx;
      minOfBlock(sharedmem);
      minx = sharedmem[0];
      __syncthreads();
      sharedmem[threadIdx.x] = tmpx;
      maxOfBlock(sharedmem);
      maxx = sharedmem[0];
      __syncthreads();
    } else {
      minx = lo[0];
      maxx = hi[0];
    }

    if(not _periodicity[1]) {
      sharedmem[threadIdx.x] = tmpy;
      minOfBlock(sharedmem);
      miny = sharedmem[0];
      __syncthreads();
      sharedmem[threadIdx.x] = tmpy;
      maxOfBlock(sharedmem);
      maxy = sharedmem[0];
      __syncthreads();
    } else {
      minx = lo[1];
      maxx = hi[1];
    }

    if(not _periodicity[2]) {
      sharedmem[threadIdx.x] = tmpz;
      minOfBlock(sharedmem);
      minz = sharedmem[0];
      __syncthreads();
      sharedmem[threadIdx.x] = tmpz;
      maxOfBlock(sharedmem);
      maxz = sharedmem[0];
      __syncthreads();
    } else {
      minz = lo[2];
      maxz = hi[2];
    }

    if(threadIdx.x == 0) {
      buf = (X_FLOAT*) _buffer;
      buf += blockIdx.x * gridDim.y + blockIdx.y;
      buf[0] = minx;
      buf += gridDim.x * gridDim.y;
      buf[0] = maxx;
      buf += gridDim.x * gridDim.y;
      buf[0] = miny;
      buf += gridDim.x * gridDim.y;
      buf[0] = maxy;
      buf += gridDim.x * gridDim.y;
      buf[0] = minz;
      buf += gridDim.x * gridDim.y;
      buf[0] = maxz;
    }
  }
}

__global__ void Domain_reduceBoxExtent(double* extent, int n)
{
  X_FLOAT* buf = (X_FLOAT*) _buffer;
  buf += blockIdx.x * n;
  copyGlobToShared(buf, sharedmem, n);

  if(blockIdx.x % 2 == 0)
    minOfData(sharedmem, n);
  else
    maxOfData(sharedmem, n);

  extent[blockIdx.x] = sharedmem[0];
}

__global__ void Domain_lamda2x_Kernel(int n)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < n) {
    X_FLOAT ytmp = _x[i + _nmax];
    X_FLOAT ztmp = _x[i + 2 * _nmax];
    _x[i] = _h[0] * _x[i] + _h[5] * ytmp + _h[4] * ztmp + _boxlo[0];
    _x[i + _nmax] = _h[1] * ytmp + _h[3] * ztmp + _boxlo[1];
    _x[i + 2 * _nmax] = _h[2] * ztmp + _boxlo[2];
  }
}

__global__ void Domain_x2lamda_Kernel(int n)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  X_FLOAT delta[3];

  if(i < n) {
    delta[0] = _x[i] - _boxlo[0];
    delta[1] = _x[i + _nmax] - _boxlo[1];
    delta[2] = _x[i + 2 * _nmax] - _boxlo[2];

    _x[i] = _h_inv[0] * delta[0] + _h_inv[5] * delta[1] + _h_inv[4] * delta[2];
    _x[i + _nmax] = _h_inv[1] * delta[1] + _h_inv[3] * delta[2];
    _x[i + 2 * _nmax] = _h_inv[2] * delta[2];
  }
}
