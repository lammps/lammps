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
#define RIMLARGER 1.000001
#define RIMSMALLER 0.999999
#define SMALL 1e-5

extern __shared__ int shared[];

template <const unsigned int data_mask>
__global__ void Cuda_AtomVecCuda_PackComm_Kernel(int* sendlist, int n, int maxlistlength, int iswap, X_CFLOAT dx, X_CFLOAT dy, X_CFLOAT dz, void* buffer)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  int* list = sendlist + iswap * maxlistlength;

  if(i < n) {
    int j = list[i];

    if(j > _nmax) _flag[0] = 1;

    int k = 0;

    if(data_mask & X_MASK) {
      ((X_CFLOAT*) buffer)[i + k * n] = _x[j] + dx;
      k++;
      ((X_CFLOAT*) buffer)[i + k * n] = _x[j + _nmax] + dy;
      k++;
      ((X_CFLOAT*) buffer)[i + k * n] = _x[j + 2 * _nmax] + dz;
      k++;
    }

    if(data_mask & V_MASK) {
      ((X_CFLOAT*) buffer)[i + k * n] = _v[j];
      k++;
      ((X_CFLOAT*) buffer)[i + k * n] = _v[j + _nmax];
      k++;
      ((X_CFLOAT*) buffer)[i + k * n] = _v[j + 2 * _nmax];
      k++;
    }

    if(data_mask & OMEGA_MASK) {
      ((X_CFLOAT*) buffer)[i + k * n] = _omega[j];
      k++;
      ((X_CFLOAT*) buffer)[i + k * n] = _omega[j + _nmax];
      k++;
      ((X_CFLOAT*) buffer)[i + k * n] = _omega[j + 2 * _nmax];
      k++;
    }

    if(data_mask & RADIUS_MASK)((X_CFLOAT*) buffer)[i + k * n] = _radius[j];

    k++;

    if(data_mask & RMASS_MASK)((X_CFLOAT*) buffer)[i + k * n] = _rmass[j];

    k++;
  }
}

template <const unsigned int data_mask>
__global__ void Cuda_AtomVecCuda_PackComm_Self_Kernel(int* sendlist, int n, int maxlistlength, int iswap, X_CFLOAT dx, X_CFLOAT dy, X_CFLOAT dz, int first)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  int* list = sendlist + iswap * maxlistlength;

  if(i < n) {
    int j = i;
    j = list[i];

    if(data_mask & X_MASK) {
      _x[i + first] = _x[j] + dx;
      _x[i + first + _nmax] = _x[j + _nmax] + dy;
      _x[i + first + 2 * _nmax] = _x[j + 2 * _nmax] + dz;
    }

    if(data_mask & V_MASK) {
      _v[i + first] = _v[j];
      _v[i + first + _nmax] = _v[j + _nmax];
      _v[i + first + 2 * _nmax] = _v[j + 2 * _nmax];
    }

    if(data_mask & OMEGA_MASK) {
      _omega[i + first] = _omega[j];
      _omega[i + first + _nmax] = _omega[j + _nmax];
      _omega[i + first + 2 * _nmax] = _omega[j + 2 * _nmax];
    }

    if(data_mask & RADIUS_MASK) _radius[i + first] = _radius[j];

    if(data_mask & RMASS_MASK) _rmass[i + first] = _rmass[j];
  }
}


template <const unsigned int data_mask>
__global__ void Cuda_AtomVecCuda_UnpackComm_Kernel(int n, int first, void* buffer)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < n) {
    int k = 0;

    if(data_mask & X_MASK) {
      _x[i + first] = ((X_CFLOAT*) buffer)[i + k * n];
      k++;
      _x[i + first + _nmax] = ((X_CFLOAT*) buffer)[i + k * n];
      k++;
      _x[i + first + 2 * _nmax] = ((X_CFLOAT*) buffer)[i + k * n];
      k++;
    }

    if(data_mask & V_MASK) {
      _v[i + first] = ((X_CFLOAT*) buffer)[i + k * n];
      k++;
      _v[i + first + _nmax] = ((X_CFLOAT*) buffer)[i + k * n];
      k++;
      _v[i + first + 2 * _nmax] = ((X_CFLOAT*) buffer)[i + k * n];
      k++;
    }

    if(data_mask & OMEGA_MASK) {
      _omega[i + first] = ((X_CFLOAT*) buffer)[i + k * n];
      k++;
      _omega[i + first + _nmax] = ((X_CFLOAT*) buffer)[i + k * n];
      k++;
      _omega[i + first + 2 * _nmax] = ((X_CFLOAT*) buffer)[i + k * n];
      k++;
    }

    if(data_mask & RADIUS_MASK) _radius[i + first] = ((X_CFLOAT*) buffer)[i + k * n];

    k++;

    if(data_mask & RMASS_MASK) _rmass[i + first] = ((X_CFLOAT*) buffer)[i + k * n];

    k++;
  }
}


__global__ void Cuda_AtomVecCuda_PackExchangeList_Kernel(int n, int dim)
{
  double* buf = (double*) _buffer;
  buf = &buf[1];

  //X_CFLOAT lo=slablo[iswap];
  //X_CFLOAT hi=slabhi[iswap];

  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  bool add = false;

  if(i < _nlocal) {
    double xdim_tmp = static_cast <double>(_x[i + dim * _nmax]);

    if(xdim_tmp < _sublo[dim] || xdim_tmp >= _subhi[dim]) {
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

  if(add && nsend + 1 < n)
    buf[nsend] = i;
}

template <const unsigned int data_mask>
__global__ void Cuda_AtomVecCuda_PackExchange_Kernel(int nsend, int* copylist)
{
  double* buf = (double*) _buffer;
  int k = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(k >= nsend) return;

  buf = &buf[1 + k];

  int i = static_cast <int>(buf[0]);
  int j = copylist[k];

  int m = 1;

  if(data_mask & X_MASK) {
    buf[(m++)*nsend] = static_cast <double>(_x[i]);
    buf[(m++)*nsend] = static_cast <double>(_x[i + _nmax]);
    buf[(m++)*nsend] = static_cast <double>(_x[i + 2 * _nmax]);
  }

  if(data_mask & V_MASK) {
    buf[(m++)*nsend] = _v[i];
    buf[(m++)*nsend] = _v[i + _nmax];
    buf[(m++)*nsend] = _v[i + 2 * _nmax];
  }

  if(data_mask & TAG_MASK) 		buf[(m++)*nsend] = _tag[i];

  if(data_mask & TYPE_MASK) 	buf[(m++)*nsend] = _type[i];

  if(data_mask & MASK_MASK) 	buf[(m++)*nsend] = _mask[i];

  if(data_mask & IMAGE_MASK) 	buf[(m++)*nsend] = _image[i];

  if(data_mask & Q_MASK) 		buf[(m++)*nsend] = _q[i];

  if(data_mask & MOLECULE_MASK) buf[(m++)*nsend] = _molecule[i];

  if(data_mask & RADIUS_MASK) 	buf[(m++)*nsend] = _radius[i];

  if(data_mask & DENSITY_MASK) 	buf[(m++)*nsend] = _density[i];

  if(data_mask & RMASS_MASK) 	buf[(m++)*nsend] = _rmass[i];

  if(data_mask & OMEGA_MASK) {
    buf[(m++)*nsend] = _omega[i];
    buf[(m++)*nsend] = _omega[i + _nmax];
    buf[(m++)*nsend] = _omega[i + 2 * _nmax];
  }

  /*  if(data_mask & NSPECIAL_MASK)
    {
    	buf[(m++)*nsend] = _nspecial[i];
    	buf[(m++)*nsend] = _nspecial[i+_nmax];
    	buf[(m++)*nsend] = _nspecial[i+2* _nmax];
    }*/

  if(i >= _nlocal) return;

  if(data_mask & X_MASK) {
    _x[i] = _x[j];
    _x[i + _nmax] = _x[j + _nmax];
    _x[i + 2 * _nmax] = _x[j + 2 * _nmax];
  }

  if(data_mask & V_MASK) {
    _v[i] = _v[j];
    _v[i + _nmax] = _v[j + _nmax];
    _v[i + 2 * _nmax] = _v[j + 2 * _nmax];
  }

  if(data_mask & TAG_MASK)		_tag[i] 	= _tag[j];

  if(data_mask & TYPE_MASK)		_type[i] 	= _type[j];

  if(data_mask & MASK_MASK)		_mask[i] 	= _mask[j];

  if(data_mask & IMAGE_MASK)	_image[i] 	= _image[j];

  if(data_mask & Q_MASK) 		_q[i] 		= _q[j];

  if(data_mask & MOLECULE_MASK) _molecule[i] = _molecule[j];

  if(data_mask & RADIUS_MASK) 	_radius[i] 	= _radius[j];

  if(data_mask & DENSITY_MASK) 	_density[i] = _density[j];

  if(data_mask & RMASS_MASK) 	_rmass[i] 	= _rmass[j];

  if(data_mask & OMEGA_MASK) {
    _omega[i] = _omega[j];
    _omega[i + _nmax] = _omega[j + _nmax];
    _omega[i + 2 * _nmax] = _omega[j + 2 * _nmax];
  }

  /* if(data_mask & NSPECIAL_MASK)
  {
  _nspecial[i] = _nspecial[j];
  _nspecial[i+_nmax] = _nspecial[j+_nmax];
  _nspecial[i+2* _nmax] = _nspecial[j+2* _nmax];
  }*/
}

template <const unsigned int data_mask>
__global__ void Cuda_AtomVecCuda_UnpackExchange_Kernel(int dim, int nsend, int* copylist)
{
  double* buf = (double*) _buffer;
  int k = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(k >= nsend) return;

  buf = &buf[1 + k];
  int i = -1;
  double xdim_tmp = buf[(1 + dim) * nsend];

  if(xdim_tmp >= _sublo[dim] - SMALL && xdim_tmp < _subhi[dim] + SMALL) {
    i = atomicAdd(_flag, 1) + _nlocal;

    int m = 1;

    if(data_mask & X_MASK) {
      _x[i] = buf[(m++) * nsend];
      _x[i + _nmax] = buf[(m++) * nsend];
      _x[i + 2 * _nmax] = buf[(m++) * nsend];
    }

    if(data_mask & V_MASK) {
      _v[i] = buf[(m++) * nsend];
      _v[i + _nmax] = buf[(m++) * nsend];
      _v[i + 2 * _nmax] = buf[(m++) * nsend];
    }

    if(data_mask & TAG_MASK) 	_tag[i] = buf[(m++) * nsend];

    if(data_mask & TYPE_MASK) 	_type[i] = buf[(m++) * nsend];

    if(data_mask & MASK_MASK) 	_mask[i] = buf[(m++) * nsend];

    if(data_mask & IMAGE_MASK) _image[i] = buf[(m++) * nsend];

    if(data_mask & Q_MASK) _q[i] = buf[(m++) * nsend];

    if(data_mask & MOLECULE_MASK) _molecule[i] = buf[(m++) * nsend];

    if(data_mask & RADIUS_MASK) _radius[i] = buf[(m++) * nsend];

    if(data_mask & DENSITY_MASK) _density[i] = buf[(m++) * nsend];

    if(data_mask & RMASS_MASK) _rmass[i] = buf[(m++) * nsend];

    if(data_mask & OMEGA_MASK) {
      _omega[i] = buf[(m++) * nsend];
      _omega[i + _nmax] = buf[(m++) * nsend];
      _omega[i + 2 * _nmax] = buf[(m++) * nsend];
    }

    /*  if(data_mask & NSPECIAL_MASK)
      {
       _nspecial[i] = buf[(m++)*nsend];
       _nspecial[i+_nmax] = buf[(m++)*nsend];
       _nspecial[i+2*_nmax] = buf[(m++)*nsend];
      }*/
  }

  copylist[k] = i;
}

template <const unsigned int data_mask>
__global__ void Cuda_AtomVecCuda_PackBorder_Kernel(int* sendlist, int n, int maxlistlength, int iswap, X_CFLOAT dx, X_CFLOAT dy, X_CFLOAT dz)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  int* list = sendlist + iswap * maxlistlength;

  if(i < n) {
    int j = list[i];
    int m = 0;

    if(data_mask & X_MASK) {
      ((X_CFLOAT*) _buffer)[i + (m++)*n] = _x[j] + dx;
      ((X_CFLOAT*) _buffer)[i + (m++)*n] = _x[j + _nmax] + dy;
      ((X_CFLOAT*) _buffer)[i + (m++)*n] = _x[j + 2 * _nmax] + dz;
    }

    if(data_mask & V_MASK) {
      ((X_CFLOAT*) _buffer)[i + (m++)*n] = _v[j];
      ((X_CFLOAT*) _buffer)[i + (m++)*n] = _v[j + _nmax];
      ((X_CFLOAT*) _buffer)[i + (m++)*n] = _v[j + 2 * _nmax];
    }

    if(data_mask & TAG_MASK)((X_CFLOAT*) _buffer)[i + (m++)*n] = _tag[j];

    if(data_mask & TYPE_MASK)((X_CFLOAT*) _buffer)[i + (m++)*n] = _type[j];

    if(data_mask & MASK_MASK)((X_CFLOAT*) _buffer)[i + (m++)*n] = _mask[j];

    if(data_mask & Q_MASK)((X_CFLOAT*) _buffer)[i + (m++)*n] = _q[j];

    if(data_mask & MOLECULE_MASK)((X_CFLOAT*) _buffer)[i + (m++)*n] = _molecule[j];

    if(data_mask & RADIUS_MASK)((X_CFLOAT*) _buffer)[i + (m++)*n] = _radius[i];

    if(data_mask & DENSITY_MASK)((X_CFLOAT*) _buffer)[i + (m++)*n] = _density[i];

    if(data_mask & RMASS_MASK)((X_CFLOAT*) _buffer)[i + (m++)*n] = _rmass[i];

    if(data_mask & OMEGA_MASK) {
      ((X_CFLOAT*) _buffer)[i + (m++)*n] = _omega[i];
      ((X_CFLOAT*) _buffer)[i + (m++)*n] = _omega[i + _nmax];
      ((X_CFLOAT*) _buffer)[i + (m++)*n] = _omega[i + 2 * _nmax];
    }
  }
}



template <const unsigned int data_mask>
__global__ void Cuda_AtomVecCuda_PackBorder_Self_Kernel(int* sendlist, int n, int maxlistlength, int iswap, X_CFLOAT dx, X_CFLOAT dy, X_CFLOAT dz, int first)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  int* list = sendlist + iswap * maxlistlength;

  if(i < n) {
    int j = list[i];

    if(data_mask & X_MASK) {
      _x[i + first] = _x[j] + dx;
      _x[i + first + _nmax] = _x[j + _nmax] + dy;
      _x[i + first + 2 * _nmax] = _x[j + 2 * _nmax] + dz;
    }

    if(data_mask & V_MASK) {
      _v[i + first] = _v[j];
      _v[i + first + _nmax] = _v[j + _nmax];
      _v[i + first + 2 * _nmax] =  _v[j + 2 * _nmax];
    }

    if(data_mask & TAG_MASK) _tag[i + first] = _tag[j];

    if(data_mask & TYPE_MASK) _type[i + first] = _type[j];

    if(data_mask & MASK_MASK) _mask[i + first] = _mask[j];

    if(data_mask & Q_MASK) _q[i + first] = _q[j];

    if(data_mask & MOLECULE_MASK) _molecule[i + first] = _molecule[j];

    if(data_mask & RADIUS_MASK) _radius[i + first] = _radius[j];

    if(data_mask & DENSITY_MASK) _density[i + first] = _density[j];

    if(data_mask & RMASS_MASK) _rmass[i + first] = _rmass[j];

    if(data_mask & OMEGA_MASK) {
      _omega[i + first] = _omega[j];
      _omega[i + first + _nmax] = _omega[j + _nmax];
      _omega[i + first + 2 * _nmax] =  _omega[j + 2 * _nmax];
    }
  }
}

template <const unsigned int data_mask>
__global__ void Cuda_AtomVecCuda_UnpackBorder_Kernel(int n, int first)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < n) {
    if(i + first < _nmax) {
      int m = 0;

      if(data_mask & X_MASK) {
        _x[i + first] = ((X_CFLOAT*) _buffer)[i + (m++) * n];
        _x[i + first + _nmax] = ((X_CFLOAT*) _buffer)[i + (m++) * n];
        _x[i + first + 2 * _nmax] = ((X_CFLOAT*) _buffer)[i + (m++) * n];
      }

      if(data_mask & V_MASK) {
        _v[i + first] = ((X_CFLOAT*) _buffer)[i + (m++) * n];
        _v[i + first + _nmax] = ((X_CFLOAT*) _buffer)[i + (m++) * n];
        _v[i + first + 2 * _nmax] = ((X_CFLOAT*) _buffer)[i + (m++) * n];
      }

      if(data_mask & TAG_MASK) _tag[i + first] = static_cast<int>(((X_CFLOAT*) _buffer)[i + (m++) * n]);

      if(data_mask & TYPE_MASK) _type[i + first] = static_cast<int>(((X_CFLOAT*) _buffer)[i + (m++) * n]);

      if(data_mask & MASK_MASK) _mask[i + first] = static_cast<int>(((X_CFLOAT*) _buffer)[i + (m++) * n]);

      if(data_mask & Q_MASK) _q[i + first] = ((X_CFLOAT*) _buffer)[i + (m++) * n];

      if(data_mask & MOLECULE_MASK) _molecule[i + first] = static_cast<int>(((X_CFLOAT*) _buffer)[i + (m++) * n]);

      if(data_mask & RADIUS_MASK) _radius[i + first] = ((X_CFLOAT*) _buffer)[i + (m++) * n];

      if(data_mask & DENSITY_MASK) _density[i + first] = ((X_CFLOAT*) _buffer)[i + (m++) * n];

      if(data_mask & RMASS_MASK) _rmass[i + first] = ((X_CFLOAT*) _buffer)[i + (m++) * n];

      if(data_mask & OMEGA_MASK) {
        _omega[i + first] = ((X_CFLOAT*) _buffer)[i + (m++) * n];
        _omega[i + first + _nmax] = ((X_CFLOAT*) _buffer)[i + (m++) * n];
        _omega[i + first + 2 * _nmax] = ((X_CFLOAT*) _buffer)[i + (m++) * n];
      }
    } else {
      _flag[0] = 1;
    }
  }
}


