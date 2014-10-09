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




static __device__ inline F_CFLOAT4 fetchRhor(int i)
{
#ifdef CUDA_USE_TEXTURE
#if F_PRECISION == 1
  return tex1Dfetch(_rhor_spline_tex, i);
#else
  return tex1Dfetch_double_f(_rhor_spline_tex, i);
#endif
#else
  return _rhor_spline[i];
#endif
}

static __device__ inline F_CFLOAT4 fetchZ2r(int i)
{
#ifdef CUDA_USE_TEXTURE
#if F_PRECISION == 1
  return tex1Dfetch(_z2r_spline_tex, i);
#else
  return tex1Dfetch_double_f(_z2r_spline_tex, i);
#endif
#else
  return _z2r_spline[i];
#endif
}

__global__ void PairEAMCuda_Kernel1(int eflag, int vflag, int eflag_atom, int vflag_atom)
{
  ENERGY_CFLOAT* sharedE;
  ENERGY_CFLOAT* sharedV = &sharedmem[threadIdx.x];


  if(eflag || eflag_atom) {
    sharedE = &sharedmem[threadIdx.x];
    sharedE[0] = ENERGY_F(0.0);
    sharedV += blockDim.x;
  }

  if(vflag || vflag_atom) {
    sharedV[0 * blockDim.x] = ENERGY_F(0.0);
    sharedV[1 * blockDim.x] = ENERGY_F(0.0);
    sharedV[2 * blockDim.x] = ENERGY_F(0.0);
    sharedV[3 * blockDim.x] = ENERGY_F(0.0);
    sharedV[4 * blockDim.x] = ENERGY_F(0.0);
    sharedV[5 * blockDim.x] = ENERGY_F(0.0);
  }

  int ii = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  X_CFLOAT xtmp, ytmp, ztmp;
  X_CFLOAT4 myxtype;
  F_CFLOAT delx, dely, delz;
  int itype;
  int i = _nlocal;
  int jnum = 0;
  int* jlist;

  if(ii < _inum) {
    i = _ilist[ii];

    myxtype = fetchXType(i);
    xtmp = myxtype.x;
    ytmp = myxtype.y;
    ztmp = myxtype.z;
    itype = static_cast <int>(myxtype.w);

    jnum = _numneigh[i];

    jlist = &_neighbors[i];

    if(i < _nlocal)
      _rho[i] = F_F(0.0);
  }

  __syncthreads();

  for(int jj = 0; jj < jnum; jj++) {
    if(ii < _inum)
      if(jj < jnum) {
        const int j = jlist[jj * _nlocal];
        myxtype = fetchXType(j);
        delx = xtmp - myxtype.x;
        dely = ytmp - myxtype.y;
        delz = ztmp - myxtype.z;
        int jtype = static_cast <int>(myxtype.w);
        const F_CFLOAT rsq = delx * delx + dely * dely + delz * delz;

        if(rsq < _cutsq_global) {
          F_CFLOAT p = sqrt(rsq) * _rdr + F_F(1.0);
          int m = static_cast<int>(p);
          m = MIN(m, _nr - 1);
          p -= m;
          p = MIN(p, F_F(1.0));

          int k = (static_cast <int>(_type2rhor[jtype * _cuda_ntypes + itype]) * (_nr + 1) + m) * 2;
          F_CFLOAT4 c = fetchRhor(k + 1);
          _rho[i] += ((c.w * p + c.x) * p + c.y) * p + c.z;
        }
      }
  }

  if(ii < _inum) {

    F_CFLOAT p = _rho[i] * _rdrho + F_F(1.0);
    int m = static_cast<int>(p);
    m = MAX(1, MIN(m, _nrho - 1));
    p -= m;
    p = MIN(p, F_F(1.0));
    F_CFLOAT* coeff = &_frho_spline[(static_cast <int>(_type2frho[itype]) * (_nrho + 1) + m) * EAM_COEFF_LENGTH];
    _fp[i] = (coeff[0] * p + coeff[1]) * p + coeff[2];

    if(eflag || eflag_atom) {
      sharedmem[threadIdx.x] += ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
    }

  }

  __syncthreads();

  if(eflag || eflag_atom) {
    if(i < _nlocal && eflag_atom)
      _eatom[i] += sharedmem[threadIdx.x];

    reduceBlock(sharedmem);
    ENERGY_CFLOAT* buffer = (ENERGY_CFLOAT*) _buffer;
    buffer[blockIdx.x * gridDim.y + blockIdx.y] = ENERGY_F(2.0) * sharedmem[0];
  }
}

__global__ void PairEAMCuda_Kernel2(int eflag, int vflag, int eflag_atom, int vflag_atom)
{
  ENERGY_CFLOAT evdwl = ENERGY_F(0.0);

  ENERGY_CFLOAT* sharedE;
  ENERGY_CFLOAT* sharedV = &sharedmem[threadIdx.x];


  if(eflag || eflag_atom) {
    sharedE = &sharedmem[threadIdx.x];
    sharedE[0] = ENERGY_F(0.0);
    sharedV += blockDim.x;
  }

  if(vflag || vflag_atom) {
    sharedV[0 * blockDim.x] = ENERGY_F(0.0);
    sharedV[1 * blockDim.x] = ENERGY_F(0.0);
    sharedV[2 * blockDim.x] = ENERGY_F(0.0);
    sharedV[3 * blockDim.x] = ENERGY_F(0.0);
    sharedV[4 * blockDim.x] = ENERGY_F(0.0);
    sharedV[5 * blockDim.x] = ENERGY_F(0.0);
  }

  int ii = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  X_CFLOAT xtmp, ytmp, ztmp;
  X_CFLOAT4 myxtype;
  F_CFLOAT fxtmp, fytmp, fztmp, fpair;
  F_CFLOAT delx, dely, delz;
  int itype, i;
  int jnum = 0;
  int* jlist;

  if(ii < _inum) {
    i = _ilist[ii];

    myxtype = fetchXType(i);
    xtmp = myxtype.x;
    ytmp = myxtype.y;
    ztmp = myxtype.z;
    itype = static_cast <int>(myxtype.w);
    fxtmp = F_F(0.0);
    fytmp = F_F(0.0);
    fztmp = F_F(0.0);

    jnum = _numneigh[i];

    jlist = &_neighbors[i];

    if(i < _nlocal)
      _rho[i] = F_F(0.0);
  }

  if(ii < gridDim.x * gridDim.y) evdwl = ((ENERGY_CFLOAT*) _buffer)[ii];

  __syncthreads();

  for(int jj = 0; jj < jnum; jj++) {
    if(ii < _inum)
      if(jj < jnum) {
        const int j = jlist[jj * _nlocal];
        myxtype = fetchXType(j);
        delx = xtmp - myxtype.x;
        dely = ytmp - myxtype.y;
        delz = ztmp - myxtype.z;
        int jtype = static_cast <int>(myxtype.w);
        const F_CFLOAT rsq = delx * delx + dely * dely + delz * delz;

        if(rsq < _cutsq_global) {
          F_CFLOAT r = _SQRT_(rsq);
          F_CFLOAT p = r * _rdr + F_F(1.0);
          int m = static_cast<int>(p);
          m = MIN(m, _nr - 1);
          p -= m;
          p = MIN(p, F_F(1.0));

          int k = (static_cast <int>(_type2rhor[itype * _cuda_ntypes + jtype]) * (_nr + 1) + m) * 2;
          F_CFLOAT4 c = fetchRhor(k);
          F_CFLOAT rhoip = (c.x * p + c.y) * p + c.z;
          k = (static_cast <int>(_type2rhor[jtype * _cuda_ntypes + itype]) * (_nr + 1) + m) * 2;
          c = fetchRhor(k);
          F_CFLOAT rhojp = (c.x * p + c.y) * p + c.z;
          k = (static_cast <int>(_type2z2r[itype * _cuda_ntypes + jtype]) * (_nr + 1) + m) * 2;
          c = fetchZ2r(k);
          F_CFLOAT z2p = (c.x * p + c.y) * p + c.z;
          c = fetchZ2r(k + 1);
          F_CFLOAT z2 = ((c.w * p + c.x) * p + c.y) * p + c.z;

          F_CFLOAT recip = F_F(1.0) / r;
          F_CFLOAT phi = z2 * recip;
          F_CFLOAT phip = z2p * recip - phi * recip;
          F_CFLOAT psip = _fp[i] * rhojp + _fp[j] * rhoip + phip;
          fpair = -psip * recip;

          F_CFLOAT dxfp, dyfp, dzfp;
          fxtmp += dxfp = delx * fpair;
          fytmp += dyfp = dely * fpair;
          fztmp += dzfp = delz * fpair;
          evdwl += phi;

          if(vflag || vflag_atom) {
            sharedV[0 * blockDim.x] += delx * dxfp;
            sharedV[1 * blockDim.x] += dely * dyfp;
            sharedV[2 * blockDim.x] += delz * dzfp;
            sharedV[3 * blockDim.x] += delx * dyfp;
            sharedV[4 * blockDim.x] += delx * dzfp;
            sharedV[5 * blockDim.x] += dely * dzfp;
          }
        }
      }
  }

  __syncthreads();

  if(ii < _inum) {
    F_CFLOAT* my_f;

    if(_collect_forces_later) {
      ENERGY_CFLOAT* buffer = (ENERGY_CFLOAT*) _buffer;

      if(eflag) {
        buffer = &buffer[1 * gridDim.x * gridDim.y];
      }

      if(vflag) {
        buffer = &buffer[6 * gridDim.x * gridDim.y];
      }

      my_f = (F_CFLOAT*) buffer;
      my_f += i;
      *my_f = fxtmp;
      my_f += _nmax;
      *my_f = fytmp;
      my_f += _nmax;
      *my_f = fztmp;
    } else {
      my_f = _f + i;
      *my_f += fxtmp;
      my_f += _nmax;
      *my_f += fytmp;
      my_f += _nmax;
      *my_f += fztmp;
    }
  }

  __syncthreads();

  if(eflag) {
    sharedE[0] = evdwl;
  }

  if(eflag_atom && i < _nlocal) {
    _eatom[i] += evdwl;
  }

  if(vflag_atom && i < _nlocal) {
    _vatom[i]         += ENERGY_F(0.5) * sharedV[0 * blockDim.x];
    _vatom[i + _nmax]   += ENERGY_F(0.5) * sharedV[1 * blockDim.x];
    _vatom[i + 2 * _nmax] += ENERGY_F(0.5) * sharedV[2 * blockDim.x];
    _vatom[i + 3 * _nmax] += ENERGY_F(0.5) * sharedV[3 * blockDim.x];
    _vatom[i + 4 * _nmax] += ENERGY_F(0.5) * sharedV[4 * blockDim.x];
    _vatom[i + 5 * _nmax] += ENERGY_F(0.5) * sharedV[5 * blockDim.x];
  }

  if(vflag || eflag) PairVirialCompute_A_Kernel(eflag, vflag, 0);
}

__global__ void PairEAMCuda_PackComm_Kernel(int* sendlist, int n, int maxlistlength, int iswap, F_CFLOAT* buffer)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  int* list = sendlist + iswap * maxlistlength;

  if(i < n) {
    int j = list[i];
    buffer[i] = _fp[j];
  }
}

__global__ void PairEAMCuda_UnpackComm_Kernel(int n, int first, F_CFLOAT* buffer)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < n) {
    _fp[i + first] = buffer[i];
  }
}
