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

static inline __device__ void check_distance(X_FLOAT &xtmp, X_FLOAT &ytmp, X_FLOAT &ztmp, int &i, int groupbit)
{
  if(_dist_check) {
    X_FLOAT tmp = xtmp - _xhold[i];
    X_FLOAT d = tmp * tmp;
    tmp = ytmp - _xhold[i + _maxhold];
    d += tmp * tmp;
    tmp = ztmp - _xhold[i + 2 * _maxhold];
    d += tmp * tmp;

    d = ((i < _nlocal) && (_mask[i] & groupbit)) ? d : X_F(0.0);

    if(not __all(d <= _triggerneighsq))
      _reneigh_flag[0] = 1;
  }
}


__global__ void FixNVECuda_InitialIntegrate_Kernel(int groupbit)
{
  X_FLOAT xtmp, ytmp, ztmp;
#ifdef CUDA_USE_BINNING

  const unsigned bin = gridDim.y * blockIdx.x + blockIdx.y;

  if(threadIdx.x < _bin_count_local[bin]) {
    const int i = 3 * blockDim.x * bin + threadIdx.x;

    if(_mask[i] & groupbit) {
      F_FLOAT* my_f = _binned_f + i;
      V_FLOAT* my_v = _binned_v + i;
      X_FLOAT* my_x = _binned_x + i;

      V_FLOAT 		dtfm = _dtf

                         if(_rmass_flag) dtfm *= V_F(1.0) / _binned_rmass[i];
      else 			dtfm *= V_F(1.0) / _mass[_binned_type[blockDim.x * bin + threadIdx.x]];

      V_FLOAT v_mem;
      v_mem = *my_v += dtfm * (*my_f);
      xtmp = *my_x += _dtv * v_mem;
      my_f += blockDim.x;
      my_v += blockDim.x;
      my_x += blockDim.x;
      v_mem = *my_v += dtfm * (*my_f);
      ytmp = *my_x += _dtv * v_mem;
      my_f += blockDim.x;
      my_v += blockDim.x;
      my_x += blockDim.x;
      v_mem = *my_v += dtfm * (*my_f);
      ztmp = *my_x += _dtv * v_mem;
    }
  }

#else

  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < _nlocal && _mask[i] & groupbit) {
    F_FLOAT* my_f = _f + i;
    V_FLOAT* my_v = _v + i;
    X_FLOAT* my_x = _x + i;

    V_FLOAT 		dtfm = _dtf;

    if(_rmass_flag) dtfm *= V_F(1.0) / _rmass[i];
    else 			dtfm *= V_F(1.0) / _mass[_type[i]];

    V_FLOAT v_mem;
    v_mem = *my_v += dtfm * (*my_f);
    xtmp = *my_x += _dtv * v_mem;
    my_f += _nmax;
    my_v += _nmax;
    my_x += _nmax;
    v_mem = *my_v += dtfm * (*my_f);
    ytmp = *my_x += _dtv * v_mem;
    my_f += _nmax;
    my_v += _nmax;
    my_x += _nmax;
    v_mem = *my_v += dtfm * (*my_f);
    ztmp = *my_x += _dtv * v_mem;
  }

#endif

  check_distance(xtmp, ytmp, ztmp, i, groupbit);
}

__global__ void FixNVECuda_FinalIntegrate_Kernel(int groupbit)
{
#ifdef CUDA_USE_BINNING

  const unsigned bin = gridDim.y * blockIdx.x + blockIdx.y;

  if(threadIdx.x < _bin_count_local[bin]) {
    const int i = 3 * blockDim.x * bin + threadIdx.x;

    if(_mask[i] & groupbit) {
      F_FLOAT* my_f = _binned_f + i;
      V_FLOAT* my_v = _binned_v + i;

      V_FLOAT 		dtfm = _dtf

                         if(_rmass_flag) dtfm *= V_F(1.0) / _binned_rmass[i];
      else 			dtfm *= V_F(1.0) / _mass[_binned_type[blockDim.x * bin + threadIdx.x]];

      *my_v += dtfm * (*my_f);
      my_f += blockDim.x;
      my_v += blockDim.x;
      *my_v += dtfm * (*my_f);
      my_f += blockDim.x;
      my_v += blockDim.x;
      *my_v += dtfm * (*my_f);
    }
  }

#else

  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < _nlocal && _mask[i] & groupbit) {
    F_FLOAT* my_f = _f + i;
    V_FLOAT* my_v = _v + i;

    V_FLOAT 		dtfm = _dtf;

    if(_rmass_flag) dtfm *= V_F(1.0) / _rmass[i];
    else 			dtfm *= V_F(1.0) / _mass[_type[i]];

    *my_v += dtfm * (*my_f);
    my_f += _nmax;
    my_v += _nmax;
    *my_v += dtfm * (*my_f);
    my_f += _nmax;
    my_v += _nmax;
    *my_v += dtfm * (*my_f);
  }

#endif
}



