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

    X_FLOAT d = X_F(0.0);

    if(i < _nlocal) {
      X_FLOAT tmp = xtmp - _xhold[i];
      d = tmp * tmp;
      tmp = ytmp - _xhold[i + _maxhold];
      d += tmp * tmp;
      tmp = ztmp - _xhold[i + 2 * _maxhold];
      d += tmp * tmp;

      d = ((_mask[i] & groupbit)) ? d : X_F(0.0);
    }

    if(not __all(d <= _triggerneighsq))
      _reneigh_flag[0] = 1;
  }
}

__global__ void FixNHCuda_nh_v_press_Kernel(int groupbit, F_FLOAT3 factor, int p_triclinic, F_FLOAT3 factor2)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < _nlocal && _mask[i] & groupbit) {
    V_FLOAT* my_v = _v + i;
    V_FLOAT vx = my_v[0];
    V_FLOAT vy = my_v[_nmax];
    V_FLOAT vz = my_v[2 * _nmax];
    vx *= factor.x;
    vy *= factor.y;
    vz *= factor.z;

    if(p_triclinic) {
      vx += vy * factor2.z + vz * factor2.y;
      vy += vz * factor2.x;
    }

    vx *= factor.x;
    vy *= factor.y;
    vz *= factor.z;
    my_v[0]       = vx;
    my_v[_nmax]   = vy;
    my_v[2 * _nmax] = vz;
  }

}

__global__ void FixNHCuda_nh_v_temp_Kernel(int groupbit, F_FLOAT factor_eta)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < _nlocal && _mask[i] & groupbit) {
    V_FLOAT* my_v = _v + i;
    my_v[0] *= factor_eta;
    my_v[_nmax] *= factor_eta;
    my_v[2 * _nmax] *= factor_eta;
  }

}

__global__ void FixNHCuda_nh_v_press_and_nve_v_NoBias_Kernel(int groupbit, F_FLOAT3 factor, int p_triclinic, F_FLOAT3 factor2)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < _nlocal && _mask[i] & groupbit) {
    F_FLOAT* my_f = _f + i;
    V_FLOAT* my_v = _v + i;

    V_FLOAT 		dtfm = _dtf;

    if(_rmass_flag) dtfm *= V_F(1.0) / _rmass[i];
    else 			dtfm *= V_F(1.0) / _mass[_type[i]];

    V_FLOAT vx = my_v[0];
    V_FLOAT vy = my_v[_nmax];
    V_FLOAT vz = my_v[2 * _nmax];
    vx *= factor.x;
    vy *= factor.y;
    vz *= factor.z;

    if(p_triclinic) {
      vx += vy * factor2.z + vz * factor2.y;
      vy += vz * factor2.x;
    }

    vx *= factor.x;
    vy *= factor.y;
    vz *= factor.z;
    my_v[0]       = vx + dtfm * my_f[0];
    my_v[_nmax]   = vy + dtfm * my_f[_nmax];
    my_v[2 * _nmax] = vz + dtfm * my_f[_nmax * 2];
  }

}

__global__ void FixNHCuda_nve_v_Kernel(int groupbit)
{

  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < _nlocal && _mask[i] & groupbit) {
    F_FLOAT* my_f = _f + i;
    V_FLOAT* my_v = _v + i;

    V_FLOAT 		dtfm = _dtf;

    if(_rmass_flag) dtfm *= V_F(1.0) / _rmass[i];
    else 			dtfm *= V_F(1.0) / _mass[_type[i]];

    *my_v = (*my_v + dtfm * (*my_f));
    my_f += _nmax;
    my_v += _nmax;
    *my_v = (*my_v + dtfm * (*my_f));
    my_f += _nmax;
    my_v += _nmax;
    *my_v = (*my_v + dtfm * (*my_f));
  }
}

__global__ void FixNHCuda_nve_x_Kernel(int groupbit)
{
  X_FLOAT xtmp, ytmp, ztmp;

  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < _nlocal && _mask[i] & groupbit) {
    V_FLOAT* my_v = _v + i;
    X_FLOAT* my_x = _x + i;

    xtmp = *my_x += _dtv * *my_v;
    my_v += _nmax;
    my_x += _nmax;
    ytmp = *my_x += _dtv * *my_v;
    my_v += _nmax;
    my_x += _nmax;
    ztmp = *my_x += _dtv * *my_v;
  }

  check_distance(xtmp, ytmp, ztmp, i, groupbit);
}


__global__ void FixNHCuda_nve_v_and_nh_v_press_NoBias_Kernel(int groupbit, F_FLOAT3 factor, int p_triclinic, F_FLOAT3 factor2)
{

  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < _nlocal && _mask[i] & groupbit) {
    F_FLOAT* my_f = _f + i;
    V_FLOAT* my_v = _v + i;

    V_FLOAT 		dtfm = _dtf;

    if(_rmass_flag) dtfm *= V_F(1.0) / _rmass[i];
    else 			dtfm *= V_F(1.0) / _mass[_type[i]];

    V_FLOAT vx = my_v[0] + dtfm * my_f[0];
    V_FLOAT vy = my_v[_nmax] + dtfm * my_f[_nmax];
    V_FLOAT vz = my_v[2 * _nmax] + dtfm * my_f[2 * _nmax];

    vx *= factor.x;
    vy *= factor.y;
    vz *= factor.z;

    if(p_triclinic) {
      vx += vy * factor2.z + vz * factor2.y;
      vy += vz * factor2.x;
    }

    vx *= factor.x;
    vy *= factor.y;
    vz *= factor.z;
    my_v[0]       = vx;
    my_v[_nmax]   = vy;
    my_v[2 * _nmax] = vz;

  }
}

