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



__global__ void Cuda_FixTempRescaleLimitCuda_EndOfStep_Kernel(int groupbit, V_FLOAT factor, V_FLOAT limit)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < _nlocal)
    if(_mask[i] & groupbit) {
      V_FLOAT vx = _v[i];
      V_FLOAT vy = _v[i + _nmax];
      V_FLOAT vz = _v[i + 2 * _nmax];
      vx *= factor;
      vy *= factor;
      vz *= factor;

      _v[i] = vx > 0 ? min(vx, limit) : max(vx, -limit);
      _v[i + _nmax] = vy > 0 ? min(vy, limit) : max(vy, -limit);
      _v[i + 2 * _nmax] = vz > 0 ? min(vz, limit) : max(vz, -limit);
    }
}

