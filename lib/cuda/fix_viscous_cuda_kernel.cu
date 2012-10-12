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

__global__ void Cuda_FixViscousCuda_PostForce_Kernel(int groupbit, F_FLOAT* gamma)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < _nlocal)
    if(_mask[i] & groupbit) {
      F_FLOAT drag = gamma[_type[i]];
      _f[i] -= drag * _v[i];
      _f[i + 1 * _nmax] -= drag * _v[i + 1 * _nmax];
      _f[i + 2 * _nmax] -= drag * _v[i + 2 * _nmax];
    }
}
