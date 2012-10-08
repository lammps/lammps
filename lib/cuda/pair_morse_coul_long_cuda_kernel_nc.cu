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
__device__ inline F_FLOAT PairMorseR6Cuda_Eval(const F_FLOAT &rsq, const int ij_type, F_FLOAT &factor_lj, int &eflag, ENERGY_FLOAT &evdwl)
{
  const F_FLOAT r2inv = F_F(1.0) / rsq;
  const F_FLOAT r = _SQRT_(rsq);
  const F_FLOAT r4inv = r2inv * r2inv;
  const F_FLOAT dr = r - _r0[ij_type];
  const F_FLOAT dexp = _EXP_(-_alpha[ij_type] * dr);

  if(eflag) evdwl += factor_lj * (_d0[ij_type] * (dexp * dexp - F_F(2.0) * dexp) + _c0[ij_type] * r4inv * r4inv * r4inv
                                    - _offset[ij_type]);

  return factor_lj * (_morse1[ij_type] * (dexp * dexp - dexp) * (F_F(1.0) / r) - F_F(12.0) * _c0[ij_type] * r4inv * r4inv * r4inv * r2inv);
}
