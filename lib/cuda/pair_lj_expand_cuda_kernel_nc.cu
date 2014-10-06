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
__device__ inline F_CFLOAT PairLJExpandCuda_Eval(const F_CFLOAT &rsq, const int ij_type, F_CFLOAT &factor_lj, int &eflag, ENERGY_CFLOAT &evdwl)
{
  const F_CFLOAT r = _SQRT_(rsq);
  const F_CFLOAT rshift = r - _shift[ij_type];
  const F_CFLOAT rshiftsq = rshift * rshift;
  const F_CFLOAT r2inv = F_F(1.0) / rshiftsq;
  const F_CFLOAT r6inv = r2inv * r2inv * r2inv;
  const F_CFLOAT forcelj = r6inv * (_lj1[ij_type] * r6inv - _lj2[ij_type]);

  if(eflag) evdwl += factor_lj * (r6inv * (_lj3[ij_type] * r6inv - _lj4[ij_type]) - _offset[ij_type]);

  return factor_lj * forcelj * (F_F(1.0) / rshift) * (F_F(1.0) / r);
}
