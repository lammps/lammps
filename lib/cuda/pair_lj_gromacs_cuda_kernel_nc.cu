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

__device__ inline F_FLOAT PairLJGromacsCuda_Eval(const F_FLOAT &rsq, const int ij_type, F_FLOAT &factor_lj, int &eflag, ENERGY_FLOAT &evdwl)
{
  F_FLOAT tlj;
  const F_FLOAT r2inv = F_F(1.0) / rsq;
  const F_FLOAT r = _RSQRT_(r2inv);
  const F_FLOAT r6inv = r2inv * r2inv * r2inv;
  F_FLOAT	forcelj = r6inv * (_lj1[ij_type] * r6inv - _lj2[ij_type]);
  const X_FLOAT cut_lj_innersq = (_cut_innersq_global > X_F(0.0) ? _cut_innersq_global : _cut_innersq[ij_type]);

  if(rsq > cut_lj_innersq) {
    tlj = r - _SQRT_(cut_lj_innersq);
    forcelj += r * tlj * tlj * (_ljsw1[ij_type] + _ljsw2[ij_type] * tlj);
  }

  if(eflag) {
    ENERGY_FLOAT evdwl_tmp = r6inv * (_lj3[ij_type] * r6inv - _lj4[ij_type]);

    if(rsq > cut_lj_innersq) {
      evdwl_tmp += tlj * tlj * tlj *
                   (_ljsw3[ij_type] + _ljsw4[ij_type] * tlj) + _ljsw5[ij_type];;
    }

    evdwl += evdwl_tmp * factor_lj;
  }

  return factor_lj * forcelj * r2inv;
}
