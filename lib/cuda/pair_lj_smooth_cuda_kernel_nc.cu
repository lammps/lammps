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

__device__ inline F_CFLOAT PairLJSmoothCuda_Eval(const F_CFLOAT &rsq, const int ij_type, F_CFLOAT &factor_lj, int &eflag, ENERGY_CFLOAT &evdwl)
{
  F_CFLOAT fskin, t, tsq, forcelj;
  const F_CFLOAT r2inv = F_F(1.0) / rsq;
  const F_CFLOAT r = _RSQRT_(r2inv);
  const F_CFLOAT r6inv = r2inv * r2inv * r2inv;


  X_CFLOAT cut_lj_innersq = (_cut_innersq_global > X_F(0.0) ? _cut_innersq_global : _cut_innersq[ij_type]);

  if(rsq < cut_lj_innersq) {
    forcelj = r6inv * (_lj1[ij_type] * r6inv - _lj2[ij_type]);
  } else {
    t = r - _SQRT_(cut_lj_innersq);
    tsq = t * t;
    fskin = _ljsw1[ij_type] +  _ljsw2[ij_type] * t +
            _ljsw3[ij_type] * tsq +  _ljsw4[ij_type] * tsq * t;
    forcelj = fskin * r;

  }

  if(eflag) {
    ENERGY_CFLOAT evdwl_tmp;

    if(rsq < cut_lj_innersq) {
      evdwl_tmp = r6inv * (_lj3[ij_type] * r6inv - _lj4[ij_type]) -
                  _offset[ij_type];
    } else {
      evdwl_tmp = _ljsw0[ij_type] - _ljsw1[ij_type] * t -
                  _ljsw2[ij_type] * tsq * F_F(0.5) - _ljsw3[ij_type] * tsq * t * (F_F(1.0) / F_F(3.0)) -
                  _ljsw4[ij_type] * tsq * tsq * (F_F(1.0) / F_F(4.0)) - _offset[ij_type];
    }

    evdwl += evdwl_tmp * factor_lj;
  }

  return factor_lj * forcelj * r2inv;
}
