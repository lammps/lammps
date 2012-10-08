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
__device__ inline F_FLOAT PairLJCharmmCuda_Eval(const F_FLOAT &rsq, const int ij_type, F_FLOAT &factor_lj, int &eflag, ENERGY_FLOAT &evdwl)
{
  const F_FLOAT r2inv = F_F(1.0) / rsq;
  const F_FLOAT r6inv = r2inv * r2inv * r2inv;
  F_FLOAT forcelj = r6inv * (_lj1[ij_type] * r6inv - _lj2[ij_type]);
  F_FLOAT philj, switch1;

  if(rsq > _cut_innersq_global) {
    switch1 = (_cutsq_global - rsq) * (_cutsq_global - rsq) *
              (_cutsq_global + F_F(2.0) * rsq - F_F(3.0) * _cut_innersq_global) * _denom_lj_inv;
    const F_FLOAT switch2 = F_F(12.0) * rsq * (_cutsq_global - rsq) *
                            (rsq - _cut_innersq_global) * _denom_lj_inv;
    philj = r6inv * (_lj3[ij_type] * r6inv - _lj4[ij_type]);
    forcelj = forcelj * switch1 + philj * switch2;
  }

  if(eflag) {
    ENERGY_FLOAT evdwl_tmp = factor_lj;

    if(rsq > _cut_innersq_global) {
      evdwl_tmp *= philj * switch1;
    } else
      evdwl_tmp *= r6inv * (_lj3[ij_type] * r6inv - _lj4[ij_type]);

    evdwl += evdwl_tmp;
  }

  return factor_lj * forcelj * r2inv;
}

__device__ inline F_FLOAT CoulCharmmCuda_Eval(const F_FLOAT &rsq, F_FLOAT &factor_coul, int &eflag, ENERGY_FLOAT &ecoul, F_FLOAT qij)
{
  F_FLOAT forcecoul;
  ENERGY_FLOAT ecoul_tmp = forcecoul = _qqrd2e * qij * _RSQRT_(rsq) * factor_coul;

  if(rsq > _cut_coul_innersq_global) {
    const F_FLOAT switch1 = (_cut_coulsq_global - rsq) * (_cut_coulsq_global - rsq) *
                            (_cut_coulsq_global + F_F(2.0) * rsq - F_F(3.0) * _cut_coul_innersq_global) * _denom_coul_inv;
    ecoul_tmp *= switch1;
    const F_FLOAT switch2 = F_F(12.0) * rsq * (_cut_coulsq_global - rsq) *
                            (rsq - _cut_coul_innersq_global) * _denom_coul_inv;
    forcecoul *= switch1 + switch2;
  }

  if(eflag) {
    ecoul += ecoul_tmp * factor_coul;
  }

  return forcecoul * (F_F(1.0) / rsq);
}

