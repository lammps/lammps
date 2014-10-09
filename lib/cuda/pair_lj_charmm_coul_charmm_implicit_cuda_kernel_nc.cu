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

__device__ inline F_CFLOAT CoulCharmmImplicitCuda_Eval(const F_CFLOAT &rsq, F_CFLOAT &factor_coul, int &eflag, ENERGY_CFLOAT &ecoul, F_CFLOAT qij)
{
  F_CFLOAT forcecoul;
  ENERGY_CFLOAT ecoul_tmp = forcecoul = _qqrd2e * qij * (F_F(1.0) / rsq) * factor_coul;

  if(rsq > _cut_coul_innersq_global) {
    const F_CFLOAT switch1 = (_cut_coulsq_global - rsq) * (_cut_coulsq_global - rsq) *
                            (_cut_coulsq_global + F_F(2.0) * rsq - F_F(3.0) * _cut_coul_innersq_global) * _denom_coul_inv;
    ecoul_tmp *= switch1;
    const F_CFLOAT switch2 = F_F(12.0) * rsq * (_cut_coulsq_global - rsq) *
                            (rsq - _cut_coul_innersq_global) * _denom_coul_inv;
    forcecoul *= (switch1 + switch2);
  }

  if(eflag) {
    ecoul += ecoul_tmp * factor_coul;
  }

  return F_F(2.0) * forcecoul * (F_F(1.0) / rsq);
}

