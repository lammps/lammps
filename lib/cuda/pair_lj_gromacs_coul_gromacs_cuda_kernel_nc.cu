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

__device__ inline F_CFLOAT CoulGromacsCuda_Eval(const F_CFLOAT &rsq, const int ij_type, F_CFLOAT &factor_coul, int &eflag, ENERGY_CFLOAT &ecoul, F_CFLOAT qij)
{
  if(qij != F_F(0.0)) {
    F_CFLOAT ecoul_tmp;
    F_CFLOAT forcecoul = _RSQRT_(rsq);

    if(eflag) ecoul_tmp = forcecoul - _coulsw5;

    if(rsq > _cut_coul_inner_global * _cut_coul_inner_global) {
      const F_CFLOAT r = F_F(1.0) / forcecoul;
      const F_CFLOAT tc = r - _cut_coul_inner_global;
      forcecoul += r * tc * tc * (_coulsw1 + _coulsw2 * tc);

      if(eflag)  ecoul_tmp -= tc * tc * tc * (_coulsw1 * (F_F(1.0) / F_F(3.0)) + _coulsw2 * tc * (F_F(1.0) / F_F(4.0)));
    }

    F_CFLOAT qprod = _qqrd2e * qij * factor_coul;
    forcecoul *= qprod;

    if(eflag) {
      ecoul += ecoul_tmp * qprod;
    }

    return forcecoul * (F_F(1.0) / rsq);
  }

  return F_F(0.0);
}
