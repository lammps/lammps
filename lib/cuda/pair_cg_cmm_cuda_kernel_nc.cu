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

__device__ inline F_FLOAT PairCGCMMCuda_Eval(const F_FLOAT &rsq, const int ij_type, F_FLOAT &factor_lj, int &eflag, ENERGY_FLOAT &evdwl) //0.11 of 0.4
{
  const F_FLOAT r2inv = F_F(1.0) / rsq;
  const int cg_type = _cg_type[ij_type];
  const F_FLOAT r4inv = r2inv * r2inv;
  const F_FLOAT rNinv_first = cg_type != CG_LJ9_6 ? r4inv : _RSQRT_(rsq);
  const F_FLOAT rNinv_second = cg_type != CG_LJ12_4 ? -r2inv : -F_F(1.0);
  const F_FLOAT forcelj = r4inv * (_lj1[ij_type] * r4inv * rNinv_first + _lj2[ij_type] * rNinv_second);

  if(eflag) evdwl += factor_lj * (r4inv * (_lj3[ij_type] * r4inv * rNinv_first + _lj4[ij_type] * rNinv_second) - _offset[ij_type]);

  return factor_lj * forcelj * r2inv;
}

/*__device__ inline F_FLOAT PairCGCMMCuda_Eval(const F_FLOAT& rsq,const int ij_type,F_FLOAT& factor_lj,int& eflag, ENERGY_FLOAT& evdwl)
{
	const int cg_type = tex1Dfetch(_coeff5_gm_tex,ij_type);
	const F_FLOAT r2inv = F_F(1.0)/rsq;
	const F_FLOAT r4inv = r2inv*r2inv;
	const F_FLOAT rNinv_first = cg_type!=CG_LJ9_6?r4inv:_RSQRT_(rsq);
	const F_FLOAT rNinv_second = cg_type!=CG_LJ12_4?r2inv:F_F(1.0);
	const F_FLOAT forcelj = r4inv * (tex1Dfetch(_coeff1_gm_tex,ij_type)*r4inv*rNinv_first - tex1Dfetch(_coeff2_gm_tex,ij_type)*rNinv_second);

    if(eflag) evdwl += factor_lj*(r4inv*(tex1Dfetch(_coeff3_gm_tex,ij_type)*r4inv*rNinv_first-tex1Dfetch(_coeff4_gm_tex,ij_type)*rNinv_second));
	return factor_lj*forcelj*r2inv;
}*/
