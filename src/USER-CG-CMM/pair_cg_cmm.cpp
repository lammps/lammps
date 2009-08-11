/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   CMM coarse grained MD potentials. Plain version w/o charges.
   Contributing author: Axel Kohlmeyer <akohlmey@gmail.com>
------------------------------------------------------------------------- */

#include "pair_cg_cmm.h"

using namespace LAMMPS_NS;
 
/* ---------------------------------------------------------------------- */

PairCGCMM::PairCGCMM(LAMMPS *lmp) : PairCMMCommon(lmp)
{
  respa_enable = 0;
  single_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairCGCMM::~PairCGCMM()
{
  /* empty */ ;
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- *
 * the real compute work is done in the PairCMMCommon::eval_XXX<>() templates
 * in the common PairCG class. Through using templates we can have one 
 * implementation for all CG varieties _and_ gain speed through having 
 * the compiler optimize away conditionals within the innerloops that 
 * can be predetermined outside the loop through instantiation of the
 * different combination of template flags.
 * ---------------------------------------------------------------------- */

void PairCGCMM::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else {
    evflag = vflag_fdotr = 0;
  }

  if (evflag) {
    if (eflag) {
      if (force->newton_pair) {
        return eval_verlet<1,1,1,CG_COUL_NONE>();
      } else {
        return eval_verlet<1,1,0,CG_COUL_NONE>();
      }
    } else {
      if (force->newton_pair) {
        return eval_verlet<1,0,1,CG_COUL_NONE>();
      } else {
        return eval_verlet<1,0,0,CG_COUL_NONE>();
      }
    }
  } else {
    if (force->newton_pair) {
      return eval_verlet<0,0,1,CG_COUL_NONE>();
    } else {
      return eval_verlet<0,0,0,CG_COUL_NONE>();
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairCGCMM::compute_inner()
{
  if (force->newton_pair) {
    return eval_inner<1,CG_COUL_NONE>();
  } else {
    return eval_inner<0,CG_COUL_NONE>();
  }
}

/* ---------------------------------------------------------------------- */

void PairCGCMM::compute_middle()
{
  if (force->newton_pair) {
    return eval_middle<1,CG_COUL_NONE>();
  } else {
    return eval_middle<0,CG_COUL_NONE>();
  }
}

/* ---------------------------------------------------------------------- */

void PairCGCMM::compute_outer(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else {
    evflag = 0;
  }

  if (evflag) {
    if (eflag) {
      if (vflag) {
        if (force->newton_pair) {
          return eval_outer<1,1,1,1,CG_COUL_NONE>();
        } else {
          return eval_outer<1,1,1,0,CG_COUL_NONE>();
        }
      } else {
        if (force->newton_pair) {
          return eval_outer<1,1,0,1,CG_COUL_NONE>();
        } else {
          return eval_outer<1,1,0,0,CG_COUL_NONE>();
        }
      }
    } else {
      if (vflag) {
        if (force->newton_pair) {
          return eval_outer<1,0,1,1,CG_COUL_NONE>();
        } else {
          return eval_outer<1,0,1,0,CG_COUL_NONE>();
        }
      } else {
        if (force->newton_pair) {
          return eval_outer<1,0,0,1,CG_COUL_NONE>();
        } else {
          return eval_outer<1,0,0,0,CG_COUL_NONE>();
        }
      }
    }
  } else {
    if (force->newton_pair) {
      return eval_outer<0,0,0,1,CG_COUL_NONE>();
    } else {
      return eval_outer<0,0,0,0,CG_COUL_NONE>();
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairCGCMM::write_restart(FILE *fp) 
{
  write_restart_settings(fp);
  PairCMMCommon::write_restart(fp);
}

/* ---------------------------------------------------------------------- */

void PairCGCMM::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();
  PairCMMCommon::read_restart(fp);
}

/* ---------------------------------------------------------------------- */

double PairCGCMM::single(int i, int j, int itype, int jtype, double rsq,
		       double factor_coul, double factor_lj, double &fforce)
{
  return eval_single(CG_COUL_NONE,i,j,itype,jtype,rsq,factor_coul,factor_lj,fforce);
}
