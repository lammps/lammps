/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   CMM coarse grained MD potentials. Coulomb with cutoff version.
   Contributing author: Axel Kohlmeyer <akohlmey@gmail.com>
------------------------------------------------------------------------- */

#include "string.h"
#include "pair_cg_cmm_coul_cut.h"
#include "memory.h"
#include "atom.h"

using namespace LAMMPS_NS;

PairCGCMMCoulCut::PairCGCMMCoulCut(LAMMPS *lmp) : PairCMMCommon(lmp)
{
  respa_enable = 0;
  single_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairCGCMMCoulCut::~PairCGCMMCoulCut()
{
  if (allocated_coul) {
    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(cut_coul);
    memory->destroy(cut_coulsq);
    allocated_coul=0;
  }
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulCut::allocate()
{
  PairCMMCommon::allocate();
  allocated_coul = 1;

  int n = atom->ntypes;

  memory->create(cut_lj,n+1,n+1,"paircg:cut_lj");
  memory->create(cut_ljsq,n+1,n+1,"paircg:cut_ljsq");
  memory->create(cut_coul,n+1,n+1,"paircg:cut_coul");
  memory->create(cut_coulsq,n+1,n+1,"paircg:cut_coulsq");
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulCut::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style cg/cut/coul/cut requires atom attribute q");

  PairCMMCommon::init_style();

  // set rRESPA cutoffs

  if (strstr(update->integrate_style,"respa") &&
      ((Respa *) update->integrate)->level_inner >= 0)
    cut_respa = ((Respa *) update->integrate)->cutoff;
  else cut_respa = NULL;
}

/* ---------------------------------------------------------------------- */

double PairCGCMMCoulCut::init_one(int i, int j)
{
  double mycut = PairCMMCommon::init_one(i,j);

  // check interior rRESPA cutoff

  if (cut_respa && MIN(cut_lj[i][j],cut_coul[i][j]) < cut_respa[3])
    error->all(FLERR,"Pair cutoff < Respa interior cutoff");

  return mycut;
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

void PairCGCMMCoulCut::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else {
    evflag = vflag_fdotr = 0;
  }

  if (evflag) {
    if (eflag) {
      if (force->newton_pair) {
        return eval_verlet<1,1,1,CG_COUL_CUT>();
      } else {
        return eval_verlet<1,1,0,CG_COUL_CUT>();
      }
    } else {
      if (force->newton_pair) {
        return eval_verlet<1,0,1,CG_COUL_CUT>();
      } else {
        return eval_verlet<1,0,0,CG_COUL_CUT>();
      }
    }
  } else {
    if (force->newton_pair) {
      return eval_verlet<0,0,1,CG_COUL_CUT>();
    } else {
      return eval_verlet<0,0,0,CG_COUL_CUT>();
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulCut::compute_inner()
{
  if (force->newton_pair) {
    return eval_inner<1,CG_COUL_CUT>();
  } else {
    return eval_inner<0,CG_COUL_CUT>();
  }
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulCut::compute_middle()
{
  if (force->newton_pair) {
    return eval_middle<1,CG_COUL_CUT>();
  } else {
    return eval_middle<0,CG_COUL_CUT>();
  }
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulCut::compute_outer(int eflag, int vflag)
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
          return eval_outer<1,1,1,1,CG_COUL_CUT>();
        } else {
          return eval_outer<1,1,1,0,CG_COUL_CUT>();
        }
      } else {
        if (force->newton_pair) {
          return eval_outer<1,1,0,1,CG_COUL_CUT>();
        } else {
          return eval_outer<1,1,0,0,CG_COUL_CUT>();
        }
      }
    } else {
      if (vflag) {
        if (force->newton_pair) {
          return eval_outer<1,0,1,1,CG_COUL_CUT>();
        } else {
          return eval_outer<1,0,1,0,CG_COUL_CUT>();
        }
      } else {
        if (force->newton_pair) {
          return eval_outer<1,0,0,1,CG_COUL_CUT>();
        } else {
          return eval_outer<1,0,0,0,CG_COUL_CUT>();
        }
      }
    }
  } else {
    if (force->newton_pair) {
      return eval_outer<0,0,0,1,CG_COUL_CUT>();
    } else {
      return eval_outer<0,0,0,0,CG_COUL_CUT>();
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulCut::write_restart(FILE *fp) 
{
  write_restart_settings(fp);
  PairCMMCommon::write_restart(fp);
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulCut::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();
  PairCMMCommon::read_restart(fp);
}

/* ---------------------------------------------------------------------- */

double PairCGCMMCoulCut::memory_usage()
{
  double bytes=PairCMMCommon::memory_usage();
  
  int n = atom->ntypes;

  // cut_coul/cut_coulsq/cut_lj/cut_ljsq;
  bytes += (n+1)*(n+1)*sizeof(double)*4; 
  
  return bytes;
}

/* ---------------------------------------------------------------------- */

double PairCGCMMCoulCut::single(int i, int j, int itype, int jtype, double rsq,
		       double factor_coul, double factor_lj, double &fforce)
{
  return eval_single(CG_COUL_CUT,i,j,itype,jtype,rsq,factor_coul,factor_lj,fforce);
}
