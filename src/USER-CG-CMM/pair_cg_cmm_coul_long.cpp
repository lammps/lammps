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
   CMM coarse grained MD potentials. Coulomb with k-space version.
   Contributing author: Axel Kohlmeyer <akohlmey@gmail.com>
------------------------------------------------------------------------- */

#include "string.h"
#include "pair_cg_cmm_coul_long.h"
#include "memory.h"
#include "atom.h"
#include "force.h"
#include "kspace.h"

using namespace LAMMPS_NS;

#define EWALD_F   1.12837917

/* ---------------------------------------------------------------------- */

PairCGCMMCoulLong::PairCGCMMCoulLong(LAMMPS *lmp) : PairCMMCommon(lmp)
{
  respa_enable = 0;
  single_enable = 0;
  ewaldflag = pppmflag = 1;
}

/* ---------------------------------------------------------------------- */

PairCGCMMCoulLong::~PairCGCMMCoulLong()
{
  if (allocated_coul) {
    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(cut_coul);
    memory->destroy(cut_coulsq);
    allocated_coul=0;
  }
  if (ftable) free_tables();
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulLong::allocate()
{
  PairCMMCommon::allocate();
  allocated_coul = 1;

  int n = atom->ntypes;

  memory->create(cut_lj,n+1,n+1,"paircg:cut_lj");
  memory->create(cut_ljsq,n+1,n+1,"paircg:cut_ljsq");
  memory->create(cut_coul,n+1,n+1,"paircg:cut_coul");
  memory->create(cut_coulsq,n+1,n+1,"paircg:cut_coulsq");
}

/* ----------------------------------------------------------------------
   free memory for tables used in pair computations
------------------------------------------------------------------------- */

void PairCGCMMCoulLong::free_tables()
{
  memory->destroy(rtable);
  memory->destroy(drtable);
  memory->destroy(ftable);
  memory->destroy(dftable);
  memory->destroy(ctable);
  memory->destroy(dctable);
  memory->destroy(etable);
  memory->destroy(detable);
  memory->destroy(vtable);
  memory->destroy(dvtable);
  memory->destroy(ptable);
  memory->destroy(dptable);
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulLong::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style cg/cut/coul/long requires atom attribute q");

  PairCMMCommon::init_style();

  // set rRESPA cutoffs

  if (strstr(update->integrate_style,"respa") &&
      ((Respa *) update->integrate)->level_inner >= 0)
    cut_respa = ((Respa *) update->integrate)->cutoff;
  else cut_respa = NULL;

  // ensure use of KSpace long-range solver, set g_ewald

  if (force->kspace == NULL)
    error->all(FLERR,"Pair style is incompatible with KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables

  if (ncoultablebits) init_tables();
}

/* ---------------------------------------------------------------------- */

double PairCGCMMCoulLong::init_one(int i, int j)
{
  double mycut = PairCMMCommon::init_one(i,j);

  // check interior rRESPA cutoff

  if (cut_respa && MIN(cut_lj[i][j],cut_coul_global) < cut_respa[3])
    error->all(FLERR,"Pair cutoff < Respa interior cutoff");

  return mycut;
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulLong::init_tables()
{
  int masklo,maskhi;
  double r,grij,expm2,derfc,rsw;
  double qqrd2e = force->qqrd2e;

  tabinnersq = tabinner*tabinner;
  init_bitmap(tabinner,cut_coul_global,ncoultablebits,
              masklo,maskhi,ncoulmask,ncoulshiftbits);

  int ntable = 1;
  for (int i = 0; i < ncoultablebits; i++) ntable *= 2;

  // linear lookup tables of length N = 2^ncoultablebits
  // stored value = value at lower edge of bin
  // d values = delta from lower edge to upper edge of bin

  if (ftable) free_tables();

  memory->create(rtable,ntable,"pair:rtable");
  memory->create(ftable,ntable,"pair:ftable");
  memory->create(ctable,ntable,"pair:ctable");
  memory->create(etable,ntable,"pair:etable");
  memory->create(drtable,ntable,"pair:drtable");
  memory->create(dftable,ntable,"pair:dftable");
  memory->create(dctable,ntable,"pair:dctable");
  memory->create(detable,ntable,"pair:detable");

  if (cut_respa == NULL) {
    vtable = ptable = dvtable = dptable = NULL;
  } else {
    memory->create(vtable,ntable,"pair:vtable");
    memory->create(ptable,ntable,"pair:ptable");
    memory->create(dvtable,ntable,"pair:dvtable");
    memory->create(dptable,ntable,"pair:dptable");
  }

  union_int_float_t rsq_lookup;
  union_int_float_t minrsq_lookup;
  int itablemin;
  minrsq_lookup.i = 0 << ncoulshiftbits;
  minrsq_lookup.i |= maskhi;
  for (int i = 0; i < ntable; i++) {
    rsq_lookup.i = i << ncoulshiftbits;
    rsq_lookup.i |= masklo;
    if (rsq_lookup.f < tabinnersq) {
      rsq_lookup.i = i << ncoulshiftbits;
      rsq_lookup.i |= maskhi;
    }
    r = sqrtf(rsq_lookup.f);
    grij = g_ewald * r;
    expm2 = exp(-grij*grij);
    derfc = erfc(grij);
    if (cut_respa == NULL) {
      rtable[i] = rsq_lookup.f;
      ftable[i] = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
      ctable[i] = qqrd2e/r;
      etable[i] = qqrd2e/r * derfc;
    } else {
      rtable[i] = rsq_lookup.f;
      ftable[i] = qqrd2e/r * (derfc + EWALD_F*grij*expm2 - 1.0);
      ctable[i] = 0.0;
      etable[i] = qqrd2e/r * derfc;
      ptable[i] = qqrd2e/r;
      vtable[i] = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
      if (rsq_lookup.f > cut_respa[2]*cut_respa[2]) {
        if (rsq_lookup.f < cut_respa[3]*cut_respa[3]) {
          rsw = (r - cut_respa[2])/(cut_respa[3] - cut_respa[2]);
          ftable[i] += qqrd2e/r * rsw*rsw*(3.0 - 2.0*rsw);
          ctable[i] = qqrd2e/r * rsw*rsw*(3.0 - 2.0*rsw);
        } else {
          ftable[i] = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
          ctable[i] = qqrd2e/r;
        }
      }
    }
    minrsq_lookup.f = MIN(minrsq_lookup.f,rsq_lookup.f);
  }
  tabinnersq = minrsq_lookup.f;

  int ntablem1 = ntable - 1;

  for (int i = 0; i < ntablem1; i++) {
    drtable[i] = 1.0/(rtable[i+1] - rtable[i]);
    dftable[i] = ftable[i+1] - ftable[i];
    dctable[i] = ctable[i+1] - ctable[i];
    detable[i] = etable[i+1] - etable[i];
  }

  if (cut_respa) {
    for (int i = 0; i < ntablem1; i++) {
      dvtable[i] = vtable[i+1] - vtable[i];
      dptable[i] = ptable[i+1] - ptable[i];
    }
  }

  // get the delta values for the last table entries
  // tables are connected periodically between 0 and ntablem1

  drtable[ntablem1] = 1.0/(rtable[0] - rtable[ntablem1]);
  dftable[ntablem1] = ftable[0] - ftable[ntablem1];
  dctable[ntablem1] = ctable[0] - ctable[ntablem1];
  detable[ntablem1] = etable[0] - etable[ntablem1];
  if (cut_respa) {
    dvtable[ntablem1] = vtable[0] - vtable[ntablem1];
    dptable[ntablem1] = ptable[0] - ptable[ntablem1];
  }

  // get the correct delta values at itablemax
  // smallest r is in bin itablemin
  // largest r is in bin itablemax, which is itablemin-1,
  //   or ntablem1 if itablemin=0
  // deltas at itablemax only needed if corresponding rsq < cut*cut
  // if so, compute deltas between rsq and cut*cut

  double f_tmp,c_tmp,e_tmp,p_tmp,v_tmp;
  itablemin = minrsq_lookup.i & ncoulmask;
  itablemin >>= ncoulshiftbits;
  int itablemax = itablemin - 1;
  if (itablemin == 0) itablemax = ntablem1;
  rsq_lookup.i = itablemax << ncoulshiftbits;
  rsq_lookup.i |= maskhi;
  if (rsq_lookup.f < cut_coulsq_global) {
    rsq_lookup.f = cut_coulsq_global;
    r = sqrtf(rsq_lookup.f);
    grij = g_ewald * r;
    expm2 = exp(-grij*grij);
    derfc = erfc(grij);

    if (cut_respa == NULL) {
      f_tmp = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
      c_tmp = qqrd2e/r;
      e_tmp = qqrd2e/r * derfc;
    } else {
      f_tmp = qqrd2e/r * (derfc + EWALD_F*grij*expm2 - 1.0);
      c_tmp = 0.0;
      e_tmp = qqrd2e/r * derfc;
      p_tmp = qqrd2e/r;
      v_tmp = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
      if (rsq_lookup.f > cut_respa[2]*cut_respa[2]) {
        if (rsq_lookup.f < cut_respa[3]*cut_respa[3]) {
          rsw = (r - cut_respa[2])/(cut_respa[3] - cut_respa[2]);
          f_tmp += qqrd2e/r * rsw*rsw*(3.0 - 2.0*rsw);
          c_tmp = qqrd2e/r * rsw*rsw*(3.0 - 2.0*rsw);
        } else {
          f_tmp = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
          c_tmp = qqrd2e/r;
        }
      }
    }

    drtable[itablemax] = 1.0/(rsq_lookup.f - rtable[itablemax]);
    dftable[itablemax] = f_tmp - ftable[itablemax];
    dctable[itablemax] = c_tmp - ctable[itablemax];
    detable[itablemax] = e_tmp - etable[itablemax];
    if (cut_respa) {
      dvtable[itablemax] = v_tmp - vtable[itablemax];
      dptable[itablemax] = p_tmp - ptable[itablemax];
    }
  }

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

void PairCGCMMCoulLong::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else {
    evflag = vflag_fdotr = 0;
  }

  if (evflag) {
    if (eflag) {
      if (force->newton_pair) {
        return eval_verlet<1,1,1,CG_COUL_LONG>();
      } else {
        return eval_verlet<1,1,0,CG_COUL_LONG>();
      }
    } else {
      if (force->newton_pair) {
        return eval_verlet<1,0,1,CG_COUL_LONG>();
      } else {
        return eval_verlet<1,0,0,CG_COUL_LONG>();
      }
    }
  } else {
    if (force->newton_pair) {
      return eval_verlet<0,0,1,CG_COUL_LONG>();
    } else {
      return eval_verlet<0,0,0,CG_COUL_LONG>();
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulLong::compute_inner()
{
  if (force->newton_pair) {
    return eval_inner<1,CG_COUL_LONG>();
  } else {
    return eval_inner<0,CG_COUL_LONG>();
  }
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulLong::compute_middle()
{
  if (force->newton_pair) {
    return eval_middle<1,CG_COUL_LONG>();
  } else {
    return eval_middle<0,CG_COUL_LONG>();
  }
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulLong::compute_outer(int eflag, int vflag)
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
          return eval_outer<1,1,1,1,CG_COUL_LONG>();
        } else {
          return eval_outer<1,1,1,0,CG_COUL_LONG>();
        }
      } else {
        if (force->newton_pair) {
          return eval_outer<1,1,0,1,CG_COUL_LONG>();
        } else {
          return eval_outer<1,1,0,0,CG_COUL_LONG>();
        }
      }
    } else {
      if (vflag) {
        if (force->newton_pair) {
          return eval_outer<1,0,1,1,CG_COUL_LONG>();
        } else {
          return eval_outer<1,0,1,0,CG_COUL_LONG>();
        }
      } else {
        if (force->newton_pair) {
          return eval_outer<1,0,0,1,CG_COUL_LONG>();
        } else {
          return eval_outer<1,0,0,0,CG_COUL_LONG>();
        }
      }
    }
  } else {
    if (force->newton_pair) {
      return eval_outer<0,0,0,1,CG_COUL_LONG>();
    } else {
      return eval_outer<0,0,0,0,CG_COUL_LONG>();
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulLong::write_restart(FILE *fp)
{
  write_restart_settings(fp);
  fwrite(&ncoultablebits,sizeof(int),1,fp);
  fwrite(&tabinner,sizeof(double),1,fp);

  PairCMMCommon::write_restart(fp);
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulLong::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  if (comm->me == 0) {
    fread(&ncoultablebits,sizeof(int),1,fp);
    fread(&tabinner,sizeof(double),1,fp);
  }
  MPI_Bcast(&ncoultablebits,1,MPI_INT,0,world);
  MPI_Bcast(&tabinner,1,MPI_DOUBLE,0,world);

  allocate();
  PairCMMCommon::read_restart(fp);
}

/* ---------------------------------------------------------------------- */

double PairCGCMMCoulLong::memory_usage()
{
  double bytes=PairCMMCommon::memory_usage();

  int n = atom->ntypes;

  // cut_coul/cut_coulsq/cut_ljsq
  bytes += (n+1)*(n+1)*sizeof(double)*4;

  return bytes;
}

/* ---------------------------------------------------------------------- */

double PairCGCMMCoulLong::single(int i, int j, int itype, int jtype, double rsq,
                                 double factor_coul, double factor_lj, double &fforce)
{
  return eval_single(CG_COUL_LONG,i,j,itype,jtype,rsq,factor_coul,factor_lj,fforce);
}

/* ---------------------------------------------------------------------- */

void *PairCGCMMCoulLong::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul_global;
  return NULL;
}
