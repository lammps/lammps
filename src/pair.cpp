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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "ctype.h"
#include "float.h"
#include "limits.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "update.h"
#include "accelerator_cuda.h"
#include "suffix.h"
#include "atom_masks.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define EWALD_F 1.12837917

enum{RLINEAR,RSQ,BMP};

/* ---------------------------------------------------------------------- */

Pair::Pair(LAMMPS *lmp) : Pointers(lmp)
{
  THIRD = 1.0/3.0;

  eng_vdwl = eng_coul = 0.0;

  comm_forward = comm_reverse = comm_reverse_off = 0;

  single_enable = 1;
  restartinfo = 1;
  respa_enable = 0;
  one_coeff = 0;
  no_virial_fdotr_compute = 0;
  writedata = 0;
  ghostneigh = 0;

  nextra = 0;
  pvector = NULL;
  single_extra = 0;
  svector = NULL;

  ewaldflag = pppmflag = msmflag = dispersionflag = tip4pflag = dipoleflag = 0;
  reinitflag = 1;

  // pair_modify settings

  compute_flag = 1;
  manybody_flag = 0;
  offset_flag = 0;
  mix_flag = GEOMETRIC;
  tail_flag = 0;
  etail = ptail = etail_ij = ptail_ij = 0.0;
  ncoultablebits = 12;
  ndisptablebits = 12;
  tabinner = sqrt(2.0);
  tabinner_disp = sqrt(2.0);

  allocated = 0;
  suffix_flag = Suffix::NONE;

  maxeatom = maxvatom = 0;
  eatom = NULL;
  vatom = NULL;

  // CUDA and KOKKOS per-fix data masks

  datamask = ALL_MASK;
  datamask_ext = ALL_MASK;

  execution_space = Host;
  datamask_read = ALL_MASK;
  datamask_modify = ALL_MASK;
}

/* ---------------------------------------------------------------------- */

Pair::~Pair()
{
  memory->destroy(eatom);
  memory->destroy(vatom);
}

/* ----------------------------------------------------------------------
   modify parameters of the pair style
   pair_hybrid has its own version of this routine for its sub-styles
------------------------------------------------------------------------- */

void Pair::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal pair_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mix") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_modify command");
      if (strcmp(arg[iarg+1],"geometric") == 0) mix_flag = GEOMETRIC;
      else if (strcmp(arg[iarg+1],"arithmetic") == 0) mix_flag = ARITHMETIC;
      else if (strcmp(arg[iarg+1],"sixthpower") == 0) mix_flag = SIXTHPOWER;
      else error->all(FLERR,"Illegal pair_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"shift") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) offset_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) offset_flag = 0;
      else error->all(FLERR,"Illegal pair_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"table") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_modify command");
      ncoultablebits = force->inumeric(FLERR,arg[iarg+1]);
      if (ncoultablebits > sizeof(float)*CHAR_BIT)
        error->all(FLERR,"Too many total bits for bitmapped lookup table");
      iarg += 2;
    } else if (strcmp(arg[iarg],"table/disp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_modify command");
      ndisptablebits = force->inumeric(FLERR,arg[iarg+1]);
      if (ndisptablebits > sizeof(float)*CHAR_BIT)
        error->all(FLERR,"Too many total bits for bitmapped lookup table");
      iarg += 2;
    } else if (strcmp(arg[iarg],"tabinner") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_modify command");
      tabinner = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"tabinner/disp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_modify command");
      tabinner_disp = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"tail") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) tail_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) tail_flag = 0;
      else error->all(FLERR,"Illegal pair_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"compute") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) compute_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) compute_flag = 0;
      else error->all(FLERR,"Illegal pair_modify command");
      iarg += 2;
    } else error->all(FLERR,"Illegal pair_modify command");
  }
}

/* ---------------------------------------------------------------------- */

void Pair::init()
{
  int i,j;

  if (offset_flag && tail_flag)
    error->all(FLERR,"Cannot have both pair_modify shift and tail set to yes");
  if (tail_flag && domain->dimension == 2)
    error->all(FLERR,"Cannot use pair tail corrections with 2d simulations");
  if (tail_flag && domain->nonperiodic && comm->me == 0)
    error->warning(FLERR,"Using pair tail corrections with nonperiodic system");

  // for manybody potentials
  // check if bonded exclusions could invalidate the neighbor list

  if (manybody_flag && atom->molecular) {
    int flag = 0;
    if (atom->nbonds > 0 && force->special_lj[1] == 0.0 && 
        force->special_coul[1] == 0.0) flag = 1;
    if (atom->nangles > 0 && force->special_lj[2] == 0.0 && 
        force->special_coul[2] == 0.0) flag = 1;
    if (atom->ndihedrals > 0 && force->special_lj[3] == 0.0 && 
        force->special_coul[3] == 0.0) flag = 1;
    if (flag && comm->me == 0)
      error->warning(FLERR,"Using a manybody potential with "
                     "bonds/angles/dihedrals and special_bond exclusions");
  }

  // I,I coeffs must be set
  // init_one() will check if I,J is set explicitly or inferred by mixing

  if (!allocated) error->all(FLERR,"All pair coeffs are not set");

  for (i = 1; i <= atom->ntypes; i++)
    if (setflag[i][i] == 0) error->all(FLERR,"All pair coeffs are not set");

  // style-specific initialization

  init_style();

  // call init_one() for each I,J
  // set cutsq for each I,J, used to neighbor
  // cutforce = max of all I,J cutoffs

  cutforce = 0.0;
  etail = ptail = 0.0;
  double cut;

  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      cut = init_one(i,j);
      cutsq[i][j] = cutsq[j][i] = cut*cut;
      cutforce = MAX(cutforce,cut);
      if (tail_flag) {
        etail += etail_ij;
        ptail += ptail_ij;
        if (i != j) {
          etail += etail_ij;
          ptail += ptail_ij;
        }
      }
    }
}

/* ----------------------------------------------------------------------
   reset all type-based params by invoking init_one() for each I,J
   called by fix adapt after it changes one or more params
------------------------------------------------------------------------- */

void Pair::reinit()
{
  // generalize this error message if reinit() is used by more than fix adapt

  if (!reinitflag)
    error->all(FLERR,"Fix adapt interface to this pair style not supported");


  etail = ptail = 0.0;

  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      init_one(i,j);
      if (tail_flag) {
        etail += etail_ij;
        ptail += ptail_ij;
        if (i != j) {
          etail += etail_ij;
          ptail += ptail_ij;
        }
      }
    }
}

/* ----------------------------------------------------------------------
   init specific to a pair style
   specific pair style can override this function
     if needs its own error checks
     if needs another kind of neighbor list
   request default neighbor list = half list
------------------------------------------------------------------------- */

void Pair::init_style()
{
  neighbor->request(this);
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   specific pair style can override this function
------------------------------------------------------------------------- */

void Pair::init_list(int which, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   setup Coulomb force tables used in compute routines
------------------------------------------------------------------------- */

void Pair::init_tables(double cut_coul, double *cut_respa)
{
  int masklo,maskhi;
  double r,grij,expm2,derfc,egamma,fgamma,rsw;
  double qqrd2e = force->qqrd2e;

  if (force->kspace == NULL)
    error->all(FLERR,"Pair style requres a KSpace style");
  double g_ewald = force->kspace->g_ewald;
  
  double cut_coulsq = cut_coul * cut_coul;
  
  tabinnersq = tabinner*tabinner;
  init_bitmap(tabinner,cut_coul,ncoultablebits,
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
    if (msmflag) {
      egamma = 1.0 - (r/cut_coul)*force->kspace->gamma(r/cut_coul);
      fgamma = 1.0 + (rsq_lookup.f/cut_coulsq)*
        force->kspace->dgamma(r/cut_coul);
    } else {
      grij = g_ewald * r;
      expm2 = exp(-grij*grij);
      derfc = erfc(grij);
    }
    if (cut_respa == NULL) {
      rtable[i] = rsq_lookup.f;
      ctable[i] = qqrd2e/r;
      if (msmflag) {
        ftable[i] = qqrd2e/r * fgamma;
        etable[i] = qqrd2e/r * egamma;
      } else {
        ftable[i] = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
        etable[i] = qqrd2e/r * derfc;
      }
    } else {
      rtable[i] = rsq_lookup.f;
      ctable[i] = 0.0;
      ptable[i] = qqrd2e/r;
      if (msmflag) {
        ftable[i] = qqrd2e/r * (fgamma - 1.0);
        etable[i] = qqrd2e/r * egamma;
        vtable[i] = qqrd2e/r * fgamma;
      } else {
        ftable[i] = qqrd2e/r * (derfc + EWALD_F*grij*expm2 - 1.0);
        etable[i] = qqrd2e/r * derfc;
        vtable[i] = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
      }
      if (rsq_lookup.f > cut_respa[2]*cut_respa[2]) {
        if (rsq_lookup.f < cut_respa[3]*cut_respa[3]) {
          rsw = (r - cut_respa[2])/(cut_respa[3] - cut_respa[2]);
          ftable[i] += qqrd2e/r * rsw*rsw*(3.0 - 2.0*rsw);
          ctable[i] = qqrd2e/r * rsw*rsw*(3.0 - 2.0*rsw);
        } else {
          if (msmflag) ftable[i] = qqrd2e/r * fgamma;
          else ftable[i] = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
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
  p_tmp = 0.0;
  v_tmp = 0.0;
  itablemin = minrsq_lookup.i & ncoulmask;
  itablemin >>= ncoulshiftbits;
  int itablemax = itablemin - 1;
  if (itablemin == 0) itablemax = ntablem1;
  rsq_lookup.i = itablemax << ncoulshiftbits;
  rsq_lookup.i |= maskhi;

  if (rsq_lookup.f < cut_coulsq) {
    rsq_lookup.f = cut_coulsq;
    r = sqrtf(rsq_lookup.f);
    if (msmflag) {
      egamma = 1.0 - (r/cut_coul)*force->kspace->gamma(r/cut_coul);
      fgamma = 1.0 + (rsq_lookup.f/cut_coulsq)*
        force->kspace->dgamma(r/cut_coul);
    } else {
      grij = g_ewald * r;
      expm2 = exp(-grij*grij);
      derfc = erfc(grij);
    }
    if (cut_respa == NULL) {
      c_tmp = qqrd2e/r;
      if (msmflag) {
        f_tmp = qqrd2e/r * fgamma;
        e_tmp = qqrd2e/r * egamma;
      } else {
        f_tmp = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
        e_tmp = qqrd2e/r * derfc;
      }
    } else {
      c_tmp = 0.0;
      p_tmp = qqrd2e/r;
      if (msmflag) {
        f_tmp = qqrd2e/r * (fgamma - 1.0);
        e_tmp = qqrd2e/r * egamma;
        v_tmp = qqrd2e/r * fgamma;
      } else {
        f_tmp = qqrd2e/r * (derfc + EWALD_F*grij*expm2 - 1.0);
        e_tmp = qqrd2e/r * derfc;
        v_tmp = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
      }
      if (rsq_lookup.f > cut_respa[2]*cut_respa[2]) {
        if (rsq_lookup.f < cut_respa[3]*cut_respa[3]) {
          rsw = (r - cut_respa[2])/(cut_respa[3] - cut_respa[2]);
          f_tmp += qqrd2e/r * rsw*rsw*(3.0 - 2.0*rsw);
          c_tmp = qqrd2e/r * rsw*rsw*(3.0 - 2.0*rsw);
        } else {
          if (msmflag) f_tmp = qqrd2e/r * fgamma;
          else f_tmp = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
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

/* ----------------------------------------------------------------------
 setup force tables for dispersion used in compute routines
 ------------------------------------------------------------------------- */

void Pair::init_tables_disp(double cut_lj_global)
{
  int masklo,maskhi;
  double rsq;
  double g_ewald_6 = force->kspace->g_ewald_6;
  double g2 = g_ewald_6*g_ewald_6, g6 = g2*g2*g2, g8 = g6*g2;
  
  tabinnerdispsq = tabinner_disp*tabinner_disp;
  init_bitmap(tabinner_disp,cut_lj_global,ndisptablebits,
              masklo,maskhi,ndispmask,ndispshiftbits);
  
  int ntable = 1;
  for (int i = 0; i < ndisptablebits; i++) ntable *= 2;
  
  // linear lookup tables of length N = 2^ndisptablebits
  // stored value = value at lower edge of bin
  // d values = delta from lower edge to upper edge of bin
  
  if (fdisptable) free_disp_tables();
  
  memory->create(rdisptable,ntable,"pair:rdisptable");
  memory->create(fdisptable,ntable,"pair:fdisptable");
  memory->create(edisptable,ntable,"pair:edisptable");
  memory->create(drdisptable,ntable,"pair:drdisptable");
  memory->create(dfdisptable,ntable,"pair:dfdisptable");
  memory->create(dedisptable,ntable,"pair:dedisptable");
  
  union_int_float_t rsq_lookup;
  union_int_float_t minrsq_lookup;
  int itablemin;
  minrsq_lookup.i = 0 << ndispshiftbits;
  minrsq_lookup.i |= maskhi;
  
  for (int i = 0; i < ntable; i++) {
    rsq_lookup.i = i << ndispshiftbits;
    rsq_lookup.i |= masklo;
    if (rsq_lookup.f < tabinnerdispsq) {
      rsq_lookup.i = i << ndispshiftbits;
      rsq_lookup.i |= maskhi;
    }
    rsq = rsq_lookup.f;
    register double x2 = g2*rsq, a2 = 1.0/x2;
    x2 = a2*exp(-x2);
    
    rdisptable[i] = rsq_lookup.f;
    fdisptable[i] = g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq;
    edisptable[i] = g6*((a2+1.0)*a2+0.5)*x2;
    
    minrsq_lookup.f = MIN(minrsq_lookup.f,rsq_lookup.f);
  }
  
  tabinnerdispsq = minrsq_lookup.f;
  
  int ntablem1 = ntable - 1;
  
  for (int i = 0; i < ntablem1; i++) {
    drdisptable[i] = 1.0/(rdisptable[i+1] - rdisptable[i]);
    dfdisptable[i] = fdisptable[i+1] - fdisptable[i];
    dedisptable[i] = edisptable[i+1] - edisptable[i];
  }
  
  // get the delta values for the last table entries
  // tables are connected periodically between 0 and ntablem1
  
  drdisptable[ntablem1] = 1.0/(rdisptable[0] - rdisptable[ntablem1]);
  dfdisptable[ntablem1] = fdisptable[0] - fdisptable[ntablem1];
  dedisptable[ntablem1] = edisptable[0] - edisptable[ntablem1];
  
  // get the correct delta values at itablemax
  // smallest r is in bin itablemin
  // largest r is in bin itablemax, which is itablemin-1,
  //   or ntablem1 if itablemin=0
  // deltas at itablemax only needed if corresponding rsq < cut*cut
  // if so, compute deltas between rsq and cut*cut
  
  double f_tmp,e_tmp;
  double cut_lj_globalsq;
  itablemin = minrsq_lookup.i & ndispmask;
  itablemin >>= ndispshiftbits;
  int itablemax = itablemin - 1;
  if (itablemin == 0) itablemax = ntablem1;
  rsq_lookup.i = itablemax << ndispshiftbits;
  rsq_lookup.i |= maskhi;
  
  if (rsq_lookup.f < (cut_lj_globalsq = cut_lj_global * cut_lj_global)) {
    rsq_lookup.f = cut_lj_globalsq;
    
    register double x2 = g2*rsq, a2 = 1.0/x2;
    x2 = a2*exp(-x2);
    f_tmp = g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq;
    e_tmp = g6*((a2+1.0)*a2+0.5)*x2;
    
    drdisptable[itablemax] = 1.0/(rsq_lookup.f - rdisptable[itablemax]);
    dfdisptable[itablemax] = f_tmp - fdisptable[itablemax];
    dedisptable[itablemax] = e_tmp - edisptable[itablemax];
  }
}

/* ----------------------------------------------------------------------
   free memory for tables used in Coulombic pair computations
------------------------------------------------------------------------- */

void Pair::free_tables()
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

/* ----------------------------------------------------------------------
  free memory for tables used in pair computations for dispersion
  ------------------------------------------------------------------------- */

void Pair::free_disp_tables()
{
  memory->destroy(rdisptable);
  memory->destroy(drdisptable);
  memory->destroy(fdisptable);
  memory->destroy(dfdisptable);
  memory->destroy(edisptable);
  memory->destroy(dedisptable);
}
/* ----------------------------------------------------------------------
   mixing of pair potential prefactors (epsilon)
------------------------------------------------------------------------- */

double Pair::mix_energy(double eps1, double eps2, double sig1, double sig2)
{
  double value;
  if (mix_flag == GEOMETRIC)
    value = sqrt(eps1*eps2);
  else if (mix_flag == ARITHMETIC)
    value = sqrt(eps1*eps2);
  else if (mix_flag == SIXTHPOWER)
    value = 2.0 * sqrt(eps1*eps2) *
      pow(sig1,3.0) * pow(sig2,3.0) / (pow(sig1,6.0) + pow(sig2,6.0));
  return value;
}

/* ----------------------------------------------------------------------
   mixing of pair potential distances (sigma, cutoff)
------------------------------------------------------------------------- */

double Pair::mix_distance(double sig1, double sig2)
{
  double value;
  if (mix_flag == GEOMETRIC)
    value = sqrt(sig1*sig2);
  else if (mix_flag == ARITHMETIC)
    value = 0.5 * (sig1+sig2);
  else if (mix_flag == SIXTHPOWER)
    value = pow((0.5 * (pow(sig1,6.0) + pow(sig2,6.0))),1.0/6.0);
  return value;
}

/* ---------------------------------------------------------------------- */

void Pair::compute_dummy(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;
}

/* ----------------------------------------------------------------------
   setup for energy, virial computation
   see integrate::ev_set() for values of eflag (0-3) and vflag (0-6)
------------------------------------------------------------------------- */

void Pair::ev_setup(int eflag, int vflag)
{
  int i,n;

  evflag = 1;

  eflag_either = eflag;
  eflag_global = eflag % 2;
  eflag_atom = eflag / 2;

  vflag_either = vflag;
  vflag_global = vflag % 4;
  vflag_atom = vflag / 4;

  // reallocate per-atom arrays if necessary

  if (eflag_atom && atom->nmax > maxeatom) {
    maxeatom = atom->nmax;
    memory->destroy(eatom);
    memory->create(eatom,comm->nthreads*maxeatom,"pair:eatom");
  }
  if (vflag_atom && atom->nmax > maxvatom) {
    maxvatom = atom->nmax;
    memory->destroy(vatom);
    memory->create(vatom,comm->nthreads*maxvatom,6,"pair:vatom");
  }

  // zero accumulators
  // use force->newton instead of newton_pair
  //   b/c some bonds/dihedrals call pair::ev_tally with pairwise info

  if (eflag_global) eng_vdwl = eng_coul = 0.0;
  if (vflag_global) for (i = 0; i < 6; i++) virial[i] = 0.0;
  if (eflag_atom) {
    n = atom->nlocal;
    if (force->newton) n += atom->nghost;
    for (i = 0; i < n; i++) eatom[i] = 0.0;
  }
  if (vflag_atom) {
    n = atom->nlocal;
    if (force->newton) n += atom->nghost;
    for (i = 0; i < n; i++) {
      vatom[i][0] = 0.0;
      vatom[i][1] = 0.0;
      vatom[i][2] = 0.0;
      vatom[i][3] = 0.0;
      vatom[i][4] = 0.0;
      vatom[i][5] = 0.0;
    }
  }

  // if vflag_global = 2 and pair::compute() calls virial_fdotr_compute()
  // compute global virial via (F dot r) instead of via pairwise summation
  // unset other flags as appropriate

  if (vflag_global == 2 && no_virial_fdotr_compute == 0) {
    vflag_fdotr = 1;
    vflag_global = 0;
    if (vflag_atom == 0) vflag_either = 0;
    if (vflag_either == 0 && eflag_either == 0) evflag = 0;
  } else vflag_fdotr = 0;

  if (lmp->cuda) lmp->cuda->evsetup_eatom_vatom(eflag_atom,vflag_atom);
}

/* ----------------------------------------------------------------------
   set all flags to zero for energy, virial computation
   called by some complicated many-body potentials that use individual flags
   to insure no holdover of flags from previous timestep
------------------------------------------------------------------------- */

void Pair::ev_unset()
{
  evflag = 0;

  eflag_either = 0;
  eflag_global = 0;
  eflag_atom = 0;

  vflag_either = 0;
  vflag_global = 0;
  vflag_atom = 0;
  vflag_fdotr = 0;
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   need i < nlocal test since called by bond_quartic and dihedral_charmm
------------------------------------------------------------------------- */

void Pair::ev_tally(int i, int j, int nlocal, int newton_pair,
                    double evdwl, double ecoul, double fpair,
                    double delx, double dely, double delz)
{
  double evdwlhalf,ecoulhalf,epairhalf,v[6];

  if (eflag_either) {
    if (eflag_global) {
      if (newton_pair) {
        eng_vdwl += evdwl;
        eng_coul += ecoul;
      } else {
        evdwlhalf = 0.5*evdwl;
        ecoulhalf = 0.5*ecoul;
        if (i < nlocal) {
          eng_vdwl += evdwlhalf;
          eng_coul += ecoulhalf;
        }
        if (j < nlocal) {
          eng_vdwl += evdwlhalf;
          eng_coul += ecoulhalf;
        }
      }
    }
    if (eflag_atom) {
      epairhalf = 0.5 * (evdwl + ecoul);
      if (newton_pair || i < nlocal) eatom[i] += epairhalf;
      if (newton_pair || j < nlocal) eatom[j] += epairhalf;
    }
  }

  if (vflag_either) {
    v[0] = delx*delx*fpair;
    v[1] = dely*dely*fpair;
    v[2] = delz*delz*fpair;
    v[3] = delx*dely*fpair;
    v[4] = delx*delz*fpair;
    v[5] = dely*delz*fpair;

    if (vflag_global) {
      if (newton_pair) {
        virial[0] += v[0];
        virial[1] += v[1];
        virial[2] += v[2];
        virial[3] += v[3];
        virial[4] += v[4];
        virial[5] += v[5];
      } else {
        if (i < nlocal) {
          virial[0] += 0.5*v[0];
          virial[1] += 0.5*v[1];
          virial[2] += 0.5*v[2];
          virial[3] += 0.5*v[3];
          virial[4] += 0.5*v[4];
          virial[5] += 0.5*v[5];
        }
        if (j < nlocal) {
          virial[0] += 0.5*v[0];
          virial[1] += 0.5*v[1];
          virial[2] += 0.5*v[2];
          virial[3] += 0.5*v[3];
          virial[4] += 0.5*v[4];
          virial[5] += 0.5*v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_pair || i < nlocal) {
        vatom[i][0] += 0.5*v[0];
        vatom[i][1] += 0.5*v[1];
        vatom[i][2] += 0.5*v[2];
        vatom[i][3] += 0.5*v[3];
        vatom[i][4] += 0.5*v[4];
        vatom[i][5] += 0.5*v[5];
      }
      if (newton_pair || j < nlocal) {
        vatom[j][0] += 0.5*v[0];
        vatom[j][1] += 0.5*v[1];
        vatom[j][2] += 0.5*v[2];
        vatom[j][3] += 0.5*v[3];
        vatom[j][4] += 0.5*v[4];
        vatom[j][5] += 0.5*v[5];
      }
    }
  }
}
 
/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   can use this version with full neighbor lists
------------------------------------------------------------------------- */

void Pair::ev_tally_full(int i, double evdwl, double ecoul, double fpair,
                         double delx, double dely, double delz)
{
  double v[6];

  if (eflag_either) {
    if (eflag_global) {
      eng_vdwl += 0.5*evdwl;
      eng_coul += 0.5*ecoul;
    }
    if (eflag_atom) eatom[i] += 0.5 * (evdwl + ecoul);
  }

  if (vflag_either) {
    v[0] = 0.5*delx*delx*fpair;
    v[1] = 0.5*dely*dely*fpair;
    v[2] = 0.5*delz*delz*fpair;
    v[3] = 0.5*delx*dely*fpair;
    v[4] = 0.5*delx*delz*fpair;
    v[5] = 0.5*dely*delz*fpair;

    if (vflag_global) {
      virial[0] += v[0];
      virial[1] += v[1];
      virial[2] += v[2];
      virial[3] += v[3];
      virial[4] += v[4];
      virial[5] += v[5];
    }

    if (vflag_atom) {
      vatom[i][0] += v[0];
      vatom[i][1] += v[1];
      vatom[i][2] += v[2];
      vatom[i][3] += v[3];
      vatom[i][4] += v[4];
      vatom[i][5] += v[5];
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   for virial, have delx,dely,delz and fx,fy,fz
------------------------------------------------------------------------- */

void Pair::ev_tally_xyz(int i, int j, int nlocal, int newton_pair,
                        double evdwl, double ecoul,
                        double fx, double fy, double fz,
                        double delx, double dely, double delz)
{
  double evdwlhalf,ecoulhalf,epairhalf,v[6];

  if (eflag_either) {
    if (eflag_global) {
      if (newton_pair) {
        eng_vdwl += evdwl;
        eng_coul += ecoul;
      } else {
        evdwlhalf = 0.5*evdwl;
        ecoulhalf = 0.5*ecoul;
        if (i < nlocal) {
          eng_vdwl += evdwlhalf;
          eng_coul += ecoulhalf;
        }
        if (j < nlocal) {
          eng_vdwl += evdwlhalf;
          eng_coul += ecoulhalf;
        }
      }
    }
    if (eflag_atom) {
      epairhalf = 0.5 * (evdwl + ecoul);
      if (newton_pair || i < nlocal) eatom[i] += epairhalf;
      if (newton_pair || j < nlocal) eatom[j] += epairhalf;
    }
  }

  if (vflag_either) {
    v[0] = delx*fx;
    v[1] = dely*fy;
    v[2] = delz*fz;
    v[3] = delx*fy;
    v[4] = delx*fz;
    v[5] = dely*fz;

    if (vflag_global) {
      if (newton_pair) {
        virial[0] += v[0];
        virial[1] += v[1];
        virial[2] += v[2];
        virial[3] += v[3];
        virial[4] += v[4];
        virial[5] += v[5];
      } else {
        if (i < nlocal) {
          virial[0] += 0.5*v[0];
          virial[1] += 0.5*v[1];
          virial[2] += 0.5*v[2];
          virial[3] += 0.5*v[3];
          virial[4] += 0.5*v[4];
          virial[5] += 0.5*v[5];
        }
        if (j < nlocal) {
          virial[0] += 0.5*v[0];
          virial[1] += 0.5*v[1];
          virial[2] += 0.5*v[2];
          virial[3] += 0.5*v[3];
          virial[4] += 0.5*v[4];
          virial[5] += 0.5*v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_pair || i < nlocal) {
        vatom[i][0] += 0.5*v[0];
        vatom[i][1] += 0.5*v[1];
        vatom[i][2] += 0.5*v[2];
        vatom[i][3] += 0.5*v[3];
        vatom[i][4] += 0.5*v[4];
        vatom[i][5] += 0.5*v[5];
      }
      if (newton_pair || j < nlocal) {
        vatom[j][0] += 0.5*v[0];
        vatom[j][1] += 0.5*v[1];
        vatom[j][2] += 0.5*v[2];
        vatom[j][3] += 0.5*v[3];
        vatom[j][4] += 0.5*v[4];
        vatom[j][5] += 0.5*v[5];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   for virial, have delx,dely,delz and fx,fy,fz
   called when using full neighbor lists
------------------------------------------------------------------------- */

void Pair::ev_tally_xyz_full(int i, double evdwl, double ecoul,
                             double fx, double fy, double fz,
                             double delx, double dely, double delz)
{
  double evdwlhalf,ecoulhalf,epairhalf,v[6];

  if (eflag_either) {
    if (eflag_global) {
      evdwlhalf = 0.5*evdwl;
      ecoulhalf = 0.5*ecoul;
      eng_vdwl += evdwlhalf;
      eng_coul += ecoulhalf;
    }
    if (eflag_atom) {
      epairhalf = 0.5 * (evdwl + ecoul);
      eatom[i] += epairhalf;
    }
  }

  if (vflag_either) {
    v[0] = 0.5*delx*fx;
    v[1] = 0.5*dely*fy;
    v[2] = 0.5*delz*fz;
    v[3] = 0.5*delx*fy;
    v[4] = 0.5*delx*fz;
    v[5] = 0.5*dely*fz;

    if (vflag_global) {
      virial[0] += v[0];
      virial[1] += v[1];
      virial[2] += v[2];
      virial[3] += v[3];
      virial[4] += v[4];
      virial[5] += v[5];
    }

    if (vflag_atom) {
      vatom[i][0] += v[0];
      vatom[i][1] += v[1];
      vatom[i][2] += v[2];
      vatom[i][3] += v[3];
      vatom[i][4] += v[4];
      vatom[i][5] += v[5];
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   called by SW and hbond potentials, newton_pair is always on
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = drji*fj + drki*fk
 ------------------------------------------------------------------------- */

void Pair::ev_tally3(int i, int j, int k, double evdwl, double ecoul,
                     double *fj, double *fk, double *drji, double *drki)
{
  double epairthird,v[6];

  if (eflag_either) {
    if (eflag_global) {
      eng_vdwl += evdwl;
      eng_coul += ecoul;
    }
    if (eflag_atom) {
      epairthird = THIRD * (evdwl + ecoul);
      eatom[i] += epairthird;
      eatom[j] += epairthird;
      eatom[k] += epairthird;
    }
  }

  if (vflag_either) {
    v[0] = drji[0]*fj[0] + drki[0]*fk[0];
    v[1] = drji[1]*fj[1] + drki[1]*fk[1];
    v[2] = drji[2]*fj[2] + drki[2]*fk[2];
    v[3] = drji[0]*fj[1] + drki[0]*fk[1];
    v[4] = drji[0]*fj[2] + drki[0]*fk[2];
    v[5] = drji[1]*fj[2] + drki[1]*fk[2];

    if (vflag_global) {
      virial[0] += v[0];
      virial[1] += v[1];
      virial[2] += v[2];
      virial[3] += v[3];
      virial[4] += v[4];
      virial[5] += v[5];
    }

    if (vflag_atom) {
      vatom[i][0] += THIRD*v[0]; vatom[i][1] += THIRD*v[1];
      vatom[i][2] += THIRD*v[2]; vatom[i][3] += THIRD*v[3];
      vatom[i][4] += THIRD*v[4]; vatom[i][5] += THIRD*v[5];

      vatom[j][0] += THIRD*v[0]; vatom[j][1] += THIRD*v[1];
      vatom[j][2] += THIRD*v[2]; vatom[j][3] += THIRD*v[3];
      vatom[j][4] += THIRD*v[4]; vatom[j][5] += THIRD*v[5];

      vatom[k][0] += THIRD*v[0]; vatom[k][1] += THIRD*v[1];
      vatom[k][2] += THIRD*v[2]; vatom[k][3] += THIRD*v[3];
      vatom[k][4] += THIRD*v[4]; vatom[k][5] += THIRD*v[5];
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   called by AIREBO potential, newton_pair is always on
 ------------------------------------------------------------------------- */

void Pair::ev_tally4(int i, int j, int k, int m, double evdwl,
                     double *fi, double *fj, double *fk,
                     double *drim, double *drjm, double *drkm)
{
  double epairfourth,v[6];

  if (eflag_either) {
    if (eflag_global) eng_vdwl += evdwl;
    if (eflag_atom) {
      epairfourth = 0.25 * evdwl;
      eatom[i] += epairfourth;
      eatom[j] += epairfourth;
      eatom[k] += epairfourth;
      eatom[m] += epairfourth;
    }
  }

  if (vflag_atom) {
    v[0] = 0.25 * (drim[0]*fi[0] + drjm[0]*fj[0] + drkm[0]*fk[0]);
    v[1] = 0.25 * (drim[1]*fi[1] + drjm[1]*fj[1] + drkm[1]*fk[1]);
    v[2] = 0.25 * (drim[2]*fi[2] + drjm[2]*fj[2] + drkm[2]*fk[2]);
    v[3] = 0.25 * (drim[0]*fi[1] + drjm[0]*fj[1] + drkm[0]*fk[1]);
    v[4] = 0.25 * (drim[0]*fi[2] + drjm[0]*fj[2] + drkm[0]*fk[2]);
    v[5] = 0.25 * (drim[1]*fi[2] + drjm[1]*fj[2] + drkm[1]*fk[2]);

    vatom[i][0] += v[0]; vatom[i][1] += v[1]; vatom[i][2] += v[2];
    vatom[i][3] += v[3]; vatom[i][4] += v[4]; vatom[i][5] += v[5];
    vatom[j][0] += v[0]; vatom[j][1] += v[1]; vatom[j][2] += v[2];
    vatom[j][3] += v[3]; vatom[j][4] += v[4]; vatom[j][5] += v[5];
    vatom[k][0] += v[0]; vatom[k][1] += v[1]; vatom[k][2] += v[2];
    vatom[k][3] += v[3]; vatom[k][4] += v[4]; vatom[k][5] += v[5];
    vatom[m][0] += v[0]; vatom[m][1] += v[1]; vatom[m][2] += v[2];
    vatom[m][3] += v[3]; vatom[m][4] += v[4]; vatom[m][5] += v[5];
  }
}

/* ----------------------------------------------------------------------
   tally ecoul and virial into each of atoms in list
   called by TIP4P potential, newton_pair is always on
   weight assignments by alpha, so contribution is all to O atom as alpha -> 0.0
   key = 0 if neither atom = water O
   key = 1 if first atom = water O
   key = 2 if second atom = water O
   key = 3 if both atoms = water O
 ------------------------------------------------------------------------- */

void Pair::ev_tally_tip4p(int key, int *list, double *v,
                          double ecoul, double alpha)
{
  int i;

  if (eflag_either) {
    if (eflag_global) eng_coul += ecoul;
    if (eflag_atom) {
      if (key == 0) {
        eatom[list[0]] += 0.5*ecoul;
        eatom[list[1]] += 0.5*ecoul;
      } else if (key == 1) {
        eatom[list[0]] += 0.5*ecoul*(1-alpha);
        eatom[list[1]] += 0.25*ecoul*alpha;
        eatom[list[2]] += 0.25*ecoul*alpha;
        eatom[list[3]] += 0.5*ecoul;
      } else if (key == 2) {
        eatom[list[0]] += 0.5*ecoul;
        eatom[list[1]] += 0.5*ecoul*(1-alpha);
        eatom[list[2]] += 0.25*ecoul*alpha;
        eatom[list[3]] += 0.25*ecoul*alpha;
      } else {
        eatom[list[0]] += 0.5*ecoul*(1-alpha);
        eatom[list[1]] += 0.25*ecoul*alpha;
        eatom[list[2]] += 0.25*ecoul*alpha;
        eatom[list[3]] += 0.5*ecoul*(1-alpha);
        eatom[list[4]] += 0.25*ecoul*alpha;
        eatom[list[5]] += 0.25*ecoul*alpha;
      }
    }
  }

  if (vflag_either) {
    if (vflag_global) {
      virial[0] += v[0];
      virial[1] += v[1];
      virial[2] += v[2];
      virial[3] += v[3];
      virial[4] += v[4];
      virial[5] += v[5];
    }

    if (vflag_atom) {
      if (key == 0) {
        for (i = 0; i <= 5; i++) {
          vatom[list[0]][i] += 0.5*v[i];
          vatom[list[1]][i] += 0.5*v[i];
        }
      } else if (key == 1) {
        for (i = 0; i <= 5; i++) {
          vatom[list[0]][i] += 0.5*v[i]*(1-alpha);
          vatom[list[1]][i] += 0.25*v[i]*alpha;
          vatom[list[2]][i] += 0.25*v[i]*alpha;
          vatom[list[3]][i] += 0.5*v[i];
        }
      } else if (key == 2) {
        for (i = 0; i <= 5; i++) {
          vatom[list[0]][i] += 0.5*v[i];
          vatom[list[1]][i] += 0.5*v[i]*(1-alpha);
          vatom[list[2]][i] += 0.25*v[i]*alpha;
          vatom[list[3]][i] += 0.25*v[i]*alpha;
        }
      } else {
        for (i = 0; i <= 5; i++) {
          vatom[list[0]][i] += 0.5*v[i]*(1-alpha);
          vatom[list[1]][i] += 0.25*v[i]*alpha;
          vatom[list[2]][i] += 0.25*v[i]*alpha;
          vatom[list[3]][i] += 0.5*v[i]*(1-alpha);
          vatom[list[4]][i] += 0.25*v[i]*alpha;
          vatom[list[5]][i] += 0.25*v[i]*alpha;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by REAX/C potential, newton_pair is always on
   fi is magnitude of force on atom i
------------------------------------------------------------------------- */

void Pair::v_tally(int i, double *fi, double *deli)
{
  double v[6];

  v[0] = 0.5*deli[0]*fi[0];
  v[1] = 0.5*deli[1]*fi[1];
  v[2] = 0.5*deli[2]*fi[2];
  v[3] = 0.5*deli[0]*fi[1];
  v[4] = 0.5*deli[0]*fi[2];
  v[5] = 0.5*deli[1]*fi[2];

  vatom[i][0] += v[0]; vatom[i][1] += v[1]; vatom[i][2] += v[2];
  vatom[i][3] += v[3]; vatom[i][4] += v[4]; vatom[i][5] += v[5];
}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by AIREBO potential, newton_pair is always on
   fpair is magnitude of force on atom I
------------------------------------------------------------------------- */

void Pair::v_tally2(int i, int j, double fpair, double *drij)
{
  double v[6];

  v[0] = 0.5 * drij[0]*drij[0]*fpair;
  v[1] = 0.5 * drij[1]*drij[1]*fpair;
  v[2] = 0.5 * drij[2]*drij[2]*fpair;
  v[3] = 0.5 * drij[0]*drij[1]*fpair;
  v[4] = 0.5 * drij[0]*drij[2]*fpair;
  v[5] = 0.5 * drij[1]*drij[2]*fpair;

  vatom[i][0] += v[0]; vatom[i][1] += v[1]; vatom[i][2] += v[2];
  vatom[i][3] += v[3]; vatom[i][4] += v[4]; vatom[i][5] += v[5];
  vatom[j][0] += v[0]; vatom[j][1] += v[1]; vatom[j][2] += v[2];
  vatom[j][3] += v[3]; vatom[j][4] += v[4]; vatom[j][5] += v[5];
}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by AIREBO and Tersoff potential, newton_pair is always on
------------------------------------------------------------------------- */

void Pair::v_tally3(int i, int j, int k,
                    double *fi, double *fj, double *drik, double *drjk)
{
  double v[6];

  v[0] = THIRD * (drik[0]*fi[0] + drjk[0]*fj[0]);
  v[1] = THIRD * (drik[1]*fi[1] + drjk[1]*fj[1]);
  v[2] = THIRD * (drik[2]*fi[2] + drjk[2]*fj[2]);
  v[3] = THIRD * (drik[0]*fi[1] + drjk[0]*fj[1]);
  v[4] = THIRD * (drik[0]*fi[2] + drjk[0]*fj[2]);
  v[5] = THIRD * (drik[1]*fi[2] + drjk[1]*fj[2]);

  vatom[i][0] += v[0]; vatom[i][1] += v[1]; vatom[i][2] += v[2];
  vatom[i][3] += v[3]; vatom[i][4] += v[4]; vatom[i][5] += v[5];
  vatom[j][0] += v[0]; vatom[j][1] += v[1]; vatom[j][2] += v[2];
  vatom[j][3] += v[3]; vatom[j][4] += v[4]; vatom[j][5] += v[5];
  vatom[k][0] += v[0]; vatom[k][1] += v[1]; vatom[k][2] += v[2];
  vatom[k][3] += v[3]; vatom[k][4] += v[4]; vatom[k][5] += v[5];
}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by AIREBO potential, newton_pair is always on
------------------------------------------------------------------------- */

void Pair::v_tally4(int i, int j, int k, int m,
                    double *fi, double *fj, double *fk,
                    double *drim, double *drjm, double *drkm)
{
  double v[6];

  v[0] = 0.25 * (drim[0]*fi[0] + drjm[0]*fj[0] + drkm[0]*fk[0]);
  v[1] = 0.25 * (drim[1]*fi[1] + drjm[1]*fj[1] + drkm[1]*fk[1]);
  v[2] = 0.25 * (drim[2]*fi[2] + drjm[2]*fj[2] + drkm[2]*fk[2]);
  v[3] = 0.25 * (drim[0]*fi[1] + drjm[0]*fj[1] + drkm[0]*fk[1]);
  v[4] = 0.25 * (drim[0]*fi[2] + drjm[0]*fj[2] + drkm[0]*fk[2]);
  v[5] = 0.25 * (drim[1]*fi[2] + drjm[1]*fj[2] + drkm[1]*fk[2]);

  vatom[i][0] += v[0]; vatom[i][1] += v[1]; vatom[i][2] += v[2];
  vatom[i][3] += v[3]; vatom[i][4] += v[4]; vatom[i][5] += v[5];
  vatom[j][0] += v[0]; vatom[j][1] += v[1]; vatom[j][2] += v[2];
  vatom[j][3] += v[3]; vatom[j][4] += v[4]; vatom[j][5] += v[5];
  vatom[k][0] += v[0]; vatom[k][1] += v[1]; vatom[k][2] += v[2];
  vatom[k][3] += v[3]; vatom[k][4] += v[4]; vatom[k][5] += v[5];
  vatom[m][0] += v[0]; vatom[m][1] += v[1]; vatom[m][2] += v[2];
  vatom[m][3] += v[3]; vatom[m][4] += v[4]; vatom[m][5] += v[5];
}

/* ----------------------------------------------------------------------
   tally virial into global and per-atom accumulators
   called by pair lubricate potential with 6 tensor components
------------------------------------------------------------------------- */

void Pair::v_tally_tensor(int i, int j, int nlocal, int newton_pair,
                          double vxx, double vyy, double vzz,
                          double vxy, double vxz, double vyz)
{
  double v[6];

  v[0] = vxx;
  v[1] = vyy;
  v[2] = vzz;
  v[3] = vxy;
  v[4] = vxz;
  v[5] = vyz;

  if (vflag_global) {
    if (newton_pair) {
      virial[0] += v[0];
      virial[1] += v[1];
      virial[2] += v[2];
      virial[3] += v[3];
      virial[4] += v[4];
      virial[5] += v[5];
    } else {
      if (i < nlocal) {
        virial[0] += 0.5*v[0];
        virial[1] += 0.5*v[1];
        virial[2] += 0.5*v[2];
        virial[3] += 0.5*v[3];
        virial[4] += 0.5*v[4];
        virial[5] += 0.5*v[5];
      }
      if (j < nlocal) {
        virial[0] += 0.5*v[0];
        virial[1] += 0.5*v[1];
        virial[2] += 0.5*v[2];
        virial[3] += 0.5*v[3];
        virial[4] += 0.5*v[4];
        virial[5] += 0.5*v[5];
      }
    }
  }

  if (vflag_atom) {
    if (newton_pair || i < nlocal) {
      vatom[i][0] += 0.5*v[0];
      vatom[i][1] += 0.5*v[1];
      vatom[i][2] += 0.5*v[2];
      vatom[i][3] += 0.5*v[3];
      vatom[i][4] += 0.5*v[4];
      vatom[i][5] += 0.5*v[5];
    }
    if (newton_pair || j < nlocal) {
      vatom[j][0] += 0.5*v[0];
      vatom[j][1] += 0.5*v[1];
      vatom[j][2] += 0.5*v[2];
      vatom[j][3] += 0.5*v[3];
      vatom[j][4] += 0.5*v[4];
      vatom[j][5] += 0.5*v[5];
    }
  }
}

/* ----------------------------------------------------------------------
   compute global pair virial via summing F dot r over own & ghost atoms
   at this point, only pairwise forces have been accumulated in atom->f
------------------------------------------------------------------------- */

void Pair::virial_fdotr_compute()
{
  double **x = atom->x;
  double **f = atom->f;

  // sum over force on all particles including ghosts

  if (neighbor->includegroup == 0) {
    int nall = atom->nlocal + atom->nghost;
    for (int i = 0; i < nall; i++) {
      virial[0] += f[i][0]*x[i][0];
      virial[1] += f[i][1]*x[i][1];
      virial[2] += f[i][2]*x[i][2];
      virial[3] += f[i][1]*x[i][0];
      virial[4] += f[i][2]*x[i][0];
      virial[5] += f[i][2]*x[i][1];
    }

  // neighbor includegroup flag is set
  // sum over force on initial nfirst particles and ghosts

  } else {
    int nall = atom->nfirst;
    for (int i = 0; i < nall; i++) {
      virial[0] += f[i][0]*x[i][0];
      virial[1] += f[i][1]*x[i][1];
      virial[2] += f[i][2]*x[i][2];
      virial[3] += f[i][1]*x[i][0];
      virial[4] += f[i][2]*x[i][0];
      virial[5] += f[i][2]*x[i][1];
    }

    nall = atom->nlocal + atom->nghost;
    for (int i = atom->nlocal; i < nall; i++) {
      virial[0] += f[i][0]*x[i][0];
      virial[1] += f[i][1]*x[i][1];
      virial[2] += f[i][2]*x[i][2];
      virial[3] += f[i][1]*x[i][0];
      virial[4] += f[i][2]*x[i][0];
      virial[5] += f[i][2]*x[i][1];
    }
  }
  
  // prevent multiple calls to update the virial
  // when a hybrid pair style uses both a gpu and non-gpu pair style
  // or when respa is used with gpu pair styles

  vflag_fdotr = 0;
}

/* ----------------------------------------------------------------------
   write a table of pair potential energy/force vs distance to a file
------------------------------------------------------------------------- */

void Pair::write_file(int narg, char **arg)
{
  if (narg < 8) error->all(FLERR,"Illegal pair_write command");
  if (single_enable == 0)
    error->all(FLERR,"Pair style does not support pair_write");

  // parse arguments

  int itype = force->inumeric(FLERR,arg[0]);
  int jtype = force->inumeric(FLERR,arg[1]);
  if (itype < 1 || itype > atom->ntypes || jtype < 1 || jtype > atom->ntypes)
    error->all(FLERR,"Invalid atom types in pair_write command");

  int n = force->inumeric(FLERR,arg[2]);

  int style;
  if (strcmp(arg[3],"r") == 0) style = RLINEAR;
  else if (strcmp(arg[3],"rsq") == 0) style = RSQ;
  else if (strcmp(arg[3],"bitmap") == 0) style = BMP;
  else error->all(FLERR,"Invalid style in pair_write command");

  double inner = force->numeric(FLERR,arg[4]);
  double outer = force->numeric(FLERR,arg[5]);
  if (inner <= 0.0 || inner >= outer)
    error->all(FLERR,"Invalid cutoffs in pair_write command");

  // open file in append mode
  // print header in format used by pair_style table

  int me;
  MPI_Comm_rank(world,&me);
  FILE *fp;
  if (me == 0) {
    fp = fopen(arg[6],"a");
    if (fp == NULL) error->one(FLERR,"Cannot open pair_write file");
    fprintf(fp,"# Pair potential %s for atom types %d %d: i,r,energy,force\n",
            force->pair_style,itype,jtype);
    if (style == RLINEAR)
      fprintf(fp,"\n%s\nN %d R %g %g\n\n",arg[7],n,inner,outer);
    if (style == RSQ)
      fprintf(fp,"\n%s\nN %d RSQ %g %g\n\n",arg[7],n,inner,outer);
  }

  // initialize potentials before evaluating pair potential
  // insures all pair coeffs are set and force constants
  // also initialize neighbor so that neighbor requests are processed
  // NOTE: might be safest to just do lmp->init()

  force->init();
  neighbor->init();

  // if pair style = any of EAM, swap in dummy fp vector

  double eamfp[2];
  eamfp[0] = eamfp[1] = 0.0;
  double *eamfp_hold;

  Pair *epair = force->pair_match("eam",0);
  if (epair) epair->swap_eam(eamfp,&eamfp_hold);

  // if atom style defines charge, swap in dummy q vec

  double q[2];
  q[0] = q[1] = 1.0;
  if (narg == 10) {
    q[0] = force->numeric(FLERR,arg[8]);
    q[1] = force->numeric(FLERR,arg[9]);
  }
  double *q_hold;

  if (atom->q) {
    q_hold = atom->q;
    atom->q = q;
  }

  // evaluate energy and force at each of N distances

  int masklo,maskhi,nmask,nshiftbits;
  if (style == BMP) {
    init_bitmap(inner,outer,n,masklo,maskhi,nmask,nshiftbits);
    int ntable = 1 << n;
    if (me == 0)
      fprintf(fp,"\n%s\nN %d BITMAP %g %g\n\n",arg[7],ntable,inner,outer);
    n = ntable;
  }

  double r,e,f,rsq;
  union_int_float_t rsq_lookup;

  for (int i = 0; i < n; i++) {
    if (style == RLINEAR) {
      r = inner + (outer-inner) * i/(n-1);
      rsq = r*r;
    } else if (style == RSQ) {
      rsq = inner*inner + (outer*outer - inner*inner) * i/(n-1);
      r = sqrt(rsq);
    } else if (style == BMP) {
      rsq_lookup.i = i << nshiftbits;
      rsq_lookup.i |= masklo;
      if (rsq_lookup.f < inner*inner) {
        rsq_lookup.i = i << nshiftbits;
        rsq_lookup.i |= maskhi;
      }
      rsq = rsq_lookup.f;
      r = sqrt(rsq);
    }

    if (rsq < cutsq[itype][jtype]) {
      e = single(0,1,itype,jtype,rsq,1.0,1.0,f);
      f *= r;
    } else e = f = 0.0;
    if (me == 0) fprintf(fp,"%d %g %g %g\n",i+1,r,e,f);
  }

  // restore original vecs that were swapped in for

  double *tmp;
  if (epair) epair->swap_eam(eamfp_hold,&tmp);
  if (atom->q) atom->q = q_hold;

  if (me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   define bitmap parameters based on inner and outer cutoffs
------------------------------------------------------------------------- */

void Pair::init_bitmap(double inner, double outer, int ntablebits,
             int &masklo, int &maskhi, int &nmask, int &nshiftbits)
{
  if (sizeof(int) != sizeof(float))
    error->all(FLERR,"Bitmapped lookup tables require int/float be same size");

  if (ntablebits > sizeof(float)*CHAR_BIT)
    error->all(FLERR,"Too many total bits for bitmapped lookup table");

  if (inner >= outer)
    error->warning(FLERR,"Table inner cutoff >= outer cutoff");

  int nlowermin = 1;
  while (!((pow(double(2),(double)nlowermin) <= inner*inner) &&
           (pow(double(2),(double)nlowermin+1.0) > inner*inner))) {
    if (pow(double(2),(double)nlowermin) <= inner*inner) nlowermin++;
    else nlowermin--;
  }

  int nexpbits = 0;
  double required_range = outer*outer / pow(double(2),(double)nlowermin);
  double available_range = 2.0;

  while (available_range < required_range) {
    nexpbits++;
    available_range = pow(double(2),pow(double(2),(double)nexpbits));
  }

  int nmantbits = ntablebits - nexpbits;

  if (nexpbits > sizeof(float)*CHAR_BIT - FLT_MANT_DIG)
    error->all(FLERR,"Too many exponent bits for lookup table");
  if (nmantbits+1 > FLT_MANT_DIG)
    error->all(FLERR,"Too many mantissa bits for lookup table");
  if (nmantbits < 3) error->all(FLERR,"Too few bits for lookup table");

  nshiftbits = FLT_MANT_DIG - (nmantbits+1);

  nmask = 1;
  for (int j = 0; j < ntablebits+nshiftbits; j++) nmask *= 2;
  nmask -= 1;

  union_int_float_t rsq_lookup;
  rsq_lookup.f = outer*outer;
  maskhi = rsq_lookup.i & ~(nmask);
  rsq_lookup.f = inner*inner;
  masklo = rsq_lookup.i & ~(nmask);
}

/* ---------------------------------------------------------------------- */

double Pair::memory_usage()
{
  double bytes = comm->nthreads*maxeatom * sizeof(double);
  bytes += comm->nthreads*maxvatom*6 * sizeof(double);
  return bytes;
}

