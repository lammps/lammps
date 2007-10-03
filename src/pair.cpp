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
#include "float.h"
#include "limits.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair.h"
#include "pair_soft.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

enum{GEOMETRIC,ARITHMETIC,SIXTHPOWER};
enum{R,RSQ,BMP};

/* ---------------------------------------------------------------------- */

Pair::Pair(LAMMPS *lmp) : Pointers(lmp)
{
  eng_vdwl = eng_coul = 0.0;

  comm_forward = comm_reverse = 0;

  single_enable = 1;
  respa_enable = 0;
  one_coeff = 0;

  // pair_modify settings

  offset_flag = 0;
  mix_flag = GEOMETRIC;
  tail_flag = 0;
  etail = ptail = etail_ij = ptail_ij = 0.0;
  ncoultablebits = 12;
  tabinner = sqrt(2.0);

  allocated = 0;
}

/* ----------------------------------------------------------------------
   modify parameters of the pair style
   pair_hybrid has its own version of this routine for its sub-styles
------------------------------------------------------------------------- */

void Pair::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all("Illegal pair_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mix") == 0) {
      if (iarg+2 > narg) error->all("Illegal pair_modify command");
      if (strcmp(arg[iarg+1],"geometric") == 0) mix_flag = GEOMETRIC;
      else if (strcmp(arg[iarg+1],"arithmetic") == 0) mix_flag = ARITHMETIC;
      else if (strcmp(arg[iarg+1],"sixthpower") == 0) mix_flag = SIXTHPOWER;
      else error->all("Illegal pair_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"shift") == 0) {
      if (iarg+2 > narg) error->all("Illegal pair_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) offset_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) offset_flag = 0;
      else error->all("Illegal pair_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"table") == 0) {
      if (iarg+2 > narg) error->all("Illegal pair_modify command");
      ncoultablebits = atoi(arg[iarg+1]);
      if (ncoultablebits > sizeof(float)*CHAR_BIT) 
        error->all("Too many total bits for bitmapped lookup table");
      iarg += 2;
    } else if (strcmp(arg[iarg],"tabinner") == 0) {
      if (iarg+2 > narg) error->all("Illegal pair_modify command");
      tabinner = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"tail") == 0) {
      if (iarg+2 > narg) error->all("Illegal pair_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) tail_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) tail_flag = 0;
      else error->all("Illegal pair_modify command");
      iarg += 2;
    } else error->all("Illegal pair_modify command");
  }
}

/* ---------------------------------------------------------------------- */

void Pair::init()
{
  int i,j;

  if (offset_flag && tail_flag)
    error->all("Cannot have both pair_modify shift and tail set to yes");
  if (tail_flag && domain->dimension == 2)
    error->all("Cannot use pair tail corrections with 2d simulations");
  if (tail_flag && domain->nonperiodic && comm->me == 0)
    error->warning("Using pair tail corrections with nonperiodic system");

  if (!allocated) error->all("All pair coeffs are not set");

  // I,I coeffs must be set
  // init_one() will check if I,J is set explicitly or inferred by mixing

  for (i = 1; i <= atom->ntypes; i++)
    if (setflag[i][i] == 0) error->all("All pair coeffs are not set");

  // style-specific initialization

  init_style();

  // call init_one() for each I,J
  // set cutsq for each I,J, used to neighbor
  // cutforce = max of all I,J cutoffs

  double cut;
  cutforce = 0.0;
  etail = ptail = 0.0;
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
   init specific to a pair style
   specific pair style can override this function
     if needs its own error checks
     if needs another kind of neighbor list
   request default neighbor list = half list
------------------------------------------------------------------------- */

void Pair::init_style()
{
  int irequest = neighbor->request(this);
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
      pow(sig1,3.0) * pow(sig2,3.0) / (pow(sig1,6.0) * pow(sig2,6.0));
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

/* ----------------------------------------------------------------------
   compute pair virial via own/ghost forces
   at this point, only pairwise forces have been accumulated in atom->f
------------------------------------------------------------------------- */

void Pair::virial_compute()
{
  double **x = atom->x;
  double **f = atom->f;
  int nall = atom->nlocal + atom->nghost;

  // sum over own & ghost atoms

  for (int i = 0; i < nall; i++) {
    virial[0] += f[i][0]*x[i][0];
    virial[1] += f[i][1]*x[i][1];
    virial[2] += f[i][2]*x[i][2];
    virial[3] += f[i][1]*x[i][0];
    virial[4] += f[i][2]*x[i][0];
    virial[5] += f[i][2]*x[i][1];
  }
}

/* ----------------------------------------------------------------------
   write a table of pair potential energy/force vs distance to a file
------------------------------------------------------------------------- */

void Pair::write_file(int narg, char **arg)
{
  if (narg < 8) error->all("Illegal pair_write command");
  if (single_enable == 0) error->all("Pair style does not support pair_write");

  // parse arguments

  int itype = atoi(arg[0]);
  int jtype = atoi(arg[1]);
  if (itype < 1 || itype > atom->ntypes || jtype < 1 || jtype > atom->ntypes)
    error->all("Invalid atom types in pair_write command");

  int n = atoi(arg[2]);

  int style;
  if (strcmp(arg[3],"r") == 0) style = R;
  else if (strcmp(arg[3],"rsq") == 0) style = RSQ;
  else if (strcmp(arg[3],"bitmap") == 0) style = BMP;
  else error->all("Invalid style in pair_write command");

  double inner = atof(arg[4]);
  double outer = atof(arg[5]);
  if (inner <= 0.0 || inner >= outer)
    error->all("Invalid cutoffs in pair_write command");

  // open file in append mode
  // print header in format used by pair_style table

  int me;
  MPI_Comm_rank(world,&me);
  FILE *fp;
  if (me == 0) {
    fp = fopen(arg[6],"a");
    if (fp == NULL) error->one("Cannot open pair_write file");
    fprintf(fp,"# Pair potential %s for atom types %d %d: i,r,energy,force\n",
	    force->pair_style,itype,jtype);
    if (style == R) 
      fprintf(fp,"\n%s\nN %d R %g %g\n\n",arg[7],n,inner,outer);
    if (style == RSQ) 
      fprintf(fp,"\n%s\nN %d RSQ %g %g\n\n",arg[7],n,inner,outer);
  }

  // initialize potentials before evaluating pair potential
  // insures all pair coeffs are set and force constants

  force->init();

  // if pair style = soft, set prefactor to final value

  Pair *spair = force->pair_match("soft");
  if (spair)
    ((PairSoft *) spair)->prefactor[itype][jtype] =
      ((PairSoft *) spair)->prestop[itype][jtype];

  // if pair style = EAM, swap in dummy fp vector

  double eamfp[2];
  eamfp[0] = eamfp[1] = 0.0;
  double *eamfp_hold;

  Pair *epair = force->pair_match("eam");
  if (epair) epair->extract_eam(eamfp,&eamfp_hold);

  // if atom style defines charge, swap in dummy q vec

  double q[2];
  q[0] = q[1] = 1.0;
  if (narg == 10) {
    q[0] = atof(arg[8]);
    q[1] = atof(arg[9]);
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
  float rsq_float;
  int *int_rsq = (int *) &rsq_float;
  One one;

  for (int i = 0; i < n; i++) {
    if (style == R) {
      r = inner + (outer-inner) * i/(n-1);
      rsq = r*r;
    } else if (style == RSQ) {
      rsq = inner*inner + (outer*outer - inner*inner) * i/(n-1);
      r = sqrt(rsq);
    } else if (style == BMP) {
      *int_rsq = i << nshiftbits;
      *int_rsq = *int_rsq | masklo;
      if (rsq_float < inner*inner) {
        *int_rsq = i << nshiftbits;
        *int_rsq = *int_rsq | maskhi;
      }
      rsq = rsq_float;
      r = sqrt(rsq);
    }

    if (rsq < cutsq[itype][jtype]) {
      single(0,1,itype,jtype,rsq,1.0,1.0,1,one);
      e = one.eng_coul + one.eng_vdwl;
      f = r * one.fforce;
    } else e = f = 0.0;
    if (me == 0) fprintf(fp,"%d %g %g %g\n",i+1,r,e,f);
  }

  // restore original vecs that were swapped in for

  double *tmp;
  if (epair) epair->extract_eam(eamfp_hold,&tmp);
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
    error->all("Bitmapped lookup tables require int/float be same size");
  
  if (ntablebits > sizeof(float)*CHAR_BIT) 
    error->all("Too many total bits for bitmapped lookup table");
          
  if (inner >= outer) error->warning("Table inner cutoff >= outer cutoff");
    
  int nlowermin = 1;
  while (!((pow(double(2),nlowermin) <= inner*inner) && 
           (pow(double(2),nlowermin+1) > inner*inner))) {
    if (pow(double(2),nlowermin) <= inner*inner) nlowermin++;
    else nlowermin--;
  }

  int nexpbits = 0;
  double required_range = outer*outer / pow(double(2),nlowermin);
  double available_range = 2.0;
  
  while (available_range < required_range) {
    nexpbits++;
    available_range = pow(double(2),pow(double(2),nexpbits));
  }
     
  int nmantbits = ntablebits - nexpbits;

  if (nexpbits > sizeof(float)*CHAR_BIT - FLT_MANT_DIG) 
    error->all("Too many exponent bits for lookup table");
  if (nmantbits+1 > FLT_MANT_DIG)
    error->all("Too many mantissa bits for lookup table");
  if (nmantbits < 3) error->all("Too few bits for lookup table");

  nshiftbits = FLT_MANT_DIG - (nmantbits+1);

  nmask = 1;
  for (int j = 0; j < ntablebits+nshiftbits; j++) nmask *= 2;
  nmask -= 1;

  float rsq;
  int *int_rsq = (int *) &rsq;
  rsq = outer*outer;
  maskhi = *int_rsq & ~(nmask);
  rsq = inner*inner;
  masklo = *int_rsq & ~(nmask);
}
