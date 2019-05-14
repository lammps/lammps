/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Stan Moore (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "pair_spin_dipolar_cut.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "fix_nve_spin.h"
#include "force.h"
#include "kspace.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "error.h"
#include "update.h"


using namespace LAMMPS_NS;
using namespace MathConst;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairSpinDipolarCut::PairSpinDipolarCut(LAMMPS *lmp) : PairSpin(lmp),
lockfixnvespin(NULL)
{
  single_enable = 0;
  spinflag = 1;
  respa_enable = 0;
  no_virial_fdotr_compute = 1;
  lattice_flag = 0;

  hbar = force->hplanck/MY_2PI;			// eV/(rad.THz)
  mub = 9.274e-4;                               // in A.Ang^2
  mu_0 = 785.15;                                // in eV/Ang/A^2
  mub2mu0 = mub * mub * mu_0 / (4.0*MY_PI);     // in eV.Ang^3
  mub2mu0 = mub * mub * mu_0 / (4.0*MY_PI);	// in eV
  mub2mu0hbinv = mub2mu0 / hbar;		// in rad.THz

}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairSpinDipolarCut::~PairSpinDipolarCut()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cut_spin_long);
    memory->destroy(cutsq);
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSpinDipolarCut::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2)
    error->all(FLERR,"Incorrect args in pair_style command");

  if (strcmp(update->unit_style,"metal") != 0)
    error->all(FLERR,"Spin simulations require metal unit style");

  if (!atom->sp) 
    error->all(FLERR,"Pair/spin style requires atom attribute sp");

  cut_spin_long_global = force->numeric(FLERR,arg[0]);
  
  // reset cutoffs that have been explicitly set
  
  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = i+1; j <= atom->ntypes; j++) {
        if (setflag[i][j]) {
          cut_spin_long[i][j] = cut_spin_long_global;
        }
      }
    }
  }

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSpinDipolarCut::coeff(int narg, char **arg)
{
  if (!allocated) allocate();
  
  // check if args correct

  if (strcmp(arg[2],"long") != 0)
    error->all(FLERR,"Incorrect args in pair_style command");
  if (narg < 1 || narg > 4)
    error->all(FLERR,"Incorrect args in pair_style command");
  
  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double spin_long_cut_one = force->numeric(FLERR,arg[3]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      setflag[i][j] = 1;
      cut_spin_long[i][j] = spin_long_cut_one;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSpinDipolarCut::init_style()
{
  if (!atom->sp_flag)
    error->all(FLERR,"Pair spin requires atom/spin style");
  
  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  
  // checking if nve/spin is a listed fix

  int ifix = 0;
  while (ifix < modify->nfix) {
    if (strcmp(modify->fix[ifix]->style,"nve/spin") == 0) break;
    ifix++;
  }
  if (ifix == modify->nfix)
    error->all(FLERR,"pair/spin style requires nve/spin");

  // get the lattice_flag from nve/spin

  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"nve/spin") == 0) {
      lockfixnvespin = (FixNVESpin *) modify->fix[i];
      lattice_flag = lockfixnvespin->lattice_flag;
    }
  }

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSpinDipolarCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  
  cut_spin_long[j][i] = cut_spin_long[i][j];
  
  return cut_spin_long_global;
}

/* ----------------------------------------------------------------------
   extract the larger cutoff if "cut" or "cut_coul"
------------------------------------------------------------------------- */

void *PairSpinDipolarCut::extract(const char *str, int &dim)
{
  if (strcmp(str,"cut") == 0) {
    dim = 0;
    return (void *) &cut_spin_long_global;
  } else if (strcmp(str,"cut_coul") == 0) {
    dim = 0;
    return (void *) &cut_spin_long_global;
  } else if (strcmp(str,"ewald_order") == 0) {
    ewald_order = 0;
    ewald_order |= 1<<1;
    ewald_order |= 1<<3;
    dim = 0;
    return (void *) &ewald_order;
  } else if (strcmp(str,"ewald_mix") == 0) {
    dim = 0;
    return (void *) &mix_flag;
  }
  return NULL;
}

/* ---------------------------------------------------------------------- */

void PairSpinDipolarCut::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;  
  double rinv,r2inv,r3inv,rsq;
  double evdwl,ecoul;
  double xi[3],rij[3],eij[3];
  double spi[4],spj[4],fi[3],fmi[3];
  double local_cut2;
  int *ilist,*jlist,*numneigh,**firstneigh;  

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double **fm = atom->fm;
  double **sp = atom->sp;	
  int *type = atom->type;  
  int nlocal = atom->nlocal;  
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // computation of the exchange interaction
  // loop over atoms and their neighbors

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xi[0] = x[i][0];
    xi[1] = x[i][1];
    xi[2] = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i]; 
    spi[0] = sp[i][0]; 
    spi[1] = sp[i][1]; 
    spi[2] = sp[i][2];
    spi[3] = sp[i][3];
    itype = type[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      spj[0] = sp[j][0]; 
      spj[1] = sp[j][1]; 
      spj[2] = sp[j][2]; 
      spj[3] = sp[j][3]; 

      evdwl = 0.0;
      fi[0] = fi[1] = fi[2] = 0.0;
      fmi[0] = fmi[1] = fmi[2] = 0.0;
     
      rij[0] = x[j][0] - xi[0];
      rij[1] = x[j][1] - xi[1];
      rij[2] = x[j][2] - xi[2];
      rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
      rinv = 1.0/sqrt(rsq);
      eij[0] = rij[0]*rinv;
      eij[1] = rij[1]*rinv;
      eij[2] = rij[2]*rinv;

      local_cut2 = cut_spin_long[itype][jtype]*cut_spin_long[itype][jtype];

      if (rsq < local_cut2) {
        r2inv = 1.0/rsq;
	r3inv = r2inv*rinv;
	
	compute_dipolar(i,j,eij,fmi,spi,spj,r3inv);
	if (lattice_flag) compute_dipolar_mech(i,j,eij,fi,spi,spj,r2inv);
      }

      // force accumulation

      f[i][0] += fi[0]; 
      f[i][1] += fi[1];  	  
      f[i][2] += fi[2];
      fm[i][0] += fmi[0];
      fm[i][1] += fmi[1];
      fm[i][2] += fmi[2];

      if (newton_pair || j < nlocal) {
	f[j][0] -= fi[0];	 
        f[j][1] -= fi[1];	  	  
        f[j][2] -= fi[2];
      }

      if (eflag) {
	if (rsq <= local_cut2) {
	  evdwl -= (spi[0]*fmi[0] + spi[1]*fmi[1] + spi[2]*fmi[2]);
	  evdwl *= hbar;
	}
      } else evdwl = 0.0;

      if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
	  evdwl,ecoul,fi[0],fi[1],fi[2],rij[0],rij[1],rij[2]);

    }
  }
}

/* ----------------------------------------------------------------------
   update the pair interaction fmi acting on the spin ii
   adding 1/r (for r in [0,rc]) contribution to the pair
   removing erf(r)/r (for r in [0,rc]) from the kspace force
------------------------------------------------------------------------- */

void PairSpinDipolarCut::compute_single_pair(int ii, double fmi[3])
{
  int i,j,jj,jnum,itype,jtype;  
  double rsq,rinv,r2inv,r3inv;
  double xi[3],rij[3],eij[3];
  double spi[4],spj[4];
  double local_cut2;
  int *ilist,*jlist,*numneigh,**firstneigh;  

  double **x = atom->x;
  double **sp = atom->sp;	
  int *type = atom->type;  

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // computation of the exchange interaction
  // loop over neighbors of atom i
    
  i = ilist[ii];
  xi[0] = x[i][0];
  xi[1] = x[i][1];
  xi[2] = x[i][2];
  spi[0] = sp[i][0]; 
  spi[1] = sp[i][1]; 
  spi[2] = sp[i][2];
  spi[3] = sp[i][3];
  jlist = firstneigh[i];
  jnum = numneigh[i]; 
  itype = type[i];
  
  for (jj = 0; jj < jnum; jj++) {
    j = jlist[jj];
    j &= NEIGHMASK;
    jtype = type[j];

    spj[0] = sp[j][0]; 
    spj[1] = sp[j][1]; 
    spj[2] = sp[j][2]; 
    spj[3] = sp[j][3]; 

    rij[0] = x[j][0] - xi[0];
    rij[1] = x[j][1] - xi[1];
    rij[2] = x[j][2] - xi[2];
    rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
    rinv = 1.0/sqrt(rsq);
    eij[0] = rij[0]*rinv;
    eij[1] = rij[1]*rinv;
    eij[2] = rij[2]*rinv;

    local_cut2 = cut_spin_long[itype][jtype]*cut_spin_long[itype][jtype];

    if (rsq < local_cut2) {
      r2inv = 1.0/rsq;
      r3inv = r2inv*rinv;
      
      // compute dipolar interaction
      
      compute_dipolar(i,j,eij,fmi,spi,spj,r3inv);
    }
  }
}

/* ----------------------------------------------------------------------
   compute dipolar interaction between spins i and j
------------------------------------------------------------------------- */

void PairSpinDipolarCut::compute_dipolar(int i, int j, double eij[3], 
    double fmi[3], double spi[4], double spj[4], double r3inv)
{
  double sjdotr;
  double gigjiri3,pre;

  sjdotr = spj[0]*eij[0] + spj[1]*eij[1] + spj[2]*eij[2];
  gigjiri3 = (spi[3] * spj[3])*r3inv;
  pre = mub2mu0hbinv * gigjiri3;

  fmi[0] += pre * (3.0 * sjdotr *eij[0] - spj[0]);
  fmi[1] += pre * (3.0 * sjdotr *eij[1] - spj[1]);
  fmi[2] += pre * (3.0 * sjdotr *eij[2] - spj[2]);
}

/* ----------------------------------------------------------------------
   compute the mechanical force due to the dipolar interaction between 
   atom i and atom j
------------------------------------------------------------------------- */

void PairSpinDipolarCut::compute_dipolar_mech(int i, int j, double eij[3],
    double fi[3], double spi[3], double spj[3], double r2inv)
{
  double sisj,sieij,sjeij;
  double gigjri4,bij,pre;

  gigjri4 = (spi[3] * spj[3])*r2inv*r2inv;
  sisj = spi[0]*spj[0] + spi[1]*spj[1] + spi[2]*spj[2];
  sieij = spi[0]*eij[0] + spi[1]*eij[1] + spi[2]*eij[2];
  sjeij = spj[0]*eij[0] + spj[1]*eij[1] + spj[2]*eij[2];
  
  bij = sisj - 5.0*sieij*sjeij;
  pre = mub2mu0*gigjri4;

  fi[0] += pre * (eij[0] * bij + (sjeij*spi[0] + sieij*spj[0]));
  fi[1] += pre * (eij[1] * bij + (sjeij*spi[1] + sieij*spj[1]));
  fi[2] += pre * (eij[2] * bij + (sjeij*spi[2] + sieij*spj[2]));
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSpinDipolarCut::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cut_spin_long,n+1,n+1,"pair/spin/long:cut_spin_long");
  memory->create(cutsq,n+1,n+1,"pair/spin/long:cutsq");
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpinDipolarCut::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
	fwrite(&cut_spin_long[i][j],sizeof(int),1,fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpinDipolarCut::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
	if (me == 0) {
	  fread(&cut_spin_long[i][j],sizeof(int),1,fp);
	}
	MPI_Bcast(&cut_spin_long[i][j],1,MPI_INT,0,world);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpinDipolarCut::write_restart_settings(FILE *fp)
{
  fwrite(&cut_spin_long_global,sizeof(double),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpinDipolarCut::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_spin_long_global,sizeof(double),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_spin_long_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}
