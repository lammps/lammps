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
   CMM coarse grained MD potentials. Plain version w/o charges.
   Contributing author: Axel Kohlmeyer <akohlmey@gmail.com>
------------------------------------------------------------------------- */

#include "pair_cg_cmm_omp.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"

#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "math.h"

using namespace LAMMPS_NS;
 
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SMALL 1.0e-6

/* ---------------------------------------------------------------------- */

PairCGCMMOMP::PairCGCMMOMP(LAMMPS *lmp) : PairOMP(lmp)
{
  respa_enable = 0;
  single_enable = 1;

  allocated = 0;
  allocated_coul = 0;

  ftable = NULL;
  kappa = 0.0;
}

/* ---------------------------------------------------------------------- */

PairCGCMMOMP::~PairCGCMMOMP()
{
  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_int_array(cg_type);

    memory->destroy_2d_double_array(cut);
    memory->destroy_2d_double_array(cutsq);
    memory->destroy_2d_double_array(epsilon);
    memory->destroy_2d_double_array(sigma);
    memory->destroy_2d_double_array(lj1);
    memory->destroy_2d_double_array(lj2);
    memory->destroy_2d_double_array(lj3);
    memory->destroy_2d_double_array(lj4);
    memory->destroy_2d_double_array(offset);

    allocated = 0;
  }
}

/* ---------------------------------------------------------------------- */

void PairCGCMMOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
    ev_setup_thr(eflag,vflag);
  } else {
    evflag = vflag_fdotr = 0;
  }

  if (evflag) {
    if (eflag) {
      if (force->newton_pair) {
        return eval<1,1,1>();
      } else {
        return eval<1,1,0>();
      }
    } else {
      if (force->newton_pair) {
        return eval<1,0,1>();
      } else {
        return eval<1,0,0>();
      }
    }
  } else {
    if (force->newton_pair) {
      return eval<0,0,1>();
    } else {
      return eval<0,0,0>();
    }
  }
}

template < int EVFLAG, int EFLAG, int NEWTON_PAIR > 
void PairCGCMMOMP::eval()  {

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {

    const int nlocal = atom->nlocal;
    const int nall = nlocal + atom->nghost;
    const int nthreads = comm->nthreads;

    const int * const type = atom->type;
    const double * const special_lj = force->special_lj;

    const int inum = list->inum;
    const int * const ilist = list->ilist;
    const int * const numneigh = list->numneigh;
    int * const * const firstneigh = list->firstneigh;
  
    // loop over neighbors of my atoms

    int ii, jj, iifrom, iito, tid;
    const double * const * const x = atom->x;
    double * const * const f = loop_setup_thr(atom->f, iifrom, iito, tid, inum, nall, nthreads);
    for (ii = iifrom; ii < iito; ++ii) {

      const int i = ilist[ii];
      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];

      const int itype = type[i];
      const int * const jlist = firstneigh[i];
      const int jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
	int j2 = jlist[jj];

	double factor_lj = 1.0;
	if (j2 >= nall) {
	  factor_lj = special_lj[j2/nall];
	  j2 %= nall;
	}
	const int j = j2;

	const double delx = xtmp - x[j][0];
	const double dely = ytmp - x[j][1];
	const double delz = ztmp - x[j][2];
	const double rsq = delx*delx + dely*dely + delz*delz;
	const int jtype = type[j];

	double evdwl = 0.0;
	double fpair = 0.0;

	const double r2inv = 1.0/rsq;
	const int cgt=cg_type[itype][jtype];

	if (rsq < cutsq[itype][jtype]) {
	
	  fpair=factor_lj;
	  if (EFLAG) evdwl=factor_lj;

	  if (cgt == CG_LJ12_4) {
	    const double r4inv=r2inv*r2inv;
	    fpair *= r4inv*(lj1[itype][jtype]*r4inv*r4inv
			    - lj2[itype][jtype]);
	    if (EFLAG) {
	      evdwl *= r4inv*(lj3[itype][jtype]*r4inv*r4inv
			      - lj4[itype][jtype]) - offset[itype][jtype];
	    }
	  } else if (cgt == CG_LJ9_6) {
	    const double r3inv = r2inv*sqrt(r2inv);
	    const double r6inv = r3inv*r3inv;
	    fpair *= r6inv*(lj1[itype][jtype]*r3inv
			    - lj2[itype][jtype]);
	    if (EFLAG) {
	      evdwl *= r6inv*(lj3[itype][jtype]*r3inv
			      - lj4[itype][jtype]) - offset[itype][jtype];
	    }
	  } else if (cgt == CG_LJ12_6) {
	    const double r6inv = r2inv*r2inv*r2inv;
	    fpair *= r6inv*(lj1[itype][jtype]*r6inv
			    - lj2[itype][jtype]);
	    if (EFLAG) {
	      evdwl *= r6inv*(lj3[itype][jtype]*r6inv
			      - lj4[itype][jtype]) - offset[itype][jtype];
	    }
	  } else {
	    /* do nothing. this is a "cannot happen(TM)" case */
	    ;
	  }
	  fpair *= r2inv;

	  f[i][0] += delx*fpair;
	  f[i][1] += dely*fpair;
	  f[i][2] += delz*fpair;
	  if (NEWTON_PAIR || j < nlocal) {
	    f[j][0] -= delx*fpair;
	    f[j][1] -= dely*fpair;
	    f[j][2] -= delz*fpair;
	  }
	  if (EVFLAG) ev_tally(i,j,nlocal,NEWTON_PAIR,
			       evdwl,0.0,fpair,delx,dely,delz);
	}
      }
    }
    // reduce per thread forces into global force array.
    force_reduce_thr(atom->f, nall, nthreads, tid);
  }
  ev_reduce_thr();
  if (vflag_fdotr) virial_compute();
}

/* ---------------------------------------------------------------------- */

void PairCGCMMOMP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  
  setflag = memory->create_2d_int_array(n+1,n+1,"paircg:setflag");
  cg_type = memory->create_2d_int_array(n+1,n+1,"paircg:cg_type");
  for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) {
      setflag[i][j] = 0;
      cg_type[i][j] = CG_NOT_SET;
    }
  }

  cut      = memory->create_2d_double_array(n+1,n+1,"paircg:cut");
  cutsq    = memory->create_2d_double_array(n+1,n+1,"paircg:cutsq");
  epsilon  = memory->create_2d_double_array(n+1,n+1,"paircg:epsilon");
  sigma    = memory->create_2d_double_array(n+1,n+1,"paircg:sigma");
  offset   = memory->create_2d_double_array(n+1,n+1,"paircg:offset"); 

  lj1      = memory->create_2d_double_array(n+1,n+1,"paircg:lj1");
  lj2      = memory->create_2d_double_array(n+1,n+1,"paircg:lj2");
  lj3      = memory->create_2d_double_array(n+1,n+1,"paircg:lj3");
  lj4      = memory->create_2d_double_array(n+1,n+1,"paircg:lj4");
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

// arguments to the pair_style command (global version)
// args = cutoff (cutoff2 (kappa))
void PairCGCMMOMP::settings(int narg, char **arg)
{
  if ((narg < 1) || (narg > 2)) error->all("Illegal pair_style command");

  cut_lj_global = force->numeric(arg[0]);
  if (narg == 1) cut_coul_global = cut_lj_global;
  else cut_coul_global = force->numeric(arg[1]);
  cut_coulsq_global = cut_coul_global*cut_coul_global;
  
  // reset cutoffs that have been explicitly set
  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = i+1; j <= atom->ntypes; j++) {
        if (setflag[i][j]) {
          cut[i][j] = cut_lj_global;
          if (allocated_coul) {
            cut[i][j] = MAX(cut_lj_global,cut_coul_global);
            cut_lj[i][j] = cut_lj_global;
            cut_coul[i][j] = cut_coul_global;
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairCGCMMOMP::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 7) error->all("Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  int cg_type_one=find_cg_type(arg[2]);
  if (cg_type_one == CG_NOT_SET) error->all("Error reading CG type flag.");
  
  double epsilon_one = force->numeric(arg[3]);
  double sigma_one = force->numeric(arg[4]);

  double cut_lj_one = cut_lj_global;
  double cut_coul_one = cut_coul_global;
  if (narg >= 6) cut_lj_one = force->numeric(arg[5]);
  if (narg == 7) cut_coul_one = force->numeric(arg[6]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cg_type[i][j] = cg_type_one;
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      setflag[i][j] = 1;

      if (allocated_coul) {
        cut_lj[i][j] = cut_lj_one;
        cut_coul[i][j] = cut_coul_one;
      } else {
        cut[i][j] = cut_lj_one;
      }

      count++;
    }
  }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCGCMMOMP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    error->all("for CG styles, epsilon and sigma need to be set explicitly for all pairs.");
  }

  const int cgt = cg_type[i][j];

  if (cgt == CG_NOT_SET)
    error->all("unrecognized LJ parameter flag");
  
  lj1[i][j] = cg_prefact[cgt] * cg_pow1[cgt] * epsilon[i][j] * pow(sigma[i][j],cg_pow1[cgt]);
  lj2[i][j] = cg_prefact[cgt] * cg_pow2[cgt] * epsilon[i][j] * pow(sigma[i][j],cg_pow2[cgt]);
  lj3[i][j] = cg_prefact[cgt] * epsilon[i][j] * pow(sigma[i][j],cg_pow1[cgt]);
  lj4[i][j] = cg_prefact[cgt] * epsilon[i][j] * pow(sigma[i][j],cg_pow2[cgt]);

  double mycut = cut[i][j];
  if (offset_flag) {
    double ratio = sigma[i][j] / mycut;
    offset[i][j] = cg_prefact[cgt] * epsilon[i][j] * (pow(ratio,cg_pow1[cgt]) - pow(ratio,cg_pow2[cgt]));
  } else offset[i][j] = 0.0;

  if (allocated_coul) {
    mycut = MAX(cut_lj[i][j],cut_coul[i][j]);
    cut[i][j] = mycut;
    cut_ljsq[i][j]=cut_lj[i][j]*cut_lj[i][j];
    cut_coulsq[i][j]=cut_coul[i][j]*cut_coul[i][j];
    if (offset_flag) {
      double ratio = sigma[i][j] / cut_lj[i][j];
      offset[i][j] = cg_prefact[cgt] * epsilon[i][j] * (pow(ratio,cg_pow1[cgt]) - pow(ratio,cg_pow2[cgt]));
    } else offset[i][j] = 0.0;
  }
  
  // make sure data is stored symmetrically
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];
  cg_type[j][i] = cg_type[i][j];
  cut[j][i] = mycut;
  
  if (allocated_coul) {
    cut_lj[j][i]=cut_lj[i][j];
    cut_ljsq[j][i]=cut_ljsq[i][j];
    cut_coul[j][i]=cut_coul[i][j];
    cut_coulsq[j][i]=cut_coulsq[i][j];
  }

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce
  if (tail_flag) {
#if 1
    error->all("tail correction not (yet) supported by CG potentials.");
#else
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);
        
    double PI = 4.0*atan(1.0);
    double sig2 = sigma[i][j]*sigma[i][j];
    double sig6 = sig2*sig2*sig2;
    double rc3 = cut[i][j]*cut[i][j]*cut[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;
    etail_ij = 8.0*PI*all[0]*all[1]*epsilon[i][j] * 
      sig6 * (sig6 - 3.0*rc6) / (9.0*rc9); 
    ptail_ij = 16.0*PI*all[0]*all[1]*epsilon[i][j] * 
      sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9); 
#endif
  } 

  return mycut;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file 
------------------------------------------------------------------------- */

void PairCGCMMOMP::write_restart(FILE *fp)
{
  int i,j;

  write_restart_settings(fp);

  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&cg_type[i][j],sizeof(int),1,fp);
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
        if (allocated_coul) {
          fwrite(&cut_lj[i][j],sizeof(double),1,fp);
          fwrite(&cut_coul[i][j],sizeof(double),1,fp);
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCGCMMOMP::read_restart(FILE *fp)
{
  int i,j;
  int me = comm->me;

  read_restart_settings(fp);
  allocate();

  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&cg_type[i][j],sizeof(int),1,fp);
          fread(&epsilon[i][j],sizeof(double),1,fp);
          fread(&sigma[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
          if(allocated_coul) {
            fread(&cut_lj[i][j],sizeof(double),1,fp);
            fread(&cut_coul[i][j],sizeof(double),1,fp);
          }
        }
        MPI_Bcast(&cg_type[i][j],1,MPI_INT,0,world);
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
        if (allocated_coul) {
          MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&cut_coul[i][j],1,MPI_DOUBLE,0,world);
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCGCMMOMP::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul_global,sizeof(double),1,fp);
  fwrite(&kappa,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCGCMMOMP::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul_global,sizeof(double),1,fp);
    fread(&kappa,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&kappa,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  cut_coulsq_global = cut_coul_global*cut_coul_global;
}

/* ---------------------------------------------------------------------- */

double PairCGCMMOMP::memory_usage()
{
  double bytes=PairOMP::memory_usage();
  
  int n = atom->ntypes;

  // setflag/cg_type
  bytes += (n+1)*(n+1)*sizeof(int)*2; 
  // cut/cutsq/epsilon/sigma/offset/lj1/lj2/lj3/lj4
  bytes += (n+1)*(n+1)*sizeof(double)*9; 
  
  return bytes;
}

/* ---------------------------------------------------------------------- */

double PairCGCMMOMP::single(int i, int j, int itype, int jtype, double rsq,
		       double factor_coul, double factor_lj, double &fforce)
{
  double lj_force, lj_erg;
  
  lj_force=lj_erg=0.0;

  if (rsq < cut_ljsq[itype][jtype]) {
      
    const int cgt = cg_type[itype][jtype];
    const double cgpow1 = cg_pow1[cgt];
    const double cgpow2 = cg_pow2[cgt];
    const double cgpref = cg_prefact[cgt];
        
    const double ratio = sigma[itype][jtype]/sqrt(rsq);
    const double eps = epsilon[itype][jtype];

    lj_force = cgpref*eps * (cgpow1*pow(ratio,cgpow1) 
                            - cgpow2*pow(ratio,cgpow2))/rsq;
    lj_erg = cgpref*eps * (pow(ratio,cgpow1) - pow(ratio,cgpow2));
  }
  
  fforce = factor_lj*lj_force;
  return factor_lj*lj_erg;
}

/* ------------------------------------------------------------------------ */
