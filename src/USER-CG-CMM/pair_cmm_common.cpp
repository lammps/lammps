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
   Common functionality for the CMM coarse grained MD potentials.
   Contributing author: Axel Kohlmeyer <akohlmey@gmail.com>
------------------------------------------------------------------------- */

#include "pair_cmm_common.h"
#include "memory.h"

#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "math.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 1.0e-6

/* ---------------------------------------------------------------------- */

PairCMMCommon::PairCMMCommon(class LAMMPS *lmp) : Pair(lmp)
{
  ftable = NULL;
  allocated_coul = 0;
  kappa = 0.0;
  respa_enable = 0;
  single_enable = 0;
}

/* ---------------------------------------------------------------------- *
 * clean up common arrays                                                 *
 * ---------------------------------------------------------------------- */

PairCMMCommon::~PairCMMCommon() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cg_type);

    memory->destroy(cut);
    memory->destroy(cutsq);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(offset);

    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);

    allocated = 0;
  }
}

/* ---------------------------------------------------------------------- *
 * allocate common arrays                                                 *
 * ---------------------------------------------------------------------- */

void PairCMMCommon::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"paircg:setflag");
  memory->create(cg_type,n+1,n+1,"paircg:cg_type");
  for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) {
      setflag[i][j] = 0;
      cg_type[i][j] = CG_NOT_SET;
    }
  }

  memory->create(cut,n+1,n+1,"paircg:cut");
  memory->create(cutsq,n+1,n+1,"paircg:cutsq");
  memory->create(epsilon,n+1,n+1,"paircg:epsilon");
  memory->create(sigma,n+1,n+1,"paircg:sigma");
  memory->create(offset,n+1,n+1,"paircg:offset");

  memory->create(lj1,n+1,n+1,"paircg:lj1");
  memory->create(lj2,n+1,n+1,"paircg:lj2");
  memory->create(lj3,n+1,n+1,"paircg:lj3");
  memory->create(lj4,n+1,n+1,"paircg:lj4");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

// arguments to the pair_style command (global version)
// args = cutoff (cutoff2 (kappa))
void PairCMMCommon::settings(int narg, char **arg)
{
  if ((narg < 1) || (narg > 3)) error->all(FLERR,"Illegal pair_style command");

  cut_lj_global = force->numeric(FLERR,arg[0]);
  if (narg == 1) cut_coul_global = cut_lj_global;
  else cut_coul_global = force->numeric(FLERR,arg[1]);
  cut_coulsq_global = cut_coul_global*cut_coul_global;

  // exponential coulomb screening (optional)
  if (narg == 3) kappa = force->numeric(FLERR,arg[2]);
  if (fabs(kappa) < SMALL) kappa=0.0;

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

void PairCMMCommon::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 7) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  int cg_type_one=find_cg_type(arg[2]);
  if (cg_type_one == CG_NOT_SET) error->all(FLERR,"Error reading CG type flag.");

  double epsilon_one = force->numeric(FLERR,arg[3]);
  double sigma_one = force->numeric(FLERR,arg[4]);

  double cut_lj_one = cut_lj_global;
  double cut_coul_one = cut_coul_global;
  if (narg >= 6) cut_lj_one = force->numeric(FLERR,arg[5]);
  if (narg == 7) cut_coul_one = force->numeric(FLERR,arg[6]);

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

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairCMMCommon::init_style()
{
  // request regular or rRESPA neighbor lists

  int irequest;

  if (update->whichflag == 0 && strstr(update->integrate_style,"respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;

    if (respa == 0) irequest = neighbor->request(this);
    else if (respa == 1) {
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respaouter = 1;
    } else {
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 2;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respamiddle = 1;
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respaouter = 1;
    }

  } else irequest = neighbor->request(this);
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   regular or rRESPA
------------------------------------------------------------------------- */

void PairCMMCommon::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  else if (id == 1) listinner = ptr;
  else if (id == 2) listmiddle = ptr;
  else if (id == 3) listouter = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCMMCommon::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    error->all(FLERR,"for CG styles, epsilon and sigma need to be set explicitly for all pairs.");
  }

  const int cgt = cg_type[i][j];

  if (cgt == CG_NOT_SET)
    error->all(FLERR,"unrecognized LJ parameter flag");

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
    error->all(FLERR,"tail correction not (yet) supported by CG potentials.");
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

    double sig2 = sigma[i][j]*sigma[i][j];
    double sig6 = sig2*sig2*sig2;
    double rc3 = cut[i][j]*cut[i][j]*cut[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;
    etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (sig6 - 3.0*rc6) / (9.0*rc9);
    ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9);
#endif
  }

  return mycut;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCMMCommon::write_restart(FILE *fp)
{
  int i,j;
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

void PairCMMCommon::read_restart(FILE *fp)
{
  int i,j;
  int me = comm->me;
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

void PairCMMCommon::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul_global,sizeof(double),1,fp);
  fwrite(&kappa,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCMMCommon::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul_global,sizeof(double),1,fp);
    fread(&kappa,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&kappa,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
  cut_coulsq_global = cut_coul_global*cut_coul_global;
}

/* ---------------------------------------------------------------------- */

double PairCMMCommon::memory_usage()
{
  double bytes=Pair::memory_usage();

  int n = atom->ntypes;

  // setflag/cg_type
  bytes += (n+1)*(n+1)*sizeof(int)*2;
  // cut/cutsq/epsilon/sigma/offset/lj1/lj2/lj3/lj4
  bytes += (n+1)*(n+1)*sizeof(double)*9;

  return bytes;
}

/* ------------------------------------------------------------------------ */

double PairCMMCommon::eval_single(int coul_type, int i, int j, int itype, int jtype,
                           double rsq, double factor_coul, double factor_lj,
                           double &fforce)
{
  double lj_force, lj_erg, coul_force, coul_erg;
  lj_force=lj_erg=coul_force=coul_erg=0.0;

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

  if (rsq < cut_coul[itype][jtype]) {
    if(coul_type == CG_COUL_LONG) {
      error->all(FLERR,"single energy computation with long-range coulomb not supported by CG potentials.");
    } else if ((coul_type == CG_COUL_CUT) || (coul_type == CG_COUL_DEBYE)) {
      const double r2inv = 1.0/rsq;
      const double rinv = sqrt(r2inv);
      const double qscreen=exp(-kappa*sqrt(rsq));
      coul_force = force->qqrd2e * atom->q[i]*atom->q[j]*rinv * qscreen * (kappa + rinv);
      coul_erg   = force->qqrd2e * atom->q[i]*atom->q[j]*rinv * qscreen;
      // error->all(FLERR,"single energy computation with coulomb not supported by CG potentials.");
    } else if (coul_type == CG_COUL_NONE) {
      ; // do nothing
    } else {
      error->all(FLERR,"unknown coulomb type with CG potentials.");
    }
  }

  fforce = factor_lj*lj_force + factor_coul*coul_force;
  return factor_lj*lj_erg + factor_coul*coul_erg;
}
