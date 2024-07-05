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
   Contributing author: Shi Liu (liushi@sas.upenn.edu)
------------------------------------------------------------------------- */
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_bvv.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

PairBVV::PairBVV(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  nmax=0;
  s0x=nullptr;
  s0y=nullptr;
  s0z=nullptr;
  Dix=nullptr;
  Diy=nullptr;
  Diz=nullptr;
  comm_forward = 3;
  comm_reverse = 3;
  single_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairBVV::~PairBVV()
{
    if(copymode) return;
    memory->destroy(s0x);
    memory->destroy(s0y);
    memory->destroy(s0z);
    memory->destroy(Dix);
    memory->destroy(Diy);
    memory->destroy(Diz);
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(alpha);
    memory->destroy(bvvsparam);
    memory->destroy(bvvv0);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairBVV::compute(int eflag, int vflag)
{
  int i,j,m,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl;
  double rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
/* BVs */
  double phi;
  double s,ss;
  double Aij,r,recip;
/* BVe */
/* BVVs */
  double recip2,Eij;
  double fx,fy,fz;
/* BVVe */
  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  if (atom->nmax > nmax) {
    memory->destroy(s0x);
    memory->destroy(s0y);
    memory->destroy(s0z);

    memory->destroy(Dix);
    memory->destroy(Diy);
    memory->destroy(Diz);


    nmax = atom->nmax;
    memory->create(s0x,nmax,"pair:s0x");
    memory->create(s0y,nmax,"pair:s0y");
    memory->create(s0z,nmax,"pair:s0z");
    memory->create(Dix,nmax,"pair:Dix");
    memory->create(Diy,nmax,"pair:Diy");
    memory->create(Diz,nmax,"pair:Diz");
  }

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero out BVVS                                                                                                                                          
  if (newton_pair) {
    for (i = 0; i < nall; i++) 
    { s0x[i] = 0.0;
      s0y[i] = 0.0;
      s0z[i] = 0.0; 
    }
  } else for (i = 0; i < nlocal; i++)
    {
      s0x[i] = 0.0;
      s0y[i] = 0.0;
      s0z[i] = 0.0; 
    }
  // BVVS

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      jtype=type[j];
      if (rsq < cutsq[itype][jtype]) {
        recip = 1.0/sqrt(rsq);
// pow must be switched by array calculated in bv
        s0x[i] += pow(r0[itype][jtype]*recip,alpha[itype][jtype])*recip*(delx);
        s0y[i] += pow(r0[itype][jtype]*recip,alpha[itype][jtype])*recip*(dely);
        s0z[i] += pow(r0[itype][jtype]*recip,alpha[itype][jtype])*recip*(delz);

        if (newton_pair || j < nlocal) {
        s0x[j] -= pow(r0[jtype][itype]*recip,alpha[jtype][itype])*recip*(delx);
        s0y[j] -= pow(r0[jtype][itype]*recip,alpha[jtype][itype])*recip*(dely);
        s0z[j] -= pow(r0[jtype][itype]*recip,alpha[jtype][itype])*recip*(delz);
        }
//
	}
    }
  }

   if (newton_pair) comm->reverse_comm(this);
    
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    s=s0x[i]*s0x[i]+s0y[i]*s0y[i]+s0z[i]*s0z[i]-bvvv0[itype][itype]*bvvv0[itype][itype];
    ss=s*s;
        Dix[i] = bvvsparam[itype][itype]*power_global*2.0*s0x[i]*s;
        Diy[i] = bvvsparam[itype][itype]*power_global*2.0*s0y[i]*s;
        Diz[i] = bvvsparam[itype][itype]*power_global*2.0*s0z[i]*s;
	    if (eflag) {
	      phi = bvvsparam[itype][itype]*ss;
	      if (eflag_global)eng_vdwl += phi;
	      if (eflag_atom) eatom[i] += phi;
	    }
  }

     comm->forward_comm(this);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];


      if (rsq < cutsq[itype][jtype]) {/*1*/
        r = sqrt(rsq);
        recip = 1.0/r;
        recip2 = recip*recip;
        Aij=pow(r0[itype][jtype]*recip,alpha[itype][jtype])*recip;
        Eij=(alpha[itype][jtype]+1.0)*recip2;
      
        fx = (Dix[j]-Dix[i])*Aij
             + (Dix[i]-Dix[j])*Eij*delx*delx*Aij
             + (Diy[i]-Diy[j])*Eij*delx*dely*Aij
             + (Diz[i]-Diz[j])*Eij*delx*delz*Aij;
        f[i][0] += fx;

        fy = (Diy[j]-Diy[i])*Aij
             + (Diy[i]-Diy[j])*Eij*dely*dely*Aij
             + (Diz[i]-Diz[j])*Eij*dely*delz*Aij
             + (Dix[i]-Dix[j])*Eij*dely*delx*Aij;

        f[i][1] += fy;

        fz = (Diz[j]-Diz[i])*Aij
             + (Diz[i]-Diz[j])*Eij*delz*delz*Aij
             + (Dix[i]-Dix[j])*Eij*delz*delx*Aij
             + (Diy[i]-Diy[j])*Eij*delz*dely*Aij;
        f[i][2] += fz;

        if (newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
        }
        if (eflag)   evdwl = 0.0;
        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
			     evdwl,0.0,fx,fy,fz,delx,dely,delz);

      } /*1*/
    }/*sum over j*/
  }/*loop over i*/
    
  if (vflag_fdotr) virial_fdotr_compute();
} /*end compute*/


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBVV::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
    
  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(r0,n+1,n+1,"pair:r0");
  memory->create(bvvsparam,n+1,n+1,"pair:bvvsparam");
  memory->create(alpha,n+1,n+1,"pair:alpha");
  memory->create(bvvv0,n+1,n+1,"pair:bvvv0");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBVV::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  power_global = utils::numeric(FLERR, arg[0], false, lmp);
  cut_global = utils::numeric(FLERR, arg[1], false, lmp);
  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBVV::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 7)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

/*  double epsilon_one = atof(arg[2]);*/
  double r0_one = utils::numeric(FLERR, arg[2], false, lmp);
  double alpha_one = utils::numeric(FLERR, arg[3], false, lmp);
  double bvvs_one = utils::numeric(FLERR, arg[4], false, lmp);
  double bvvv0_one = utils::numeric(FLERR, arg[5], false, lmp);
  double cut_one = cut_global;

  if (narg == 7) double cut_one = utils::numeric(FLERR, arg[4], false, lmp);
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
/* BVVs */
      r0[i][j]=r0_one;
      alpha[i][j]=alpha_one;
      bvvsparam[i][j]=bvvs_one;
      bvvv0[i][j]=bvvv0_one;
/* BVVs */
      setflag[i][j] = 1;
      count++;

    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}
/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBVV::init_style()
{

  neighbor->add_request(this);

}

/* ----------------------------------------------------------------------
   ----------------------------------------------------------------------*/


double PairBVV::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");


  offset[i][j] = 0.0;
  alpha[j][i] = alpha[i][j];
  r0[j][i] = r0[i][j];
  bvvsparam[j][i] = bvvsparam[i][j];
  bvvv0[j][i] = bvvv0[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* BVs */

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBVV::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
/* BVVs */
        fwrite(&r0[i][j],sizeof(double),1,fp);
        fwrite(&alpha[i][j],sizeof(double),1,fp);
        fwrite(&bvvsparam[i][j],sizeof(double),1,fp);
        fwrite(&bvvv0[i][j],sizeof(double),1,fp);
/* BVVs */
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBVV::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
/*          fread(&epsilon[i][j],sizeof(double),1,fp);
          fread(&sigma[i][j],sizeof(double),1,fp);*/
          utils::sfread(FLERR, &r0[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &alpha[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &bvvsparam[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &bvvv0[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
/*        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);*/
/* BVVs */
        MPI_Bcast(&r0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&alpha[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&bvvsparam[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&bvvv0[i][j],1,MPI_DOUBLE,0,world);
/* BVVs */
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBVV::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&power_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBVV::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &power_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&power_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

int PairBVV::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;

    m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = Dix[j];
    buf[m++] = Diy[j];
    buf[m++] = Diz[j];
   }

 return m; 
 }

/* ---------------------------------------------------------------------- */

void PairBVV::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  //  for (i = first; i < last; i++) fp[i] = buf[m++];

  for (i = first; i < last; i++){
    Dix[i] = buf[m++];
    Diy[i] = buf[m++];
    Diz[i] = buf[m++];  
  }
}

/* ---------------------------------------------------------------------- */

int PairBVV::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
   m = 0;
   last = first + n;
  for (i = first; i < last; i++){
      buf[m++] = s0x[i];
      buf[m++] = s0y[i];
      buf[m++] = s0z[i];
  }
return m;
}

/* ---------------------------------------------------------------------- */

void PairBVV::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  m = 0;
     for (i = 0; i < n; i++) {
       j = list[i];
       s0x[j] += buf[m++];
       s0y[j] += buf[m++];
       s0z[j] += buf[m++];
    }
}

