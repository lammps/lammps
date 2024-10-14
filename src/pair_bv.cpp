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
   Contributing author: Chuanfu Luo (luochuanfu@gmail.com)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_bv.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

PairBV::PairBV(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  nmax = 0;
  s0 = nullptr;
  fp = nullptr;
  comm_forward = 1;
  comm_reverse = 1;
  single_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairBV::~PairBV()
{
    if (copymode) return;
    memory->destroy(s0);
    memory->destroy(fp);
    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
        memory->destroy(cut);
        memory->destroy(alpha);
        memory->destroy(sparam);
        memory->destroy(v0);
        memory->destroy(offset);
    }
}

/* ---------------------------------------------------------------------- */

void PairBV::compute(int eflag, int vflag)
{
  int i,j,m,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r3inv,r6inv;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rinv,phi;
  double s,ss;
  double Aij,r,recip,psip;
  evdwl = 0.0;
    
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

  if (atom->nmax > nmax) {
    memory->destroy(s0);
    memory->destroy(fp);
    nmax = atom->nmax;
    memory->create(s0,nmax,"pair:s0");
    memory->create(fp,nmax,"pair:fp");
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

  // zero out BVS                                                               
    if (newton_pair) {
    for (i = 0; i < nall; i++) s0[i] = 0.0;
    } else for (i = 0; i < nlocal; i++) s0[i] = 0.0;

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

      jtype = type[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (alpha[itype][jtype]!=0.0) { /*1*/
        if (rsq < cutsq[itype][jtype]) {/*2*/
          recip = 1.0/sqrt(rsq);
	      r=sqrt(rsq);
          s0[i] += pow(r0[itype][jtype]/r,alpha[itype][jtype]);
          if (newton_pair || j < nlocal) {
            s0[j] += pow(r0[jtype][itype]/r,alpha[jtype][itype]);
           }
         }/*2*/
      }/*1*/
     }
  }

  if (newton_pair) comm->reverse_comm(this);


  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    s=s0[i]-v0[itype][itype];
    ss=s*s;
    fp[i] = sparam[itype][itype]*power_global*s;
    if (eflag) {
      phi = sparam[itype][itype]*ss+energy0[itype];
      if (eflag_global) eng_vdwl += phi;     
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
      if (alpha[itype][jtype]!=0.0) { /*1*/
        if (rsq < cutsq[itype][jtype]) { /*2*/
          r = sqrt(rsq);
          recip = 1.0/sqrt(rsq);
          Aij=alpha[itype][jtype]*pow(r0[itype][jtype]*recip,alpha[itype][jtype])*recip;
          psip = (fp[i]+fp[j])*Aij;
          fpair = psip*recip;

          f[i][0] += delx*fpair;
          f[i][1] += dely*fpair;
          f[i][2] += delz*fpair;

          if (newton_pair || j < nlocal) { /*3*/
            f[j][0] -= delx*fpair;
            f[j][1] -= dely*fpair;
            f[j][2] -= delz*fpair;
          } /*3*/

         if (eflag){ /*4*/
	       evdwl = 0.0;
	     }/*4*/
         if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
	   }/*2*/
      } /*1*/
    }
  }
  
  if (vflag_fdotr) virial_fdotr_compute();
}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBV::allocate()
{
  allocated = 1;
  int n = atom->ntypes+1;


  memory->create(setflag,n,n,"pair:setflag");
  for (int i = 1; i < n; i++)
    for (int j = i; j < n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n,n,"pair:cutsq");
  memory->create(cut,n,n,"pair:cut");
  memory->create(r0,n,n,"pair:r0");
  memory->create(alpha,n,n,"pair:alpha");
  memory->create(sparam,n,n,"pair:sparam");
  memory->create(v0,n,n,"pair:v0");
  memory->create(energy0,n,"pair:energy0");
  memory->create(offset,n,n,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBV::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");
  power_global = utils::numeric(FLERR, arg[0], false, lmp);
  cut_global = utils::numeric(FLERR, arg[1], false, lmp);

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
void PairBV::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 7)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double r0_one = utils::numeric(FLERR, arg[2], false, lmp);
  double alpha_one = utils::numeric(FLERR, arg[3], false, lmp);
  double s_one = utils::numeric(FLERR, arg[4], false, lmp);
  double v0_one = utils::numeric(FLERR, arg[5], false, lmp);

  double cut_one = cut_global;
  if (narg == 7) double cut_one = utils::numeric(FLERR, arg[6], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      r0[i][j]=r0_one;
      alpha[i][j]=alpha_one;
      sparam[i][j]=s_one;
      v0[i][j]=v0_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}
/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBV::init_style()
{

  neighbor->add_request(this);

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */
/* BVs */
double PairBV::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
    
  offset[i][j] = 0.0;
  alpha[j][i] = alpha[i][j];
  r0[j][i] = r0[i][j];
  sparam[j][i] = sparam[i][j];
  v0[j][i] = v0[i][j];
  offset[j][i] = offset[i][j];
    
  if(i == j){
    energy0[i]=0.0;
  }

  return cut[i][j];
}
/* BVs */

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBV::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
/* BVs */
        fwrite(&r0[i][j],sizeof(double),1,fp);
        fwrite(&alpha[i][j],sizeof(double),1,fp);
        fwrite(&sparam[i][j],sizeof(double),1,fp);
        fwrite(&v0[i][j],sizeof(double),1,fp);
/* BVs */
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBV::read_restart(FILE *fp)
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
          utils::sfread(FLERR, &r0[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &alpha[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &sparam[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &v0[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
/* BVs */
        MPI_Bcast(&r0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&alpha[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sparam[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&v0[i][j],1,MPI_DOUBLE,0,world);
/* BVs */
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBV::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&power_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBV::read_restart_settings(FILE *fp)
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

int PairBV::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int */*pbc*/)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = fp[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairBV::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) fp[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

int PairBV::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = s0[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void PairBV::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    s0[j] += buf[m++];
  }
}

