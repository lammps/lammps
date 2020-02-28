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

/*-----------------------------------------------------------------------

  pair_activity.cpp file for simplest case active brownian particles,
  where the forces are just self-propelled with some activity coefficient.

  File initially started by Sam Cameron, 28/02/20.

  -----------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "pair_activity.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "compute_coord_atom.h"


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairActivity::PairActivity(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  //cvec(NULL), id_coord(NULL)
}

/* ---------------------------------------------------------------------- */

PairActivity::~PairActivity()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(offset);

  }
}

/* ---------------------------------------------------------------------- */

void PairActivity::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl;
  double delxsq,delysq,r2inv,rinv, fpair;
  double rsq,forceyukawa,factor;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  double **mu = atom->mu;

  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;



  // loop over neighbors of my atoms 
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    fpair = activity; 
    f[i][0] += fpair * mu[i][0];
    f[i][1] += fpair * mu[i][1];
    /*

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor = special_lj[sbmask(j)];
      j &= NEIGHMASK;
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        delxsq = delx * delx;
        delysq = dely * dely;
        r2inv = 1.0/rsq;
        rinv  = sqrt(r2inv);

        if (newton_pair || j < nlocal) {
          g2  = a2[itype][jtype]*(-delx*mu[j][1] + dely*mu[j][0])*rinv;
          tpairj = density_factor*(-g1 + g2 + g3);          
          torque[j][2] += factor*tpairj;
        }

        if (eflag) {
          g1  = a1[itype][jtype]*(g11+g22);
          g2  = a2[itype][jtype]*0.5*( delx*mu[i][0] + dely*mu[i][1] - delx*mu[j][0] - dely*mu[j][1] )*rinv;    
          g3  = a3[itype][jtype]*( (delx-dely)*(delx+dely)*(g11-g22) + 2*delx*dely*(g21+g12)  )*r2inv;          
          evdwl = - density_factor*(g1 + g2 + g3);
          evdwl *= factor;
        }


        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                              evdwl,0.0,0.0,0.0,0.0, delx,dely,delz);
  

      }

    }
    */
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairActivity::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  int natom = atom->nmax;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(offset,n+1,n+1,"pair:offset");
 
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairActivity::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);  
  activity = force->numeric(FLERR,arg[1]);  

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut[i][j] = cut_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairActivity::coeff(int narg, char **arg)
{
  if (narg > 3)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);


  double cut_one = cut_global;

  if (narg >= 3) cut_one = force->numeric(FLERR,arg[2]);


  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairActivity::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }


  if (offset_flag) {
    offset[i][j] = 0.0;
  } else offset[i][j] = 0.0;

  cut[j][i] = cut[i][j];

  offset[j][i] = offset[i][j];

  return cut[i][j];
}


void PairActivity::init_style()
{
  int irequest = neighbor->request(this,instance_me);
  int newton_pair = force->newton_pair;
  if (newton_pair){ 
  neighbor->requests[irequest]->half = 1;
  neighbor->requests[irequest]->full = 0;
  //printf("build half list.\n");
  } else{
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;  
  //printf("build full list.\n");
  }
}


/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairActivity::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairActivity::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairActivity::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&activity,sizeof(double),1,fp);  
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairActivity::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&activity,sizeof(double),1,fp);  
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&activity,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairActivity::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g\n",i,cut[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairActivity::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g\n",i,j,cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairActivity::single(int i, int j, int itype, int jtype, double rsq,
                          double factor_coul, double factor_lj,
                          double &fforce)
{
  

  double phi;
  fforce = 0.0;

  phi = 0.0;
  return factor_lj*phi;
}
