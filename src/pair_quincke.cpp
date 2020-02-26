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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "pair_quincke.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "compute_coord_atom.h"


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairQuincke::PairQuincke(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  //cvec(NULL), id_coord(NULL)
}

/* ---------------------------------------------------------------------- */

PairQuincke::~PairQuincke()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(nb_cutsq);

    memory->destroy(cut);
    memory->destroy(nb_cut);
    memory->destroy(a1);
    memory->destroy(a2);
    memory->destroy(a3);
    memory->destroy(offset);

    memory->destroy(nbvec);
  }
}

/* ---------------------------------------------------------------------- */

void PairQuincke::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl;
  double g11,g22,g12,g21,g1,g2,g3,delxsq,delysq,r2inv,rinv,tpairi,tpairj, fpair;
  double rsq,forceyukawa,factor;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double density_factor=1.0;
  int cnt;

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

  // loop over neighbors of my atoms to calculate number of nb's of each atom 
  if (density_flag){
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i]; 
    cnt = 0;
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      if (rsq < nb_cutsq[itype][jtype]) cnt ++;
    }
  nbvec[i] = cnt;  
  } 
  } 


  // loop over neighbors of my atoms 
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (density_flag) density_factor = 2.0/(pow(nbvec[i], nb_alpha) + 1.0);

    fpair = density_factor * activity; 
    f[i][0] += fpair * mu[i][0];
    f[i][1] += fpair * mu[i][1];
    // printf ("nbvec[%d]=%d, density_factor=%g, alpha=%g, nbcut=%g\n",i, nbvec[i], density_factor, nb_alpha, nb_cut[1][1]);
    // printf ("%d : fpair=%g density_factor=%g activity=%g nb=%d\n",i, fpair, density_factor, activity, nbvec[i]);


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
        g11 = mu[i][0]*mu[j][0];
        g22 = mu[i][1]*mu[j][1];
        g21 = mu[i][1]*mu[j][0];
        g12 = mu[i][0]*mu[j][1];
        delxsq = delx * delx;
        delysq = dely * dely;
        r2inv = 1.0/rsq;
        rinv  = sqrt(r2inv);

        g1  = a1[itype][jtype]*(g21 - g12);
        g2  = a2[itype][jtype]*(delx*mu[i][1] - dely*mu[i][0])*rinv;
        g3  = a3[itype][jtype]*( (delx-dely)*(delx+dely)*(g21+g12) - 2*delx*dely*(g11-g22)  )*r2inv;

        tpairi = density_factor*(g1 + g2 + g3);
        torque[i][2] += factor*tpairi;

        //printf ("%d %d: tij=%g, tji=%g\n",i,j, tpairi, tpairj);

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
          //printf("evdwl=%g\n", evdwl);
        }

        //if (force_flag){
        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                              evdwl,0.0,0.0,0.0,0.0, delx,dely,delz);
        //}
        // else{
        // if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
        //                       evdwl,0.0,0.0,0.0,0.0, delx,dely,delz);          

        // } 



        // if ((i==0 && j==1) || (i==1 && j==0) ){
        // double norm[2];
        // norm[0] = mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1];
        // norm[1] = mu[j][0]*mu[j][0] + mu[j][1]*mu[j][1];
        // printf ( "%d %d xi,j={%g %g}{%g %g} pi,j={%g %g}{%g %g},tij=%g, tji=%g, e=%g, eflag=%d\n"
        //         ,i,j, x[i][0], x[i][1], x[j][0], x[j][1], mu[i][0], mu[i][1], mu[j][0], mu[j][1], tpairi, tpairj, evdwl, eflag );
         
        // printf ( "{x,y}={%g, %g}; {pi1,pi2}={%g, %g}; {pj1,pj2}={%g, %g}; a1=%g; a2=%g;\nnewton=%d\n\n"
        //         ,delx, dely, mu[i][0], mu[i][1], mu[j][0], mu[j][1], a1[itype][jtype],a2[itype][jtype], newton_pair );
        // }
        // //printf ( "    %d: torque=(%.8f %.8f %.8f) tpair=%g \n",i ,torque[i][0], torque[i][1], torque[i][2], tpair );

        

      }

    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairQuincke::allocate()
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
  memory->create(nb_cutsq,n+1,n+1,"pair:nb_cutsq");
  memory->create(nb_cut,n+1,n+1,"pair:nb_cut");
  memory->create(a1,n+1,n+1,"pair:a1");
  memory->create(a2,n+1,n+1,"pair:a2");
  memory->create(a3,n+1,n+1,"pair:a3");
  memory->create(offset,n+1,n+1,"pair:offset");
 
  memory->create(nbvec,natom+1,"pair:nbvec");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairQuincke::settings(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal pair_style command");
  if (narg == 1) {
      cut_global = force->numeric(FLERR,arg[0]);  
      activity = force->numeric(FLERR,arg[1]);  
      nb_alpha = 0;
      density_flag = 0;
  }
  else if (narg == 4 ) {
      cut_global = force->numeric(FLERR,arg[0]);  
      activity = force->numeric(FLERR,arg[1]);  
      nb_alpha = force->numeric(FLERR,arg[2]); 
      nb_cut_global = force->numeric(FLERR,arg[3]);  
      density_flag = 1;
      if ( abs(nb_alpha-0.0) < 0.000001) density_flag = 0;
      //if (cut_global<nb_cut_global)  error->all(FLERR,"Illegal pair_style: rcut must be larger than nb_rcut.");
  }
  else error->all(FLERR,"Illegal pair_style command");


  // else if (narg == 3 ) {
  //     int n = strlen(arg[1]) + 1;
  //     id_coord = new char[n];
  //     strcpy(id_coord,arg[1]);
  //     int icoord = modify->find_compute(id_coord);
     
  //     if (icoord < 0) 
  //         error->all(FLERR,"Could not find compute coord/atom compute ID");
  //     else if (strcmp(modify->compute[icoord]->style,"coord/atom") != 0)
  //         error->all(FLERR,"Pair Quincke compute ID is not coord/atom");
  //     else {
  //         exponent = force->numeric(FLERR,arg[2]); 
  //         density_flag = 1;
  //     }    
  // }
  // else error->all(FLERR,"Illegal pair_style command");

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut[i][j] = cut_global;
          nb_cut[i][j] = nb_cut_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairQuincke::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 8)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double a1_one = force->numeric(FLERR,arg[2]);
  double a2_one = force->numeric(FLERR,arg[3]);
  double a3_one = force->numeric(FLERR,arg[4]);

  double cut_one = cut_global;
  double nb_cut_one = nb_cut_global;
  if (narg >= 6) cut_one = force->numeric(FLERR,arg[5]);
  if (narg >= 7) nb_cut_one = force->numeric(FLERR,arg[6]);
  //if (cut_one<nb_cut_one)  error->all(FLERR,"rcut must be larger than nb_rcut.");


  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a1[i][j] = a1_one;
      a2[i][j] = a2_one;
      a3[i][j] = a3_one;
      cut[i][j] = cut_one;
      nb_cut[i][j] = nb_cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairQuincke::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    a1[i][j] = mix_energy(a1[i][i],a1[j][j],1.0,1.0);
    a2[i][j] = mix_energy(a2[i][i],a2[j][j],1.0,1.0);
    a3[i][j] = mix_energy(a3[i][i],a3[j][j],1.0,1.0);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
    nb_cut[i][j] = mix_distance(nb_cut[i][i],nb_cut[j][j]);
  }

  nb_cutsq[i][j] = nb_cut[i][j]*nb_cut[i][j];


  if (offset_flag) {
    //double screening = exp(-kappa * cut[i][j]);
    offset[i][j] = 0.0;
  } else offset[i][j] = 0.0;

  a1[j][i] = a1[i][j];
  a2[j][i] = a2[i][j];
  a3[j][i] = a3[i][j];
  cut[j][i] = cut[i][j];
  nb_cut[j][i] = nb_cut[i][j];

  offset[j][i] = offset[i][j];

  return cut[i][j];
}


void PairQuincke::init_style()
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

void PairQuincke::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&a1[i][j],sizeof(double),1,fp);
        fwrite(&a2[i][j],sizeof(double),1,fp);
        fwrite(&a3[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
        fwrite(&nb_cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairQuincke::read_restart(FILE *fp)
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
          fread(&a1[i][j],sizeof(double),1,fp);
          fread(&a2[i][j],sizeof(double),1,fp);
          fread(&a3[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
          fread(&nb_cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&a1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&nb_cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairQuincke::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&activity,sizeof(double),1,fp);  
  fwrite(&nb_alpha,sizeof(double),1,fp);
  fwrite(&nb_cut_global,sizeof(double),1,fp);
  fwrite(&density_flag,sizeof(int),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairQuincke::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&activity,sizeof(double),1,fp);  
    fread(&nb_alpha,sizeof(double),1,fp);
    fread(&nb_cut_global,sizeof(double),1,fp);
    fread(&density_flag,sizeof(int),1,fp);    
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&activity,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&nb_alpha,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&nb_cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&density_flag,1,MPI_INT,0,world);  
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairQuincke::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g %g\n",i,a1[i][i],a2[i][i],a3[i][i],cut[i][i],nb_cut[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairQuincke::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g\n",i,j,a1[i][j],a2[i][j],a3[i][i],cut[i][j],nb_cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairQuincke::single(int i, int j, int itype, int jtype, double rsq,
                          double factor_coul, double factor_lj,
                          double &fforce)
{
  
  //double r2inv,r,rinv,screening,forceyukawa,phi;

  //r2inv = 1.0/rsq;
  //r = sqrt(rsq);
  //rinv = 1.0/r;
  //screening = exp(-kappa*r);
  //forceyukawa = a[itype][jtype] * screening * (kappa + rinv);
  //fforce = factor_lj*forceyukawa * r2inv;

  double phi;
  fforce = 0.0;

  //printf("this is a test*********************\n");

  phi = 0.0;
  return factor_lj*phi;
}
