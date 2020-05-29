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

#include "pair_morse_soft.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "math_special.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace MathSpecial;

/* ----------------------------------------------------------------------
   Contributing author: Stefan Paquay (TU/e)
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

PairMorseSoft::~PairMorseSoft()
{
  if(allocated){
    memory->destroy(lambda);
  }
}

/* ---------------------------------------------------------------------- */

void PairMorseSoft::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,dr,dexp,dexp2,dexp3,factor_lj;
  double ea,phi,V0,iea2;
  double D, a, x0, l, B, s1, llf;

  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
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

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        dr = r - r0[itype][jtype];

        D  = d0[itype][jtype];
        a  = alpha[itype][jtype];
        x0 = r0[itype][jtype];
        dexp = exp( -a * dr );
        dexp2 = dexp*dexp;
        dexp3 = dexp2*dexp;

        l = lambda[itype][jtype];

        ea  = exp( a * x0 );
        iea2 = exp( -2.*a*x0 );

        V0 = D * dexp * ( dexp - 2.0 );
        B = -2.0 * D * iea2 * ( ea - 1.0 ) / 3.0;

        if (l >= shift_range){
          s1  = (l - 1.0) / (shift_range - 1.0);
          phi = V0 + B*dexp3 * s1;

          // Force computation:
          fpair = 3.0*a*B*dexp3*s1 + 2.0*a*D*(dexp2 - dexp);
          fpair /= r;
        } else {
          llf = MathSpecial::powint( l / shift_range, nlambda );
          phi = V0 + B*dexp3;
          phi *= llf;

          // Force computation:
          if (r == 0.0){
            fpair = 0.0;
          } else {
            fpair = 3.0*a*B*dexp3 + 2.0*a*D*(dexp2 - dexp);
            fpair *= llf / r;
          }
        }

        fpair *= factor_lj;


        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;

        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = phi*factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMorseSoft::allocate()
{
  PairMorse::allocate();
  int n = atom->ntypes;
  memory->create(lambda,n+1,n+1,"pair:lambda");

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMorseSoft::coeff(int narg, char **arg)
{
  if (narg < 6 || narg > 7) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double d0_one     = force->numeric(FLERR,arg[2]);
  double alpha_one  = force->numeric(FLERR,arg[3]);
  double r0_one     = force->numeric(FLERR,arg[4]);
  double lambda_one = force->numeric(FLERR,arg[5]);

  double cut_one = cut_global;
  if (narg == 7) cut_one = force->numeric(FLERR,arg[6]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      d0[i][j]      = d0_one;
      alpha[i][j]   = alpha_one;
      r0[i][j]      = r0_one;
      lambda[i][j]  = lambda_one;
      cut[i][j]     = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   Set global stuff.
------------------------------------------------------------------------- */

void PairMorseSoft::settings(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal pair_style command");

  nlambda     = force->inumeric(FLERR,arg[0]);
  shift_range = force->numeric(FLERR,arg[1]);
  cut_global  = force->numeric(FLERR,arg[2]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMorseSoft::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    error->all(FLERR,"All pair coeffs are not set");

  morse1[i][j] = 2.0*d0[i][j]*alpha[i][j];

  if (offset_flag) {
    double l, s1, V0, B, llf;
    double alpha_dr = -alpha[i][j] * (cut[i][j] - r0[i][j]);
    double D  = d0[i][j];
    double a  = alpha[i][j];
    double x0 = r0[i][j];
    double dexp  = exp( alpha_dr );
    double dexp2 = dexp*dexp;
    double dexp3 = dexp2*dexp;

    l = lambda[i][j];

    double ea  = exp( a*x0 );
    double iea2 = exp( -2.*a*x0 );

    V0 = D * dexp * ( dexp - 2.0 );
    B = -2.0 * D * iea2 * ( ea - 1.0 ) / 3.0;

    if (l >= shift_range){
      s1  = (l - 1.0) / (shift_range - 1.0);
      offset[i][j] = V0 + B*dexp3 * s1;
    } else {
      llf = MathSpecial::powint( l / shift_range, nlambda );
      offset[i][j] = V0 + B*dexp3;
      offset[i][j] *= llf;
    }
  } else offset[i][j] = 0.0;

  d0[j][i]     = d0[i][j];
  alpha[j][i]  = alpha[i][j];
  r0[j][i]     = r0[i][j];
  morse1[j][i] = morse1[i][j];
  lambda[j][i] = lambda[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMorseSoft::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&d0[i][j],sizeof(double),1,fp);
        fwrite(&alpha[i][j],sizeof(double),1,fp);
        fwrite(&r0[i][j],sizeof(double),1,fp);
        fwrite(&lambda[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMorseSoft::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,NULL,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&d0[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&alpha[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&r0[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&lambda[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,NULL,error);
        }
        MPI_Bcast(&d0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&alpha[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&r0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&lambda[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
  }
}


/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairMorseSoft::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,d0[i][i],alpha[i][i],r0[i][i],
            lambda[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairMorseSoft::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %g %g %g %g\n",i,d0[i][j],alpha[i][j],r0[i][j],
              lambda[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairMorseSoft::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                             double /*factor_coul*/, double factor_lj,
                             double &fforce)
{
  double r, dr, dexp, dexp2, dexp3, phi;
  double B, D, a, ea, iea2;
  double x0, V0, s1, l, llf;

  D  = d0[itype][jtype];
  a  = alpha[itype][jtype];
  x0 = r0[itype][jtype];
  r = sqrt(rsq);
  dr = r - r0[itype][jtype];
  dexp = exp( -a * dr );
  dexp2 = dexp*dexp;
  dexp3 = dexp2*dexp;

  l = lambda[itype][jtype];

  ea  = exp( a * x0 );
  iea2 = exp( -2.*a*x0 );

  V0 = D * dexp * ( dexp - 2.0 );
  B = -2.0 * D * iea2 * ( ea - 1.0 ) / 3.0;

  if (l >= shift_range){
    s1  = (l - 1.0) / (shift_range - 1.0);
    phi = V0 + B*dexp3 * s1;

    // Force computation:
    fforce = 3.0*a*B*dexp3*s1 + 2.0*a*D*(dexp2 - dexp);
    fforce /= r;
  } else {
    llf = MathSpecial::powint( l / shift_range, nlambda );
    phi = V0 + B*dexp3;
    phi *= llf;

    // Force computation:
    if (r == 0.0){
      fforce = 0.0;
    } else {
      fforce = 3.0*a*B*dexp3 + 2.0*a*D*(dexp2 - dexp);
      fforce *= llf / r;
    }
  }

  fforce *= factor_lj;
  phi -= offset[itype][jtype];
  return factor_lj*phi;
}

/* ---------------------------------------------------------------------- */

void *PairMorseSoft::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"d0") == 0) return (void *) d0;
  if (strcmp(str,"r0") == 0) return (void *) r0;
  if (strcmp(str,"alpha") == 0) return (void *) alpha;
  if (strcmp(str,"lambda") == 0) return (void *) lambda;
  return NULL;
}
