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
   Contributing authors: Amalie Frischknecht and Ahmed Ismail (SNL)
   simpler force assignment added by Rolf Isele-Holder (Aachen University)
   Soft-core version: Agilio Padua (Univ Blaise Pascal & CNRS)
------------------------------------------------------------------------- */

#include "pair_lj_cut_tip4p_long_soft.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairLJCutTIP4PLongSoft::PairLJCutTIP4PLongSoft(LAMMPS *lmp) :
  PairLJCutCoulLongSoft(lmp)
{
  tip4pflag = 1;
  single_enable = 0;
  respa_enable = 0;

  nmax = 0;
  hneigh = NULL;
  newsite = NULL;

  // TIP4P cannot compute virial as F dot r
  // due to finding bonded H atoms which are not near O atom

  no_virial_fdotr_compute = 1;
}

/* ---------------------------------------------------------------------- */

PairLJCutTIP4PLongSoft::~PairLJCutTIP4PLongSoft()
{
  memory->destroy(hneigh);
  memory->destroy(newsite);
}

/* ---------------------------------------------------------------------- */

void PairLJCutTIP4PLongSoft::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,key;
  int n,vlist[6];
  int iH1,iH2,jH1,jH2;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul;
  double r,forcecoul,forcelj,cforce;
  double factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  double denc, denlj, r4sig6;
  double fO[3],fH[3],fd[3],v[6];
  double *x1,*x2,*xH1,*xH2;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  // reallocate hneigh & newsite if necessary
  // initialize hneigh[0] to -1 on steps when reneighboring occurred
  // initialize hneigh[2] to 0 every step

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(hneigh);
    memory->create(hneigh,nmax,3,"pair:hneigh");
    memory->destroy(newsite);
    memory->create(newsite,nmax,3,"pair:newsite");
  }
  if (neighbor->ago == 0)
    for (i = 0; i < nall; i++) hneigh[i][0] = -1;
  for (i = 0; i < nall; i++) hneigh[i][2] = 0;

  double **f = atom->f;
  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;
  double cut_coulsqplus = (cut_coul+2.0*qdist) * (cut_coul+2.0*qdist);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    // if atom I = water O, set x1 = offset charge site
    // else x1 = x of atom I

    if (itype == typeO) {
      if (hneigh[i][0] < 0) {
        iH1 = atom->map(atom->tag[i] + 1);
        iH2 = atom->map(atom->tag[i] + 2);
        if (iH1 == -1 || iH2 == -1)
          error->one(FLERR,"TIP4P hydrogen is missing");
        if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
          error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
        // set iH1,iH2 to closest image to O
        iH1 = domain->closest_image(i,iH1);
        iH2 = domain->closest_image(i,iH2);
        compute_newsite(x[i],x[iH1],x[iH2],newsite[i]);
        hneigh[i][0] = iH1;
        hneigh[i][1] = iH2;
        hneigh[i][2] = 1;

      } else {
        iH1 = hneigh[i][0];
        iH2 = hneigh[i][1];
        if (hneigh[i][2] == 0) {
          hneigh[i][2] = 1;
          compute_newsite(x[i],x[iH1],x[iH2],newsite[i]);
        }
      }
      x1 = newsite[i];
    } else x1 = x[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      // LJ interaction based on true rsq

      if (rsq < cut_ljsq[itype][jtype]) {

        r4sig6 = rsq*rsq / lj2[itype][jtype];
        denlj = lj3[itype][jtype] + rsq*r4sig6;
        forcelj = lj1[itype][jtype] * epsilon[itype][jtype] *
          (48.0*r4sig6/(denlj*denlj*denlj) - 24.0*r4sig6/(denlj*denlj));

        forcelj *= factor_lj;

        f[i][0] += delx*forcelj;
        f[i][1] += dely*forcelj;
        f[i][2] += delz*forcelj;
        f[j][0] -= delx*forcelj;
        f[j][1] -= dely*forcelj;
        f[j][2] -= delz*forcelj;

        if (eflag) {
          evdwl = lj1[itype][jtype] * 4.0 * epsilon[itype][jtype] *
            (1.0/(denlj*denlj) - 1.0/denlj) - offset[itype][jtype];
          evdwl *= factor_lj;
        } else evdwl = 0.0;

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,forcelj,delx,dely,delz);
      }

      // adjust rsq and delxyz for off-site O charge(s) if necessary
      // but only if they are within reach

      if (rsq < cut_coulsqplus) {
        if (itype == typeO || jtype == typeO) {

          // if atom J = water O, set x2 = offset charge site
          // else x2 = x of atom J

          if (jtype == typeO) {
            if (hneigh[j][0] < 0) {
              jH1 = atom->map(atom->tag[j] + 1);
              jH2 = atom->map(atom->tag[j] + 2);
              if (jH1 == -1 || jH2 == -1)
                error->one(FLERR,"TIP4P hydrogen is missing");
              if (atom->type[jH1] != typeH || atom->type[jH2] != typeH)
                error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
              // set jH1,jH2 to closest image to O
              jH1 = domain->closest_image(j,jH1);
              jH2 = domain->closest_image(j,jH2);
              compute_newsite(x[j],x[jH1],x[jH2],newsite[j]);
              hneigh[j][0] = jH1;
              hneigh[j][1] = jH2;
              hneigh[j][2] = 1;

            } else {
              jH1 = hneigh[j][0];
              jH2 = hneigh[j][1];
              if (hneigh[j][2] == 0) {
                hneigh[j][2] = 1;
                compute_newsite(x[j],x[jH1],x[jH2],newsite[j]);
              }
            }
            x2 = newsite[j];
          } else x2 = x[j];

          delx = x1[0] - x2[0];
          dely = x1[1] - x2[1];
          delz = x1[2] - x2[2];
          rsq = delx*delx + dely*dely + delz*delz;
        }

        // Coulombic interaction based on modified rsq

        if (rsq < cut_coulsq) {
          r = sqrt(rsq);
          grij = g_ewald * r;
          expm2 = exp(-grij*grij);
          t = 1.0 / (1.0 + EWALD_P*grij);
          erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;

          denc = sqrt(lj4[itype][jtype] + rsq);
          prefactor = qqrd2e * lj1[itype][jtype] * qtmp*q[j] / (denc*denc*denc);

          forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
          if (factor_coul < 1.0) {
            forcecoul -= (1.0-factor_coul)*prefactor;
          }

          cforce = forcecoul;

          // if i,j are not O atoms, force is applied directly
          // if i or j are O atoms, force is on fictitious atom & partitioned
          // force partitioning due to Feenstra, J Comp Chem, 20, 786 (1999)
          // f_f = fictitious force, fO = f_f (1 - 2 alpha), fH = alpha f_f
          // preserves total force and torque on water molecule
          // virial = sum(r x F) where each water's atoms are near xi and xj
          // vlist stores 2,4,6 atoms whose forces contribute to virial

          n = 0;
          key = 0;

          if (itype != typeO) {
            f[i][0] += delx * cforce;
            f[i][1] += dely * cforce;
            f[i][2] += delz * cforce;

            if (vflag) {
              v[0] = x[i][0] * delx * cforce;
              v[1] = x[i][1] * dely * cforce;
              v[2] = x[i][2] * delz * cforce;
              v[3] = x[i][0] * dely * cforce;
              v[4] = x[i][0] * delz * cforce;
              v[5] = x[i][1] * delz * cforce;
            }
            vlist[n++] = i;

          } else {
            key++;

            fd[0] = delx*cforce;
            fd[1] = dely*cforce;
            fd[2] = delz*cforce;

            fO[0] = fd[0]*(1 - alpha);
            fO[1] = fd[1]*(1 - alpha);
            fO[2] = fd[2]*(1 - alpha);

            fH[0] = 0.5 * alpha * fd[0];
            fH[1] = 0.5 * alpha * fd[1];
            fH[2] = 0.5 * alpha * fd[2];

            f[i][0] += fO[0];
            f[i][1] += fO[1];
            f[i][2] += fO[2];

            f[iH1][0] += fH[0];
            f[iH1][1] += fH[1];
            f[iH1][2] += fH[2];

            f[iH2][0] += fH[0];
            f[iH2][1] += fH[1];
            f[iH2][2] += fH[2];

            if (vflag) {
              xH1 = x[iH1];
              xH2 = x[iH2];
              v[0] = x[i][0]*fO[0] + xH1[0]*fH[0] + xH2[0]*fH[0];
              v[1] = x[i][1]*fO[1] + xH1[1]*fH[1] + xH2[1]*fH[1];
              v[2] = x[i][2]*fO[2] + xH1[2]*fH[2] + xH2[2]*fH[2];
              v[3] = x[i][0]*fO[1] + xH1[0]*fH[1] + xH2[0]*fH[1];
              v[4] = x[i][0]*fO[2] + xH1[0]*fH[2] + xH2[0]*fH[2];
              v[5] = x[i][1]*fO[2] + xH1[1]*fH[2] + xH2[1]*fH[2];
            }
            vlist[n++] = i;
            vlist[n++] = iH1;
            vlist[n++] = iH2;
          }

          if (jtype != typeO) {
            f[j][0] -= delx * cforce;
            f[j][1] -= dely * cforce;
            f[j][2] -= delz * cforce;

            if (vflag) {
              v[0] -= x[j][0] * delx * cforce;
              v[1] -= x[j][1] * dely * cforce;
              v[2] -= x[j][2] * delz * cforce;
              v[3] -= x[j][0] * dely * cforce;
              v[4] -= x[j][0] * delz * cforce;
              v[5] -= x[j][1] * delz * cforce;
            }
            vlist[n++] = j;

          } else {
            key += 2;

            fd[0] = -delx*cforce;
            fd[1] = -dely*cforce;
            fd[2] = -delz*cforce;

            fO[0] = fd[0]*(1 - alpha);
            fO[1] = fd[1]*(1 - alpha);
            fO[2] = fd[2]*(1 - alpha);

            fH[0] = 0.5 * alpha * fd[0];
            fH[1] = 0.5 * alpha * fd[1];
            fH[2] = 0.5 * alpha * fd[2];

            f[j][0] += fO[0];
            f[j][1] += fO[1];
            f[j][2] += fO[2];

            f[jH1][0] += fH[0];
            f[jH1][1] += fH[1];
            f[jH1][2] += fH[2];

            f[jH2][0] += fH[0];
            f[jH2][1] += fH[1];
            f[jH2][2] += fH[2];

            if (vflag) {
              xH1 = x[jH1];
              xH2 = x[jH2];
              v[0] += x[j][0]*fO[0] + xH1[0]*fH[0] + xH2[0]*fH[0];
              v[1] += x[j][1]*fO[1] + xH1[1]*fH[1] + xH2[1]*fH[1];
              v[2] += x[j][2]*fO[2] + xH1[2]*fH[2] + xH2[2]*fH[2];
              v[3] += x[j][0]*fO[1] + xH1[0]*fH[1] + xH2[0]*fH[1];
              v[4] += x[j][0]*fO[2] + xH1[0]*fH[2] + xH2[0]*fH[2];
              v[5] += x[j][1]*fO[2] + xH1[1]*fH[2] + xH2[1]*fH[2];
            }
            vlist[n++] = j;
            vlist[n++] = jH1;
            vlist[n++] = jH2;
          }

          if (eflag) {
            prefactor = qqrd2e * lj1[itype][jtype] * qtmp*q[j] / denc;
            ecoul = prefactor*erfc;
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;

          if (evflag) ev_tally_tip4p(key,vlist,v,ecoul,alpha);
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJCutTIP4PLongSoft::settings(int narg, char **arg)
{
  if (narg < 9 || narg > 10) error->all(FLERR,"Illegal pair_style command");

  typeO = force->inumeric(FLERR,arg[0]);
  typeH = force->inumeric(FLERR,arg[1]);
  typeB = force->inumeric(FLERR,arg[2]);
  typeA = force->inumeric(FLERR,arg[3]);
  qdist = force->numeric(FLERR,arg[4]);
  nlambda = force->numeric(FLERR,arg[5]);
  alphalj = force->numeric(FLERR,arg[6]);
  alphac = force->numeric(FLERR,arg[7]);

  cut_lj_global = force->numeric(FLERR,arg[8]);
  if (narg == 9) cut_coul = cut_lj_global;
  else cut_coul = force->numeric(FLERR,arg[9]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut_lj[i][j] = cut_lj_global;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutTIP4PLongSoft::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style lj/cut/tip4p/long/soft requires atom IDs");
  if (!force->newton_pair)
    error->all(FLERR,
               "Pair style lj/cut/tip4p/long/soft requires newton pair on");
  if (!atom->q_flag)
    error->all(FLERR,
               "Pair style lj/cut/tip4p/long/soft requires atom attribute q");
  if (force->bond == NULL)
    error->all(FLERR,"Must use a bond style with TIP4P potential");
  if (force->angle == NULL)
    error->all(FLERR,"Must use an angle style with TIP4P potential");

  PairLJCutCoulLongSoft::init_style();

  // set alpha parameter

  double theta = force->angle->equilibrium_angle(typeA);
  double blen = force->bond->equilibrium_distance(typeB);
  alpha = qdist / (cos(0.5*theta) * blen);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCutTIP4PLongSoft::init_one(int i, int j)
{
  double cut = PairLJCutCoulLongSoft::init_one(i,j);

  // check that LJ epsilon = 0.0 for water H
  // set LJ cutoff to 0.0 for any interaction involving water H
  // so LJ term isn't calculated in compute()

  if ((i == typeH && epsilon[i][i] != 0.0) ||
      (j == typeH && epsilon[j][j] != 0.0))
    error->all(FLERR,"Water H epsilon must be 0.0 for "
               "pair style lj/cut/tip4p/long/soft");

  if (i == typeH || j == typeH)
    cut_ljsq[j][i] = cut_ljsq[i][j] = 0.0;

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutTIP4PLongSoft::write_restart_settings(FILE *fp)
{
  PairLJCutCoulLongSoft::write_restart_settings(fp);

  fwrite(&typeO,sizeof(int),1,fp);
  fwrite(&typeH,sizeof(int),1,fp);
  fwrite(&typeB,sizeof(int),1,fp);
  fwrite(&typeA,sizeof(int),1,fp);
  fwrite(&qdist,sizeof(double),1,fp);

  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutTIP4PLongSoft::read_restart_settings(FILE *fp)
{
  PairLJCutCoulLongSoft::read_restart_settings(fp);

  if (comm->me == 0) {
    utils::sfread(FLERR,&typeO,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&typeH,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&typeB,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&typeA,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&qdist,sizeof(double),1,fp,NULL,error);

    utils::sfread(FLERR,&cut_lj_global,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&cut_coul,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&tail_flag,sizeof(int),1,fp,NULL,error);
  }

  MPI_Bcast(&typeO,1,MPI_INT,0,world);
  MPI_Bcast(&typeH,1,MPI_INT,0,world);
  MPI_Bcast(&typeB,1,MPI_INT,0,world);
  MPI_Bcast(&typeA,1,MPI_INT,0,world);
  MPI_Bcast(&qdist,1,MPI_DOUBLE,0,world);

  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
  compute position xM of fictitious charge site for O atom and 2 H atoms
  return it as xM
------------------------------------------------------------------------- */

void PairLJCutTIP4PLongSoft::compute_newsite(double *xO, double *xH1,
                                         double *xH2, double *xM)
{
  double delx1 = xH1[0] - xO[0];
  double dely1 = xH1[1] - xO[1];
  double delz1 = xH1[2] - xO[2];

  double delx2 = xH2[0] - xO[0];
  double dely2 = xH2[1] - xO[1];
  double delz2 = xH2[2] - xO[2];

  xM[0] = xO[0] + alpha * 0.5 * (delx1 + delx2);
  xM[1] = xO[1] + alpha * 0.5 * (dely1 + dely2);
  xM[2] = xO[2] + alpha * 0.5 * (delz1 + delz2);
}

/* ---------------------------------------------------------------------- */

void *PairLJCutTIP4PLongSoft::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"qdist") == 0) return (void *) &qdist;
  if (strcmp(str,"typeO") == 0) return (void *) &typeO;
  if (strcmp(str,"typeH") == 0) return (void *) &typeH;
  if (strcmp(str,"typeA") == 0) return (void *) &typeA;
  if (strcmp(str,"typeB") == 0) return (void *) &typeB;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul;
  dim = 2;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  if (strcmp(str,"lambda") == 0) return (void *) lambda;
  return NULL;
}

/* ----------------------------------------------------------------------
   memory usage of hneigh
------------------------------------------------------------------------- */

double PairLJCutTIP4PLongSoft::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  bytes += 2 * nmax * sizeof(double);
  return bytes;
}
