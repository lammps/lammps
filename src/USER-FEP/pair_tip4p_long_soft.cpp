// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
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
   Soft-core version: Agilio Padua (ENS de Lyon & CNRS)
------------------------------------------------------------------------- */

#include "pair_tip4p_long_soft.h"

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


using namespace LAMMPS_NS;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairTIP4PLongSoft::PairTIP4PLongSoft(LAMMPS *lmp) : PairCoulLongSoft(lmp)
{
  tip4pflag = 1;
  single_enable = 0;
  respa_enable = 0;

  nmax = 0;
  hneigh = nullptr;
  newsite = nullptr;

  // TIP4P cannot compute virial as F dot r
  // due to finding bonded H atoms which are not near O atom

  no_virial_fdotr_compute = 1;
}

/* ---------------------------------------------------------------------- */

PairTIP4PLongSoft::~PairTIP4PLongSoft()
{
  memory->destroy(hneigh);
  memory->destroy(newsite);
}

/* ---------------------------------------------------------------------- */

void PairTIP4PLongSoft::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,key;
  int n,vlist[6];
  int iH1,iH2,jH1,jH2;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul;
  double r,forcecoul,cforce;
  double factor_coul;
  double grij,expm2,prefactor,t,erfc;
  double denc;
  double fO[3],fH[3],fd[3],v[6];
  double *x1,*x2,*xH1,*xH2;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq;

  ecoul = 0.0;
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
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

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

          denc = sqrt(lam2[itype][jtype] + rsq);
          prefactor = qqrd2e * lam1[itype][jtype] * qtmp*q[j] / (denc*denc*denc);

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
            prefactor = qqrd2e * lam1[itype][jtype] * qtmp*q[j] / denc;
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

void PairTIP4PLongSoft::settings(int narg, char **arg)
{
  if (narg != 8) error->all(FLERR,"Illegal pair_style command");

  typeO = utils::inumeric(FLERR,arg[0],false,lmp);
  typeH = utils::inumeric(FLERR,arg[1],false,lmp);
  typeB = utils::inumeric(FLERR,arg[2],false,lmp);
  typeA = utils::inumeric(FLERR,arg[3],false,lmp);
  qdist = utils::numeric(FLERR,arg[4],false,lmp);

  nlambda = utils::numeric(FLERR,arg[5],false,lmp);
  alphac  = utils::numeric(FLERR,arg[6],false,lmp);

  cut_coul = utils::numeric(FLERR,arg[7],false,lmp);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairTIP4PLongSoft::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style tip4p/long requires atom IDs");
  if (!force->newton_pair)
    error->all(FLERR,
               "Pair style tip4p/long requires newton pair on");
  if (!atom->q_flag)
    error->all(FLERR,
               "Pair style tip4p/long requires atom attribute q");
  if (force->bond == nullptr)
    error->all(FLERR,"Must use a bond style with TIP4P potential");
  if (force->angle == nullptr)
    error->all(FLERR,"Must use an angle style with TIP4P potential");

  PairCoulLongSoft::init_style();

  // set alpha parameter

  double theta = force->angle->equilibrium_angle(typeA);
  double blen = force->bond->equilibrium_distance(typeB);
  alpha = qdist / (cos(0.5*theta) * blen);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairTIP4PLongSoft::init_one(int i, int j)
{
  return PairCoulLongSoft::init_one(i,j);
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTIP4PLongSoft::write_restart_settings(FILE *fp)
{
  PairCoulLongSoft::write_restart_settings(fp);

  fwrite(&typeO,sizeof(int),1,fp);
  fwrite(&typeH,sizeof(int),1,fp);
  fwrite(&typeB,sizeof(int),1,fp);
  fwrite(&typeA,sizeof(int),1,fp);
  fwrite(&qdist,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTIP4PLongSoft::read_restart_settings(FILE *fp)
{
  PairCoulLongSoft::read_restart_settings(fp);

  if (comm->me == 0) {
    utils::sfread(FLERR,&typeO,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&typeH,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&typeB,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&typeA,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&qdist,sizeof(double),1,fp,nullptr,error);
  }
  MPI_Bcast(&typeO,1,MPI_INT,0,world);
  MPI_Bcast(&typeH,1,MPI_INT,0,world);
  MPI_Bcast(&typeB,1,MPI_INT,0,world);
  MPI_Bcast(&typeA,1,MPI_INT,0,world);
  MPI_Bcast(&qdist,1,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
  compute position xM of fictitious charge site for O atom and 2 H atoms
  return it as xM
------------------------------------------------------------------------- */

void PairTIP4PLongSoft::compute_newsite(double *xO, double *xH1,
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

void *PairTIP4PLongSoft::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"qdist") == 0) return (void *) &qdist;
  if (strcmp(str,"typeO") == 0) return (void *) &typeO;
  if (strcmp(str,"typeH") == 0) return (void *) &typeH;
  if (strcmp(str,"typeA") == 0) return (void *) &typeA;
  if (strcmp(str,"typeB") == 0) return (void *) &typeB;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul;
  dim = 2;
  if (strcmp(str,"lambda") == 0) return (void *) lambda;
  return nullptr;
}

/* ----------------------------------------------------------------------
   memory usage of hneigh
------------------------------------------------------------------------- */

double PairTIP4PLongSoft::memory_usage()
{
  double bytes = (double)maxeatom * sizeof(double);
  bytes += (double)maxvatom*6 * sizeof(double);
  bytes += (double)2 * nmax * sizeof(double);
  return bytes;
}
