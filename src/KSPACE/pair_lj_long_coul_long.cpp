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
   Contributing author: Pieter J. in 't Veld (SNL)
   Tabulation for long-range dispersion added by Wayne Mitchell (Loyola
   University New Orleans)
------------------------------------------------------------------------- */

#include "pair_lj_long_coul_long.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathExtra;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairLJLongCoulLong::PairLJLongCoulLong(LAMMPS *lmp) : Pair(lmp)
{
  dispersionflag = ewaldflag = pppmflag = 1;
  respa_enable = 1;
  writedata = 1;
  ftable = nullptr;
  fdisptable = nullptr;
  qdist = 0.0;
  cut_respa = nullptr;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJLongCoulLong::options(char **arg, int order)
{
  const char *option[] = {"long", "cut", "off", nullptr};
  int i;

  if (!*arg) error->all(FLERR,"Illegal pair_style lj/long/coul/long command");
  for (i=0; option[i]&&strcmp(arg[0], option[i]); ++i);
  switch (i) {
    case 0: ewald_order |= 1<<order; break;
    case 2: ewald_off |= 1<<order; break;
    case 1: break;
    default: error->all(FLERR,"Illegal pair_style lj/long/coul/long command");
  }
}

void PairLJLongCoulLong::settings(int narg, char **arg)
{
  if (narg != 3 && narg != 4) error->all(FLERR,"Illegal pair_style command");

  ewald_order = 0;
  ewald_off = 0;

  options(arg,6);
  options(++arg,1);

  if (!comm->me && ewald_order == ((1<<1) | (1<<6)))
    error->warning(FLERR,"Using largest cutoff for lj/long/coul/long");
  if (!*(++arg))
    error->all(FLERR,"Cutoffs missing in pair_style lj/long/coul/long");
  if (!((ewald_order^ewald_off) & (1<<6)))
    dispersionflag = 0;
  if (!((ewald_order^ewald_off) & (1<<1)))
    error->all(FLERR,"Coulomb cut not supported in pair_style lj/long/coul/long");
  cut_lj_global = utils::numeric(FLERR,*(arg++),false,lmp);
  if (narg == 4 && ((ewald_order & 0x42) == 0x42))
    error->all(FLERR,"Only one cutoff allowed when requesting all long");
  if (narg == 4) cut_coul = utils::numeric(FLERR,*arg,false,lmp);
  else cut_coul = cut_lj_global;

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut_lj[i][j] = cut_lj_global;
  }
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairLJLongCoulLong::~PairLJLongCoulLong()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj_read);
    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(epsilon_read);
    memory->destroy(epsilon);
    memory->destroy(sigma_read);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
  }
  if (ftable) free_tables();
  if (fdisptable) free_disp_tables();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJLongCoulLong::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut_lj_read,n+1,n+1,"pair:cut_lj_read");
  memory->create(cut_lj,n+1,n+1,"pair:cut_lj");
  memory->create(cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  memory->create(epsilon_read,n+1,n+1,"pair:epsilon_read");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma_read,n+1,n+1,"pair:sigma_read");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   extract protected data from object
------------------------------------------------------------------------- */

void *PairLJLongCoulLong::extract(const char *id, int &dim)
{
  const char *ids[] = {
    "B", "sigma", "epsilon", "ewald_order", "ewald_cut", "ewald_mix",
    "cut_coul", "cut_LJ", nullptr};
  void *ptrs[] = {
    lj4, sigma, epsilon, &ewald_order, &cut_coul, &mix_flag,
    &cut_coul, &cut_lj_global, nullptr};
  int i;

  for (i=0; ids[i]&&strcmp(ids[i], id); ++i);
  if (i <= 2) dim = 2;
  else dim = 0;
  return ptrs[i];
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJLongCoulLong::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double epsilon_one = utils::numeric(FLERR,arg[2],false,lmp);
  double sigma_one = utils::numeric(FLERR,arg[3],false,lmp);

  double cut_lj_one = cut_lj_global;
  if (narg == 5) cut_lj_one = utils::numeric(FLERR,arg[4],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon_read[i][j] = epsilon_one;
      sigma_read[i][j] = sigma_one;
      cut_lj_read[i][j] = cut_lj_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJLongCoulLong::init_style()
{
  // require an atom style with charge defined

  if (!atom->q_flag && (ewald_order&(1<<1)))
    error->all(FLERR,
        "Invoking coulombic in pair style lj/long/coul/long requires atom attribute q");

  // ensure use of KSpace long-range solver, set two g_ewalds

  if (force->kspace == nullptr)
    error->all(FLERR,"Pair style requires a KSpace style");
  if (ewald_order&(1<<1)) g_ewald = force->kspace->g_ewald;
  if (ewald_order&(1<<6)) g_ewald_6 = force->kspace->g_ewald_6;

  // set rRESPA cutoffs

  if (utils::strmatch(update->integrate_style,"^respa") &&
      ((Respa *) update->integrate)->level_inner >= 0)
    cut_respa = ((Respa *) update->integrate)->cutoff;
  else cut_respa = nullptr;

  // setup force tables

  if (ncoultablebits && (ewald_order&(1<<1))) init_tables(cut_coul,cut_respa);
  if (ndisptablebits && (ewald_order&(1<<6))) init_tables_disp(cut_lj_global);

  // request regular or rRESPA neighbor lists if neighrequest_flag != 0

  if (force->kspace->neighrequest_flag) {
    int irequest;
    int respa = 0;

    if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
      if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
      if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    }

    irequest = neighbor->request(this,instance_me);

    if (respa >= 1) {
      neighbor->requests[irequest]->respaouter = 1;
      neighbor->requests[irequest]->respainner = 1;
    }
    if (respa == 2) neighbor->requests[irequest]->respamiddle = 1;
  }

  cut_coulsq = cut_coul * cut_coul;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJLongCoulLong::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon_read[i][i],epsilon_read[j][j],
                               sigma_read[i][i],sigma_read[j][j]);
    sigma[i][j] = mix_distance(sigma_read[i][i],sigma_read[j][j]);
    if (ewald_order&(1<<6))
      cut_lj[i][j] = cut_lj_global;
    else
      cut_lj[i][j] = mix_distance(cut_lj_read[i][i],cut_lj_read[j][j]);
  }
  else {
    sigma[i][j] = sigma_read[i][j];
    epsilon[i][j] = epsilon_read[i][j];
    cut_lj[i][j] = cut_lj_read[i][j];
  }

  double cut = MAX(cut_lj[i][j], cut_coul + 2.0*qdist);
  cutsq[i][j] = cut*cut;
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);

  // check interior rRESPA cutoff

  if (cut_respa && MIN(cut_lj[i][j],cut_coul) < cut_respa[3])
    error->all(FLERR,"Pair cutoff < Respa interior cutoff");

  if (offset_flag && (cut_lj[i][j] > 0.0)) {
    double ratio = sigma[i][j] / cut_lj[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  cutsq[j][i] = cutsq[i][j];
  cut_ljsq[j][i] = cut_ljsq[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJLongCoulLong::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon_read[i][j],sizeof(double),1,fp);
        fwrite(&sigma_read[i][j],sizeof(double),1,fp);
        fwrite(&cut_lj_read[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJLongCoulLong::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&epsilon_read[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&sigma_read[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_lj_read[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&epsilon_read[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma_read[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_lj_read[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJLongCoulLong::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&ncoultablebits,sizeof(int),1,fp);
  fwrite(&tabinner,sizeof(double),1,fp);
  fwrite(&ewald_order,sizeof(int),1,fp);
  fwrite(&dispersionflag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJLongCoulLong::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_lj_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_coul,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&ncoultablebits,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&tabinner,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&ewald_order,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&dispersionflag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&ncoultablebits,1,MPI_INT,0,world);
  MPI_Bcast(&tabinner,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&ewald_order,1,MPI_INT,0,world);
  MPI_Bcast(&dispersionflag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJLongCoulLong::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fmt::print(fp,"{} {} {}\n",i,epsilon_read[i][i],sigma_read[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file. must use the "mixed" parameters.
   also must not write out cutoff for lj = long
------------------------------------------------------------------------- */

void PairLJLongCoulLong::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (ewald_order & (1<<6)) {
        fmt::print(fp,"{} {} {} {}\n",i,j,
                   epsilon[i][j],sigma[i][j]);
      } else {
        fmt::print(fp,"{} {} {} {} {}\n",i,j,
                   epsilon[i][j],sigma[i][j],cut_lj[i][j]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute pair interactions
------------------------------------------------------------------------- */

void PairLJLongCoulLong::compute(int eflag, int vflag)
{
  double evdwl,ecoul,fpair;
  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x, *x0 = x[0];
  double **f = atom->f, *f0 = f[0], *fi = f0;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  int i, j, order1 = ewald_order&(1<<1), order6 = ewald_order&(1<<6);
  int *ineigh, *ineighn, *jneigh, *jneighn, typei, typej, ni;
  double qi = 0.0, qri = 0.0;
  double *cutsqi, *cut_ljsqi, *lj1i, *lj2i, *lj3i, *lj4i, *offseti;
  double rsq, r2inv, force_coul, force_lj;
  double g2 = g_ewald_6*g_ewald_6, g6 = g2*g2*g2, g8 = g6*g2;
  double xi[3], d[3];

  ineighn = (ineigh = list->ilist)+list->inum;

  for (; ineigh<ineighn; ++ineigh) {                        // loop over my atoms
    i = *ineigh; fi = f0+3*i;
    if (order1) qri = (qi = q[i])*qqrd2e;                // initialize constants
    offseti = offset[typei = type[i]];
    lj1i = lj1[typei]; lj2i = lj2[typei]; lj3i = lj3[typei]; lj4i = lj4[typei];
    cutsqi = cutsq[typei]; cut_ljsqi = cut_ljsq[typei];
    memcpy(xi, x0+(i+(i<<1)), 3*sizeof(double));
    jneighn = (jneigh = list->firstneigh[i])+list->numneigh[i];

    for (; jneigh<jneighn; ++jneigh) {                       // loop over neighbors
      j = *jneigh;
      ni = sbmask(j);
      j &= NEIGHMASK;

      { double *xj = x0+(j+(j<<1));
        d[0] = xi[0] - xj[0];                               // pair vector
        d[1] = xi[1] - xj[1];
        d[2] = xi[2] - xj[2]; }

      if ((rsq = dot3(d, d)) >= cutsqi[typej = type[j]]) continue;
      r2inv = 1.0/rsq;

      if (order1 && (rsq < cut_coulsq)) {                   // coulombic
        if (!ncoultablebits || rsq <= tabinnersq) {         // series real space
          double r = sqrt(rsq), x = g_ewald*r;
          double s = qri*q[j], t = 1.0/(1.0+EWALD_P*x);
          if (ni == 0) {
            s *= g_ewald*exp(-x*x);
            force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s;
            if (eflag) ecoul = t;
          } else {                                          // special case
            r = s*(1.0-special_coul[ni])/r; s *= g_ewald*exp(-x*x);
            force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-r;
            if (eflag) ecoul = t-r;
          }
        } else {                                            // table real space
          union_int_float_t t;
          t.f = rsq;
          const int k = (t.i & ncoulmask)>>ncoulshiftbits;
          double f = (rsq-rtable[k])*drtable[k], qiqj = qi*q[j];
          if (ni == 0) {
            force_coul = qiqj*(ftable[k]+f*dftable[k]);
            if (eflag) ecoul = qiqj*(etable[k]+f*detable[k]);
          } else {                                          // special case
            t.f = (1.0-special_coul[ni])*(ctable[k]+f*dctable[k]);
            force_coul = qiqj*(ftable[k]+f*dftable[k]-t.f);
            if (eflag) ecoul = qiqj*(etable[k]+f*detable[k]-t.f);
          }
        }
      } else force_coul = ecoul = 0.0;

      if (rsq < cut_ljsqi[typej]) {                         // lj
        if (order6) {                                       // long-range lj
          if (!ndisptablebits || rsq <= tabinnerdispsq) {    // series real space
            double rn = r2inv*r2inv*r2inv;
            double x2 = g2*rsq, a2 = 1.0/x2;
            x2 = a2*exp(-x2)*lj4i[typej];
            if (ni == 0) {
              force_lj =
              (rn*=rn)*lj1i[typej]-g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq;
              if (eflag)
                evdwl = rn*lj3i[typej]-g6*((a2+1.0)*a2+0.5)*x2;
            } else {                                        // special case
              double f = special_lj[ni], t = rn*(1.0-f);
              force_lj = f*(rn *= rn)*lj1i[typej]-
              g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq+t*lj2i[typej];
              if (eflag)
                evdwl = f*rn*lj3i[typej]-g6*((a2+1.0)*a2+0.5)*x2+t*lj4i[typej];
            }
          } else {                                          // table real space
            union_int_float_t disp_t;
            disp_t.f = rsq;
            const int disp_k = (disp_t.i & ndispmask)>>ndispshiftbits;
            double f_disp = (rsq-rdisptable[disp_k])*drdisptable[disp_k];
            double rn = r2inv*r2inv*r2inv;
            if (ni == 0) {
              force_lj = (rn*=rn)*lj1i[typej]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*lj4i[typej];
              if (eflag) evdwl = rn*lj3i[typej]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*lj4i[typej];
            } else {                                        // special case
              double f = special_lj[ni], t = rn*(1.0-f);
              force_lj = f*(rn *= rn)*lj1i[typej]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*lj4i[typej]+t*lj2i[typej];
              if (eflag) evdwl = f*rn*lj3i[typej]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*lj4i[typej]+t*lj4i[typej];
            }
          }
        }
        else {                                                // cut lj
          double rn = r2inv*r2inv*r2inv;
          if (ni == 0) {
            force_lj = rn*(rn*lj1i[typej]-lj2i[typej]);
            if (eflag) evdwl = rn*(rn*lj3i[typej]-lj4i[typej])-offseti[typej];
          }
          else {                                        // special case
            double f = special_lj[ni];
            force_lj = f*rn*(rn*lj1i[typej]-lj2i[typej]);
            if (eflag)
              evdwl = f * (rn*(rn*lj3i[typej]-lj4i[typej])-offseti[typej]);
          }
        }
      }

      else force_lj = evdwl = 0.0;

      fpair = (force_coul+force_lj)*r2inv;

      if (newton_pair || j < nlocal) {
        double *fj = f0+(j+(j<<1)), f;
        fi[0] += f = d[0]*fpair; fj[0] -= f;
        fi[1] += f = d[1]*fpair; fj[1] -= f;
        fi[2] += f = d[2]*fpair; fj[2] -= f;
      }
      else {
        fi[0] += d[0]*fpair;
        fi[1] += d[1]*fpair;
        fi[2] += d[2]*fpair;
      }

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,ecoul,fpair,d[0],d[1],d[2]);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairLJLongCoulLong::compute_inner()
{
  double rsq, r2inv, force_coul = 0.0, force_lj, fpair;

  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *x0 = atom->x[0], *f0 = atom->f[0], *fi = f0, *q = atom->q;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;


  double cut_out_on = cut_respa[0];
  double cut_out_off = cut_respa[1];


  double cut_out_diff = cut_out_off - cut_out_on;
  double cut_out_on_sq = cut_out_on*cut_out_on;
  double cut_out_off_sq = cut_out_off*cut_out_off;

  int *ineigh, *ineighn, *jneigh, *jneighn, typei, typej, ni;
  int i, j, order1 = (ewald_order|(ewald_off^-1))&(1<<1);
  double qri, *cut_ljsqi, *lj1i, *lj2i;
  double xi[3], d[3];

  ineighn = (ineigh = list->ilist_inner)+list->inum_inner;
  for (; ineigh<ineighn; ++ineigh) {                        // loop over my atoms
    i = *ineigh; fi = f0+3*i;
    memcpy(xi, x0+(i+(i<<1)), 3*sizeof(double));
    cut_ljsqi = cut_ljsq[typei = type[i]];
    lj1i = lj1[typei]; lj2i = lj2[typei];
    jneighn = (jneigh = list->firstneigh_inner[i])+list->numneigh_inner[i];
    for (; jneigh<jneighn; ++jneigh) {                        // loop over neighbors
      j = *jneigh;
      ni = sbmask(j);
      j &= NEIGHMASK;

      { double *xj = x0+(j+(j<<1));
        d[0] = xi[0] - xj[0];                                // pair vector
        d[1] = xi[1] - xj[1];
        d[2] = xi[2] - xj[2]; }

      if ((rsq = dot3(d, d)) >= cut_out_off_sq) continue;
      r2inv = 1.0/rsq;

      if (order1 && (rsq < cut_coulsq)) {                       // coulombic
        qri = qqrd2e*q[i];
        force_coul = ni == 0 ?
          qri*q[j]*sqrt(r2inv) : qri*q[j]*sqrt(r2inv)*special_coul[ni];
      }

      if (rsq < cut_ljsqi[typej = type[j]]) {                // lennard-jones
        double rn = r2inv*r2inv*r2inv;
        force_lj = ni == 0 ?
          rn*(rn*lj1i[typej]-lj2i[typej]) :
          rn*(rn*lj1i[typej]-lj2i[typej])*special_lj[ni];
      }
      else force_lj = 0.0;

      fpair = (force_coul + force_lj) * r2inv;

      if (rsq > cut_out_on_sq) {                        // switching
        double rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
        fpair  *= 1.0 + rsw*rsw*(2.0*rsw-3.0);
      }

      if (newton_pair || j < nlocal) {                        // force update
        double *fj = f0+(j+(j<<1)), f;
        fi[0] += f = d[0]*fpair; fj[0] -= f;
        fi[1] += f = d[1]*fpair; fj[1] -= f;
        fi[2] += f = d[2]*fpair; fj[2] -= f;
      }
      else {
        fi[0] += d[0]*fpair;
        fi[1] += d[1]*fpair;
        fi[2] += d[2]*fpair;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLJLongCoulLong::compute_middle()
{
  double rsq, r2inv, force_coul = 0.0, force_lj, fpair;

  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *x0 = atom->x[0], *f0 = atom->f[0], *fi = f0, *q = atom->q;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  double cut_in_off = cut_respa[0];
  double cut_in_on = cut_respa[1];
  double cut_out_on = cut_respa[2];
  double cut_out_off = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_out_diff = cut_out_off - cut_out_on;
  double cut_in_off_sq = cut_in_off*cut_in_off;
  double cut_in_on_sq = cut_in_on*cut_in_on;
  double cut_out_on_sq = cut_out_on*cut_out_on;
  double cut_out_off_sq = cut_out_off*cut_out_off;

  int *ineigh, *ineighn, *jneigh, *jneighn, typei, typej, ni;
  int i, j, order1 = (ewald_order|(ewald_off^-1))&(1<<1);
  double qri, *cut_ljsqi, *lj1i, *lj2i;
  double xi[3], d[3];

  ineighn = (ineigh = list->ilist_middle)+list->inum_middle;

  for (; ineigh<ineighn; ++ineigh) {                        // loop over my atoms
    i = *ineigh; fi = f0+3*i;
    if (order1) qri = qqrd2e*q[i];
    memcpy(xi, x0+(i+(i<<1)), 3*sizeof(double));
    cut_ljsqi = cut_ljsq[typei = type[i]];
    lj1i = lj1[typei]; lj2i = lj2[typei];
    jneighn = (jneigh = list->firstneigh_middle[i])+list->numneigh_middle[i];

    for (; jneigh<jneighn; ++jneigh) {
      j = *jneigh;
      ni = sbmask(j);
      j &= NEIGHMASK;

      { double *xj = x0+(j+(j<<1));
        d[0] = xi[0] - xj[0];                                // pair vector
        d[1] = xi[1] - xj[1];
        d[2] = xi[2] - xj[2]; }

      if ((rsq = dot3(d, d)) >= cut_out_off_sq) continue;
      if (rsq <= cut_in_off_sq) continue;
      r2inv = 1.0/rsq;

      if (order1 && (rsq < cut_coulsq))                        // coulombic
        force_coul = ni == 0 ?
          qri*q[j]*sqrt(r2inv) : qri*q[j]*sqrt(r2inv)*special_coul[ni];

      if (rsq < cut_ljsqi[typej = type[j]]) {                // lennard-jones
        double rn = r2inv*r2inv*r2inv;
        force_lj = ni == 0 ?
          rn*(rn*lj1i[typej]-lj2i[typej]) :
          rn*(rn*lj1i[typej]-lj2i[typej])*special_lj[ni];
      }
      else force_lj = 0.0;

      fpair = (force_coul + force_lj) * r2inv;

      if (rsq < cut_in_on_sq) {                                // switching
        double rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff;
        fpair  *= rsw*rsw*(3.0 - 2.0*rsw);
      }
      if (rsq > cut_out_on_sq) {
        double rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
        fpair  *= 1.0 + rsw*rsw*(2.0*rsw-3.0);
      }

      if (newton_pair || j < nlocal) {                        // force update
        double *fj = f0+(j+(j<<1)), f;
        fi[0] += f = d[0]*fpair; fj[0] -= f;
        fi[1] += f = d[1]*fpair; fj[1] -= f;
        fi[2] += f = d[2]*fpair; fj[2] -= f;
      }
      else {
        fi[0] += d[0]*fpair;
        fi[1] += d[1]*fpair;
        fi[2] += d[2]*fpair;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLJLongCoulLong::compute_outer(int eflag, int vflag)
{
  double evdwl,ecoul,fvirial,fpair;
  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x, *x0 = x[0];
  double **f = atom->f, *f0 = f[0], *fi = f0;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  int i, j, order1 = ewald_order&(1<<1), order6 = ewald_order&(1<<6);
  int *ineigh, *ineighn, *jneigh, *jneighn, typei, typej, ni, respa_flag;
  double qi = 0.0, qri = 0.0;
  double *cutsqi, *cut_ljsqi, *lj1i, *lj2i, *lj3i, *lj4i, *offseti;
  double rsq, r2inv, force_coul, force_lj;
  double g2 = g_ewald_6*g_ewald_6, g6 = g2*g2*g2, g8 = g6*g2;
  double respa_lj = 0.0, respa_coul = 0.0, frespa = 0.0;
  double xi[3], d[3];

  double cut_in_off = cut_respa[2];
  double cut_in_on = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_in_off_sq = cut_in_off*cut_in_off;
  double cut_in_on_sq = cut_in_on*cut_in_on;

  ineighn = (ineigh = list->ilist)+list->inum;

  for (; ineigh<ineighn; ++ineigh) {                        // loop over my atoms
    i = *ineigh; fi = f0+3*i;
    if (order1) qri = (qi = q[i])*qqrd2e;                // initialize constants
    offseti = offset[typei = type[i]];
    lj1i = lj1[typei]; lj2i = lj2[typei]; lj3i = lj3[typei]; lj4i = lj4[typei];
    cutsqi = cutsq[typei]; cut_ljsqi = cut_ljsq[typei];
    memcpy(xi, x0+(i+(i<<1)), 3*sizeof(double));
    jneighn = (jneigh = list->firstneigh[i])+list->numneigh[i];

    for (; jneigh<jneighn; ++jneigh) {                        // loop over neighbors
      j = *jneigh;
      ni = sbmask(j);
      j &= NEIGHMASK;

      { double *xj = x0+(j+(j<<1));
        d[0] = xi[0] - xj[0];                                // pair vector
        d[1] = xi[1] - xj[1];
        d[2] = xi[2] - xj[2]; }

      if ((rsq = dot3(d, d)) >= cutsqi[typej = type[j]]) continue;
      r2inv = 1.0/rsq;

      frespa = 1.0;                                       // check whether and how to compute respa corrections
      respa_coul = 0;
      respa_lj = 0;
      respa_flag = rsq < cut_in_on_sq ? 1 : 0;
      if (respa_flag && (rsq > cut_in_off_sq)) {
        double rsw = (sqrt(rsq)-cut_in_off)/cut_in_diff;
        frespa = 1-rsw*rsw*(3.0-2.0*rsw);
      }

      if (order1 && (rsq < cut_coulsq)) {                // coulombic
        if (!ncoultablebits || rsq <= tabinnersq) {        // series real space
          double r = sqrt(rsq), s = qri*q[j];
          if (respa_flag)                                // correct for respa
            respa_coul = ni == 0 ? frespa*s/r : frespa*s/r*special_coul[ni];
          double x = g_ewald*r, t = 1.0/(1.0+EWALD_P*x);
          if (ni == 0) {
            s *= g_ewald*exp(-x*x);
            force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-respa_coul;
            if (eflag) ecoul = t;
          }
          else {                                        // correct for special
            r = s*(1.0-special_coul[ni])/r; s *= g_ewald*exp(-x*x);
            force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-r-respa_coul;
            if (eflag) ecoul = t-r;
          }
        }                                                // table real space
        else {
          if (respa_flag) {
            double r = sqrt(rsq), s = qri*q[j];
            respa_coul = ni == 0 ? frespa*s/r : frespa*s/r*special_coul[ni];
          }
          union_int_float_t t;
          t.f = rsq;
          const int k = (t.i & ncoulmask) >> ncoulshiftbits;
          double f = (rsq-rtable[k])*drtable[k], qiqj = qi*q[j];
          if (ni == 0) {
            force_coul = qiqj*(ftable[k]+f*dftable[k]);
            if (eflag) ecoul = qiqj*(etable[k]+f*detable[k]);
          }
          else {                                        // correct for special
            t.f = (1.0-special_coul[ni])*(ctable[k]+f*dctable[k]);
            force_coul = qiqj*(ftable[k]+f*dftable[k]-t.f);
            if (eflag) {
              t.f = (1.0-special_coul[ni])*(ptable[k]+f*dptable[k]);
              ecoul = qiqj*(etable[k]+f*detable[k]-t.f);
            }
          }
        }
      }

      else force_coul = respa_coul = ecoul = 0.0;

      if (rsq < cut_ljsqi[typej]) {                        // lennard-jones
        double rn = r2inv*r2inv*r2inv;
        if (respa_flag) respa_lj = ni == 0 ?                 // correct for respa
            frespa*rn*(rn*lj1i[typej]-lj2i[typej]) :
            frespa*rn*(rn*lj1i[typej]-lj2i[typej])*special_lj[ni];
        if (order6) {                                        // long-range form
          if (!ndisptablebits || rsq <= tabinnerdispsq) {
            double x2 = g2*rsq, a2 = 1.0/x2;
            x2 = a2*exp(-x2)*lj4i[typej];
            if (ni == 0) {
              force_lj =
                (rn*=rn)*lj1i[typej]-g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq-respa_lj;
              if (eflag) evdwl = rn*lj3i[typej]-g6*((a2+1.0)*a2+0.5)*x2;
            }
            else {                                        // correct for special
              double f = special_lj[ni], t = rn*(1.0-f);
              force_lj = f*(rn *= rn)*lj1i[typej]-
                g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq+t*lj2i[typej]-respa_lj;
              if (eflag)
                evdwl = f*rn*lj3i[typej]-g6*((a2+1.0)*a2+0.5)*x2+t*lj4i[typej];
            }
          }
          else {                        // table real space
            union_int_float_t disp_t;
            disp_t.f = rsq;
            const int disp_k = (disp_t.i & ndispmask)>>ndispshiftbits;
            double f_disp = (rsq-rdisptable[disp_k])*drdisptable[disp_k];
            double rn = r2inv*r2inv*r2inv;
            if (ni == 0) {
              force_lj = (rn*=rn)*lj1i[typej]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*lj4i[typej]-respa_lj;
              if (eflag) evdwl = rn*lj3i[typej]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*lj4i[typej];
            }
            else {                  // special case
              double f = special_lj[ni], t = rn*(1.0-f);
              force_lj = f*(rn *= rn)*lj1i[typej]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*lj4i[typej]+t*lj2i[typej]-respa_lj;
              if (eflag) evdwl = f*rn*lj3i[typej]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*lj4i[typej]+t*lj4i[typej];
            }
          }
        }
        else {                                                // cut form
          if (ni == 0) {
            force_lj = rn*(rn*lj1i[typej]-lj2i[typej])-respa_lj;
            if (eflag) evdwl = rn*(rn*lj3i[typej]-lj4i[typej])-offseti[typej];
          }
          else {                                        // correct for special
            double f = special_lj[ni];
            force_lj = f*rn*(rn*lj1i[typej]-lj2i[typej])-respa_lj;
            if (eflag)
              evdwl = f*(rn*(rn*lj3i[typej]-lj4i[typej])-offseti[typej]);
          }
        }
      }
      else force_lj = respa_lj = evdwl = 0.0;

      fpair = (force_coul+force_lj)*r2inv;

      if (newton_pair || j < nlocal) {
        double *fj = f0+(j+(j<<1)), f;
        fi[0] += f = d[0]*fpair; fj[0] -= f;
        fi[1] += f = d[1]*fpair; fj[1] -= f;
        fi[2] += f = d[2]*fpair; fj[2] -= f;
      }
      else {
        fi[0] += d[0]*fpair;
        fi[1] += d[1]*fpair;
        fi[2] += d[2]*fpair;
      }

      if (evflag) {
        fvirial = (force_coul + force_lj + respa_coul + respa_lj)*r2inv;
        ev_tally(i,j,nlocal,newton_pair,
                 evdwl,ecoul,fvirial,d[0],d[1],d[2]);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairLJLongCoulLong::single(int i, int j, int itype, int jtype,
                          double rsq, double factor_coul, double factor_lj,
                          double &fforce)
{
  double r2inv, r6inv, force_coul, force_lj;
  double g2 = g_ewald_6*g_ewald_6, g6 = g2*g2*g2, g8 = g6*g2, *q = atom->q;

  double eng = 0.0;

  r2inv = 1.0/rsq;
  if ((ewald_order&2) && (rsq < cut_coulsq)) {                // coulombic
    if (!ncoultablebits || rsq <= tabinnersq) {                // series real space
      double r = sqrt(rsq), x = g_ewald*r;
      double s = force->qqrd2e*q[i]*q[j], t = 1.0/(1.0+EWALD_P*x);
      r = s*(1.0-factor_coul)/r; s *= g_ewald*exp(-x*x);
      force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-r;
      eng += t-r;
    } else {                                                // table real space
      union_int_float_t t;
      t.f = rsq;
      const int k = (t.i & ncoulmask) >> ncoulshiftbits;
      double f = (rsq-rtable[k])*drtable[k], qiqj = q[i]*q[j];
      t.f = (1.0-factor_coul)*(ctable[k]+f*dctable[k]);
      force_coul = qiqj*(ftable[k]+f*dftable[k]-t.f);
      eng += qiqj*(etable[k]+f*detable[k]-t.f);
    }
  } else force_coul = 0.0;

  if (rsq < cut_ljsq[itype][jtype]) {                       // lennard-jones
    r6inv = r2inv*r2inv*r2inv;
    if (ewald_order&64) {                                   // long-range
      double x2 = g2*rsq, a2 = 1.0/x2, t = r6inv*(1.0-factor_lj);
      x2 = a2*exp(-x2)*lj4[itype][jtype];
      force_lj = factor_lj*(r6inv *= r6inv)*lj1[itype][jtype]-
               g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq+t*lj2[itype][jtype];
      eng += factor_lj*r6inv*lj3[itype][jtype]-
        g6*((a2+1.0)*a2+0.5)*x2+t*lj4[itype][jtype];
    } else {                                                // cut
      force_lj = factor_lj*r6inv*(lj1[itype][jtype]*r6inv-lj2[itype][jtype]);
      eng += factor_lj*(r6inv*(r6inv*lj3[itype][jtype]-
                               lj4[itype][jtype])-offset[itype][jtype]);
    }
  } else force_lj = 0.0;

  fforce = (force_coul+force_lj)*r2inv;
  return eng;
}
