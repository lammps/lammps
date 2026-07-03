// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Pieter J. in 't Veld and Stan Moore (Sandia)
------------------------------------------------------------------------- */

#include "pair_lj_long_dipole_long.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kspace.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathExtra;
using namespace EwaldConst;

// ----------------------------------------------------------------------

PairLJLongDipoleLong::PairLJLongDipoleLong(LAMMPS *lmp) :
    Pair(lmp), cut_lj(nullptr), cut_lj_read(nullptr), cut_ljsq(nullptr), epsilon_read(nullptr),
    epsilon(nullptr), sigma_read(nullptr), sigma(nullptr), lj1(nullptr), lj2(nullptr),
    lj3(nullptr), lj4(nullptr), offset(nullptr), cut_respa(nullptr)
{
  dispersionflag = ewaldflag = dipoleflag = 1;
  respa_enable = 0;
  single_enable = 0;
}

// ----------------------------------------------------------------------
// global settings
// ----------------------------------------------------------------------

void PairLJLongDipoleLong::options(char **arg, int mask)
{
  if (!*arg) error->all(FLERR,"Illegal pair_style lj/long/dipole/long command");

  // "long" : treat this interaction with the long-range (k-space) solver
  // "cut"  : plain real-space cutoff (default, leaves both flags clear)
  // "off"  : interaction disabled

  if (strcmp(*arg,"long") == 0) ewald_order |= mask;
  else if (strcmp(*arg,"off") == 0) ewald_off |= mask;
  else if (strcmp(*arg,"cut") != 0)
    error->all(FLERR,"Illegal pair_style lj/long/dipole/long command");
}

void PairLJLongDipoleLong::settings(int narg, char **arg)
{
  if (narg != 3 && narg != 4) error->all(FLERR,"Illegal pair_style command");

  ewald_order = 0;
  ewald_off = 0;

  options(arg, EWALD_DISP);
  options(++arg, EWALD_DIPOLE);
  options(arg, EWALD_COUL);

  if (!comm->me && ewald_order & EWALD_DISP)
    error->warning(FLERR,"Geometric mixing assumed for 1/r^6 coefficients" + utils::errorurl(21));
  if (!comm->me && ewald_order == (EWALD_DIPOLE | EWALD_DISP))
    error->warning(FLERR,"Using largest cut-off for lj/long/dipole/long long long");
  if (!*(++arg))
    error->all(FLERR,"Cutoffs missing in pair_style lj/long/dipole/long");
  if (!((ewald_order^ewald_off) & EWALD_DISP))
    dispersionflag = 0;
  if (!((ewald_order^ewald_off) & EWALD_DIPOLE))
    error->all(FLERR,"Coulombic cut not supported in pair_style lj/long/dipole/long");
  cut_lj_global = utils::numeric(FLERR,*(arg++),false,lmp);
  if (narg == 4 && (ewald_order == (EWALD_COUL | EWALD_DIPOLE | EWALD_DISP)))
    error->all(FLERR,"Only one cut-off allowed when requesting all long");
  if (narg == 4) cut_coul = utils::numeric(FLERR,*(arg++),false,lmp);
  else cut_coul = cut_lj_global;

  if (allocated) {                                      // reset explicit cuts
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut_lj[i][j] = cut_lj_global;
  }
}

// ----------------------------------------------------------------------
// free all arrays
// ----------------------------------------------------------------------

PairLJLongDipoleLong::~PairLJLongDipoleLong()
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
  //if (ftable) free_tables();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJLongDipoleLong::allocate()
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

void *PairLJLongDipoleLong::extract(const char *id, int &dim)
{
  // per-type-pair arrays are 2d, scalars are 0d

  dim = 2;
  if (strcmp(id,"B") == 0) return (void *) lj4;
  if (strcmp(id,"sigma") == 0) return (void *) sigma;
  if (strcmp(id,"epsilon") == 0) return (void *) epsilon;

  dim = 0;
  if (strcmp(id,"ewald_order") == 0) return (void *) &ewald_order;
  if (strcmp(id,"ewald_cut") == 0) return (void *) &cut_coul;
  if (strcmp(id,"ewald_mix") == 0) return (void *) &mix_flag;
  if (strcmp(id,"cut_coul") == 0) return (void *) &cut_coul;
  if (strcmp(id,"cut_vdwl") == 0) return (void *) &cut_lj_global;
  return nullptr;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJLongDipoleLong::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));
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

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJLongDipoleLong::init_style()
{
  if (strcmp(update->unit_style,"electron") == 0)
    error->all(FLERR,"Cannot (yet) use 'electron' units with dipoles");

  // require an atom style with charge defined

  if (!atom->q_flag && (ewald_order & EWALD_COUL))
    error->all(FLERR,
               "Invoking coulombic in pair style lj/long/dipole/long requires atom attribute q");
  if (!atom->mu_flag || !atom->torque_flag)
    error->all(FLERR,"Pair lj/long/dipole/long requires atom attributes mu, torque");

  neighbor->add_request(this);

  cut_coulsq = cut_coul * cut_coul;

  // ensure use of KSpace long-range solver, set g_ewald

  if (ewald_order & EWALD_DIPOLE) {                     // r^-1 kspace
    if (force->kspace == nullptr)
      error->all(FLERR,"Pair style requires a KSpace style");
    if (!force->kspace->dipoleflag)
      error->all(FLERR,"Pair style requires use of kspace_style with dipole support");
  }
  if (ewald_order & EWALD_DISP) {                       // r^-6 kspace
    if (force->kspace == nullptr)
      error->all(FLERR,"Pair style requires a KSpace style");
    if (!force->kspace->dispersionflag)
      error->all(FLERR,"Pair style requires use of kspace_style with dispersion support");
  }
  if (force->kspace) g_ewald = force->kspace->g_ewald;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJLongDipoleLong::init_one(int i, int j)
{
  if ((ewald_order & EWALD_DISP)||(setflag[i][j] == 0)) {
    epsilon[i][j] = mix_energy(epsilon_read[i][i],epsilon_read[j][j],
                               sigma_read[i][i],sigma_read[j][j]);
    sigma[i][j] = mix_distance(sigma_read[i][i],sigma_read[j][j]);
    if (ewald_order & EWALD_DISP)
      cut_lj[i][j] = cut_lj_global;
    else
      cut_lj[i][j] = mix_distance(cut_lj_read[i][i],cut_lj_read[j][j]);
  }
  else {
    sigma[i][j] = sigma_read[i][j];
    epsilon[i][j] = epsilon_read[i][j];
    cut_lj[i][j] = cut_lj_read[i][j];
  }

  double cut = MAX(cut_lj[i][j], cut_coul);
  cutsq[i][j] = cut*cut;
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);

  // check interior rRESPA cutoff

  //if (cut_respa && MIN(cut_lj[i][j],cut_coul) < cut_respa[3])
    //error->all(FLERR,"Pair cutoff < Respa interior cutoff");

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

void PairLJLongDipoleLong::write_restart(FILE *fp)
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

void PairLJLongDipoleLong::read_restart(FILE *fp)
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

void PairLJLongDipoleLong::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&ewald_order,sizeof(int),1,fp);
  fwrite(&ewald_off,sizeof(int),1,fp);
  fwrite(&dispersionflag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJLongDipoleLong::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_lj_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_coul,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&ewald_order,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&ewald_off,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&dispersionflag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&ewald_order,1,MPI_INT,0,world);
  MPI_Bcast(&ewald_off,1,MPI_INT,0,world);
  MPI_Bcast(&dispersionflag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   compute pair interactions
------------------------------------------------------------------------- */

void PairLJLongDipoleLong::compute(int eflag, int vflag)
{
  double evdwl,ecoul,fpair;
  evdwl = ecoul = 0.0;

  ev_init(eflag,vflag);

  double **x = atom->x;
  double **mu = atom->mu;
  double **tq = atom->torque;
  double **f = atom->f;
  double fx, fy, fz;
  double *q = atom->q, qi = 0, qj;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  int i, j, ii, jj, typei, typej, ni;
  int order3 = ewald_order & EWALD_DIPOLE, order6 = ewald_order & EWALD_DISP;
  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  double *cutsqi, *cut_ljsqi, *lj1i, *lj2i, *lj3i, *lj4i, *offseti, *tqi;
  double rsq, r2inv, force_coul, force_lj;
  double g2 = g_ewald*g_ewald, g6 = g2*g2*g2, g8 = g6*g2;
  double B0, B1, B2, B3, G0, G1, G2, mudi, mudj, muij;
  double mui[3], muj[3], xi[3], d[3];

  double C1 = 2.0 * g_ewald / MY_PIS;
  double C2 = 2.0 * g2 * C1;
  double C3 = 2.0 * g2 * C2;

  for (ii = 0; ii < inum; ii++) {                       // loop over all neighs
    i = ilist[ii];
    tqi = tq[i];
    qi = q[i];                          // initialize constants
    typei = type[i];
    offseti = offset[typei];
    lj1i = lj1[typei]; lj2i = lj2[typei]; lj3i = lj3[typei]; lj4i = lj4[typei];
    cutsqi = cutsq[typei]; cut_ljsqi = cut_ljsq[typei];
    xi[0] = x[i][0]; xi[1] = x[i][1]; xi[2] = x[i][2];
    mui[0] = mu[i][0]; mui[1] = mu[i][1]; mui[2] = mu[i][2];

    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {                     // loop over neighbors
      j = jlist[jj];
      ni = sbmask(j);                                   // special index
      j &= NEIGHMASK;

      d[0] = xi[0] - x[j][0];                           // pair vector
      d[1] = xi[1] - x[j][1];
      d[2] = xi[2] - x[j][2];

      rsq = dot3(d, d);
      typej = type[j];
      if (rsq >= cutsqi[typej]) continue;
      r2inv = 1.0/rsq;

      double force_d[3] = {0.0,0.0,0.0}, ti[3] = {0.0,0.0,0.0}, tj[3] = {0.0,0.0,0.0};

      if (order3 && (rsq < cut_coulsq)) {               // dipole
        muj[0] = mu[j][0]; muj[1] = mu[j][1]; muj[2] = mu[j][2];
        {                                               // series real space
          double r = sqrt(rsq);
          double grij = g_ewald*r;
          double f = exp(-grij*grij)*qqrd2e;

          B0 = 1.0/(1.0+EWALD_P*grij);                  // eqn 2.8
          B0 *= ((((A5*B0+A4)*B0+A3)*B0+A2)*B0+A1)*f/r;
          B1 = (B0 + C1 * f) * r2inv;
          B2 = (3.0*B1 + C2 * f) * r2inv;
          B3 = (5.0*B2 + C3 * f) * r2inv;

          mudi = mui[0]*d[0]+mui[1]*d[1]+mui[2]*d[2];
          mudj = muj[0]*d[0]+muj[1]*d[1]+muj[2]*d[2];
          muij = mui[0]*muj[0]+mui[1]*muj[1]+mui[2]*muj[2];
          G0 = qi*(qj = q[j]);                          // eqn 2.10
          G1 = qi*mudj-qj*mudi+muij;
          G2 = -mudi*mudj;
          force_coul = G0*B1+G1*B2+G2*B3;

          mudi *= B2; mudj *= B2;                       // torque contribs
          ti[0] = mudj*d[0]+(qj*d[0]-muj[0])*B1;
          ti[1] = mudj*d[1]+(qj*d[1]-muj[1])*B1;
          ti[2] = mudj*d[2]+(qj*d[2]-muj[2])*B1;

          if (newton_pair || j < nlocal) {
            tj[0] = mudi*d[0]-(qi*d[0]+mui[0])*B1;
            tj[1] = mudi*d[1]-(qi*d[1]+mui[1])*B1;
            tj[2] = mudi*d[2]-(qi*d[2]+mui[2])*B1;
          }

          if (eflag) ecoul = G0*B0+G1*B1+G2*B2;
          if (ni > 0) {                                 // adj part, eqn 2.13
            force_coul -= (f = qqrd2e*(1.0-special_coul[ni])/r)*(
              (3.0*G1+15.0*G2*r2inv)*r2inv+G0)*r2inv;
            if (eflag)
              ecoul -= f*((G1+3.0*G2*r2inv)*r2inv+G0);
            B1 -= f*r2inv;
          }
          B0 = mudj+qj*B1; B3 = -qi*B1+mudi;            // position independent
          if (ni > 0) B0 -= f*3.0*mudj*r2inv*r2inv/B2;
          if (ni > 0) B3 -= f*3.0*mudi*r2inv*r2inv/B2;
          force_d[0] = B0*mui[0]+B3*muj[0];             // force contribs
          force_d[1] = B0*mui[1]+B3*muj[1];
          force_d[2] = B0*mui[2]+B3*muj[2];
          if (ni > 0) {
            ti[0] -= f*(3.0*mudj*r2inv*r2inv*d[0]/B2+(qj*r2inv*d[0]-muj[0]*r2inv));
            ti[1] -= f*(3.0*mudj*r2inv*r2inv*d[1]/B2+(qj*r2inv*d[1]-muj[1]*r2inv));
            ti[2] -= f*(3.0*mudj*r2inv*r2inv*d[2]/B2+(qj*r2inv*d[2]-muj[2]*r2inv));
            if (newton_pair || j < nlocal) {
              tj[0] -= f*(3.0*mudi*r2inv*r2inv*d[0]/B2-(qi*r2inv*d[0]+mui[0]*r2inv));
              tj[1] -= f*(3.0*mudi*r2inv*r2inv*d[1]/B2-(qi*r2inv*d[1]+mui[1]*r2inv));
              tj[2] -= f*(3.0*mudi*r2inv*r2inv*d[2]/B2-(qi*r2inv*d[2]+mui[2]*r2inv));
            }
          }
        }                                               // table real space
      } else {
        force_coul = ecoul = 0.0;
        memset(force_d, 0, 3*sizeof(double));
      }

      if (rsq < cut_ljsqi[typej]) {                     // lj
        if (order6) {                                   // long-range lj
          double r6inv = r2inv*r2inv*r2inv;
          double r12inv = r6inv*r6inv;
          double gr2 = g2*rsq, a2 = 1.0/gr2;
          double expterm = a2*exp(-gr2)*lj4i[typej];      // damped 1/r^6 reciprocal term
          double g6term = g6*((a2+1.0)*a2+0.5)*expterm;
          double g8term = g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*expterm*rsq;
          if (ni == 0) {
            force_lj = r12inv*lj1i[typej]-g8term;
            if (eflag) evdwl = r12inv*lj3i[typej]-g6term;
          } else {                                        // special case
            double factor = special_lj[ni], t = r6inv*(1.0-factor);
            force_lj = factor*r12inv*lj1i[typej]-g8term+t*lj2i[typej];
            if (eflag) evdwl = factor*r12inv*lj3i[typej]-g6term+t*lj4i[typej];
          }
        } else {                                          // cut lj
          double r6inv = r2inv*r2inv*r2inv;
          if (ni == 0) {
            force_lj = r6inv*(r6inv*lj1i[typej]-lj2i[typej]);
            if (eflag) evdwl = r6inv*(r6inv*lj3i[typej]-lj4i[typej])-offseti[typej];
          } else {                                        // special case
            double factor = special_lj[ni];
            force_lj = factor*r6inv*(r6inv*lj1i[typej]-lj2i[typej]);
            if (eflag) evdwl = factor*(
              r6inv*(r6inv*lj3i[typej]-lj4i[typej])-offseti[typej]);
          }
        }
        force_lj *= r2inv;
      } else force_lj = evdwl = 0.0;

      fpair = force_coul+force_lj;                      // force
      fx = d[0]*fpair+force_d[0];
      fy = d[1]*fpair+force_d[1];
      fz = d[2]*fpair+force_d[2];
      f[i][0] += fx;
      f[i][1] += fy;
      f[i][2] += fz;
      tqi[0] += mui[1]*ti[2]-mui[2]*ti[1];              // torque
      tqi[1] += mui[2]*ti[0]-mui[0]*ti[2];
      tqi[2] += mui[0]*ti[1]-mui[1]*ti[0];
      if (newton_pair || j < nlocal) {
        f[j][0] -= fx;
        f[j][1] -= fy;
        f[j][2] -= fz;
        double *tqj = tq[j];
        tqj[0] += muj[1]*tj[2]-muj[2]*tj[1];
        tqj[1] += muj[2]*tj[0]-muj[0]*tj[2];
        tqj[2] += muj[0]*tj[1]-muj[1]*tj[0];
      }

      if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                               evdwl,ecoul,fx,fy,fz,d[0],d[1],d[2]);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

