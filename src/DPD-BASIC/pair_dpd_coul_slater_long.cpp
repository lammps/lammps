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
   Contributing author: Eddy BARRAUD (IFPEN/Sorbonne)
------------------------------------------------------------------------- */

#include "pair_dpd_coul_slater_long.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "random_mars.h"
#include "update.h"

#include "ewald_const.h"
#include "kspace.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace EwaldConst;

static constexpr double EPSILON = 1.0e-10;

/* ---------------------------------------------------------------------- */

PairDPDCoulSlaterLong::PairDPDCoulSlaterLong(LAMMPS *lmp) :
    Pair(lmp), cut_dpd(nullptr), cut_dpdsq(nullptr), cut_slatersq(nullptr),
    a0(nullptr), gamma(nullptr), sigma(nullptr), random(nullptr)
{
  writedata = 1;
  ewaldflag = pppmflag = 1;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairDPDCoulSlaterLong::~PairDPDCoulSlaterLong()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut_dpd);
    memory->destroy(cut_dpdsq);
    memory->destroy(cut_slatersq);

    memory->destroy(cut);
    memory->destroy(a0);
    memory->destroy(gamma);
    memory->destroy(sigma);
  }

  if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void PairDPDCoulSlaterLong::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double r2inv,forcedpd,forcecoul,factor_coul;
  double grij,expm2,prefactor,t,erfc;
  double rsq,r,rinv,dot,wd,randnum,factor_dpd,factor_sqrt;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double slater_term;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double dtinvsqrt = 1.0/sqrt(update->dt);

  double *q = atom->q;
  double *special_coul = force->special_coul;
  double qqrd2e = force->qqrd2e;

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
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_dpd = special_lj[sbmask(j)];
      factor_sqrt = special_sqrt[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      // forces if below maximum cutoff
      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        if (evflag) evdwl = ecoul = 0.0;

        // apply DPD force if distance below DPD cutoff
        if (rsq < cut_dpdsq[itype][jtype] && r > EPSILON) {
          rinv = 1.0/r;
          delvx = vxtmp - v[j][0];
          delvy = vytmp - v[j][1];
          delvz = vztmp - v[j][2];
          dot = delx*delvx + dely*delvy + delz*delvz;
          wd = 1.0 - r/cut_dpd[itype][jtype];
          randnum = random->gaussian();

          // conservative force = a0 * wd
          // drag force = -gamma * wd^2 * (delx dot delv) / r
          // random force = sigma * wd * rnd * dtinvsqrt;
          // random force must be scaled by sqrt(factor_dpd)

          forcedpd = a0[itype][jtype]*wd;
          forcedpd -= gamma[itype][jtype]*wd*wd*dot*rinv;
          forcedpd *= factor_dpd;
          forcedpd += factor_sqrt*sigma[itype][jtype]*wd*randnum*dtinvsqrt;
          forcedpd *= rinv;

          if (eflag) {
            // eng shifted to 0.0 at cutoff
            evdwl = 0.5*a0[itype][jtype]*cut_dpd[itype][jtype] * wd*wd;
            evdwl *= factor_dpd;
          }

        } else forcedpd = 0.0;

        // apply Slater electrostatic force if distance below Slater cutoff
        // and the two species are charged
        if (rsq < cut_slatersq[itype][jtype]){
          r2inv = 1.0/rsq;
          grij = g_ewald * r;
          expm2 = exp(-grij*grij);
          t = 1.0 / (1.0 + EWALD_P*grij);
          erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
          slater_term = exp(-2*r/lamda)*(1 + (2*r/lamda*(1+r/lamda)));
          prefactor = qqrd2e * qtmp*q[j]/r;
          forcecoul = prefactor * (erfc + EWALD_F*grij*expm2 - slater_term);
          if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor*(1-slater_term);
          forcecoul *= r2inv;

          if (eflag) {
            ecoul = prefactor*(erfc - (1 + r/lamda)*exp(-2*r/lamda));
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor*(1.0-(1 + r/lamda)*exp(-2*r/lamda));
          }

        } else forcecoul = 0.0;

        fpair = forcedpd + forcecoul;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }


        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairDPDCoulSlaterLong::allocate()
{
  int i,j;
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(cut_dpd,n+1,n+1,"pair:cut_dpd");
  memory->create(cut_dpdsq,n+1,n+1,"pair:cut_dpdsq");
  memory->create(cut_slatersq,n+1,n+1,"pair:cut_slatersq");
  memory->create(a0,n+1,n+1,"pair:a0");
  memory->create(gamma,n+1,n+1,"pair:gamma");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  for (i = 0; i <= atom->ntypes; i++)
    for (j = 0; j <= atom->ntypes; j++)
      sigma[i][j] = gamma[i][j] = 0.0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairDPDCoulSlaterLong::settings(int narg, char **arg)
{
  // params : T cut_dpd seed lambda cut_coul
  if (narg != 5) error->all(FLERR,"Illegal pair_style command");

  temperature = utils::numeric(FLERR,arg[0],false,lmp);
  cut_global = utils::numeric(FLERR,arg[1],false,lmp);
  seed = utils::inumeric(FLERR,arg[2],false,lmp);
  lamda = utils::numeric(FLERR,arg[3],false,lmp);
  cut_coul = utils::numeric(FLERR,arg[4],false,lmp);

  // initialize Marsaglia RNG with processor-unique seed

  if (seed <= 0)
    error->all(FLERR,"Invalid random seed {} for pair_style dpd/coul/slater/long command", seed);
  delete random;
  random = new RanMars(lmp,seed + comm->me);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut_dpd[i][j] = MAX(cut_global,cut_coul);
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairDPDCoulSlaterLong::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 6)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double a0_one = utils::numeric(FLERR,arg[2],false,lmp);
  double gamma_one = utils::numeric(FLERR,arg[3],false,lmp);

  double cut_one = cut_global;
  double cut_two = 0.0;

  if (narg > 4) {
    bool do_slater = utils::logical(FLERR,arg[4],false,lmp);
    if (do_slater) cut_two = cut_coul;
  }

  if (narg > 5) cut_one = utils::numeric(FLERR,arg[5],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a0[i][j] = a0_one;
      gamma[i][j] = gamma_one;
      cut_dpd[i][j] = cut_one;
      cut_slatersq[i][j] = cut_two * cut_two;
      cut[i][j] = MAX(cut_one, cut_two);
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairDPDCoulSlaterLong::init_style()
{
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair dpd requires ghost atoms store velocity");
  if (!atom->q_flag)
    error->all(FLERR,"Pair style coul/slater/long requires atom attribute q");

  // if newton off, forces between atoms ij will be double computed
  // using different random numbers

  if (force->newton_pair == 0 && comm->me == 0)
    error->warning(FLERR, "Pair dpd needs newton pair on for momentum conservation");

  neighbor->add_request(this);

  // precompute random force scaling factors

  for (int i = 0; i < 4; ++i) special_sqrt[i] = sqrt(force->special_lj[i]);


  // ensure use of KSpace long-range solver, set g_ewald

 if (force->kspace == nullptr)
    error->all(FLERR,"Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
   return the maximum cutoff between Slater or DPD cutoff if charged
   return the DPD cutoff for uncharged
------------------------------------------------------------------------- */

double PairDPDCoulSlaterLong::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  sigma[i][j] = sqrt(2.0*force->boltz*temperature*gamma[i][j]);

  cut_dpdsq[i][j] = cut_dpd[i][j] * cut_dpd[i][j];

  a0[j][i] = a0[i][j];
  gamma[j][i] = gamma[i][j];
  sigma[j][i] = sigma[i][j];
  cut_dpd[j][i] = cut_dpd[i][j];
  cut[j][i] = cut[i][j];
  cut_dpdsq[j][i] = cut_dpdsq[i][j];
  cut_slatersq[j][i] = cut_slatersq[i][j];

  return MAX(cut_dpd[i][j], sqrt(cut_slatersq[i][j]));
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPDCoulSlaterLong::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&a0[i][j],sizeof(double),1,fp);
        fwrite(&gamma[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
        fwrite(&cut_dpd[i][j],sizeof(double),1,fp);
        fwrite(&cut_slatersq[i][j],sizeof(double),1,fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPDCoulSlaterLong::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&a0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&gamma[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_dpd[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_slatersq[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&a0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gamma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_dpd[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_slatersq[i][j],1,MPI_DOUBLE,0,world);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPDCoulSlaterLong::write_restart_settings(FILE *fp)
{
  fwrite(&temperature,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&lamda,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPDCoulSlaterLong::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&temperature,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&seed,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&lamda,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_coul,sizeof(double),1,fp,nullptr,error);
  }
  MPI_Bcast(&temperature,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&seed,1,MPI_INT,0,world);
  MPI_Bcast(&lamda,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);

  // initialize Marsaglia RNG with processor-unique seed
  // same seed that pair_style command initially specified

  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairDPDCoulSlaterLong::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,a0[i][i],gamma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairDPDCoulSlaterLong::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %s %g\n",i,j,a0[i][j],gamma[i][j],
              (cut_slatersq[i][j] == 0.0) ? "yes" : "no", cut_dpd[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairDPDCoulSlaterLong::single(int i, int j, int itype, int jtype, double rsq,
                                     double factor_coul, double factor_dpd, double &fforce)
{
  double r,rinv,wd,phi;
  double r2inv,grij,expm2,t,erfc,prefactor;
  double slater_term;
  double forcecoul,phicoul;

  double energy = 0.0;
  fforce = 0.0;

  r = sqrt(rsq);

  // compute DPD force and energy
  if (rsq < cut_dpdsq[itype][jtype] && r > EPSILON) {
    rinv = 1.0/r;
    wd = 1.0 - r/cut_dpd[itype][jtype];
    fforce += a0[itype][jtype]*wd * factor_dpd*rinv;

    phi = 0.5*a0[itype][jtype]*cut_dpd[itype][jtype] * wd*wd;
    energy += factor_dpd*phi;
  }

  // compute Slater coulombic force and energy
  if (atom->q[i]*atom->q[j] != 0.0 && rsq < cut_slatersq[itype][jtype]) {
    r2inv = 1.0/rsq;
    grij = g_ewald * r;
    expm2 = exp(-grij*grij);
    t = 1.0 / (1.0 + EWALD_P*grij);
    erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
    slater_term = exp(-2*r/lamda)*(1 + (2*r/lamda*(1+r/lamda)));
    prefactor = force->qqrd2e * atom->q[i]*atom->q[j]/r;
    forcecoul = prefactor * (erfc + EWALD_F*grij*expm2 - slater_term);
    if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
    fforce += forcecoul * r2inv;
    phicoul = prefactor*(erfc - (1 + r/lamda)*exp(-2*r/lamda));
    if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor;
    energy += phicoul;
  }

  return energy;
}

void *PairDPDCoulSlaterLong::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul;
  if (strcmp(str,"lamda") == 0) return (void *) &lamda;
  dim = 2;
  if (strcmp(str,"a0") == 0) return (void *) a0;
  if (strcmp(str,"gamma") == 0) return (void *) gamma;
  return nullptr;
}
