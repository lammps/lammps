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
   Contributing author: Zhen Li (Clemson University)
   Email: zli7@clemson.edu
------------------------------------------------------------------------- */

#include "pair_tdpd.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "random_mars.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

static constexpr double EPSILON = 1.0e-10;

static const char cite_pair_tdpd[] =
  "pair tdpd command: doi:10.1063/1.4923254\n\n"
  "@Article{ZLi2015_JCP,\n"
  " author = {Li, Z. and Yazdani, A. and Tartakovsky, A. and Karniadakis, G. E.},\n"
  " title = {Transport Dissipative Particle Dynamics Model for Mesoscopic Advection-Diffusion-Reaction Problems},\n"
  " journal = {The Journal of Chemical Physics},\n"
  " year = {2015},\n"
  " volume = {143},\n"
  " number = {1},\n"
  " pages = {014101}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

PairTDPD::PairTDPD(LAMMPS *lmp) : Pair(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_pair_tdpd);
  cc_species = atom->cc_species;

  writedata = 1;
  random = nullptr;
}

/* ---------------------------------------------------------------------- */

PairTDPD::~PairTDPD()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(cutcc);

    memory->destroy(a0);
    memory->destroy(gamma);
    memory->destroy(sigma);

    memory->destroy(power);
    memory->destroy(kappa);
    memory->destroy(epsilon);
    memory->destroy(powercc);
  }

  if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void PairTDPD::compute(int eflag, int vflag)
{
  double evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **cc = atom->cc;
  double **cc_flux = atom->cc_flux;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double dtinvsqrt = 1.0/sqrt(update->dt);

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    double xtmp = x[i][0];
    double ytmp = x[i][1];
    double ztmp = x[i][2];
    double vxtmp = v[i][0];
    double vytmp = v[i][1];
    double vztmp = v[i][2];
    int itype = type[i];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      double factor_dpd = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      double delx = xtmp - x[j][0];
      double dely = ytmp - x[j][1];
      double delz = ztmp - x[j][2];
      double rsq = delx*delx + dely*dely + delz*delz;
      int jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        double r = sqrt(rsq);
        if (r < EPSILON) continue;     // r can be 0.0 in DPD systems
        double rinv = 1.0/r;
        double delvx = vxtmp - v[j][0];
        double delvy = vytmp - v[j][1];
        double delvz = vztmp - v[j][2];
        double dot = delx*delvx + dely*delvy + delz*delvz;
        double wc = 1.0 - r/cut[itype][jtype];
        wc = MAX(0,MIN(1.0,wc));
        double wr = pow(wc, 0.5*power[itype][jtype]);
        double randnum = random->gaussian();

        // conservative force = a0 * wc
        // drag force = -gamma * wr^2 * (delx dot delv) / r
        // random force = sigma * wr^(power/2) * rnd * dtinvsqrt;

        double fpair = a0[itype][jtype]*wc;
        fpair -= gamma[itype][jtype]*wr*wr*dot*rinv;
        fpair += sigma[itype][jtype]*wr*randnum*dtinvsqrt;
        fpair *= factor_dpd*rinv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;

        // chemical concentration transport
        if (r < cutcc[itype][jtype]) {
          for (int k=0; k<cc_species; k++) {
            double wcr = 1.0 - r/cutcc[itype][jtype];
            wcr = MAX(0,wcr);
            wcr = pow(wcr, 0.5*powercc[itype][jtype][k]);
            double randnum = random->gaussian();
            randnum = MAX(-5.0,MIN(randnum,5.0));
            double dQc = -kappa[itype][jtype][k] * wcr*wcr *(cc[i][k]-cc[j][k]);
            double dQr = epsilon[itype][jtype][k] *wcr* (cc[i][k]+cc[j][k]) *randnum*dtinvsqrt;
            cc_flux[i][k] += (dQc + dQr);
            if (newton_pair || j < nlocal)
              cc_flux[j][k] -= ( dQc + dQr );
          }
        }
        //-----------------------------------------------------------

        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = 0.5*a0[itype][jtype]*cut[itype][jtype] * wc*wc;
          evdwl *= factor_dpd;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairTDPD::allocate()
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
  memory->create(cutcc,n+1,n+1,"pair:cutcc");
  memory->create(a0,n+1,n+1,"pair:a0");
  memory->create(gamma,n+1,n+1,"pair:gamma");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(power,n+1,n+1,"pair:power");
  memory->create(kappa,n+1,n+1,cc_species,"pair:kappa");
  memory->create(epsilon,n+1,n+1,cc_species,"pair:epsilon");
  memory->create(powercc,n+1,n+1,cc_species,"pair:powercc");

  for (i = 0; i <= atom->ntypes; i++)
  for (j = 0; j <= atom->ntypes; j++)
    sigma[i][j] = gamma[i][j] = 0.0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairTDPD::settings(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal pair_style command");

  temperature = utils::numeric(FLERR,arg[0],false,lmp);
  cut_global = utils::numeric(FLERR,arg[1],false,lmp);
  seed = utils::inumeric(FLERR,arg[2],false,lmp);

  // initialize Marsaglia RNG with processor-unique seed
  // create a positive seed based on the system clock, if requested.

  if (seed <= 0) {
    constexpr double LARGE_NUM = 2<<30;
    seed = int(fmod(platform::walltime() * LARGE_NUM, LARGE_NUM)) + 1;
  }

  delete random;
  random = new RanMars(lmp,(seed + comm->me) % 900000000);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
    for (j = i+1; j <= atom->ntypes; j++)
    if (setflag[i][j])
      cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairTDPD::coeff(int narg, char **arg)
{
  if (narg != 7 + 3*cc_species)
    error->all(FLERR,"Incorrect args for pair tdpd coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double a0_one = utils::numeric(FLERR,arg[2],false,lmp);
  double gamma_one = utils::numeric(FLERR,arg[3],false,lmp);
  double power_one = utils::numeric(FLERR,arg[4],false,lmp);
  double cut_one   = utils::numeric(FLERR,arg[5],false,lmp);
  double cutcc_one = utils::numeric(FLERR,arg[6],false,lmp);
  auto kappa_one = new double[cc_species];
  auto epsilon_one = new double[cc_species];
  auto powercc_one = new double[cc_species];
  for (int k=0; k<cc_species; k++) {
    kappa_one[k]   = utils::numeric(FLERR,arg[7+3*k],false,lmp);
    epsilon_one[k] = utils::numeric(FLERR,arg[8+3*k],false,lmp);
    powercc_one[k] = utils::numeric(FLERR,arg[9+3*k],false,lmp);
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++)
  for (int j = MAX(jlo,i); j <= jhi; j++) {
    a0[i][j] = a0_one;
    gamma[i][j] = gamma_one;
    power[i][j] = power_one;
    cut[i][j]   = cut_one;
    cutcc[i][j] = cutcc_one;
    for (int k=0; k<cc_species; k++)
    {
      kappa [i][j][k] = kappa_one[k];
      epsilon[i][j][k]= epsilon_one[k];
      powercc[i][j][k]= powercc_one[k];
    }
    setflag[i][j] = 1;
    count++;
  }
  delete[] kappa_one;
  delete[] epsilon_one;
  delete[] powercc_one;

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairTDPD::init_style()
{
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair style tdpd requires ghost atoms store velocity");

  if (!atom->tdpd_flag)
    error->all(FLERR,"Pair style tdpd requires atom properties cc/cc_flux");

  // if newton off, forces between atoms ij will be double computed
  // using different random numbers

  if (force->newton_pair == 0 && comm->me == 0)
    error->warning(FLERR,"Pair tdpd needs newton pair on for momentum conservation");

  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairTDPD::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  sigma[i][j] = sqrt(2.0*force->boltz*temperature*gamma[i][j]);

  cut[j][i]   = cut[i][j];
  cutcc[j][i] = cutcc[i][j];
  a0[j][i]    = a0[i][j];
  gamma[j][i] = gamma[i][j];
  sigma[j][i] = sigma[i][j];
  power[j][i] = power[i][j];
  for (int k=0; k<cc_species; k++) {
    kappa[j][i][k] = kappa[i][j][k];
    epsilon[j][i][k] = epsilon[i][j][k];
    powercc[j][i][k] = powercc[i][j][k];
  }
  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTDPD::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  for (int i = 1; i <= atom->ntypes; i++)
  for (int j = i; j <= atom->ntypes; j++) {
    fwrite(&setflag[i][j],sizeof(int),1,fp);
    if (setflag[i][j]) {
      fwrite(&a0[i][j],sizeof(double),1,fp);
      fwrite(&gamma[i][j],sizeof(double),1,fp);
      fwrite(&power[i][j],sizeof(double),1,fp);
      fwrite(&cut[i][j],sizeof(double),1,fp);
      fwrite(&cutcc[i][j],sizeof(double),1,fp);
      for (int k=0; k<cc_species; k++) {
        fwrite(&kappa[i][j][k],sizeof(double),1,fp);
        fwrite(&epsilon[i][j][k],sizeof(double),1,fp);
        fwrite(&powercc[i][j][k],sizeof(double),1,fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTDPD::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int me = comm->me;
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&a0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&gamma[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&power[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cutcc[i][j],sizeof(double),1,fp,nullptr,error);
          for (int k=0; k<cc_species; k++) {
            utils::sfread(FLERR,&kappa[i][j][k],sizeof(double),1,fp,nullptr,error);
            utils::sfread(FLERR,&epsilon[i][j][k],sizeof(double),1,fp,nullptr,error);
            utils::sfread(FLERR,&powercc[i][j][k],sizeof(double),1,fp,nullptr,error);
          }
        }
        MPI_Bcast(&a0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gamma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&power[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cutcc[i][j],1,MPI_DOUBLE,0,world);
        for (int k=0; k<cc_species; k++) {
          MPI_Bcast(&kappa[i][j][k],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&epsilon[i][j][k],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&powercc[i][j][k],1,MPI_DOUBLE,0,world);
        }
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTDPD::write_restart_settings(FILE *fp)
{
  fwrite(&temperature,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTDPD::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&temperature,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&seed,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&temperature,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&seed,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);

  // initialize Marsaglia RNG with processor-unique seed
  // same seed that pair_style command initially specified

  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);
}

/* ---------------------------------------------------------------------- */

double PairTDPD::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                       double /*factor_coul*/, double factor_dpd, double &fforce)
{
  double r,rinv,wc,phi;

  r = sqrt(rsq);
  if (r < EPSILON) {
    fforce = 0.0;
    return 0.0;
  }

  rinv = 1.0/r;
  wc = 1.0 - r/cut[itype][jtype];
  fforce = a0[itype][jtype]*wc*factor_dpd*rinv;

  phi = 0.5*a0[itype][jtype]*cut[itype][jtype]*wc*wc;
  return factor_dpd*phi;
}
