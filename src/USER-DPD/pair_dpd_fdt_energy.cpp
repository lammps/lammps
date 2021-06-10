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
   Contributing author: James Larentzos (U.S. Army Research Laboratory)
------------------------------------------------------------------------- */

#include "pair_dpd_fdt_energy.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "random_mars.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

PairDPDfdtEnergy::PairDPDfdtEnergy(LAMMPS *lmp) : Pair(lmp)
{
  random = nullptr;
  duCond = nullptr;
  duMech = nullptr;
  splitFDT_flag = false;
  a0_is_zero = false;

  comm_reverse = 2;
}

/* ---------------------------------------------------------------------- */

PairDPDfdtEnergy::~PairDPDfdtEnergy()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(a0);
    memory->destroy(sigma);
    memory->destroy(kappa);
    memory->destroy(alpha);
    memory->destroy(duCond);
    memory->destroy(duMech);
  }

  if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void PairDPDfdtEnergy::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double rsq,r,rinv,wd,wr,factor_dpd,uTmp;
  double dot,randnum;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double dtinvsqrt = 1.0/sqrt(update->dt);

  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double *dpdTheta = atom->dpdTheta;
  double kappa_ij, alpha_ij, theta_ij, gamma_ij;
  double mass_i, mass_j;
  double massinv_i, massinv_j;
  double randPair, mu_ij;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  if (splitFDT_flag) {
    if (!a0_is_zero) for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      itype = type[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        factor_dpd = special_lj[sbmask(j)];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        jtype = type[j];

        if (rsq < cutsq[itype][jtype]) {
          r = sqrt(rsq);
          if (r < EPSILON) continue;     // r can be 0.0 in DPD systems
          rinv = 1.0/r;
          wr = 1.0 - r/cut[itype][jtype];
          wd = wr*wr;

          // conservative force = a0 * wr
          fpair = a0[itype][jtype]*wr;
          fpair *= factor_dpd*rinv;

          f[i][0] += delx*fpair;
          f[i][1] += dely*fpair;
          f[i][2] += delz*fpair;
          if (newton_pair || j < nlocal) {
            f[j][0] -= delx*fpair;
            f[j][1] -= dely*fpair;
            f[j][2] -= delz*fpair;
          }

          if (eflag) {
            // unshifted eng of conservative term:
            // evdwl = -a0[itype][jtype]*r * (1.0-0.5*r/cut[itype][jtype]);
            // eng shifted to 0.0 at cutoff
            evdwl = 0.5*a0[itype][jtype]*cut[itype][jtype] * wd;
            evdwl *= factor_dpd;
          }

          if (evflag) ev_tally(i,j,nlocal,newton_pair,
                               evdwl,0.0,fpair,delx,dely,delz);
        }
      }
    }
  } else {

    // Allocate memory for duCond and duMech
    if (allocated) {
      memory->destroy(duCond);
      memory->destroy(duMech);
    }
    memory->create(duCond,nlocal+nghost,"pair:duCond");
    memory->create(duMech,nlocal+nghost,"pair:duMech");
    for (int ii = 0; ii < nlocal+nghost; ii++) {
      duCond[ii] = 0.0;
      duMech[ii] = 0.0;
    }

    // loop over neighbors of my atoms
    for (int ii = 0; ii < inum; ii++) {
      i = ilist[ii];
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
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        jtype = type[j];

        if (rsq < cutsq[itype][jtype]) {
          r = sqrt(rsq);
          if (r < EPSILON) continue;     // r can be 0.0 in DPD systems
          rinv = 1.0/r;
          wr = 1.0 - r/cut[itype][jtype];
          wd = wr*wr;

          delvx = vxtmp - v[j][0];
          delvy = vytmp - v[j][1];
          delvz = vztmp - v[j][2];
          dot = delx*delvx + dely*delvy + delz*delvz;
          randnum = random->gaussian();

          // Compute the current temperature
          theta_ij = 0.5*(1.0/dpdTheta[i] + 1.0/dpdTheta[j]);
          theta_ij = 1.0/theta_ij;

          gamma_ij = sigma[itype][jtype]*sigma[itype][jtype]
                     / (2.0*force->boltz*theta_ij);

          // conservative force = a0 * wr
          // drag force = -gamma * wr^2 * (delx dot delv) / r
          // random force = sigma * wr * rnd * dtinvsqrt;

          fpair = a0[itype][jtype]*wr;
          fpair -= gamma_ij*wd*dot*rinv;
          fpair += sigma[itype][jtype]*wr*randnum*dtinvsqrt;
          fpair *= factor_dpd*rinv;

          f[i][0] += delx*fpair;
          f[i][1] += dely*fpair;
          f[i][2] += delz*fpair;
          if (newton_pair || j < nlocal) {
            f[j][0] -= delx*fpair;
            f[j][1] -= dely*fpair;
            f[j][2] -= delz*fpair;
          }

          if (rmass) {
            mass_i = rmass[i];
            mass_j = rmass[j];
          } else {
            mass_i = mass[itype];
            mass_j = mass[jtype];
          }
          massinv_i = 1.0 / mass_i;
          massinv_j = 1.0 / mass_j;

          // Compute the mechanical and conductive energy, uMech and uCond
          mu_ij = massinv_i + massinv_j;
          mu_ij *= force->ftm2v;

          uTmp = gamma_ij*wd*rinv*rinv*dot*dot
                 - 0.5*sigma[itype][jtype]*sigma[itype][jtype]*mu_ij*wd;
          uTmp -= sigma[itype][jtype]*wr*rinv*dot*randnum*dtinvsqrt;
          uTmp *= 0.5;

          duMech[i] += uTmp;
          if (newton_pair || j < nlocal) {
            duMech[j] += uTmp;
          }

          // Compute uCond
          randnum = random->gaussian();
          kappa_ij = kappa[itype][jtype];
          alpha_ij = alpha[itype][jtype];
          randPair = alpha_ij*wr*randnum*dtinvsqrt;

          uTmp = kappa_ij*(1.0/dpdTheta[i] - 1.0/dpdTheta[j])*wd;
          uTmp += randPair;

          duCond[i] += uTmp;
          if (newton_pair || j < nlocal) {
            duCond[j] -= uTmp;
          }

          if (eflag) {
            // unshifted eng of conservative term:
            // evdwl = -a0[itype][jtype]*r * (1.0-0.5*r/cut[itype][jtype]);
            // eng shifted to 0.0 at cutoff
            evdwl = 0.5*a0[itype][jtype]*cut[itype][jtype] * wd;
            evdwl *= factor_dpd;
          }

          if (evflag) ev_tally(i,j,nlocal,newton_pair,
                               evdwl,0.0,fpair,delx,dely,delz);
        }
      }
    }
    // Communicate the ghost delta energies to the locally owned atoms
    comm->reverse_comm_pair(this);
  }
  if (vflag_fdotr) virial_fdotr_compute();

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairDPDfdtEnergy::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(a0,n+1,n+1,"pair:a0");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(kappa,n+1,n+1,"pair:kappa");
  memory->create(alpha,n+1,n+1,"pair:alpha");
  if (!splitFDT_flag) {
    memory->create(duCond,nlocal+nghost+1,"pair:duCond");
    memory->create(duMech,nlocal+nghost+1,"pair:duMech");
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairDPDfdtEnergy::settings(int narg, char **arg)
{
  // process keywords
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);
  seed = utils::inumeric(FLERR,arg[1],false,lmp);
  if (atom->dpd_flag != 1)
    error->all(FLERR,"pair_style dpd/fdt/energy requires atom_style with internal temperature and energies (e.g. dpd)");

  // initialize Marsaglia RNG with processor-unique seed

  if (seed <= 0) error->all(FLERR,"Illegal pair_style command");
  delete random;
  random = new RanMars(lmp,seed + comm->me);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairDPDfdtEnergy::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double a0_one = utils::numeric(FLERR,arg[2],false,lmp);
  double sigma_one = utils::numeric(FLERR,arg[3],false,lmp);
  double cut_one = cut_global;
  double kappa_one, alpha_one;

  a0_is_zero = (a0_one == 0.0); // Typical use with SSA is to set a0 to zero

  kappa_one = utils::numeric(FLERR,arg[4],false,lmp);
  alpha_one = sqrt(2.0*force->boltz*kappa_one);
  if (narg == 6) cut_one = utils::numeric(FLERR,arg[5],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a0[i][j] = a0_one;
      sigma[i][j] = sigma_one;
      kappa[i][j] = kappa_one;
      alpha[i][j] = alpha_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairDPDfdtEnergy::init_style()
{
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair dpd/fdt/energy requires ghost atoms store velocity");

  splitFDT_flag = false;
  neighbor->request(this,instance_me);
  for (int i = 0; i < modify->nfix; i++)
    if (utils::strmatch(modify->fix[i]->style,"^shardlow")) {
      splitFDT_flag = true;
    }

  // if newton off, forces between atoms ij will be double computed
  // using different random numbers if splitFDT_flag is false
  if (!splitFDT_flag && (force->newton_pair == 0) && (comm->me == 0)) error->warning(FLERR,
      "Pair dpd/fdt/energy requires newton pair on if not also using fix shardlow");

  bool eos_flag = false;
  for (int i = 0; i < modify->nfix; i++)
    if (utils::strmatch(modify->fix[i]->style,"^eos")) eos_flag = true;
  if (!eos_flag) error->all(FLERR,"pair_style dpd/fdt/energy requires an EOS fix to be specified");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairDPDfdtEnergy::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  cut[j][i] = cut[i][j];
  a0[j][i] = a0[i][j];
  sigma[j][i] = sigma[i][j];
  kappa[j][i] = kappa[i][j];
  alpha[j][i] = alpha[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPDfdtEnergy::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&a0[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&kappa[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPDfdtEnergy::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  a0_is_zero = true; // start with assumption that a0 is zero
  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&a0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&sigma[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&kappa[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&a0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&kappa[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
        alpha[i][j] = sqrt(2.0*force->boltz*kappa[i][j]);
        a0_is_zero = a0_is_zero && (a0[i][j] == 0.0); // verify the zero assumption
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPDfdtEnergy::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPDfdtEnergy::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&seed,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&seed,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);

  // initialize Marsaglia RNG with processor-unique seed
  // same seed that pair_style command initially specified

  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);
}

/* ---------------------------------------------------------------------- */

double PairDPDfdtEnergy::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                       double /*factor_coul*/, double factor_dpd, double &fforce)
{
  double r,rinv,wr,wd,phi;

  r = sqrt(rsq);
  if (r < EPSILON) {
    fforce = 0.0;
    return 0.0;
  }

  rinv = 1.0/r;
  wr = 1.0 - r/cut[itype][jtype];
  wd = wr*wr;
  fforce = a0[itype][jtype]*wr * factor_dpd*rinv;

  phi = 0.5*a0[itype][jtype]*cut[itype][jtype] * wd;
  return factor_dpd*phi;
}

/* ---------------------------------------------------------------------- */

int PairDPDfdtEnergy::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = duCond[i];
    buf[m++] = duMech[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairDPDfdtEnergy::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];

    duCond[j] += buf[m++];
    duMech[j] += buf[m++];
  }
}
