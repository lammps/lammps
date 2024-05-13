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
   Contributing author: Sergey Lishchuk
------------------------------------------------------------------------- */

#include "pair_atm.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;

static const char cite_atm_package[] =
  "ATM package: doi:10.1063/1.4704930\n\n"
  "@Article{Lishchuk:2012:164501,\n"
  " author = {S. V. Lishchuk},\n"
  " title = {Role of Three-Body Interactions in Formation of Bulk Viscosity in Liquid Argon},\n"
  " journal = {J.~Chem.\\ Phys.},\n"
  " year =    2012,\n"
  " volume =  136,\n"
  " number =  16,\n"
  " pages =   {164501}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

PairATM::PairATM(LAMMPS *lmp) : Pair(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_atm_package);

  single_enable = 0;
  restartinfo = 1;
  one_coeff = 0;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairATM::~PairATM()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(nu);
  }
}

/* ----------------------------------------------------------------------
   workhorse routine that computes pairwise interactions
------------------------------------------------------------------------- */

void PairATM::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum,jnumm1;
  double xi,yi,zi,evdwl;
  double rij2,rik2,rjk2;
  double rij[3],rik[3],rjk[3],fj[3],fk[3];
  double nu_local;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;

  double cutoff_squared = cut_global*cut_global;
  double triple = cut_triple*cut_triple*cut_triple;
  double cutoff_triple_sixth = triple*triple;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // triple loop over local atoms and neighbors twice
  // must compute each IJK triplet interaction exactly once
  // by proc that owns the triplet atom with smallest x coord
  //   special logic to break ties if multiple atoms have same x or y coords
  // inner two loops for jj=1,Jnum and kk=jj+1,Jnum ensure
  //   the pair of other 2 non-minimum-x atoms is only considered once
  // triplet geometry criteria for calculation:
  //   each pair distance <= cutoff
  //   produce of 3 pair distances <= cutoff_triple^3

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xi = x[i][0];
    yi = x[i][1];
    zi = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    jnumm1 = jnum - 1;

    for (jj = 0; jj < jnumm1; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      rij[0] = x[j][0] - xi;
      if (rij[0] < 0.0) continue;
      rij[1] = x[j][1] - yi;
      if (rij[0] == 0.0 && rij[1] < 0.0) continue;
      rij[2] = x[j][2] - zi;
      if (rij[0] == 0.0 && rij[1] == 0.0 && rij[2] < 0.0) continue;
      rij2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
      if (rij2 > cutoff_squared) continue;

      for (kk = jj+1; kk < jnum; kk++) {
        k = jlist[kk];
        k &= NEIGHMASK;

        rik[0] = x[k][0] - xi;
        if (rik[0] < 0.0) continue;
        rik[1] = x[k][1] - yi;
        if (rik[0] == 0.0 && rik[1] < 0.0) continue;
        rik[2] = x[k][2] - zi;
        if (rik[0] == 0.0 && rik[1] == 0.0 && rik[2] < 0.0) continue;
        rik2 = rik[0]*rik[0] + rik[1]*rik[1] + rik[2]*rik[2];
        if (rik2 > cutoff_squared) continue;

        rjk[0] = x[k][0] - x[j][0];
        rjk[1] = x[k][1] - x[j][1];
        rjk[2] = x[k][2] - x[j][2];
        rjk2 = rjk[0]*rjk[0] + rjk[1]*rjk[1] + rjk[2]*rjk[2];
        if (rjk2 > cutoff_squared) continue;

        double r6 = rij2*rjk2*rik2;
        if (r6 > cutoff_triple_sixth) continue;

        nu_local = nu[type[i]][type[j]][type[k]];
        if (nu_local == 0.0) continue;

        interaction_ddd(nu_local,
                        r6,rij2,rik2,rjk2,rij,rik,rjk,fj,fk,eflag,evdwl);

        f[i][0] -= fj[0] + fk[0];
        f[i][1] -= fj[1] + fk[1];
        f[i][2] -= fj[2] + fk[2];
        f[j][0] += fj[0];
        f[j][1] += fj[1];
        f[j][2] += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];

        if (evflag) ev_tally3(i,j,k,evdwl,0.0,fj,fk,rij,rik);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairATM::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(nu,n+1,n+1,n+1,"pair:nu");

  // initialize all nu values to 0.0

  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      for (int k = 1; k <= n; k++)
        nu[i][j][k] = 0.0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairATM::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);
  cut_triple = utils::numeric(FLERR,arg[1],false,lmp);
}

/* ----------------------------------------------------------------------
   set coefficients for one I,J,K type triplet
------------------------------------------------------------------------- */

void PairATM::coeff(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi,klo,khi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);
  utils::bounds(FLERR,arg[2],1,atom->ntypes,klo,khi,error);

  double nu_one = utils::numeric(FLERR,arg[3],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j<=jhi; j++) {
      for (int k = MAX(klo,j); k<=khi; k++) {
        nu[i][j][k] = nu_one;
        count++;
      }
      setflag[i][j] = 1;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairATM::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style ATM requires newton pair on");

  // need a full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one i,j type pair and corresponding j,i
   also for all k type permutations
------------------------------------------------------------------------- */

double PairATM::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  // set all 6 symmetric permutations of I,J,K types to same nu value

  int ntypes = atom->ntypes;
  for (int k = j; k <= ntypes; k++)
    nu[i][k][j] = nu[j][i][k] = nu[j][k][i] = nu[k][i][j] = nu[k][j][i] =
      nu[i][j][k];

  return cut_global;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairATM::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j,k;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j])
        for (k = j; k <= atom->ntypes; k++)
          fwrite(&nu[i][j][k],sizeof(double),1,fp);
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairATM::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j,k;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) for (k = j; k <= atom->ntypes; k++) {
          if (me == 0) utils::sfread(FLERR,&nu[i][j][k],sizeof(double),1,fp,nullptr,error);
        MPI_Bcast(&nu[i][j][k],1,MPI_DOUBLE,0,world);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairATM::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&cut_triple,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairATM::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_triple,sizeof(double),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_triple,1,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   Axilrod-Teller-Muto (dipole-dipole-dipole) potential
------------------------------------------------------------------------- */

void PairATM::interaction_ddd(double nu, double r6,
                              double rij2, double rik2, double rjk2,
                              double *rij, double *rik, double *rjk,
                              double *fj, double *fk, int eflag, double &eng)
{
  double r5inv,rri,rrj,rrk,rrr;
  r5inv = nu / (r6*r6*sqrt(r6));
  rri = rik[0]*rij[0] + rik[1]*rij[1] + rik[2]*rij[2];
  rrj = rij[0]*rjk[0] + rij[1]*rjk[1] + rij[2]*rjk[2];
  rrk = rjk[0]*rik[0] + rjk[1]*rik[1] + rjk[2]*rik[2];
  rrr = 5.0*rri*rrj*rrk;
  for (int i = 0; i < 3; i++) {
    fj[i] = rrj*(rrk - rri)*rik[i] -
      (rrk*rri - rjk2*rik2 + rrr/rij2) * rij[i] +
      (rrk*rri - rik2*rij2 + rrr/rjk2) * rjk[i];
    fj[i] *= 3.0*r5inv;
    fk[i] = rrk*(rri + rrj)*rij[i] +
      (rri*rrj + rik2*rij2 - rrr/rjk2) * rjk[i] +
      (rri*rrj + rij2*rjk2 - rrr/rik2) * rik[i];
    fk[i] *= 3.0*r5inv;
  }
  if (eflag) eng = (r6 - 0.6*rrr)*r5inv;
}
