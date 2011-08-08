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
   CMM coarse grained MD potentials. Coulomb with MSM version.
   Contributing author: Mike Brown <brownw@ornl.gov>
------------------------------------------------------------------------- */

#include "pair_cg_cmm_coul_msm.h"
#include "memory.h"
#include "atom.h"
#include "force.h"
#include "kspace.h"

#include "string.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
enum {C3=0,C4=1};

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairCGCMMCoulMSM::PairCGCMMCoulMSM(LAMMPS *lmp) : PairCMMCommon(lmp)
{
  respa_enable = 0;
  single_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairCGCMMCoulMSM::~PairCGCMMCoulMSM()
{
  if (allocated_coul) {
    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(cut_coul);
    memory->destroy(cut_coulsq);
    allocated_coul=0;
  }
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulMSM::allocate()
{
  PairCMMCommon::allocate();
  allocated_coul = 1;

  int n = atom->ntypes;

  memory->create(cut_lj,n+1,n+1,"paircg:cut_lj");
  memory->create(cut_ljsq,n+1,n+1,"paircg:cut_ljsq");
  memory->create(cut_coul,n+1,n+1,"paircg:cut_coul");
  memory->create(cut_coulsq,n+1,n+1,"paircg:cut_coulsq");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCGCMMCoulMSM::settings(int narg, char **arg)
{
  // strip off smoothing type and send args to parent

  if (narg < 1) error->all("Illegal pair_style command");

  if (strcmp(arg[0],"C3") == 0)
    _smooth = C3;
  else if (strcmp(arg[0],"C4") == 0)
    _smooth = C4;
  else error->all("Illegal pair_style command");

  PairCMMCommon::settings(narg-1,&arg[1]);
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulMSM::init_style()
{
  if (!atom->q_flag)
    error->all("Pair style cg/cut/coul/msm requires atom attribute q");
  
  PairCMMCommon::init_style();
  _ia=-1.0/cut_coul_global;
  _ia2=_ia*_ia;
  _ia3=_ia2*_ia;

  cut_respa = NULL;
}

/* ---------------------------------------------------------------------- */

double PairCGCMMCoulMSM::init_one(int i, int j)
{
  double mycut = PairCMMCommon::init_one(i,j);

  return mycut;
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulMSM::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,itable;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz;
  double fraction,table;
  double r,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  double qqrd2e = force->qqrd2e;
  int newton_pair = force->newton_pair;

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
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq = delx*delx + dely*dely + delz*delz;
      const int jtype = type[j];

      double evdwl = 0.0;
      double ecoul = 0.0;
      double fpair = 0.0;

      if (rsq < cutsq[itype][jtype]) {
        const double r2inv = 1.0/rsq;
        const int cgt=cg_type[itype][jtype];

        double forcelj  = 0.0;
        double forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
          forcelj=factor_lj;
          if (eflag) evdwl=factor_lj;

          if (cgt == CG_LJ12_4) {
            const double r4inv=r2inv*r2inv;
            forcelj *= r4inv*(lj1[itype][jtype]*r4inv*r4inv
                       - lj2[itype][jtype]);
            if (eflag) {
              evdwl *= r4inv*(lj3[itype][jtype]*r4inv*r4inv
                       - lj4[itype][jtype]) - offset[itype][jtype];
            }
          } else if (cgt == CG_LJ9_6) {
            const double r3inv = r2inv*sqrt(r2inv);
            const double r6inv = r3inv*r3inv;
            forcelj *= r6inv*(lj1[itype][jtype]*r3inv
                       - lj2[itype][jtype]);
            if (eflag) {
              evdwl *= r6inv*(lj3[itype][jtype]*r3inv
                        - lj4[itype][jtype]) - offset[itype][jtype];
            }
          } else {
            const double r6inv = r2inv*r2inv*r2inv;
            forcelj *= r6inv*(lj1[itype][jtype]*r6inv
                       - lj2[itype][jtype]);
            if (eflag) {
              evdwl *= r6inv*(lj3[itype][jtype]*r6inv
                       - lj4[itype][jtype]) - offset[itype][jtype];
            }
          }
        }

        if (rsq < cut_coulsq_global) {
          const double ir = 1.0/sqrt(rsq);
          const double prefactor = qqrd2e * qtmp*q[j];
          const double r2_ia2 = rsq*_ia2;
          const double r4_ia4 = r2_ia2*r2_ia2;
          if (_smooth==C3) {
            forcecoul = prefactor*(_ia3*(-4.375+5.25*r2_ia2-1.875*r4_ia4)-
                                   ir/rsq);
            if (eflag)
              ecoul = prefactor*(ir+_ia*(2.1875-2.1875*r2_ia2+
                                         1.3125*r4_ia4-0.3125*r2_ia2*r4_ia4));
          } else {
            const double r6_ia6 = r2_ia2*r4_ia4;
            forcecoul = prefactor*(_ia3*(-6.5625+11.8125*r2_ia2-8.4375*r4_ia4+
                                         2.1875*r6_ia6)-ir/rsq);
            if (eflag)
              ecoul = prefactor*(ir+_ia*(2.4609375-3.28125*r2_ia2+
                                         2.953125*r4_ia4-1.40625*r6_ia6+
                                         0.2734375*r4_ia4*r4_ia4));
          }
          if (factor_coul < 1.0) {
            forcecoul -= (1.0-factor_coul)*prefactor*ir;
            if (eflag) ecoul -= (1.0-factor_coul)*prefactor*ir;
          }
        }
        fpair = forcecoul + forcelj * r2inv;

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

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulMSM::write_restart(FILE *fp) 
{
  write_restart_settings(fp);
  PairCMMCommon::write_restart(fp);
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulMSM::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();
  PairCMMCommon::read_restart(fp);
}

/* ---------------------------------------------------------------------- */

double PairCGCMMCoulMSM::memory_usage()
{
  double bytes=PairCMMCommon::memory_usage();
  
  int n = atom->ntypes;

  // cut_coul/cut_coulsq/cut_ljsq
  bytes += (n+1)*(n+1)*sizeof(double)*4; 
  
  return bytes;
}

/* ---------------------------------------------------------------------- */

void *PairCGCMMCoulMSM::extract(char *str)
{
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul_global;
  return NULL;
}

/* ---------------------------------------------------------------------- */
