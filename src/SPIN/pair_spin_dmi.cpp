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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "fix.h"
#include "fix_nve_spin.h"
#include "pair_hybrid.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "pair_spin_dmi.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairSpinDmi::PairSpinDmi(LAMMPS *lmp) : PairSpin(lmp),
lockfixnvespin(NULL)
{
  single_enable = 0;
  no_virial_fdotr_compute = 1;
  lattice_flag = 0;
}

/* ---------------------------------------------------------------------- */

PairSpinDmi::~PairSpinDmi()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cut_spin_dmi);
    memory->destroy(DM);
    memory->destroy(v_dmx);
    memory->destroy(v_dmy);
    memory->destroy(v_dmz);
    memory->destroy(vmech_dmx);
    memory->destroy(vmech_dmy);
    memory->destroy(vmech_dmz);
    memory->destroy(cutsq);
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSpinDmi::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2)
    error->all(FLERR,"Incorrect number of args in pair/spin/dmi command");

  if (strcmp(update->unit_style,"metal") != 0)
    error->all(FLERR,"Spin simulations require metal unit style");

  cut_spin_dmi_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = i+1; j <= atom->ntypes; j++) {
        if (setflag[i][j]) {
          cut_spin_dmi[i][j] = cut_spin_dmi_global;
        }
      }
    }
  }

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type spin pairs (only one for now)
------------------------------------------------------------------------- */

void PairSpinDmi::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  // check if args correct

  if (strcmp(arg[2],"dmi") != 0)
    error->all(FLERR,"Incorrect args in pair_style command");
  if (narg != 8)
    error->all(FLERR,"Incorrect args in pair_style command");

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  const double rij = force->numeric(FLERR,arg[3]);
  const double dm = (force->numeric(FLERR,arg[4]));
  double dmx = force->numeric(FLERR,arg[5]);
  double dmy = force->numeric(FLERR,arg[6]);
  double dmz = force->numeric(FLERR,arg[7]);

  double inorm = 1.0/(dmx*dmx+dmy*dmy+dmz*dmz);
  dmx *= inorm;
  dmy *= inorm;
  dmz *= inorm;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut_spin_dmi[i][j] = rij;
      DM[i][j] = dm;
      v_dmx[i][j] = dmx * dm / hbar;
      v_dmy[i][j] = dmy * dm / hbar;
      v_dmz[i][j] = dmz * dm / hbar;
      vmech_dmx[i][j] = dmx * dm;
      vmech_dmy[i][j] = dmy * dm;
      vmech_dmz[i][j] = dmz * dm;
      setflag[i][j] = 1;
      count++;
    }
  }
  if (count == 0)
    error->all(FLERR,"Incorrect args in pair_style command");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSpinDmi::init_style()
{
  if (!atom->sp_flag)
    error->all(FLERR,"Pair spin requires atom/spin style");

  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  // checking if nve/spin is a listed fix

  int ifix = 0;
  while (ifix < modify->nfix) {
    if (strcmp(modify->fix[ifix]->style,"nve/spin") == 0) break;
    ifix++;
  }
  // test remove if test
  //if (ifix == modify->nfix)
  //  error->all(FLERR,"pair/spin style requires nve/spin");

  // get the lattice_flag from nve/spin

  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"nve/spin") == 0) {
      lockfixnvespin = (FixNVESpin *) modify->fix[i];
      lattice_flag = lockfixnvespin->lattice_flag;
    }
  }

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSpinDmi::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  DM[j][i] = DM[i][j];
  v_dmx[j][i] = v_dmx[i][j];
  v_dmy[j][i] = v_dmy[i][j];
  v_dmz[j][i] = v_dmz[i][j];
  vmech_dmx[j][i] = vmech_dmx[i][j];
  vmech_dmy[j][i] = vmech_dmy[i][j];
  vmech_dmz[j][i] = vmech_dmz[i][j];
  cut_spin_dmi[j][i] = cut_spin_dmi[i][j];

  return cut_spin_dmi_global;
}

/* ----------------------------------------------------------------------
   extract the larger cutoff
------------------------------------------------------------------------- */

void *PairSpinDmi::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"cut") == 0) return (void *) &cut_spin_dmi_global;
  return NULL;
}

/* ---------------------------------------------------------------------- */

void PairSpinDmi::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double evdwl, ecoul;
  double xi[3], eij[3];
  double delx,dely,delz;
  double spi[3], spj[3];
  double fi[3], fmi[3];
  double local_cut2;
  double rsq, inorm;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double **fm = atom->fm;
  double **sp = atom->sp;	
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // dmi computation
  // loop over all atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    xi[0] = x[i][0];
    xi[1] = x[i][1];
    xi[2] = x[i][2];
    spi[0] = sp[i][0];
    spi[1] = sp[i][1];
    spi[2] = sp[i][2];


    // loop on neighbors

    for (jj = 0; jj < jnum; jj++) {

      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      spj[0] = sp[j][0];
      spj[1] = sp[j][1];
      spj[2] = sp[j][2];

      evdwl = 0.0;
      fi[0] = fi[1] = fi[2] = 0.0;
      fmi[0] = fmi[1] = fmi[2] = 0.0;

      delx = xi[0] - x[j][0];
      dely = xi[1] - x[j][1];
      delz = xi[2] - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      inorm = 1.0/sqrt(rsq);
      eij[0] = -inorm*delx;
      eij[1] = -inorm*dely;
      eij[2] = -inorm*delz;

      local_cut2 = cut_spin_dmi[itype][jtype]*cut_spin_dmi[itype][jtype];

      // compute dmi interaction

      if (rsq <= local_cut2) {
	compute_dmi(i,j,eij,fmi,spj);
	if (lattice_flag) {
	  compute_dmi_mech(i,j,rsq,eij,fi,spi,spj);
	}
      }

      f[i][0] += fi[0];	
      f[i][1] += fi[1];	  	
      f[i][2] += fi[2];
      fm[i][0] += fmi[0];	
      fm[i][1] += fmi[1];	  	
      fm[i][2] += fmi[2];

      if (newton_pair || j < nlocal) {
	f[j][0] -= fi[0];	
        f[j][1] -= fi[1];	  	
        f[j][2] -= fi[2];
      }

      if (eflag) {
	evdwl -= (spi[0]*fmi[0] + spi[1]*fmi[1] + spi[2]*fmi[2]);
	evdwl *= hbar;
      } else evdwl = 0.0;

      if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
	  evdwl,ecoul,fi[0],fi[1],fi[2],delx,dely,delz);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();

}

/* ---------------------------------------------------------------------- */

void PairSpinDmi::compute_single_pair(int ii, double fmi[3])
{
  int *type = atom->type;
  double **x = atom->x;
  double **sp = atom->sp;
  double local_cut2;
  double xi[3], eij[3];
  double delx,dely,delz;
  double spj[3];

  int i,j,jnum,itype,jtype;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double rsq, inorm;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  i = ilist[ii];
  itype = type[i];

  xi[0] = x[i][0];
  xi[1] = x[i][1];
  xi[2] = x[i][2];

  jlist = firstneigh[i];
  jnum = numneigh[i];

  for (int jj = 0; jj < jnum; jj++) {

    j = jlist[jj];
    j &= NEIGHMASK;
    jtype = type[j];

    spj[0] = sp[j][0];
    spj[1] = sp[j][1];
    spj[2] = sp[j][2];

    delx = xi[0] - x[j][0];
    dely = xi[1] - x[j][1];
    delz = xi[2] - x[j][2];
    rsq = delx*delx + dely*dely + delz*delz;
    inorm = 1.0/sqrt(rsq);
    eij[0] = -inorm*delx;
    eij[1] = -inorm*dely;
    eij[2] = -inorm*delz;

    local_cut2 = cut_spin_dmi[itype][jtype]*cut_spin_dmi[itype][jtype];

    if (rsq <= local_cut2) {
      compute_dmi(i,j,eij,fmi,spj);
    }

  }

}

/* ----------------------------------------------------------------------
   compute the dmi interaction between spin i and spin j
------------------------------------------------------------------------- */

void PairSpinDmi::compute_dmi(int i, int j, double eij[3], double fmi[3], double spj[3])
{
  int *type = atom->type;
  int itype, jtype;
  double dmix, dmiy, dmiz;	
  itype = type[i];
  jtype = type[j];

  dmix = eij[1]*v_dmz[itype][jtype] - eij[2]*v_dmy[itype][jtype];
  dmiy = eij[2]*v_dmx[itype][jtype] - eij[0]*v_dmz[itype][jtype];
  dmiz = eij[0]*v_dmy[itype][jtype] - eij[1]*v_dmx[itype][jtype];

  fmi[0] -= (spj[1]*dmiz - spj[2]*dmiy);
  fmi[1] -= (spj[2]*dmix - spj[0]*dmiz);
  fmi[2] -= (spj[0]*dmiy - spj[1]*dmix);
}

/* ----------------------------------------------------------------------
   compute the mechanical force due to the dmi interaction between atom i and atom j
------------------------------------------------------------------------- */

void PairSpinDmi::compute_dmi_mech(int i, int j, double rsq, double /*eij*/[3],
    double fi[3],  double spi[3], double spj[3])
{
  int *type = atom->type;
  int itype, jtype;
  double dmix,dmiy,dmiz;	
  itype = type[i];
  jtype = type[j];
  double csx,csy,csz,cdmx,cdmy,cdmz,irij;

  irij = 1.0/sqrt(rsq);

  dmix = vmech_dmx[itype][jtype];
  dmiy = vmech_dmy[itype][jtype];
  dmiz = vmech_dmz[itype][jtype];

  csx = (spi[1]*spj[2] - spi[2]*spj[1]);
  csy = (spi[2]*spj[0] - spi[0]*spj[2]);
  csz = (spi[0]*spj[1] - spi[1]*spj[0]);

  cdmx = (dmiy*csz - dmiz*csy);
  cdmy = (dmiz*csx - dmix*csz);
  cdmz = (dmix*csy - dmiy*csz);

  fi[0] += irij*cdmx;
  fi[1] += irij*cdmy;
  fi[2] += irij*cdmz;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSpinDmi::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cut_spin_dmi,n+1,n+1,"pair:cut_spin_dmi");
  memory->create(DM,n+1,n+1,"pair:DM");
  memory->create(v_dmx,n+1,n+1,"pair:DM_vector_x");
  memory->create(v_dmy,n+1,n+1,"pair:DM_vector_y");
  memory->create(v_dmz,n+1,n+1,"pair:DM_vector_z");
  memory->create(vmech_dmx,n+1,n+1,"pair:DMmech_vector_x");
  memory->create(vmech_dmy,n+1,n+1,"pair:DMmech_vector_y");
  memory->create(vmech_dmz,n+1,n+1,"pair:DMmech_vector_z");

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

}


/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpinDmi::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
	fwrite(&DM[i][j],sizeof(double),1,fp);
        fwrite(&v_dmx[i][j],sizeof(double),1,fp);
        fwrite(&v_dmy[i][j],sizeof(double),1,fp);
        fwrite(&v_dmz[i][j],sizeof(double),1,fp);
        fwrite(&vmech_dmx[i][j],sizeof(double),1,fp);
        fwrite(&vmech_dmy[i][j],sizeof(double),1,fp);
        fwrite(&vmech_dmz[i][j],sizeof(double),1,fp);
        fwrite(&cut_spin_dmi[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpinDmi::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&DM[i][j],sizeof(double),1,fp);
          fread(&v_dmx[i][j],sizeof(double),1,fp);
          fread(&v_dmy[i][j],sizeof(double),1,fp);
          fread(&v_dmz[i][j],sizeof(double),1,fp);
          fread(&vmech_dmx[i][j],sizeof(double),1,fp);
          fread(&vmech_dmy[i][j],sizeof(double),1,fp);
          fread(&vmech_dmz[i][j],sizeof(double),1,fp);
          fread(&cut_spin_dmi[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&DM[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&v_dmx[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&v_dmy[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&v_dmz[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&vmech_dmx[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&vmech_dmy[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&vmech_dmz[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_spin_dmi[i][j],1,MPI_DOUBLE,0,world);
      }
    }
  }
}


/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpinDmi::write_restart_settings(FILE *fp)
{
  fwrite(&cut_spin_dmi_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpinDmi::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_spin_dmi_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_spin_dmi_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}
