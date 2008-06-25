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
   Contributing authors: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_gran_history.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_pour.h"
#include "fix_shear_history.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PairGranHistory::PairGranHistory(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  history = 1;
  fix_history = NULL;
}

/* ---------------------------------------------------------------------- */

PairGranHistory::~PairGranHistory()
{
  if (fix_history) modify->delete_fix("SHEAR_HISTORY");

  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairGranHistory::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double xmeff,damp,ccel,tor1,tor2,tor3;
  double fn,fs,fs1,fs2,fs3;
  double shrmag,rsht;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = listgranhistory->firstneigh;
  firstshear = listgranhistory->firstdouble;
 
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    touch = firsttouch[i];
    allshear = firstshear[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      if (rsq >= radsum*radsum) {

	// unset non-touching neighbors

        touch[jj] = 0;
	shear = &allshear[3*jj];
        shear[0] = 0.0;
        shear[1] = 0.0;
        shear[2] = 0.0;

      } else {
	r = sqrt(rsq);

	// relative translational velocity

	vr1 = v[i][0] - v[j][0];
	vr2 = v[i][1] - v[j][1];
	vr3 = v[i][2] - v[j][2];

 	vr1 *= dt;
 	vr2 *= dt;
 	vr3 *= dt;

	// normal component

	vnnr = vr1*delx + vr2*dely + vr3*delz;
	vn1 = delx*vnnr / rsq;
	vn2 = dely*vnnr / rsq;
	vn3 = delz*vnnr / rsq;

	// tangential component

	vt1 = vr1 - vn1;
	vt2 = vr2 - vn2;
	vt3 = vr3 - vn3;

	// relative rotational velocity

	wr1 = radi*omega[i][0] + radj*omega[j][0];
	wr2 = radi*omega[i][1] + radj*omega[j][1];
	wr3 = radi*omega[i][2] + radj*omega[j][2];

	wr1 *= dt/r;
	wr2 *= dt/r;
	wr3 *= dt/r;

	// normal damping term
	// this definition of DAMP includes the extra 1/r term

	xmeff = rmass[i]*rmass[j] / (rmass[i]+rmass[j]);
	if (mask[i] & freeze_group_bit) xmeff = rmass[j];
	if (mask[j] & freeze_group_bit) xmeff = rmass[i];
	damp = xmeff*gamman_dl*vnnr/rsq;
	ccel = xkk*(radsum-r)/r - damp;

	// relative velocities

	vtr1 = vt1 - (delz*wr2-dely*wr3);
	vtr2 = vt2 - (delx*wr3-delz*wr1);
	vtr3 = vt3 - (dely*wr1-delx*wr2);
	vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
	vrel = sqrt(vrel);

	// shear history effects
	// shrmag = magnitude of shear

	touch[jj] = 1;
	shear = &allshear[3*jj];
        shear[0] += vtr1;
        shear[1] += vtr2;
        shear[2] += vtr3;
        shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] +
		      shear[2]*shear[2]);

	// rotate shear displacements correctly

	rsht = shear[0]*delx + shear[1]*dely + shear[2]*delz;
	rsht /= rsq;
        shear[0] -= rsht*delx;
        shear[1] -= rsht*dely;
        shear[2] -= rsht*delz;

	// tangential forces

	fs1 = - (xkkt*shear[0] + xmeff*gammas_dl*vtr1);
	fs2 = - (xkkt*shear[1] + xmeff*gammas_dl*vtr2);
	fs3 = - (xkkt*shear[2] + xmeff*gammas_dl*vtr3);

	// force normalization
	// rescale frictional displacements and forces if needed

	fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
	fn = xmu * fabs(ccel*r);

	if (fs > fn) {
	  if (shrmag != 0.0) {
	    shear[0] = (fn/fs) * (shear[0] + xmeff*gammas_dl*vtr1/xkkt) -
	      xmeff*gammas_dl*vtr1/xkkt;
	    shear[1] = (fn/fs) * (shear[1] + xmeff*gammas_dl*vtr2/xkkt) -
	      xmeff*gammas_dl*vtr2/xkkt;
	    shear[2] = (fn/fs) * (shear[2] + xmeff*gammas_dl*vtr3/xkkt) -
	      xmeff*gammas_dl*vtr3/xkkt;
	    fs1 *= fn/fs;
	    fs2 *= fn/fs;
	    fs3 *= fn/fs;
	  } else {
	    fs1 = 0.0;
	    fs2 = 0.0;
	    fs3 = 0.0;
	  }
	}

	// forces & torques

	fx = delx*ccel + fs1;
	fy = dely*ccel + fs2;
	fz = delz*ccel + fs3;
	f[i][0] += fx;
	f[i][1] += fy;
	f[i][2] += fz;

	rinv = 1/r;
	tor1 = rinv * (dely*fs3 - delz*fs2);
	tor2 = rinv * (delz*fs1 - delx*fs3);
	tor3 = rinv * (delx*fs2 - dely*fs1);
	torque[i][0] -= radi*tor1;
	torque[i][1] -= radi*tor2;
	torque[i][2] -= radi*tor3;

	if (j < nlocal) {
	  f[j][0] -= fx;
	  f[j][1] -= fy;
	  f[j][2] -= fz;
	  torque[j][0] -= radj*tor1;
	  torque[j][1] -= radj*tor2;
	  torque[j][2] -= radj*tor3;
	}

	if (evflag) ev_tally_xyz(i,j,nlocal,0,
				 0.0,0.0,fx,fy,fz,delx,dely,delz);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairGranHistory::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranHistory::settings(int narg, char **arg)
{
  if (narg != 4) error->all("Illegal pair_style command");

  xkk = atof(arg[0]);
  gamman = atof(arg[1]);
  xmu = atof(arg[2]);
  dampflag = atoi(arg[3]);

  // granular styles do not use pair_coeff, so set setflag for everything now

  if (!allocated) allocate();

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++)
      setflag[i][j] = 1;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGranHistory::coeff(int narg, char **arg)
{
  error->all("Granular pair styles do not use pair_coeff settings");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGranHistory::init_style()
{
  int i;

  if (!atom->radius_flag || !atom->omega_flag || !atom->torque_flag)
    error->all("Pair granular requires atom attributes radius, omega, torque");
  if (atom->avec->ghost_velocity == 0)
    error->all("Pair granular requires ghost atoms store velocity");

  // need a half neigh list and optionally a granular history neigh list

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->gran = 1;
  if (history) {
    irequest = neighbor->request(this);
    neighbor->requests[irequest]->id = 1;
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->granhistory = 1;
    neighbor->requests[irequest]->dnum = 3;
  }

  xkkt = xkk * 2.0/7.0;
  dt = update->dt;
  double gammas = 0.5*gamman;
  if (dampflag == 0) gammas = 0.0;
  gamman_dl = gamman/dt;
  gammas_dl = gammas/dt;

  // if shear history is stored:
  // check if newton flag is valid
  // if first init, create Fix needed for storing shear history

  if (history && force->newton_pair == 1)
    error->all("Potential with shear history requires newton pair off");

  if (history && fix_history == NULL) {
    char **fixarg = new char*[3];
    fixarg[0] = (char *) "SHEAR_HISTORY";
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "SHEAR_HISTORY";
    modify->add_fix(3,fixarg);
    delete [] fixarg;
    fix_history = (FixShearHistory *) modify->fix[modify->nfix-1];
    fix_history->pair = this;
  }

  // check for freeze Fix and set freeze_group_bit

  for (i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"freeze") == 0) break;
  if (i < modify->nfix) freeze_group_bit = modify->fix[i]->groupbit;
  else freeze_group_bit = 0;

  // set cutoff by largest particles
  // maxrad_dynamic = radius of largest dynamic particle, including inserted
  // maxrad_frozen = radius of largest dynamic particle
  // include frozen-dynamic interactions
  // do not include frozen-frozen interactions
  // include future inserted particles as dynamic
  // cutforce was already set in pair::init(), but this sets it correctly

  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double maxrad_dynamic = 0.0;
  for (i = 0; i < nlocal; i++)
    if (!(mask[i] & freeze_group_bit))
      maxrad_dynamic = MAX(maxrad_dynamic,radius[i]);
  double mine = maxrad_dynamic;
  MPI_Allreduce(&mine,&maxrad_dynamic,1,MPI_DOUBLE,MPI_MAX,world);

  for (i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"pour") == 0)
      maxrad_dynamic =
	MAX(maxrad_dynamic,((FixPour *) modify->fix[i])->radius_hi);

  double maxrad_frozen = 0.0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & freeze_group_bit)
      maxrad_frozen = MAX(maxrad_frozen,radius[i]);
  mine = maxrad_frozen;
  MPI_Allreduce(&mine,&maxrad_frozen,1,MPI_DOUBLE,MPI_MAX,world);

  cutforce = maxrad_dynamic + MAX(maxrad_dynamic,maxrad_frozen);
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   optional granular history list
------------------------------------------------------------------------- */

void PairGranHistory::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  else if (id == 1) listgranhistory = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGranHistory::init_one(int i, int j)
{
  if (!allocated) allocate();

  // return max diameter of any particle
  // not used in granular neighbor calculation since particles have radius
  // but will insure cutoff is big enough for any other neighbor lists built
  // if no particles, return 1.0

  double *radius = atom->radius;
  int nlocal = atom->nlocal;

  double maxrad = 0.0;
  for (int m = 0; m < nlocal; m++) maxrad = MAX(maxrad,radius[m]);
  double mine = 2.0*maxrad;

  double maxdiam;
  MPI_Allreduce(&mine,&maxdiam,1,MPI_DOUBLE,MPI_MAX,world);
  if (maxdiam == 0.0) maxdiam = 1.0;

  return maxdiam;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGranHistory::write_restart(FILE *fp)
{
  write_restart_settings(fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGranHistory::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGranHistory::write_restart_settings(FILE *fp)
{
  fwrite(&xkk,sizeof(double),1,fp);
  fwrite(&gamman,sizeof(double),1,fp);
  fwrite(&xmu,sizeof(double),1,fp);
  fwrite(&dampflag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGranHistory::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&xkk,sizeof(double),1,fp);
    fread(&gamman,sizeof(double),1,fp);
    fread(&xmu,sizeof(double),1,fp);
    fread(&dampflag,sizeof(int),1,fp);
  }
  MPI_Bcast(&xkk,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&gamman,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&xmu,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&dampflag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

void *PairGranHistory::extract(char *str)
{
  if (strcmp(str,"xkk") == 0) return (void *) &xkk;
  else if (strcmp(str,"gamman") == 0) return (void *) &gamman;
  else if (strcmp(str,"xmu") == 0) return (void *) &xmu;
  else if (strcmp(str,"dampflag") == 0) return (void *) &dampflag;
  return NULL;
}

/* ---------------------------------------------------------------------- */

void PairGranHistory::reset_dt()
{
  dt = update->dt;
  double gammas = 0.5*gamman;
  if (dampflag == 0) gammas = 0.0;
  gamman_dl = gamman/dt;
  gammas_dl = gammas/dt;
}
