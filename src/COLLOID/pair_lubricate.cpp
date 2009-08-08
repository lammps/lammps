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
   Contributing author: Randy Schunk (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_lubricate.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "memory.h"
#include "random_mars.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PairLubricate::PairLubricate(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;

  random = NULL;
}

/* ---------------------------------------------------------------------- */

PairLubricate::~PairLubricate()
{
  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);

    memory->destroy_2d_double_array(cut);
    memory->destroy_2d_double_array(cut_inner);
  }

  delete random;
}

/* ---------------------------------------------------------------------- */

void PairLubricate::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fpair,fx,fy,fz,tx,ty,tz;
  double rsq,r,h_sep,radi,tfmag;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3;
  double vt1,vt2,vt3,w1,w2,w3,v_shear1,v_shear2,v_shear3;
  double omega_t_1,omega_t_2,omega_t_3;
  double n_cross_omega_t_1,n_cross_omega_t_2,n_cross_omega_t_3;
  double wr1,wr2,wr3,wnnr,wn1,wn2,wn3,inv_inertia;
  double P_dot_wrel_1,P_dot_wrel_2,P_dot_wrel_3;
  double a_squeeze,a_shear,a_pump,a_twist;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double PI = 4.0*atan(1.0);

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double vxmu2f = force->vxmu2f;

  double prethermostat = sqrt(2.0 * force->boltz * t_target / update->dt);
  prethermostat *= sqrt(force->vxmu2f/force->ftm2v/force->mvv2e);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  a_squeeze = a_shear = a_pump = a_twist = 0.0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = atom->shape[itype][0];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {

	r = sqrt(rsq);

        // relative translational velocity 

        vr1 = v[i][0] - v[j][0];
        vr2 = v[i][1] - v[j][1];
        vr3 = v[i][2] - v[j][2];

        // normal component N.(v1-v2) = nn.(v1-v2)

        vnnr = vr1*delx + vr2*dely + vr3*delz;
	vnnr /= r;
	vn1 = delx*vnnr / r;
        vn2 = dely*vnnr / r;
        vn3 = delz*vnnr / r;

        // tangential component -P.(v1-v2)
	// P = (I - nn) where n is vector between centers
     
        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // additive rotational velocity = omega_1 + omega_2
	// use omega directly if it exists, else angmom
	// angular momentum = I*omega = 2/5 * M*R^2 * omega

	if (omega_flag) {
	  w1 = omega[i][0] + omega[j][0];
	  w2 = omega[i][1] + omega[j][1];
	  w3 = omega[i][2] + omega[j][2];
	} else {
	  inv_inertia = 1.0 / (0.4*atom->mass[itype]*radi*radi);
	  w1 = inv_inertia * (angmom[i][0] + angmom[j][0]);
	  w2 = inv_inertia * (angmom[i][1] + angmom[j][1]);
	  w3 = inv_inertia * (angmom[i][2] + angmom[j][2]);
	}

        // relative velocities n X P . (v1-v2) = n X (I-nn) . (v1-v2)

        v_shear1 = (dely*vt3 - delz*vt2) / r;
        v_shear2 = -(delx*vt3 - delz*vt1) / r;
        v_shear3 = (delx*vt2 - dely*vt1) / r;

        // relative rotation rate P.(omega1 + omega2)

	omega_t_1 = w1 - delx*(delx*w1) / rsq;
	omega_t_2 = w2 - dely*(dely*w2) / rsq;
	omega_t_3 = w3 - delz*(delz*w3) / rsq;

        // n X omega_t

        n_cross_omega_t_1 =  (dely*omega_t_3 - delz*omega_t_2) / r;
        n_cross_omega_t_2 =  -(delx*omega_t_3 - delz*omega_t_1) / r;
        n_cross_omega_t_3 =  (delx*omega_t_2 - dely*omega_t_1) / r;

        // N.(w1-w2) and P.(w1-w2)

	if (omega_flag) {
	  wr1 = omega[i][0] - omega[j][0];
	  wr2 = omega[i][1] - omega[j][1];
	  wr3 = omega[i][2] - omega[j][2];
	} else {
	  wr1 = inv_inertia * (angmom[i][0] - angmom[j][0]);
	  wr2 = inv_inertia * (angmom[i][1] - angmom[j][1]);
	  wr3 = inv_inertia * (angmom[i][2] - angmom[j][2]);
	}
 
	wnnr = wr1*delx + wr2*dely + wr3*delz;
	wn1 = delx*wnnr / rsq;
	wn2 = dely*wnnr / rsq;
	wn3 = delz*wnnr / rsq;

        P_dot_wrel_1 = wr1 - delx*(delx*wr1)/rsq; 
        P_dot_wrel_2 = wr2 - dely*(dely*wr2)/rsq; 
        P_dot_wrel_3 = wr3 - delz*(delz*wr3)/rsq; 

        // compute components of pair-hydro

        h_sep = r - 2.0*radi;

	if (flag1)
	  a_squeeze = (3.0*PI*mu*2.0*radi/2.0) * (2.0*radi/4.0/h_sep);
	if (flag2) 
	  a_shear = (PI*mu*2.*radi/2.0) *
	    log(2.0*radi/2.0/h_sep)*(2.0*radi+h_sep)*(2.0*radi+h_sep)/4.0;
	if (flag3) 
	  a_pump = (PI*mu*pow(2.0*radi,4)/8.0) *
	    ((3.0/20.0) * log(2.0*radi/2.0/h_sep) + 
	     (63.0/250.0) * (h_sep/2.0/radi) * log(2.0*radi/2.0/h_sep));
	if (flag4)
	  a_twist = (PI*mu*pow(2.0*radi,4)/4.0) *
	    (h_sep/2.0/radi) * log(2.0/(2.0*h_sep));

        if (h_sep >= cut_inner[itype][jtype]) {
          fx = -a_squeeze*vn1 - a_shear*(2.0/r)*(2.0/r)*vt1 + 
	    (2.0/r)*a_shear*n_cross_omega_t_1;
          fy = -a_squeeze*vn2 - a_shear*(2.0/r)*(2.0/r)*vt2 + 
	    (2.0/r)*a_shear*n_cross_omega_t_2;
          fz = -a_squeeze*vn3 - a_shear*(2.0/r)*(2.0/r)*vt3 +
	    (2.0/r)*a_shear*n_cross_omega_t_3;
	  fx *= vxmu2f;
	  fy *= vxmu2f;
	  fz *= vxmu2f;

	  // add in thermostat force

	  tfmag = prethermostat*sqrt(a_squeeze)*(random->uniform()-0.5);
	  fx -= tfmag * delx/r;
	  fy -= tfmag * dely/r;
	  fz -= tfmag * delz/r;
	  
	  tx = -(2.0/r)*a_shear*v_shear1 - a_shear*omega_t_1 - 
	    a_pump*P_dot_wrel_1 - a_twist*wn1;
	  ty = -(2.0/r)*a_shear*v_shear2 - a_shear*omega_t_2 - 
	    a_pump*P_dot_wrel_2 - a_twist*wn2;
	  tz = -(2.0/r)*a_shear*v_shear3 - a_shear*omega_t_3 - 
	    a_pump*P_dot_wrel_3 - a_twist*wn3;
	  torque[i][0] += vxmu2f * tx;
	  torque[i][1] += vxmu2f * ty;
	  torque[i][2] += vxmu2f * tz;

        } else {
	  a_squeeze = (3.0*PI*mu*2.0*radi/2.0) * 
	    (2.0*radi/4.0/cut_inner[itype][jtype]);
	  fpair = -a_squeeze*vnnr;
	  fpair *= vxmu2f;

	  // add in thermostat force

	  fpair -= prethermostat*sqrt(a_squeeze)*(random->uniform()-0.5);

	  fx = fpair * delx/r;
	  fy = fpair * dely/r;
	  fz = fpair * delz/r;
	}

	f[i][0] += fx;
	f[i][1] += fy;
	f[i][2] += fz;

	if (newton_pair || j < nlocal) {
	  f[j][0] -= fx;
	  f[j][1] -= fy;
	  f[j][2] -= fz;

	  if (h_sep >= cut_inner[itype][jtype]) {
	    tx = -(2.0/r)*a_shear*v_shear1 - a_shear*omega_t_1 + 
	      a_pump*P_dot_wrel_1 + a_twist*wn1;
	    ty = -(2.0/r)*a_shear*v_shear2 - a_shear*omega_t_2 + 
	      a_pump*P_dot_wrel_2 + a_twist*wn2;
	    tz = -(2.0/r)*a_shear*v_shear3 - a_shear*omega_t_3 + 
	      a_pump*P_dot_wrel_3 + a_twist*wn3;
	    torque[j][0] += vxmu2f * tx;
	    torque[j][1] += vxmu2f * ty;
	    torque[j][2] += vxmu2f * tz;
	  }
	}

	if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
				 0.0,0.0,fx,fy,fz,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays 
------------------------------------------------------------------------- */

void PairLubricate::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  cut = memory->create_2d_double_array(n+1,n+1,"pair:cut");
  cut_inner = memory->create_2d_double_array(n+1,n+1,"pair:cut_inner");
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairLubricate::settings(int narg, char **arg)
{
  if (narg != 9) error->all("Illegal pair_style command");

  mu = force->numeric(arg[0]);
  flag1 = force->inumeric(arg[1]);
  flag2 = force->inumeric(arg[2]);
  flag3 = force->inumeric(arg[3]);
  flag4 = force->inumeric(arg[4]);
  cut_inner_global = force->numeric(arg[5]);
  cut_global = force->numeric(arg[6]);
  t_target = force->numeric(arg[7]);
  seed = force->inumeric(arg[8]);

  // initialize Marsaglia RNG with processor-unique seed

  if (seed <= 0) error->all("Illegal pair_style command");
  delete random;
  random = new RanMars(lmp,seed + comm->me);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	if (setflag[i][j]) {
	  cut_inner[i][j] = cut_inner_global;
	  cut[i][j] = cut_global;
	}
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLubricate::coeff(int narg, char **arg)
{
  if (narg != 2 && narg != 4)
    error->all("Incorrect args for pair coefficients");

  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double cut_inner_one = cut_inner_global;
  double cut_one = cut_global;
  if (narg == 4) {
    cut_inner_one = force->numeric(arg[2]);
    cut_one = force->numeric(arg[3]);
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut_inner[i][j] = cut_inner_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLubricate::init_style()
{
  if (!atom->torque_flag || !atom->avec->shape_type)
    error->all("Pair lubricate requires atom attributes torque and shape");
  if (atom->radius_flag || atom->rmass_flag)
    error->all("Pair lubricate cannot be used with atom attributes "
	       "diameter or rmass");

  if (atom->omega_flag) omega_flag = 1;
  else if (atom->angmom_flag) omega_flag = 0;
  else error->all("Pair lubricate requires atom attribute omega or angmom");

  // insure all particle shapes are finite-size, spherical, and monodisperse

  double value = atom->shape[1][0];
  if (value == 0.0) error->all("Pair lubricate requires extended particles");
  for (int i = 1; i <= atom->ntypes; i++)
    if (atom->shape[i][0] != value || atom->shape[i][0] != value || 
	atom->shape[i][0] != value)
      error->all("Pair lubricate requires spherical, mono-disperse particles");
  
  int irequest = neighbor->request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLubricate::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    cut_inner[i][j] = mix_distance(cut_inner[i][i],cut_inner[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  cut_inner[j][i] = cut_inner[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file 
------------------------------------------------------------------------- */

void PairLubricate::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
	fwrite(&cut_inner[i][j],sizeof(double),1,fp);
	fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLubricate::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
	if (me == 0) {
	  fread(&cut_inner[i][j],sizeof(double),1,fp);
	  fread(&cut[i][j],sizeof(double),1,fp);
	}
	MPI_Bcast(&cut_inner[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLubricate::write_restart_settings(FILE *fp)
{
  fwrite(&mu,sizeof(double),1,fp);
  fwrite(&flag1,sizeof(int),1,fp);
  fwrite(&flag2,sizeof(int),1,fp);
  fwrite(&flag3,sizeof(int),1,fp);
  fwrite(&flag4,sizeof(int),1,fp);
  fwrite(&cut_inner_global,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&t_target,sizeof(double),1,fp);
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLubricate::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&mu,sizeof(double),1,fp);
    fread(&flag1,sizeof(int),1,fp);
    fread(&flag2,sizeof(int),1,fp);
    fread(&flag3,sizeof(int),1,fp);
    fread(&flag4,sizeof(int),1,fp);
    fread(&cut_inner_global,sizeof(double),1,fp);
    fread(&cut_global,sizeof(double),1,fp);
    fread(&t_target, sizeof(double),1,fp);
    fread(&seed, sizeof(int),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&mu,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&flag1,1,MPI_INT,0,world);
  MPI_Bcast(&flag2,1,MPI_INT,0,world);
  MPI_Bcast(&flag3,1,MPI_INT,0,world);
  MPI_Bcast(&flag4,1,MPI_INT,0,world);
  MPI_Bcast(&cut_inner_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&t_target,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&seed,1,MPI_INT,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);

  // additional setup based on restart parameters

  delete random;
  random = new RanMars(lmp,seed + comm->me);
}
