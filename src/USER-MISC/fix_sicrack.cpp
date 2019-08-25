/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   //***CETTINA
   //Mode 1 crack
/*------------------------------------------------------------------------- */

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "fix_sicrack.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "lattice.h"
#include "math.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixSICrack::FixSICrack(LAMMPS *lmp, int narg, char **arg) : 
  Fix(lmp, narg, arg)
{
  if (narg < 10) error->all(FLERR,"Illegal fix sicrack command");

  K1 = 0.0;
  dK1 = 0.0;
  K2 = 0.0;
  dK2 = 0.0;
  K3 = 0.0;
  dK3 = 0.0;
  T = 0.0;
  dT = 0.0;

  // read arguments from input and check validity

  xtip = utils::numeric(FLERR,arg[3],true,lmp);  // give coords in (rotated) lattice units, not in angstrøm
  ytip = utils::numeric(FLERR,arg[4],true,lmp); 
  mu = utils::numeric(FLERR,arg[5],true,lmp);  // shear modulus, in GPa
  vu = utils::numeric(FLERR,arg[6],true,lmp);  // poisson ratio

  int iarg = 7;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"M1") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix sicrack command");
      else {
	K1 = utils::numeric(FLERR,arg[iarg+1],true,lmp);  // Mode I stress-intensity factor, in MPa*sqrt(m)
	dK1 = utils::numeric(FLERR,arg[iarg+2],true,lmp);  // Increase in Mode I stress-intensity factor, in MPa*sqrt(m)
	iarg += 3;
      }
    } else if (strcmp(arg[iarg],"M2") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix sicrack command");
      else {
	K2 = utils::numeric(FLERR,arg[iarg+1],true,lmp);  // Mode II stress-intensity factor, in MPa*sqrt(m)
	dK2 = utils::numeric(FLERR,arg[iarg+2],true,lmp);  // Increase in Mode II stress-intensity factor, in MPa*sqrt(m)
	iarg += 3;
      }
    } else if (strcmp(arg[iarg],"M3") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix sicrack command");
      else {
	K3 = utils::numeric(FLERR,arg[iarg+1],true,lmp);  // Mode III stress-intensity factor, in MPa*sqrt(m)
	dK3 = utils::numeric(FLERR,arg[iarg+2],true,lmp);  // Increase in Mode III stress-intensity factor, in MPa*sqrt(m)
	iarg += 3;
      }
    } else if (strcmp(arg[iarg],"T") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix sicrack command");
      else {
	T = utils::numeric(FLERR,arg[iarg+1],true,lmp);  // T-stress, in MPa
	dT = utils::numeric(FLERR,arg[iarg+2],true,lmp);  // Increase in T-stress, in MPa
	iarg += 3;
      }
    } else error->all(FLERR,"Illegal fix sicrack command");
  }

  PI = 4.0*atan(1.0);
  Xk = 3.0-4.0*vu;

  double xlattice = domain->lattice->xlattice;  // Get lattice scale factors
  double ylattice = domain->lattice->ylattice;
  xtip *= xlattice;  // Set the true coords of crack front
  ytip *= ylattice;
  //  printf("K1: %f, dK1: %f, K2: %f, dK2: %f, K3: %f, dK3: %f, T: %f, dT: %f\n", K1, dK1, K2, dK2, K3, dK3, T, dT);

  // perform initial allocation of atom-based array
  // register with Atom class

  xoriginal = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  // xoriginal = initial unwrapped positions of atoms

  double **x = atom->x;
  int *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) domain->unmap(x[i],image[i],xoriginal[i]);
    else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
  }

}

/* ---------------------------------------------------------------------- */

FixSICrack::~FixSICrack()
{
  // unregister callbacks to this fix from Atom class
 
  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored arrays

  memory->destroy(xoriginal);
}

/* ---------------------------------------------------------------------- */


int FixSICrack::setmask()
{
  int mask = 0;
    mask |= FixConst::INITIAL_INTEGRATE;
  mask |= FixConst::FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSICrack::init()
{

  // int ifix = modify->find_fix(id_fix);
  // if (ifix < 0) error->all("Could not find store/state fix ID");
  // fix = modify->fix[ifix];
}

/* ---------------------------------------------------------------------- */

void FixSICrack::setup(int vflag)
{
  //  double **xorig = fix->array_atom;
  double xold[3];
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      xold[0] = x[i][0];
      xold[1] = x[i][1];
      xold[2] = x[i][2];

      xtemp = x[i][0] - xtip;
      ytemp = x[i][1] - ytip;
      rdist = sqrt(xtemp*xtemp + ytemp*ytemp);
      
      // theta is the angle made by vector (x,y) with the crack plane
      
      // a domain error occurs if both xtemp and ytemp are 0
      if (rdist != 0) {
      	theta2 = atan2(ytemp,xtemp);
	theta = 0.5*theta2;
	cost = cos(theta);
	sint = sin(theta);

	ux = 0.5/mu*sqrt(0.5*rdist/MathConst::MY_PI)*(K1*cost*(Xk - 1.0 + 2.0*sint*sint) + K2*sint*(Xk + 1.0 + 2.0*cost*cost))*100;
	uy = 0.5/mu*sqrt(0.5*rdist/MathConst::MY_PI)*(K1*sint*(Xk + 1.0 - 2.0*cost*cost) - K2*cost*(Xk - 1.0 - 2.0*sint*sint))*100;
	uz = K3/mu*sqrt(0.5*rdist/MathConst::MY_PI)*sint*100;
	// factor 100 is to make the units right, so we get a displacement in Ångstrøm

	uTx = T/8/mu*(Xk+1)*rdist*cos(theta2)/1000;
	uTy = T/8/mu*(3-Xk)*rdist*sin(theta2)/1000;
	// factor 1000 is to make the units right, since [T]=MPa and [mu]=GPa
	
	// if (xold[0] == 0.0) {
	//   printf("xold0: %f, xold1: %f, xtemp: %f, ytemp: %f, rdist: %f, ux: %f, uy: %f, theta: %f\n", xold[0], xold[1], xtemp, ytemp, rdist, ux, uy, theta);
	// }
	

	x[i][0] += ux + uTx;
	x[i][1] += uy + uTy;
	x[i][2] += uz;

	domain->remap_near(x[i],xold);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSICrack::initial_integrate(int vflag)
{
  double xold[3];
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      xold[0] = x[i][0];
      xold[1] = x[i][1];
      xold[2] = x[i][2];

      xtemp = x[i][0] - xtip;
      ytemp = x[i][1] - ytip;
      rdist = sqrt(xtemp*xtemp + ytemp*ytemp);
      
      // theta is the angle made by vector (x,y) with the crack plane
      
      // a domain error occurs if both xtemp and ytemp are 0
      if (rdist != 0) {
      	theta2 = atan2(ytemp,xtemp);
	theta = 0.5*theta2;
	cost = cos(theta);
	sint = sin(theta);

	ux = 0.5/mu*sqrt(0.5*rdist/PI)*(dK1*cost*(Xk - 1.0 + 2.0*sint*sint) + dK2*sint*(Xk + 1.0 + 2.0*cost*cost))*100;
	uy = 0.5/mu*sqrt(0.5*rdist/PI)*(dK1*sint*(Xk + 1.0 - 2.0*cost*cost) - dK2*cost*(Xk - 1.0 - 2.0*sint*sint))*100;
	uz = dK3/mu*sqrt(0.5*rdist/PI)*sint*100;
	// factor 100 is to make the units right, so we get a displacement in Ångstrøm

	uTx = dT/8/mu*(Xk+1)*rdist*cos(theta2)/1000;
	uTy = dT/8/mu*(3-Xk)*rdist*sin(theta2)/1000;
	// factor 1000 is to make the units right, since [T]=MPa and [mu]=GPa
	
	x[i][0] += ux + uTx;
	x[i][1] += uy + uTy;
	x[i][2] += uz;

	domain->remap_near(x[i],xold);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSICrack::final_integrate()
{
  //KI = dKI;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixSICrack::memory_usage()
{
  double bytes = atom->nmax*3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixSICrack::grow_arrays(int nmax)
{
  xoriginal =
    memory->grow(xoriginal,nmax,3,"move:xoriginal");
  array_atom = xoriginal;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixSICrack::copy_arrays(int i, int j)
{
  xoriginal[j][0] = xoriginal[i][0];
  xoriginal[j][1] = xoriginal[i][1];
  xoriginal[j][2] = xoriginal[i][2];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixSICrack::pack_exchange(int i, double *buf)
{
  buf[0] = xoriginal[i][0];
  buf[1] = xoriginal[i][1];
  buf[2] = xoriginal[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixSICrack::unpack_exchange(int nlocal, double *buf)
{
  xoriginal[nlocal][0] = buf[0];
  xoriginal[nlocal][1] = buf[1];
  xoriginal[nlocal][2] = buf[2];
  return 3;
}
