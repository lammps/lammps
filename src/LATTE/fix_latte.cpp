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

#include <stdio.h>
#include <string.h>
#include "fix_latte.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "modify.h"
#include "compute.h"
#include "comm.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

 extern "C" {
   void latte(int*, double*, int*, double*);
  }

#define INVOKED_PERATOM 8

/* ---------------------------------------------------------------------- */

FixLatte::FixLatte(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal fix latte command");

  // store pe/atom ID used for input of Coulomb potential to LATTE
  // insure it is valid for these computations

  int n = strlen(arg[3]) + 1;
  id_pe = new char[n];
  strcpy(id_pe,arg[3]);

  int ipe = modify->find_compute(id_pe);
  if (ipe < 0) error->all(FLERR,"Could not find fix latte compute ID");
  if (modify->compute[ipe]->peatomflag == 0)
    error->all(FLERR,"Fix latte compute ID does not compute pe/atom");

//   latte(arg[4]);

  // initialize LATTE with LAMMPS info about box, atoms, atom types, etc ?
  // may need to be done in init() ??

  // any per-atom quantities to allocate/initialize for LAMMPS?
  // i.e. quantities carried with atoms across timesteps ??
}

/* ---------------------------------------------------------------------- */

FixLatte::~FixLatte()
{
  delete [] id_pe;
}

/* ---------------------------------------------------------------------- */

int FixLatte::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLatte::init()
{
  // error checks

  int ipe = modify->find_compute(id_pe);
  if (ipe < 0) error->all(FLERR,"Could not find fix latte compute ID");
  c_pe = modify->compute[ipe];

  // warn if any integrate fix comes after this one
  // is it actually necessary for q(n) update to come after x,v update ??

  int after = 0;
  int flag = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(id,modify->fix[i]->id) == 0) after = 1;
    else if ((modify->fmask[i] & INITIAL_INTEGRATE) && after) flag = 1;
  }
  if (flag && comm->me == 0)
    error->warning(FLERR,"Fix latte should come after all other "
                   "integration fixes");

  // need a full neighbor list
  // perpetual list, built whenever re-neighboring occurs

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  // integrator timesteps

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  // any more LATTE initialization to do ?
}

/* ---------------------------------------------------------------------- */

void FixLatte::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixLatte::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixLatte::min_setup(int vflag)
{
  // for minimize, what about charge DOFs ??

  post_force(vflag);
}

/* ----------------------------------------------------------------------
   integrate electronic degrees of freedom
------------------------------------------------------------------------- */

void FixLatte::initial_integrate(int vflag)
{
  // what do I do here for q(n) update?
  // is it the same variable as LAMMPS q, I think so
  // is it an Euler update or VV update ??

  /*
  double dtfm;

  // update v and x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
  */
}

/* ---------------------------------------------------------------------- */

void FixLatte::post_force(int vflag)
{
  // what should cutoffs be for passing neighlist info to LATTE ??
  // do cutoffs include many self image atoms for tiny periodic system ??

/*  int i,j,ii,jj,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      // enforce a different cutoff than pair style?
      // what are optimal cutoffs for pair/Kspace for tiny system?
      // operations on I/J pairs if necessary
    } 
  }

  // invoke compute pe/atom
  // wrap with clear/add and trigger pe/atom calculation every step
  // Coulomb potential should just be pe[i]/q[i] ??

  modify->clearstep_compute();

  if (!(c_pe->invoked_flag & INVOKED_PERATOM)) {
    c_pe->compute_peratom();
    c_pe->invoked_flag |= INVOKED_PERATOM;
  }

  double *pe = c_pe->vector_atom;

  modify->addstep_compute(update->ntimestep+1);
*/
  int natoms = (int) atom->natoms;
  latte(&natoms,&atom->x[0][0],atom->type,&atom->f[0][0]);  
//   latte(&natoms,&atom->x[0][0],atom->type);
  
  // construct H0,S,Z
  // setup full Hamiltonian H(R,n)
  // calculate density matrix D and charge q(n)
  // calculate DFTB forces via LATTE = F(R,H0,S,Z,D,q,n)
  // how to pass neighbor list and Coulomb potential to LATTE ??

  // HERE is where main call to LATTE goes to get forces

  // simply add the returned forces to atom->f
  // how to request/get global or per-atom energy back from LATTE ??
  // how to request/get global or per-atom virial back from LATTE ??
}

/* ----------------------------------------------------------------------
   integrate electronic degrees of freedom
------------------------------------------------------------------------- */

void FixLatte::final_integrate()
{
  // possibly nothing to do here if Euler step of q(n) ??

  /*
  double dtfm;

  // update v of atoms in group

  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
    }
  */
}

/* ---------------------------------------------------------------------- */

void FixLatte::reset_dt()
{
  // will we ever auto-adjust the timestep ??

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ---------------------------------------------------------------------- */

double FixLatte::compute_scalar()
{
  // return DFTB global energy
  return 0.0;
}


