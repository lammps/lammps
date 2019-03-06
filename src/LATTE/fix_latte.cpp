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
   Contributing author: Christian Negre (LANL)
------------------------------------------------------------------------- */

#include <cstdio>
#include <cstring>
#include "fix_latte.h"
#include "atom.h"
#include "comm.h"
#include "update.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "modify.h"
#include "compute.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

extern "C" {
  void latte(int *, int *, double *, int *, int *,
             double *, double *, double *, double *,
             double *, double *, double *, int *,
             double *, double *, double *, double *, int * , bool *);
  int latte_abiversion();
}

// the ABIVERSION number here must be kept consistent
// with its counterpart in the LATTE library and the
// prototype above. We want to catch mismatches with
// a meaningful error messages, as they can cause
// difficult to debug crashes or memory corruption.

#define LATTE_ABIVERSION 20180622
#define INVOKED_PERATOM 8

/* ---------------------------------------------------------------------- */

FixLatte::FixLatte(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(update->unit_style,"metal") != 0)
    error->all(FLERR,"Must use units metal with fix latte command");

  if (comm->nprocs != 1)
    error->all(FLERR,"Fix latte currently runs only in serial");

  if (LATTE_ABIVERSION != latte_abiversion())
    error->all(FLERR,"LAMMPS is linked against incompatible LATTE library");

  if (narg != 4) error->all(FLERR,"Illegal fix latte command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  virial_flag = 1;
  thermo_virial = 1;

  // store ID of compute pe/atom used to generate Coulomb potential for LATTE
  // NULL means LATTE will compute Coulombic potential

  coulomb = 0;
  id_pe = NULL;

  if (strcmp(arg[3],"NULL") != 0) {
    coulomb = 1;
    error->all(FLERR,"Fix latte does not yet support a LAMMPS calculation "
               "of a Coulomb potential");

    int n = strlen(arg[3]) + 1;
    id_pe = new char[n];
    strcpy(id_pe,arg[3]);

    int ipe = modify->find_compute(id_pe);
    if (ipe < 0) error->all(FLERR,"Could not find fix latte compute ID");
    if (modify->compute[ipe]->peatomflag == 0)
      error->all(FLERR,"Fix latte compute ID does not compute pe/atom");
  }

  // initializations

  nmax = 0;
  qpotential = NULL;
  flatte = NULL;

  latte_energy = 0.0;
}

/* ---------------------------------------------------------------------- */

FixLatte::~FixLatte()
{
  delete [] id_pe;
  memory->destroy(qpotential);
  memory->destroy(flatte);
}

/* ---------------------------------------------------------------------- */

int FixLatte::setmask()
{
  int mask = 0;
  //mask |= INITIAL_INTEGRATE;
  //mask |= FINAL_INTEGRATE;
  mask |= PRE_REVERSE;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLatte::init()
{
  // error checks

  if (domain->dimension == 2)
    error->all(FLERR,"Fix latte requires 3d problem");

  if (coulomb) {
    if (atom->q_flag == 0 || force->pair == NULL || force->kspace == NULL)
      error->all(FLERR,"Fix latte cannot compute Coulomb potential");

    int ipe = modify->find_compute(id_pe);
    if (ipe < 0) error->all(FLERR,"Could not find fix latte compute ID");
    c_pe = modify->compute[ipe];
  }

  // must be fully periodic or fully non-periodic

  if (domain->nonperiodic == 0) pbcflag = 1;
  else if (!domain->xperiodic && !domain->yperiodic && !domain->zperiodic)
    pbcflag = 0;
  else error->all(FLERR,"Fix latte requires 3d simulation");

  // create qpotential & flatte if needed
  // for now, assume nlocal will never change

  if (coulomb && qpotential == NULL) {
    memory->create(qpotential,atom->nlocal,"latte:qpotential");
    memory->create(flatte,atom->nlocal,3,"latte:flatte");
  }

  /*
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
  */

  /*
  // need a full neighbor list
  // could we use a half list?
  // perpetual list, built whenever re-neighboring occurs

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  */
}

/* ---------------------------------------------------------------------- */

void FixLatte::init_list(int /*id*/, NeighList * /*ptr*/)
{
  // list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixLatte::setup(int vflag)
{
  newsystem = 1;
  post_force(vflag);
  newsystem = 0;
}

/* ---------------------------------------------------------------------- */

void FixLatte::min_setup(int vflag)
{
  newsystem = 1;
  post_force(vflag);
  newsystem = 0;
}

/* ---------------------------------------------------------------------- */

void FixLatte::setup_pre_reverse(int eflag, int vflag)
{
  pre_reverse(eflag,vflag);
}

/* ----------------------------------------------------------------------
   integrate electronic degrees of freedom
------------------------------------------------------------------------- */

void FixLatte::initial_integrate(int /*vflag*/) {}

/* ----------------------------------------------------------------------
   store eflag, so can use it in post_force to tally per-atom energies
------------------------------------------------------------------------- */

void FixLatte::pre_reverse(int eflag, int /*vflag*/)
{
  eflag_caller = eflag;
}

/* ---------------------------------------------------------------------- */

void FixLatte::post_force(int vflag)
{
  int eflag = eflag_caller;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = eflag_global = vflag_global = eflag_atom = vflag_atom = 0;

  // compute Coulombic potential = pe[i]/q[i]
  // invoke compute pe/atom
  // wrap with clear/add and trigger pe/atom calculation every step

  if (coulomb) {
    modify->clearstep_compute();

    if (!(c_pe->invoked_flag & INVOKED_PERATOM)) {
      c_pe->compute_peratom();
      c_pe->invoked_flag |= INVOKED_PERATOM;
    }

    modify->addstep_compute(update->ntimestep+1);

    double *pe = c_pe->vector_atom;
    double *q = atom->q;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (q[i]) qpotential[i] = pe[i]/q[i];
      else qpotential[i] = 0.0;
  }

  // hardwire these unsupported flags for now

  int coulombflag = 0;
  // pe_peratom = 0;
  // virial_global = 1;              // set via vflag_global at some point
  // virial_peratom = 0;
  neighflag = 0;

  // set flags used by LATTE
  // NOTE: LATTE does not compute per-atom energies or virials

  int flags[6];

  flags[0] = pbcflag;         // 1 for fully periodic, 0 for fully non-periodic
  flags[1] = coulombflag;     // 1 for LAMMPS computes Coulombics, 0 for LATTE
  flags[2] = eflag_atom;      // 1 to return per-atom energies, 0 for no
  flags[3] = vflag_global && thermo_virial;    // 1 to return global/per-atom
  flags[4] = vflag_atom && thermo_virial;      //   virial, 0 for no
  flags[5] = neighflag;       // 1 to pass neighbor list to LATTE, 0 for no

  // setup LATTE arguments

  int natoms = atom->nlocal;
  double *coords = &atom->x[0][0];
  int *type = atom->type;
  int ntypes = atom->ntypes;
  double *mass = &atom->mass[1];
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *forces;
  bool latteerror = 0;
  if (coulomb) forces = &flatte[0][0];
  else forces = &atom->f[0][0];
  int maxiter = -1;

  latte(flags,&natoms,coords,type,&ntypes,mass,boxlo,boxhi,&domain->xy,
        &domain->xz,&domain->yz,forces,&maxiter,&latte_energy,
        &atom->v[0][0],&update->dt,virial,&newsystem,&latteerror);

  if (latteerror) error->all(FLERR,"Internal LATTE problem");

  // sum LATTE forces to LAMMPS forces
  // e.g. LAMMPS may compute Coulombics at some point

  if (coulomb) {
    double **f = atom->f;
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++) {
      f[i][0] += flatte[i][0];
      f[i][1] += flatte[i][1];
      f[i][2] += flatte[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixLatte::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   integrate electronic degrees of freedom
------------------------------------------------------------------------- */

void FixLatte::final_integrate() {}

/* ---------------------------------------------------------------------- */

void FixLatte::reset_dt()
{
  //dtv = update->dt;
  //dtf = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
   DFTB energy from LATTE
------------------------------------------------------------------------- */

double FixLatte::compute_scalar()
{
  return latte_energy;
}

/* ----------------------------------------------------------------------
   memory usage of local arrays
------------------------------------------------------------------------- */

double FixLatte::memory_usage()
{
  double bytes = 0.0;
  if (coulomb) bytes += nmax * sizeof(double);
  if (coulomb) bytes += nmax*3 * sizeof(double);
  return bytes;
}
