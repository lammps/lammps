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

#include "stdlib.h"
#include "string.h"
#include "kspace.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "suffix.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

KSpace::KSpace(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  energy = 0.0;
  virial[0] = virial[1] = virial[2] = virial[3] = virial[4] = virial[5] = 0.0;

  compute_flag = 1;
  order = 5;
  gridflag = 0;
  gewaldflag = 0;
  slabflag = 0;
  differentiation_flag = 0;
  slab_volfactor = 1;
  suffix_flag = Suffix::NONE;

  accuracy_absolute = -1.0;
  two_charge_force = force->qqr2e * 
    (force->qelectron * force->qelectron) / 
    (force->angstrom * force->angstrom);

  maxeatom = maxvatom = 0;
  eatom = NULL;
  vatom = NULL;
}

/* ---------------------------------------------------------------------- */

KSpace::~KSpace()
{
  memory->destroy(eatom);
  memory->destroy(vatom);
}

/* ---------------------------------------------------------------------- */

void KSpace::compute_dummy(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = evflag_atom = eflag_global = vflag_global = 
	 eflag_atom = vflag_atom = 0;
}

/* ----------------------------------------------------------------------
   setup for energy, virial computation
   see integrate::ev_set() for values of eflag (0-3) and vflag (0-6)
------------------------------------------------------------------------- */

void KSpace::ev_setup(int eflag, int vflag)
{
  int i,n;

  evflag = 1;

  eflag_either = eflag;
  eflag_global = eflag % 2;
  eflag_atom = eflag / 2;

  vflag_either = vflag;
  vflag_global = vflag % 4;
  vflag_atom = vflag / 4;
  
  if (eflag_atom || vflag_atom) evflag_atom = 1;
  else evflag_atom = 0;

  // reallocate per-atom arrays if necessary

  if (eflag_atom && atom->nlocal > maxeatom) {
    maxeatom = atom->nmax;
    memory->destroy(eatom);
    memory->create(eatom,maxeatom,"kspace:eatom");
  }
  if (vflag_atom && atom->nlocal > maxvatom) {
    maxvatom = atom->nmax;
    memory->destroy(vatom);
    memory->create(vatom,maxvatom,6,"kspace:vatom");
  }

  // zero accumulators

  if (eflag_global) energy = 0.0;
  if (vflag_global) for (i = 0; i < 6; i++) virial[i] = 0.0;
  if (eflag_atom) {
    n = atom->nlocal;
    for (i = 0; i < n; i++) eatom[i] = 0.0;
  }
  if (vflag_atom) {
    n = atom->nlocal;
    for (i = 0; i < n; i++) {
      vatom[i][0] = 0.0;
      vatom[i][1] = 0.0;
      vatom[i][2] = 0.0;
      vatom[i][3] = 0.0;
      vatom[i][4] = 0.0;
      vatom[i][5] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   modify parameters of the KSpace style 
------------------------------------------------------------------------- */

void KSpace::modify_params(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mesh") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal kspace_modify command");
      nx_pppm = atoi(arg[iarg+1]);
      ny_pppm = atoi(arg[iarg+2]);
      nz_pppm = atoi(arg[iarg+3]);
      if (nx_pppm == 0 && ny_pppm == 0 && nz_pppm == 0) gridflag = 0;
      else gridflag = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg],"order") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      order = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"force") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      accuracy_absolute = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"gewald") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      g_ewald = atof(arg[iarg+1]);
      if (g_ewald == 0.0) gewaldflag = 0;
      else gewaldflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"slab") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      slab_volfactor = atof(arg[iarg+1]);
      iarg += 2;
      if (slab_volfactor <= 1.0)
	error->all(FLERR,"Bad kspace_modify slab parameter");
      if (slab_volfactor < 2.0 && comm->me == 0) 
	error->warning(FLERR,"Kspace_modify slab param < 2.0 may "
		       "cause unphysical behavior");
      slabflag = 1;
    } else if (strcmp(arg[iarg],"compute") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) compute_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) compute_flag = 0;
      else error->all(FLERR,"Illegal kspace_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"diff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      if (strcmp(arg[iarg+1],"ad") == 0) differentiation_flag = 1;
      else if (strcmp(arg[iarg+1],"ik") == 0) differentiation_flag = 0;
      else error->all(FLERR, "Illegal kspace_modify command");
      iarg += 2;
    } else error->all(FLERR,"Illegal kspace_modify command");
  }
}

/* ---------------------------------------------------------------------- */

void *KSpace::extract(const char *str)
{
  if (strcmp(str,"scale") == 0) return (void *) &scale;
  return NULL;
}
