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
   Contributing author: Charles Sievers (UC Davis)
------------------------------------------------------------------------- */

#include "fix_numdiff.h"
#include <mpi.h>
#include <memory.h>
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "respa.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNumDiff::FixNumDiff(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), numdiff_forces(NULL), temp_f(NULL), 
  id_pe(NULL)
{
  if (narg < 5) error->all(FLERR,"Illegal fix numdiff command");

  peratom_flag = 1;
  peratom_freq = nevery;
  size_peratom_cols = 3;
  respa_level_support = 1;

  nevery = force->inumeric(FLERR,arg[3]);
  delta = force->numeric(FLERR,arg[4]);
  if (nevery <= 0 || delta <= 0.0) 
    error->all(FLERR,"Illegal fix numdiff command");

  int n = strlen(id) + 6;
  id_pe = new char[n];
  strcpy(id_pe,id);
  strcat(id_pe,"_pe");

  char **newarg = new char*[3];
  newarg[0] = id_pe;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pe";
  modify->add_compute(3,newarg);
  delete [] newarg;

  maxatom = 0;
  numdiff_forces = NULL;
  temp_f = NULL;
  array_atom = NULL;

}

/* ---------------------------------------------------------------------- */

FixNumDiff::~FixNumDiff()
{
  memory->destroy(numdiff_forces);
  memory->destroy(temp_f);

  modify->delete_compute(id_pe);
  delete [] id_pe;
}

/* ---------------------------------------------------------------------- */

int FixNumDiff::setmask()
{
  datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNumDiff::init()
{
  // require consecutive atom IDs

  if (!atom->tag_enable || !atom->tag_consecutive()) 
    error->all(FLERR,"Fix numdiff requires consecutive atom IDs");

  // check for PE compute

  int icompute = modify->find_compute(id_pe);
  if (icompute < 0) error->all(FLERR,"Compute ID for fix numdiff does not exist");
  pe = modify->compute[icompute];

  if (force->pair && force->pair->compute_flag) pair_compute_flag = 1;
  else pair_compute_flag = 0;
  if (force->kspace && force->kspace->compute_flag) kspace_compute_flag = 1;
  else kspace_compute_flag = 0;

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }

  if (atom->nmax > maxatom) {
    memory->destroy(numdiff_forces);
    memory->destroy(temp_f);
    maxatom = atom->nmax;
    memory->create(numdiff_forces,maxatom,3,"numdiff:numdiff_force");
    memory->create(temp_f,maxatom,3,"numdiff:temp_f");
    array_atom = numdiff_forces;
  }
}

/* ---------------------------------------------------------------------- */

void FixNumDiff::post_force(int vflag)
{
  if (update->ntimestep % nevery) return;

  calculate_forces(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNumDiff::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNumDiff::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   create dynamical matrix
------------------------------------------------------------------------- */

void FixNumDiff::calculate_forces(int vflag)
{
  int i,j,ilocal;
  double energy;

  if (atom->nmax > maxatom) {
    memory->destroy(numdiff_forces);
    memory->destroy(temp_f);
    maxatom = atom->nmax;
    memory->create(numdiff_forces,maxatom,3,"numdiff:numdiff_force");
    memory->create(temp_f,maxatom,3,"numdiff:temp_f");
    array_atom = numdiff_forces;
  }

  // store copy of current forces for owned and ghost atoms

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double position;

  for (i = 0; i < nall; i++)
    for (j = 0; j < 3; j++) {
      temp_f[i][j] = f[i][j];
    }

  // initialize numerical forces to zero

  force_clear(numdiff_forces);

  int flag,allflag;
  double denominator = 0.5 / delta;

  int *mask = atom->mask;
  int ntotal = static_cast<tagint> (atom->natoms);
  int dimension = domain->dimension;

  for (tagint m = 1; m <= ntotal; m++) {
    ilocal = atom->map(m);
    flag = 0;
    if ((ilocal >= 0 && ilocal < nlocal) && (mask[ilocal] & groupbit)) flag = 1;
    MPI_Allreduce(&flag,&allflag,1,MPI_INT,MPI_SUM,world);
    if (!allflag) continue;

    for (int idim = 0; idim < dimension; idim++) {
      position = x[ilocal][idim];
      displace_atom(ilocal,idim,1);
      energy = update_energy(vflag);
      if (ilocal >= 0 && ilocal < nlocal)
        numdiff_forces[ilocal][idim] -= energy;

      displace_atom(ilocal,idim,-2);
      energy = update_energy(vflag);
      if (ilocal >= 0 && ilocal < nlocal) {
        numdiff_forces[ilocal][idim] += energy;
        numdiff_forces[ilocal][idim] *= denominator;
      }

      reset_atom_position(ilocal, idim, position);
    }
  }

  // restore original forces for owned and ghost atoms

  for (i = 0; i < nall; i++)
    for (j = 0; j < 3; j++) {
      f[i][j] = temp_f[i][j];
    }
}

/* ----------------------------------------------------------------------
   displace coord of all owned and ghost copies of ilocal
---------------------------------------------------------------------- */

void FixNumDiff::displace_atom(int ilocal, int idim, int magnitude)
{
  if (ilocal < 0) return;

  double **x = atom->x;
  int *sametag = atom->sametag;
  int j = ilocal;
  x[ilocal][idim] += delta*magnitude;

  while (sametag[j] >= 0) {
    j = sametag[j];
    x[j][idim] += delta*magnitude;
  }
}

/* ----------------------------------------------------------------------
   reset coord of all owned and ghost copies of ilocal (moved atom)
---------------------------------------------------------------------- */

void FixNumDiff::reset_atom_position(int ilocal, int idim, double position)
{
  if (ilocal < 0) return;

  double **x = atom->x;
  int *sametag = atom->sametag;
  int j = ilocal;
  x[ilocal][idim] = position;

  while (sametag[j] >= 0) {
    j = sametag[j];
    x[j][idim] = position;
  }
}

/* ----------------------------------------------------------------------
   evaluate potential energy and forces
   same logic as in Verlet
------------------------------------------------------------------------- */

double FixNumDiff::update_energy(int vflag)
{
  force_clear(atom->f);

  int eflag = 1;
  int vvflag = 0;

  if (pair_compute_flag) force->pair->compute(eflag,vvflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vvflag);
    if (force->angle) force->angle->compute(eflag,vvflag);
    if (force->dihedral) force->dihedral->compute(eflag,vvflag);
    if (force->improper) force->improper->compute(eflag,vvflag);
  }

  if (kspace_compute_flag) force->kspace->compute(eflag,vvflag);

  double energy = pe->compute_scalar();
  return energy;
}

/* ----------------------------------------------------------------------
   clear forces needed
------------------------------------------------------------------------- */

void FixNumDiff::force_clear(double **forces)
{
  size_t nbytes = sizeof(double) * atom->nlocal;
  if (force->newton) nbytes += sizeof(double) * atom->nghost;
  if (nbytes) memset(&forces[0][0],0,3*nbytes);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixNumDiff::memory_usage()
{
  bigint bytes = 0.0;
  bytes += 2 * maxatom*3 * sizeof(double);
  return bytes;
}
