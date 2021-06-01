// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
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

#include <cstring>
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
  Fix(lmp, narg, arg), id_pe(nullptr), numdiff_forces(nullptr),
  temp_x(nullptr), temp_f(nullptr)
{
  if (narg < 5) error->all(FLERR,"Illegal fix numdiff command");

  peratom_flag = 1;
  peratom_freq = nevery;
  size_peratom_cols = 3;
  respa_level_support = 1;

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  delta = utils::numeric(FLERR,arg[4],false,lmp);
  if (nevery <= 0 || delta <= 0.0)
    error->all(FLERR,"Illegal fix numdiff command");

  std::string cmd = id + std::string("_pe");
  id_pe = new char[cmd.size()+1];
  strcpy(id_pe,cmd.c_str());

  cmd += " all pe";
  modify->add_compute(cmd);

  maxatom = 0;

  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR,"Fix numdiff requires an atom map, see atom_modify");

  // perform initial allocation of atom-based arrays
  // zero numdiff_forces since dump may access it on timestep 0
  // zero numdiff_forces since a variable may access it before first run

  reallocate();
  force_clear(numdiff_forces);
}

/* ---------------------------------------------------------------------- */

FixNumDiff::~FixNumDiff()
{
  memory->destroy(numdiff_forces);
  memory->destroy(temp_x);
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

  if (utils::strmatch(update->integrate_style,"^respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixNumDiff::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixNumDiff::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNumDiff::post_force(int /* vflag */)
{
  if (update->ntimestep % nevery) return;

  calculate_forces();
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
   compute finite difference forces
------------------------------------------------------------------------- */

void FixNumDiff::calculate_forces()
{
  int i,j,ilocal;
  double energy;

  // grow arrays if necessary

  if (atom->nlocal + atom->nghost > maxatom) reallocate();

  // store copy of current forces for owned and ghost atoms

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  for (i = 0; i < nall; i++)
    for (j = 0; j < 3; j++) {
      temp_x[i][j] = x[i][j];
      temp_f[i][j] = f[i][j];
    }

  // initialize numerical forces to zero

  force_clear(numdiff_forces);

  // loop over all atoms in system
  // compute a finite difference force in each dimension

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
      displace_atoms(ilocal,idim,1);
      energy = update_energy();
      if (ilocal >= 0 && ilocal < nlocal)
        numdiff_forces[ilocal][idim] -= energy;

      displace_atoms(ilocal,idim,-2);
      energy = update_energy();
      if (ilocal >= 0 && ilocal < nlocal) {
        numdiff_forces[ilocal][idim] += energy;
        numdiff_forces[ilocal][idim] *= denominator;
      }

      restore_atoms(ilocal,idim);
    }
  }

  // restore original forces for owned and ghost atoms

  for (i = 0; i < nall; i++)
    for (j = 0; j < 3; j++) {
      f[i][j] = temp_f[i][j];
    }
}

/* ----------------------------------------------------------------------
   displace position of all owned and ghost copies of ilocal
---------------------------------------------------------------------- */

void FixNumDiff::displace_atoms(int ilocal, int idim, int magnitude)
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
   restore position of all owned and ghost copies of ilocal
---------------------------------------------------------------------- */

void FixNumDiff::restore_atoms(int ilocal, int idim)
{
  if (ilocal < 0) return;

  double **x = atom->x;
  int *sametag = atom->sametag;
  int j = ilocal;
  x[ilocal][idim] = temp_x[ilocal][idim];

  while (sametag[j] >= 0) {
    j = sametag[j];
    x[j][idim] = temp_x[j][idim];
  }
}

/* ----------------------------------------------------------------------
   evaluate potential energy and forces
   same logic as in Verlet
------------------------------------------------------------------------- */

double FixNumDiff::update_energy()
{
  force_clear(atom->f);

  int eflag = 1;

  if (pair_compute_flag) force->pair->compute(eflag,0);

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) force->bond->compute(eflag,0);
    if (force->angle) force->angle->compute(eflag,0);
    if (force->dihedral) force->dihedral->compute(eflag,0);
    if (force->improper) force->improper->compute(eflag,0);
  }

  if (kspace_compute_flag) force->kspace->compute(eflag,0);

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
   reallocated local per-atoms arrays
------------------------------------------------------------------------- */

void FixNumDiff::reallocate()
{
  memory->destroy(numdiff_forces);
  memory->destroy(temp_x);
  memory->destroy(temp_f);
  maxatom = atom->nmax;
  memory->create(numdiff_forces,maxatom,3,"numdiff:numdiff_force");
  memory->create(temp_x,maxatom,3,"numdiff:temp_x");
  memory->create(temp_f,maxatom,3,"numdiff:temp_f");
  array_atom = numdiff_forces;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixNumDiff::memory_usage()
{
  double bytes = 0.0;
  bytes += (double)3 * maxatom*3 * sizeof(double);
  return bytes;
}
