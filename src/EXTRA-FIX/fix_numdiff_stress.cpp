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
   Contributing author: Aidan Thompson (Sandia)
------------------------------------------------------------------------- */

#include "fix_numdiff_stress.h"

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

FixNumDiffStress::FixNumDiffStress(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), id_pe(nullptr),
  temp_x(nullptr), temp_f(nullptr)
{
  if (narg < 5) error->all(FLERR,"Illegal fix numdiff/stress command");

  peratom_freq = nevery;
  respa_level_support = 1;
  vector_flag = 1;
  size_vector = NDIR_STRESS;
  extvector = 0;
  
  maxatom = 0;
  
  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  delta = utils::numeric(FLERR,arg[4],false,lmp);
  if (nevery <= 0 || delta <= 0.0)
    error->all(FLERR,"Illegal fix numdiff command");

  std::string cmd = id + std::string("_pe");
  id_pe = utils::strdup(cmd);

  cmd += " all pe";
  modify->add_compute(cmd);

  // perform initial allocation of atom-based arrays
  // zero numdiff_forces since dump may access it on timestep 0
  // zero numdiff_forces since a variable may access it before first run

  reallocate();
  stress_clear();


  // set fixed-point to default = center of cell

  fixedpoint[0] = 0.5*(domain->boxlo[0]+domain->boxhi[0]);
  fixedpoint[1] = 0.5*(domain->boxlo[1]+domain->boxhi[1]);
  fixedpoint[2] = 0.5*(domain->boxlo[2]+domain->boxhi[2]);

  // define the cartesian indices for each strain (Voigt order)

  dirlist[0][0] = 0;
  dirlist[0][1] = 0;
  dirlist[1][0] = 1;
  dirlist[1][1] = 1;
  dirlist[2][0] = 2;
  dirlist[2][1] = 2;
  
  dirlist[3][0] = 1;
  dirlist[3][1] = 2;
  dirlist[4][0] = 0;
  dirlist[4][1] = 2;
  dirlist[5][0] = 0;
  dirlist[5][1] = 1;
}

/* ---------------------------------------------------------------------- */

FixNumDiffStress::~FixNumDiffStress()
{
  memory->destroy(temp_x);
  memory->destroy(temp_f);

  modify->delete_compute(id_pe);
  delete[] id_pe;
}

/* ---------------------------------------------------------------------- */

int FixNumDiffStress::setmask()
{
  datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNumDiffStress::init()
{
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

void FixNumDiffStress::setup(int vflag)
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

void FixNumDiffStress::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNumDiffStress::post_force(int /* vflag */)
{
  if (update->ntimestep % nevery) return;

  calculate_stress();
}

/* ---------------------------------------------------------------------- */

void FixNumDiffStress::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNumDiffStress::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   compute finite difference stress tensor
------------------------------------------------------------------------- */

void FixNumDiffStress::calculate_stress()
{
  double energy;

  // grow arrays if necessary

  if (atom->nlocal + atom->nghost > maxatom) reallocate();

  // store copy of current forces for owned and ghost atoms

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  for (int i = 0; i < nall; i++)
    for (int j = 0; j < 3; j++) {
      temp_x[i][j] = x[i][j];
      temp_f[i][j] = f[i][j];
    }
  
  // initialize numerical stress to zero

  stress_clear();

  // loop over 6 strain directions
  // compute a finite difference force in each dimension

  int flag,allflag;
  double denominator = 0.5 / delta;

  for (int idir = 0; idir < NDIR_STRESS; idir++) {
    displace_atoms(nall, idir, 1.0);
    energy = update_energy();
    stress[idir] += energy;
    displace_atoms(nall, idir,-1.0);
    energy = update_energy();
    stress[idir] -= energy;
    stress[idir] *= denominator;
    restore_atoms(nall, idir);
  }

  // restore original forces for owned and ghost atoms
  
  for (int i = 0; i < nall; i++)
    for (int j = 0; j < 3; j++)
      f[i][j] = temp_f[i][j];
  
}

/* ----------------------------------------------------------------------
   displace position of all owned and ghost atoms
---------------------------------------------------------------------- */

void FixNumDiffStress::displace_atoms(int nall, int idir, double magnitude)
{
  double **x = atom->x;
  int k = dirlist[idir][0]; 
  int l = dirlist[idir][1]; 
  for (int i = 0; i < nall; i++)
    x[i][k] += temp_x[i][k] + delta*magnitude*
      (temp_x[i][l]-fixedpoint[l]);
}

/* ----------------------------------------------------------------------
   restore position of all owned and ghost atoms
---------------------------------------------------------------------- */

void FixNumDiffStress::restore_atoms(int nall, int idir)
{
  double **x = atom->x;
  int j = dirlist[idir][0]; 
  for (int i = 0; i < nall; i++) {
    x[i][j] = temp_x[i][j];
  }
}

/* ----------------------------------------------------------------------
   evaluate potential energy and forces
   same logic as in Verlet
------------------------------------------------------------------------- */

double FixNumDiffStress::update_energy()
{
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

void FixNumDiffStress::stress_clear()
{
  size_t nbytes = sizeof(double) * size_vector;
  memset(&stress[0],0,nbytes);
}

/* ----------------------------------------------------------------------
   return Ith vector value, assume in range of size_vector
------------------------------------------------------------------------- */

double FixNumDiffStress::compute_vector(int i)
{
  return stress[i];
}

/* ----------------------------------------------------------------------
   reallocated local per-atoms arrays
------------------------------------------------------------------------- */

void FixNumDiffStress::reallocate()
{
  memory->destroy(temp_x);
  memory->destroy(temp_f);
  maxatom = atom->nmax;
  memory->create(temp_x,maxatom,3,"numdiff/stress:temp_x");
  memory->create(temp_f,maxatom,3,"numdiff/stress:temp_f");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixNumDiffStress::memory_usage()
{
  double bytes = 0.0;
  bytes += (double)2 * maxatom*3 * sizeof(double);
  return bytes;
}
