/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Aidan Thompson (Sandia)
------------------------------------------------------------------------- */

#include "fix_numdiff_virial.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "improper.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "respa.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNumDiffVirial::FixNumDiffVirial(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), id_pe(nullptr), pe(nullptr), temp_x(nullptr), temp_f(nullptr)
{
  if (narg < 5) error->all(FLERR, "Illegal fix numdiff/virial command");
  if (igroup) error->all(FLERR, "Compute numdiff/virial must use group all");

  peratom_freq = nevery;
  respa_level_support = 1;
  vector_flag = 1;
  size_vector = NDIR_VIRIAL;
  extvector = 0;
  maxatom = 0;

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  delta = utils::numeric(FLERR, arg[4], false, lmp);
  if (nevery <= 0 || delta <= 0.0) error->all(FLERR, "Illegal fix numdiff command");

  std::string cmd = id + std::string("_pe");
  id_pe = utils::strdup(cmd);
  cmd += " all pe";
  modify->add_compute(cmd);

  // perform initial allocation of atom-based arrays
  // zero numdiff_forces since dump may access it on timestep 0
  // zero numdiff_forces since a variable may access it before first run

  reallocate();

  // set fixed-point to default = center of cell

  fixedpoint[0] = 0.5 * (domain->boxlo[0] + domain->boxhi[0]);
  fixedpoint[1] = 0.5 * (domain->boxlo[1] + domain->boxhi[1]);
  fixedpoint[2] = 0.5 * (domain->boxlo[2] + domain->boxhi[2]);

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

FixNumDiffVirial::~FixNumDiffVirial()
{
  memory->destroy(temp_x);
  memory->destroy(temp_f);

  modify->delete_compute(id_pe);
  delete[] id_pe;
}

/* ---------------------------------------------------------------------- */

int FixNumDiffVirial::setmask()
{
  datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNumDiffVirial::init()
{
  // check for PE compute

  pe = modify->get_compute_by_id(id_pe);
  if (!pe) error->all(FLERR, "PE compute ID for fix numdiff/virial does not exist");

  if (force->pair && force->pair->compute_flag)
    pair_compute_flag = 1;
  else
    pair_compute_flag = 0;
  if (force->kspace && force->kspace->compute_flag)
    kspace_compute_flag = 1;
  else
    kspace_compute_flag = 0;

  if (utils::strmatch(update->integrate_style, "^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixNumDiffVirial::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet"))
    post_force(vflag);
  else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag, ilevel_respa, 0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixNumDiffVirial::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNumDiffVirial::post_force(int /* vflag */)
{
  if (update->ntimestep % nevery) return;

  calculate_virial();
}

/* ---------------------------------------------------------------------- */

void FixNumDiffVirial::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNumDiffVirial::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   compute finite difference virial stress tensor
------------------------------------------------------------------------- */

void FixNumDiffVirial::calculate_virial()
{
  double energy;

  // grow arrays if necessary

  if (atom->nlocal + atom->nghost > maxatom) reallocate();

  // store copy of current forces for owned and ghost atoms

  double **x = atom->x;
  double **f = atom->f;
  int nall = atom->nlocal + atom->nghost;

  for (int i = 0; i < nall; i++)
    for (int k = 0; k < 3; k++) {
      temp_x[i][k] = x[i][k];
      temp_f[i][k] = f[i][k];
    }

  // loop over 6 strain directions
  // compute a finite difference force in each dimension

  double nktv2p = force->nktv2p;
  double inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);

  double denominator = -0.5 / delta * inv_volume * nktv2p;

  for (int idir = 0; idir < NDIR_VIRIAL; idir++) {
    displace_atoms(nall, idir, 1.0);
    energy = update_energy();
    virial[idir] = energy;
    restore_atoms(nall, idir);
    displace_atoms(nall, idir, -1.0);
    energy = update_energy();
    virial[idir] -= energy;
    virial[idir] *= denominator;
    restore_atoms(nall, idir);
  }

  // recompute energy so all contributions are as before

  energy = update_energy();

  // restore original forces for owned and ghost atoms

  for (int i = 0; i < nall; i++)
    for (int k = 0; k < 3; k++) f[i][k] = temp_f[i][k];
}

/* ----------------------------------------------------------------------
   displace position of all owned and ghost atoms
---------------------------------------------------------------------- */

void FixNumDiffVirial::displace_atoms(int nall, int idir, double magnitude)
{
  double **x = atom->x;
  int k = dirlist[idir][0];
  int l = dirlist[idir][1];
  for (int i = 0; i < nall; i++)
    x[i][k] = temp_x[i][k] + delta * magnitude * (temp_x[i][l] - fixedpoint[l]);
}

/* ----------------------------------------------------------------------
   restore position of all owned and ghost atoms
---------------------------------------------------------------------- */

void FixNumDiffVirial::restore_atoms(int nall, int idir)
{
  double **x = atom->x;
  int k = dirlist[idir][0];
  for (int i = 0; i < nall; i++) { x[i][k] = temp_x[i][k]; }
}

/* ----------------------------------------------------------------------
   evaluate potential energy and forces
   same logic as in Verlet
------------------------------------------------------------------------- */

double FixNumDiffVirial::update_energy()
{
  int eflag = 1;

  if (pair_compute_flag) force->pair->compute(eflag, 0);

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) force->bond->compute(eflag, 0);
    if (force->angle) force->angle->compute(eflag, 0);
    if (force->dihedral) force->dihedral->compute(eflag, 0);
    if (force->improper) force->improper->compute(eflag, 0);
  }

  if (kspace_compute_flag) force->kspace->compute(eflag, 0);

  double energy = pe->compute_scalar();
  return energy;
}

/* ----------------------------------------------------------------------
   return Ith vector value, assume in range of size_vector
------------------------------------------------------------------------- */

double FixNumDiffVirial::compute_vector(int i)
{
  return virial[i];
}

/* ----------------------------------------------------------------------
   reallocated local per-atoms arrays
------------------------------------------------------------------------- */

void FixNumDiffVirial::reallocate()
{
  memory->destroy(temp_x);
  memory->destroy(temp_f);
  maxatom = atom->nmax;
  memory->create(temp_x, maxatom, 3, "numdiff/virial:temp_x");
  memory->create(temp_f, maxatom, 3, "numdiff/virial:temp_f");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixNumDiffVirial::memory_usage()
{
  double bytes = 0.0;
  bytes += (double) 2 * maxatom * 3 * sizeof(double);
  return bytes;
}
