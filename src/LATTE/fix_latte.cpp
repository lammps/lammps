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
   Contributing author: Christian Negre (LANL)
------------------------------------------------------------------------- */

#include "fix_latte.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>

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

/* ---------------------------------------------------------------------- */

FixLatte::FixLatte(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 3) utils::missing_cmd_args(FLERR, "fix latte", error);

  if (strcmp(update->unit_style,"metal") != 0)
    error->all(FLERR,"Must use units metal with fix latte command");

  if (comm->nprocs != 1)
    error->all(FLERR,"Fix latte currently runs only in serial");

  if (LATTE_ABIVERSION != latte_abiversion())
    error->all(FLERR,"LAMMPS is linked against incompatible LATTE library");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  virial_global_flag = 1;
  thermo_energy = thermo_virial = 1;

  // process optional args

  coulomb = 0;
  id_pe = nullptr;
  exclude = 0;
  id_exclude = nullptr;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"coulomb") == 0) {
      if (iarg+2 > narg)
        utils::missing_cmd_args(FLERR, "fix latte coulomb", error);
      coulomb = 1;
      error->all(FLERR,"Fix latte does not yet support LAMMPS calculation of a Coulomb potential");
      delete[] id_pe;
      id_pe = utils::strdup(arg[iarg+1]);
      c_pe = modify->get_compute_by_id(id_pe);
      if (!c_pe) error->all(FLERR,"Could not find fix latte compute ID {}", id_pe);
      if (c_pe->peatomflag == 0) error->all(FLERR,"Fix latte compute ID does not compute pe/atom");
      iarg += 2;

    } else if (strcmp(arg[iarg],"exclude") == 0) {
      if (iarg+2 > narg)
        utils::missing_cmd_args(FLERR, "fix latte exclude", error);
      exclude = 1;
      delete[] id_exclude;
      id_exclude = utils::strdup(arg[iarg+1]);
      iarg += 2;

    } else
      error->all(FLERR, "Unknown fix latte keyword: {}", arg[iarg]);
  }

  // initializations

  nmax = 0;
  qpotential = nullptr;
  flatte = nullptr;

  latte_energy = 0.0;
}

/* ---------------------------------------------------------------------- */

FixLatte::~FixLatte()
{
  delete[] id_pe;
  delete[] id_exclude;
  memory->destroy(qpotential);
  memory->destroy(flatte);
}

/* ---------------------------------------------------------------------- */

int FixLatte::setmask()
{
  int mask = 0;
  mask |= PRE_REVERSE;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLatte::init()
{
  // error checks

  if (domain->dimension == 2)
    error->all(FLERR,"Fix latte requires 3d problem");

  if (coulomb) {
    if (atom->q_flag == 0 || force->pair == nullptr || force->kspace == nullptr)
      error->all(FLERR,"Fix latte cannot compute Coulomb potential");
    c_pe = modify->get_compute_by_id(id_pe);
    if (!c_pe) error->all(FLERR,"Fix latte could not find Coulomb compute ID {}",id_pe);
  }

  // must be fully periodic or fully non-periodic

  if (domain->nonperiodic == 0) pbcflag = 1;
  else if (!domain->xperiodic && !domain->yperiodic && !domain->zperiodic)
    pbcflag = 0;
  else error->all(FLERR,"Fix latte requires 3d simulation");

  // create qpotential & flatte if needed
  // for now, assume nlocal will never change

  if (coulomb && qpotential == nullptr) {
    memory->create(qpotential,atom->nlocal,"latte:qpotential");
    memory->create(flatte,atom->nlocal,3,"latte:flatte");
  }

  // extract pointer to exclusion_group variable from id_exclude
  // exclusion_group is index of a group the Fix defines

  if (exclude) {
    Fix *f_exclude = modify->get_fix_by_id(id_exclude);
    if (!f_exclude) error->all(FLERR,"Fix latte could not find exclude fix ID {}", id_exclude);
    int dim;
    exclusion_group_ptr = (int *) f_exclude->extract("exclusion_group", dim);
    if (!exclusion_group_ptr || dim != 0)
      error->all(FLERR,"Fix latte could not query exclude_group of fix ID {}", id_exclude);
  }
}

/* ---------------------------------------------------------------------- */

void FixLatte::init_list(int /*id*/, NeighList * /*ptr*/)
{
  // list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixLatte::setup(int vflag)
{
  natoms_last = -1;
  setupflag = 1;
  post_force(vflag);
  setupflag = 0;
}

/* ---------------------------------------------------------------------- */

void FixLatte::min_setup(int vflag)
{
  natoms_last = -1;
  setupflag = 1;
  post_force(vflag);
  setupflag = 0;
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
  ev_init(eflag,vflag);

  // compute Coulombic potential = pe[i]/q[i]
  // invoke compute pe/atom
  // wrap with clear/add and trigger pe/atom calculation every step

  if (coulomb) {
    modify->clearstep_compute();

    if (!(c_pe->invoked_flag & Compute::INVOKED_PERATOM)) {
      c_pe->compute_peratom();
      c_pe->invoked_flag |= Compute::INVOKED_PERATOM;
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
  neighflag = 0;

  // set flags used by LATTE
  // note that LATTE does not compute per-atom energies or virials

  flags_latte[0] = pbcflag;   // 1 for fully periodic, 0 for fully non-periodic
  flags_latte[1] = coulombflag;  // 1 for LAMMPS computes Coulombics, 0 for LATTE
  flags_latte[2] = eflag_atom;   // 1 to return per-atom energies, 0 for no
  flags_latte[3] = vflag_global && thermo_virial; // 1 to return global/per-atom
  flags_latte[4] = vflag_atom && thermo_virial;   //   virial, 0 for no
  flags_latte[5] = neighflag;    // 1 to pass neighbor list to LATTE, 0 for no

  // newsystem flag determines whether LATTE treats snapshot
  //   as new system (more work) or increment to last system
  // if setup or atom count changed then newsystem = 1
  // else newsystem = 0

  if (setupflag || atom->natoms != natoms_last) newsystem = 1;
  else newsystem = 0;

  // setup arguments for latte() function and invoke it
  // either for all atoms or excluding some atoms
  // in latter case, need to construct reduced-size per-atom vectors/arrays

  if (!exclude) latte_wrapper_all();
  else {
    int anyexclude = 0;

    int exclusion_group = *exclusion_group_ptr;
    if (exclusion_group) {
      int excludebit = group->bitmask[exclusion_group];

      int *mask = atom->mask;
      int nlocal = atom->nlocal;

      for (int i = 0; i < nlocal; i++)
        if (mask[i] & excludebit) anyexclude = 1;
    }

    if (!anyexclude) latte_wrapper_all();
    else latte_wrapper_exclude();
  }

  newsystem = 0;
  natoms_last = atom->natoms;

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

/* ----------------------------------------------------------------------
   invoke LATTE on all LAMMPS atoms
------------------------------------------------------------------------- */

void FixLatte::latte_wrapper_all()
{
  int natoms = atom->nlocal;
  double *coords = &atom->x[0][0];
  int *types = atom->type;
  int ntypes = atom->ntypes;
  double *mass = &atom->mass[1];
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *forces;
  bool latteerror = false;
  if (coulomb) forces = &flatte[0][0];
  else forces = &atom->f[0][0];
  int maxiter = -1;

  latte(flags_latte,&natoms,coords,types,&ntypes,mass,boxlo,boxhi,
        &domain->xy,&domain->xz,&domain->yz,forces,&maxiter,&latte_energy,
        &atom->v[0][0],&update->dt,virial,&newsystem,&latteerror);

  if (latteerror) error->all(FLERR,"Internal LATTE problem");
}

/* ----------------------------------------------------------------------
   invoke LATTE on only LAMMPS atoms not in exclude group
------------------------------------------------------------------------- */

void FixLatte::latte_wrapper_exclude()
{
  int m;

  int exclusion_group = *exclusion_group_ptr;
  int excludebit = group->bitmask[exclusion_group];

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // nlatte = number of non-excluded atoms to pass to LATTE

  int nlatte = 0;
  for (int i = 0; i < nlocal; i++)
    if (!(mask[i] & excludebit)) nlatte++;

  // created compressed type vector and coords array

  int *typeinclude;
  double **xinclude,**finclude;
  memory->create(typeinclude,nlatte,"latte:typeinclude");
  memory->create(xinclude,nlatte,3,"latte:xinclude");
  memory->create(finclude,nlatte,3,"latte:finclude");

  double *coords = &xinclude[0][0];
  int *types = typeinclude;
  double *forces = &finclude[0][0];

  double **x = atom->x;
  int *type = atom->type;

  nlatte = 0;
  m = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & excludebit) continue;
    types[nlatte] = type[i];
    nlatte++;
    coords[m+0] = x[i][0];
    coords[m+1] = x[i][1];
    coords[m+2] = x[i][2];
    m += 3;
  }

  int ntypes = atom->ntypes;
  double *mass = &atom->mass[1];
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  bool latteerror = false;
  int maxiter = -1;

  latte(flags_latte,&nlatte,coords,types,&ntypes,mass,boxlo,boxhi,
        &domain->xy,&domain->xz,&domain->yz,forces,&maxiter,&latte_energy,
        &atom->v[0][0],&update->dt,virial,&newsystem,&latteerror);

  if (latteerror) error->all(FLERR,"Internal LATTE problem");

  // expand compressed forces array

  double **f = atom->f;

  m = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & excludebit) continue;
    f[i][0] = forces[m+0];
    f[i][1] = forces[m+1];
    f[i][2] = forces[m+2];
    m += 3;
  }

  memory->destroy(typeinclude);
  memory->destroy(xinclude);
  memory->destroy(finclude);
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
  if (coulomb) bytes += (double)nmax * sizeof(double);
  if (coulomb) bytes += (double)nmax*3 * sizeof(double);
  return bytes;
}
