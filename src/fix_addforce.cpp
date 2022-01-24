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

#include "fix_addforce.h"

#include "atom.h"
#include "atom_masks.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixAddForce::FixAddForce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  xstr(nullptr), ystr(nullptr), zstr(nullptr), estr(nullptr),
  idregion(nullptr), sforce(nullptr)
{
  if (narg < 6) error->all(FLERR,"Illegal fix addforce command");

  dynamic_group_allow = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  energy_global_flag = 1;
  virial_global_flag = virial_peratom_flag = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  xstr = ystr = zstr = nullptr;

  if (utils::strmatch(arg[3],"^v_")) {
    xstr = utils::strdup(arg[3]+2);
  } else {
    xvalue = utils::numeric(FLERR,arg[3],false,lmp);
    xstyle = CONSTANT;
  }
  if (utils::strmatch(arg[4],"^v_")) {
    ystr = utils::strdup(arg[4]+2);
  } else {
    yvalue = utils::numeric(FLERR,arg[4],false,lmp);
    ystyle = CONSTANT;
  }
  if (utils::strmatch(arg[5],"^v_")) {
    zstr = utils::strdup(arg[5]+2);
  } else {
    zvalue = utils::numeric(FLERR,arg[5],false,lmp);
    zstyle = CONSTANT;
  }

  // optional args

  nevery = 1;
  iregion = -1;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix addforce command");
      nevery = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nevery <= 0) error->all(FLERR,"Illegal fix addforce command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix addforce command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix addforce does not exist");
      idregion = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"energy") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix addforce command");
      if (utils::strmatch(arg[iarg+1],"^v_")) {
        estr = utils::strdup(arg[iarg+1]+2);
      } else error->all(FLERR,"Illegal fix addforce command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix addforce command");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

  maxatom = 1;
  memory->create(sforce,maxatom,4,"addforce:sforce");
}

/* ---------------------------------------------------------------------- */

FixAddForce::~FixAddForce()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] estr;
  delete [] idregion;
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixAddForce::setmask()
{
  datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAddForce::init()
{
  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix addforce does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for fix addforce is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix addforce does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR,"Variable for fix addforce is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix addforce does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for fix addforce is invalid style");
  }
  if (estr) {
    evar = input->variable->find(estr);
    if (evar < 0)
      error->all(FLERR,"Variable name for fix addforce does not exist");
    if (input->variable->atomstyle(evar)) estyle = ATOM;
    else error->all(FLERR,"Variable for fix addforce is invalid style");
  } else estyle = NONE;

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix addforce does not exist");
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  if (varflag == CONSTANT && estyle != NONE)
    error->all(FLERR,"Cannot use variable energy with "
               "constant force in fix addforce");
  if ((varflag == EQUAL || varflag == ATOM) &&
      update->whichflag == 2 && estyle == NONE)
    error->all(FLERR,"Must use variable energy with fix addforce");

  if (utils::strmatch(update->integrate_style,"^respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixAddForce::setup(int vflag)
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

void FixAddForce::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddForce::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  double v[6];
  int nlocal = atom->nlocal;

  if (update->ntimestep % nevery) return;

  // virial setup

  v_init(vflag);

  if (lmp->kokkos)
    atom->sync_modify(Host, (unsigned int) (F_MASK | MASK_MASK),
                      (unsigned int) F_MASK);

  // update region if necessary

  Region *region = nullptr;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // reallocate sforce array if necessary

  if ((varflag == ATOM || estyle == ATOM) && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce,maxatom,4,"addforce:sforce");
  }

  // foriginal[0] = "potential energy" for added force
  // foriginal[123] = force on atoms before extra force added

  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;

  // constant force
  // potential energy = - x dot f in unwrapped coords

  if (varflag == CONSTANT) {
    double unwrap[3];
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        domain->unmap(x[i],image[i],unwrap);
        foriginal[0] -= xvalue*unwrap[0] + yvalue*unwrap[1] + zvalue*unwrap[2];
        foriginal[1] += f[i][0];
        foriginal[2] += f[i][1];
        foriginal[3] += f[i][2];
        f[i][0] += xvalue;
        f[i][1] += yvalue;
        f[i][2] += zvalue;
        if (evflag) {
          v[0] = xvalue * unwrap[0];
          v[1] = yvalue * unwrap[1];
          v[2] = zvalue * unwrap[2];
          v[3] = xvalue * unwrap[1];
          v[4] = xvalue * unwrap[2];
          v[5] = yvalue * unwrap[2];
          v_tally(i,v);
        }
      }

  // variable force, wrap with clear/add
  // potential energy = evar if defined, else 0.0
  // wrap with clear/add

  } else {
    double unwrap[3];

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar,igroup,&sforce[0][0],4,0);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM)
      input->variable->compute_atom(yvar,igroup,&sforce[0][1],4,0);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM)
      input->variable->compute_atom(zvar,igroup,&sforce[0][2],4,0);
    if (estyle == ATOM)
      input->variable->compute_atom(evar,igroup,&sforce[0][3],4,0);

    modify->addstep_compute(update->ntimestep + 1);

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        domain->unmap(x[i],image[i],unwrap);
        if (xstyle == ATOM) xvalue = sforce[i][0];
        if (ystyle == ATOM) yvalue = sforce[i][1];
        if (zstyle == ATOM) zvalue = sforce[i][2];

        if (estyle == ATOM) {
          foriginal[0] += sforce[i][3];
        } else {
          if (xstyle) foriginal[0] -= xvalue*unwrap[0];
          if (ystyle) foriginal[0] -= yvalue*unwrap[1];
          if (zstyle) foriginal[0] -= zvalue*unwrap[2];
        }
        foriginal[1] += f[i][0];
        foriginal[2] += f[i][1];
        foriginal[3] += f[i][2];

        if (xstyle) f[i][0] += xvalue;
        if (ystyle) f[i][1] += yvalue;
        if (zstyle) f[i][2] += zvalue;
        if (evflag) {
          v[0] = xstyle ? xvalue*unwrap[0] : 0.0;
          v[1] = ystyle ? yvalue*unwrap[1] : 0.0;
          v[2] = zstyle ? zvalue*unwrap[2] : 0.0;
          v[3] = xstyle ? xvalue*unwrap[1] : 0.0;
          v[4] = xstyle ? xvalue*unwrap[2] : 0.0;
          v[5] = ystyle ? yvalue*unwrap[2] : 0.0;
          v_tally(i,v);
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixAddForce::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddForce::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixAddForce::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[0];
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixAddForce::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n+1];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixAddForce::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = maxatom*4 * sizeof(double);
  return bytes;
}
