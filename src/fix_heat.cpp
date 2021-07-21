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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "fix_heat.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixHeat::FixHeat(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg),
idregion(nullptr), hstr(nullptr), vheat(nullptr), vscale(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal fix heat command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix heat command");

  hstr = nullptr;

  if (utils::strmatch(arg[4],"^v_")) {
    hstr = utils::strdup(arg[4]+2);
  } else {
    heat_input = utils::numeric(FLERR,arg[4],false,lmp);
    hstyle = CONSTANT;
  }

  // optional args

  iregion = -1;

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix heat command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix heat does not exist");
      idregion = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix heat command");
  }

  scale = 1.0;

  maxatom = 0;
}

/* ---------------------------------------------------------------------- */

FixHeat::~FixHeat()
{
  delete [] hstr;
  delete [] idregion;
  memory->destroy(vheat);
  memory->destroy(vscale);
}

/* ---------------------------------------------------------------------- */

int FixHeat::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHeat::init()
{
  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix heat does not exist");
  }

  // check variable

  if (hstr) {
    hvar = input->variable->find(hstr);
    if (hvar < 0)
      error->all(FLERR,"Variable name for fix heat does not exist");
    if (input->variable->equalstyle(hvar)) hstyle = EQUAL;
    else if (input->variable->atomstyle(hvar)) hstyle = ATOM;
    else error->all(FLERR,"Variable for fix heat is invalid style");
  }

  // check for rigid bodies in region (done here for performance reasons)
  if (iregion >= 0 && modify->check_rigid_region_overlap(groupbit,domain->regions[iregion]))
    error->warning(FLERR,"Cannot apply fix heat to atoms in rigid bodies");

  // cannot have 0 atoms in group

  if (group->count(igroup) == 0)
    error->all(FLERR,"Fix heat group has no atoms");
  masstotal = group->mass(igroup);
  if (masstotal <= 0.0)
    error->all(FLERR,"Fix heat group has invalid mass");
}

/* ---------------------------------------------------------------------- */

void FixHeat::end_of_step()
{
  int i;
  double heat,ke,massone;
  double vsub[3],vcm[3];

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;

  // reallocate per-atom arrays if necessary

  if (hstyle == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(vheat);
    memory->destroy(vscale);
    memory->create(vheat,maxatom,"heat:vheat");
    memory->create(vscale,maxatom,"heat:vscale");
  }

  // evaluate variable

  if (hstyle != CONSTANT) {
    modify->clearstep_compute();
    if (hstyle == EQUAL) heat_input = input->variable->compute_equal(hvar);
    else input->variable->compute_atom(hvar,igroup,vheat,1,0);
    modify->addstep_compute(update->ntimestep + nevery);
  }

  // vcm = center-of-mass velocity of scaled atoms

  if (iregion < 0) {
    ke = group->ke(igroup)*force->ftm2v;
    group->vcm(igroup,masstotal,vcm);
  } else {
    masstotal = group->mass(igroup,iregion);
    if (masstotal == 0.0) error->all(FLERR,"Fix heat group has no atoms");
    ke = group->ke(igroup,iregion)*force->ftm2v;
    group->vcm(igroup,masstotal,vcm,iregion);
  }
  double vcmsq = vcm[0]*vcm[0] + vcm[1]*vcm[1] + vcm[2]*vcm[2];

  // add heat via scale factor on velocities for CONSTANT and EQUAL cases
  // scale = velocity scale factor to accomplish eflux change in energy
  // vsub = velocity subtracted from each atom to preserve momentum
  // overall KE cannot go negative

  Region *region = nullptr;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  if (hstyle != ATOM) {
    heat = heat_input*nevery*update->dt*force->ftm2v;
    double escale =
      (ke + heat - 0.5*vcmsq*masstotal)/(ke - 0.5*vcmsq*masstotal);
    if (escale < 0.0) error->all(FLERR,"Fix heat kinetic energy went negative");
    scale = sqrt(escale);
    vsub[0] = (scale-1.0) * vcm[0];
    vsub[1] = (scale-1.0) * vcm[1];
    vsub[2] = (scale-1.0) * vcm[2];

    if (iregion < 0) {
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          v[i][0] = scale*v[i][0] - vsub[0];
          v[i][1] = scale*v[i][1] - vsub[1];
          v[i][2] = scale*v[i][2] - vsub[2];
        }
    } else {
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
          v[i][0] = scale*v[i][0] - vsub[0];
          v[i][1] = scale*v[i][1] - vsub[1];
          v[i][2] = scale*v[i][2] - vsub[2];
        }
    }

  // add heat via per-atom scale factor on velocities for ATOM case
  // vscale = velocity scale factor to accomplish eflux change in energy
  // vsub = velocity subtracted from each atom to preserve momentum
  // KE of an atom cannot go negative

  } else {
    vsub[0] = vsub[1] = vsub[2] = 0.0;
    if (iregion < 0) {
      for (i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          heat = vheat[i]*nevery*update->dt*force->ftm2v;
          vscale[i] =
            (ke + heat - 0.5*vcmsq*masstotal)/(ke - 0.5*vcmsq*masstotal);
          if (vscale[i] < 0.0)
            error->all(FLERR,
                       "Fix heat kinetic energy of an atom went negative");
          scale = sqrt(vscale[i]);
          if (rmass) massone = rmass[i];
          else massone = mass[type[i]];
          vsub[0] += (scale-1.0) * v[i][0]*massone;
          vsub[1] += (scale-1.0) * v[i][1]*massone;
          vsub[2] += (scale-1.0) * v[i][2]*massone;
        }
      }

      vsub[0] /= masstotal;
      vsub[1] /= masstotal;
      vsub[2] /= masstotal;

      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          scale = sqrt(vscale[i]);
          v[i][0] = scale*v[i][0] - vsub[0];
          v[i][1] = scale*v[i][1] - vsub[1];
          v[i][2] = scale*v[i][2] - vsub[2];
        }

    } else {
      for (i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
          heat = vheat[i]*nevery*update->dt*force->ftm2v;
          vscale[i] =
            (ke + heat - 0.5*vcmsq*masstotal)/(ke - 0.5*vcmsq*masstotal);
          if (vscale[i] < 0.0)
            error->all(FLERR,
                       "Fix heat kinetic energy of an atom went negative");
          scale = sqrt(vscale[i]);
          if (rmass) massone = rmass[i];
          else massone = mass[type[i]];
          vsub[0] += (scale-1.0) * v[i][0]*massone;
          vsub[1] += (scale-1.0) * v[i][1]*massone;
          vsub[2] += (scale-1.0) * v[i][2]*massone;
        }
      }

      vsub[0] /= masstotal;
      vsub[1] /= masstotal;
      vsub[2] /= masstotal;

      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
          scale = sqrt(vscale[i]);
          v[i][0] = scale*v[i][0] - vsub[0];
          v[i][1] = scale*v[i][1] - vsub[1];
          v[i][2] = scale*v[i][2] - vsub[2];
        }
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixHeat::compute_scalar()
{
  double average_scale = scale;
  if (hstyle == ATOM) {
    if (!vscale) return 1.0;
    double scale_sum = 0.0;
    int ncount = 0;
    int *mask = atom->mask;
    double **x = atom->x;
    int nlocal = atom->nlocal;
    if (iregion < 0) {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          scale_sum += sqrt(vscale[i]);
          ncount++;
        }
      }
    } else {
      Region *region = domain->regions[iregion];
      region->prematch();
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
          scale_sum += sqrt(vscale[i]);
          ncount++;
        }
      }
    }
    double scale_sum_all = 0.0;
    int ncount_all = 0;
    MPI_Allreduce(&scale_sum,&scale_sum_all,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&ncount,&ncount_all,1,MPI_INT,MPI_SUM,world);
    if (ncount_all == 0) average_scale = 0.0;
    else average_scale = scale_sum_all/static_cast<double>(ncount_all);
  }
  return average_scale;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixHeat::memory_usage()
{
  double bytes = 0.0;
  if (hstyle == ATOM) bytes = atom->nmax*2 * sizeof(double);
  return bytes;
}
