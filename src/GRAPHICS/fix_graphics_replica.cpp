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

#include "fix_graphics_replica.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "graphics.h"
#include "group.h"
#include "math_special.h"
#include "memory.h"
#include "universe.h"
#include "update.h"

#include <algorithm>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathSpecial::square;

/* ---------------------------------------------------------------------- */

FixGraphicsReplica::FixGraphicsReplica(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 4) error->universe_all(FLERR, "Too few arguments for fix graphics/replica");

  // parse mandatory arg

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->universe_all(FLERR, "Illegal fix graphics/replica nevery value");
  global_freq = nevery;
  dynamic_group_allow = 1;

  dflag = false;
  aflag = false;
  dradius = 1.0;
  aradius = 1.0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "display") == 0) {
      if (iarg + 2 > narg)
        error->universe_all(FLERR, "Too few arguments for fix graphics/replica display");
      dflag = true;
      dradius = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "average") == 0) {
      if (iarg + 2 > narg)
        error->universe_all(FLERR, "Too few arguments for fix graphics/replica average");
      aflag = true;
      aradius = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->universe_all(FLERR, std::string("Unknown fix graphics/replica keyword: ") + arg[iarg]);
    }
  }

  // require atom map to sort atoms by ID

  if (atom->map_style == Atom::MAP_NONE)
    error->universe_all(FLERR, "Fix graphics/replica requires an atom map, see atom_modify");

  imggroup = fmt::format("GRAPHICS/REPLICA:{}", id);
  group->find_or_create(imggroup);
  numobjs = 0;
}

/* ---------------------------------------------------------------------- */

FixGraphicsReplica::~FixGraphicsReplica()
{
  memory->destroy(imgobjs);
  memory->destroy(imgparms);

  // delete group created by the constructor,
  // unless the Group class has already been deleted at program exit

  if (group) {
    try {
      group->assign(imggroup + " delete");
    } catch (std::exception &e) {
      if (comm->me == 0)
        fprintf(stderr, "Error deleting group %s: %s\n", imggroup.c_str(), e.what());
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixGraphicsReplica::setmask()
{
  return END_OF_STEP | MIN_POST_FORCE;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsReplica::setup(int /*vflag*/)
{
  // must come after any dynamic group is updated for the first time
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixGraphicsReplica::min_post_force(int)
{
  if (update->ntimestep % nevery) return;
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixGraphicsReplica::end_of_step()
{
  const int me = comm->me;
  const int nprocs = comm->nprocs;

  memory->destroy(imgobjs);
  memory->destroy(imgparms);

  // per-atom data of local atoms

  const auto *const *const x = atom->x;
  const auto *const tag = atom->tag;
  const auto *const type = atom->type;
  const auto *const image = atom->image;
  int *const mask = atom->mask;
  const auto nlocal = atom->nlocal;

  // clear assignment of atoms to the custom group on all ranks

  const int repgroup = group->find_or_create(imggroup);
  const int repbit = group->bitmask[repgroup];
  group->assign(imggroup + " clear");

  // create sorted list of the IDs of all fix group atoms on replica 0

  bigint nper = 0;
  std::vector<tagint> taglist;
  if (universe->iworld == 0) {
    std::vector<tagint> tagme;
    for (int i = 0; i < nlocal; ++i)
      if (mask[i] & groupbit) tagme.emplace_back(tag[i]);
    int ngroup = (int) tagme.size();
    std::vector<int> recvcounts(nprocs, 0);
    std::vector<int> displs(nprocs, 0);
    MPI_Allgather(&ngroup, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, world);
    for (int i = 1; i < nprocs; ++i) displs[i] = displs[i - 1] + recvcounts[i - 1];
    taglist.resize(displs[nprocs - 1] + recvcounts[nprocs - 1], 0);
    MPI_Allgatherv(tagme.data(), ngroup, MPI_LMP_TAGINT, taglist.data(), recvcounts.data(),
                   displs.data(), MPI_LMP_TAGINT, world);
    std::sort(taglist.begin(), taglist.end());
    nper = (bigint) taglist.size();
  }
  MPI_Bcast(&nper, 1, MPI_LMP_BIGINT, 0, universe->uworld);

  // determine number of spheres to draw and check for overflow

  bigint numtotal = 0;
  if (dflag) numtotal += nper * universe->nworlds;
  if (aflag) numtotal += nper;
  if (numtotal >= MAXSMALLINT) error->universe_all(FLERR, "Too many graphics objects");
  numobjs = (int) numtotal;

  // broadcast ID list to all replicas and assign matching local atoms to the
  // custom group.  the custom group now contains the same atoms on all replicas,
  // even if a dynamic fix group selects different atoms on different replicas.

  taglist.resize(nper, 0);
  MPI_Bcast(taglist.data(), (int) nper, MPI_LMP_TAGINT, 0, universe->uworld);
  for (const auto id : taglist) {
    const int i = atom->map(id);
    if ((i >= 0) && (i < nlocal)) mask[i] |= repbit;
  }

  // consistency check: all atoms in the ID list must exist on all replicas

  int missing = (group->count(repgroup) != nper) ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &missing, 1, MPI_INT, MPI_MAX, universe->uworld);
  if (missing)
    error->universe_all(FLERR, "Atoms in fix group on replica 0 must exist on all replicas");

  // collect unwrapped positions and atom type data
  // for selected atoms sorted by ID on the root of each replica

  std::vector<double> pos(3 * nper, 0.0);
  std::vector<int> types(nper, 0);
  int n = 0;
  for (const auto id : taglist) {
    int i = atom->map(id);
    if ((i >= 0) && (i < nlocal)) {
      double tmp[3] = {x[i][0], x[i][1], x[i][2]};
      domain->unmap(tmp, image[i]);
      pos[3 * n] = tmp[0];
      pos[3 * n + 1] = tmp[1];
      pos[3 * n + 2] = tmp[2];
      types[n] = type[i];
    }
    ++n;
  }

  if (me == 0) {
    MPI_Reduce(MPI_IN_PLACE, types.data(), nper, MPI_INT, MPI_SUM, 0, world);
    MPI_Reduce(MPI_IN_PLACE, pos.data(), 3 * nper, MPI_DOUBLE, MPI_SUM, 0, world);
  } else {
    MPI_Reduce(types.data(), nullptr, nper, MPI_INT, MPI_SUM, 0, world);
    MPI_Reduce(pos.data(), nullptr, 3 * nper, MPI_DOUBLE, MPI_SUM, 0, world);
  }

  // now we are ready to create the graphics items
  // only universe root creates the objects

  if (universe->me == 0) {
    memory->create(imgobjs, numobjs, "fix_graphics:imgobjs");
    memory->create(imgparms, numobjs, 5, "fix_graphics:imgparms");

    // reset counter for total graphics objects
    int n = 0;

    // first process our own data;
    std::vector<double> avg(pos);
    std::vector<double> var(3 * nper, 0.0);

    if (dflag) {
      for (int i = 0; i < nper; ++i) {
        imgobjs[n] = Graphics::SPHERE;
        imgparms[n][0] = types[i];
        domain->remap(pos.data() + 3 * i);
        imgparms[n][1] = pos[3 * i];
        imgparms[n][2] = pos[3 * i + 1];
        imgparms[n][3] = pos[3 * i + 2];
        imgparms[n][4] = dradius;
        ++n;
      }
    }

    // now get data from other replicas

    for (int j = 1; j < universe->nworlds; ++j) {
      MPI_Recv(pos.data(), 3 * nper, MPI_DOUBLE, universe->root_proc[j], 0, universe->uworld,
               MPI_STATUS_IGNORE);
      for (int i = 0; i < nper; ++i) {
        if (aflag) {
          var[3 * i] += (double) j / ((double) j + 1) * square(pos[3 * i] - avg[3 * i]);
          var[3 * i + 1] += (double) j / ((double) j + 1) * square(pos[3 * i + 1] - avg[3 * i + 1]);
          var[3 * i + 2] += (double) j / ((double) j + 1) * square(pos[3 * i + 2] - avg[3 * i + 2]);
          avg[3 * i] += (pos[3 * i] - avg[3 * i]) / ((double) j + 1.0);
          avg[3 * i + 1] += (pos[3 * i + 1] - avg[3 * i + 1]) / ((double) j + 1.0);
          avg[3 * i + 2] += (pos[3 * i + 2] - avg[3 * i + 2]) / ((double) j + 1.0);
        }
        if (dflag) {
          imgobjs[n] = Graphics::SPHERE;
          imgparms[n][0] = types[i];
          domain->remap(pos.data() + 3 * i);
          imgparms[n][1] = pos[3 * i];
          imgparms[n][2] = pos[3 * i + 1];
          imgparms[n][3] = pos[3 * i + 2];
          imgparms[n][4] = dradius;
          ++n;
        }
      }
    }
    if (aflag) {
      double norm = 1.0 / (double) universe->nworlds;
      for (int i = 0; i < nper; ++i) {
        imgobjs[n] = Graphics::SPHERE;
        imgparms[n][0] = types[i];
        var[3 * i] *= norm;
        var[3 * i + 1] *= norm;
        var[3 * i + 2] *= norm;
        domain->remap(avg.data() + 3 * i);
        imgparms[n][1] = avg[3 * i];
        imgparms[n][2] = avg[3 * i + 1];
        imgparms[n][3] = avg[3 * i + 2];
        imgparms[n][4] =
            (aradius == 0.0) ? 1.0 : aradius * sqrt(var[3 * i] + var[3 * i + 1] + var[3 * i + 2]);
        ++n;
      }
    }
  } else {
    if (me == 0)
      MPI_Send(pos.data(), 3 * nper, MPI_DOUBLE, universe->root_proc[0], 0, universe->uworld);
  }
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image on universe root only
------------------------------------------------------------------------- */

int FixGraphicsReplica::image(int *&objs, double **&parms)
{
  if (universe->me == 0) {
    objs = imgobjs;
    parms = imgparms;
    return numobjs;
  }
  return 0;
}
