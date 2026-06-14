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

#include "fix_graphics_lines.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_store_atom.h"
#include "graphics.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGraphicsLines::FixGraphicsLines(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), imgobjs(nullptr), imgparms(nullptr), id_cprop(nullptr), cprop(nullptr),
    id_fave(nullptr), fave(nullptr), id_fstore(nullptr), fstore(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix graphics/lines", error);

  // parse mandatory args

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  nrepeat = utils::inumeric(FLERR, arg[4], false, lmp);
  nfreq = utils::inumeric(FLERR, arg[5], false, lmp);
  nlength = utils::inumeric(FLERR, arg[6], false, lmp);

  if (nevery <= 0) error->all(FLERR, 3, "Illegal fix graphics/lines nevery value: {}", nevery);
  if (nrepeat <= 0) error->all(FLERR, 4, "Illegal fix graphics/lines nrepeat value: {}", nrepeat);
  if (nfreq <= 0) error->all(FLERR, 5, "Illegal fix graphics/lines nfreq value: {}", nfreq);
  if (nlength <= 0) error->all(FLERR, 6, "Illegal fix graphics/lines nlength value: {}", nlength);
  if ((nfreq % nevery) || (nrepeat * nevery > nfreq))
    error->all(FLERR, Error::NOPOINTER,
               "Inconsistent fix graphics/lines nevery/nrepeat/nfreq values");

  // fix settings
  scalar_flag = 1;
  extscalar = 0;
  time_depend = 1;
  restart_global = 1;
  dynamic_group_allow = 0;    // there is no clean way to use dynamic groups for this
  firstflag = 1;

  nvalues = ivalue = numobjs = 0;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsLines::post_constructor()
{
  const auto *gname = group->names[igroup];
  id_cprop = utils::strdup(id + std::string("_COMPUTE_PROPERTY_ATOM"));
  cprop = modify->add_compute(fmt::format("{0} {1} property/atom xu yu zu", id_cprop, gname));
  if (!cprop)
    error->all(FLERR, Error::NOLASTLINE,
               "Error creating internal unwrapped coodinate compute for fix graphics/lines");
  id_fave = utils::strdup(id + std::string("_FIX_AVE_ATOM"));
  fave = modify->add_fix(fmt::format("{0} {1} ave/atom {2} {3} {4} c_{5}[1] c_{5}[2] c_{5}[3]",
                                     id_fave, gname, nevery, nrepeat, nfreq, id_cprop));
  if (!fave)
    error->all(FLERR, Error::NOLASTLINE,
               "Error creating internal position averaging for fix graphics/lines");
  id_fstore = utils::strdup(id + std::string("_FIX_STORE_ATOM"));
  auto cmd = fmt::format("{0} {1} STORE/ATOM {2} 3 0 1", id_fstore, gname, nlength);
  fstore = dynamic_cast<FixStoreAtom *>(modify->add_fix(cmd));
  if (!fstore)
    error->all(FLERR, Error::NOLASTLINE, "Error creating internal storage for fix graphics/lines");

  // turn off automatic end_of_step() processing for fix ave/atom.
  // We call it manually instead to ensure it is called *before* we access its data
  int ifave = modify->find_fix(id_fave);
  if (ifave < 0) error->all(FLERR, Error::NOLASTLINE, "Internal fix information corrupted");
  modify->fmask[ifave] &= ~END_OF_STEP;
}

/* ---------------------------------------------------------------------- */

FixGraphicsLines::~FixGraphicsLines()
{
  if (modify->ncompute && modify->get_compute_by_id(id_cprop)) modify->delete_compute(id_cprop);
  if (modify->nfix && modify->get_fix_by_id(id_fave)) modify->delete_fix(id_fave);
  if (modify->nfix && modify->get_fix_by_id(id_fstore)) modify->delete_fix(id_fstore);
  delete[] id_cprop;
  delete[] id_fave;
  delete[] id_fstore;
  memory->destroy(imgobjs);
  memory->destroy(imgparms);
}

/* ---------------------------------------------------------------------- */

int FixGraphicsLines::setmask()
{
  return END_OF_STEP;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsLines::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixGraphicsLines::write_restart(FILE *fp)
{
  if (comm->me == 0) {
    int size = 6 * sizeof(double);
    double list[6];
    list[0] = nevery;
    list[1] = nrepeat;
    list[2] = nfreq;
    list[3] = nlength;
    list[4] = nvalues;
    list[5] = ivalue;
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(list, sizeof(double), 6, fp);
  }
}

/* ---------------------------------------------------------------------- */

void FixGraphicsLines::restart(char *buf)
{
  auto *list = (double *) buf;
  if ((list[0] != nevery) || (list[1] != nrepeat) || (list[2] != nfreq) || (list[3] != nlength))
    error->all(FLERR, Error::NOLASTLINE, "Cannot restart fix graphics/lines: settings don't match");

  nvalues = (int) list[4];
  ivalue = (int) list[5];
}

/* ---------------------------------------------------------------------- */

void FixGraphicsLines::end_of_step()
{
  const auto *const *const x = atom->x;
  const auto *const *const xave = fave->array_atom;
  auto *const *const *const xstore = fstore->tstore;
  const auto *const mask = atom->mask;
  const auto *const type = atom->type;
  const auto &prd_half = domain->prd_half;
  const auto nlocal = atom->nlocal;

  int n = 0;

  if (restart_reset == 0) {
    // when we are not restarting, we have to call end_of_step() for fix ave/atom explicitly
    fave->end_of_step();
    if (firstflag) {
      firstflag = 0;
      return;
    }
    // we only continue when one block of averaging is complete
    if (update->ntimestep % nfreq) return;

    // update history data storage
    n = 0;
    for (int i = 0; i < nlocal; ++i) {
      if (mask[i] & groupbit) {
        xstore[i][ivalue][0] = xave[i][0];
        xstore[i][ivalue][1] = xave[i][1];
        xstore[i][ivalue][2] = xave[i][2];
        domain->remap(xstore[i][ivalue]);
        ++n;
      }
    }
    if (nvalues < nlength) ++nvalues;
    ++ivalue;
    if (ivalue == nlength) ivalue = 0;
  } else {
    firstflag = 0;
    // when we are restarting, we also get called from setup()
    // we only need to know the number of local atoms in the fix group
    // then we use the restarted data to fill the image data arrays
    n = 0;
    for (int i = 0; i < nlocal; ++i) {
      if (mask[i] & groupbit) ++n;
    }
    restart_reset = 0;
    // we do not have any stored data, so we cannot generate any graphics
    if (nvalues == 0) return;
  }

  // allocate storage for the largest possible number of graphics objects

  memory->destroy(imgobjs);
  memory->destroy(imgparms);
  numobjs = n * nvalues;
  memory->create(imgobjs, numobjs, "fix_graphics_lines:imgobjs");
  memory->create(imgparms, numobjs, 8, "fix_graphics_lines:imgparms");

  n = 0;
  double zoffs = 0.0;
  if (domain->dimension == 2) zoffs = prd_half[2];
  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      int j1, j2 = (ivalue - 1 + nvalues) % nvalues;

      // draw cylinder from current position to first averaged history point
      // skip if there is a jump across periodic boundaries in any direction
      if ((fabs(x[i][0] - xstore[i][j2][0]) < prd_half[0]) &&
          (fabs(x[i][1] - xstore[i][j2][1]) < prd_half[1]) &&
          (fabs(x[i][2] - xstore[i][j2][2]) < prd_half[2])) {
        imgobjs[n] = Graphics::CYLINDER;
        imgparms[n][0] = type[i];
        imgparms[n][1] = xstore[i][j2][0];
        imgparms[n][2] = xstore[i][j2][1];
        imgparms[n][3] = xstore[i][j2][2] + zoffs;
        imgparms[n][4] = x[i][0];
        imgparms[n][5] = x[i][1];
        imgparms[n][6] = x[i][2] + zoffs;
        imgparms[n][7] = 0.0;
        ++n;
      }
      // draw cylinders for the available stored averaged history
      // skip if there is a jump across periodic boundaries in any direction
      for (int jj = 1; jj < nvalues; ++jj) {
        j1 = (ivalue + jj - 1) % nvalues;
        j2 = (ivalue + jj + 0) % nvalues;
        if ((fabs(xstore[i][j1][0] - xstore[i][j2][0]) < prd_half[0]) &&
            (fabs(xstore[i][j1][1] - xstore[i][j2][1]) < prd_half[1]) &&
            (fabs(xstore[i][j1][2] - xstore[i][j2][2]) < prd_half[2])) {
          imgobjs[n] = Graphics::CYLINDER;
          imgparms[n][0] = type[i];
          imgparms[n][1] = xstore[i][j1][0];
          imgparms[n][2] = xstore[i][j1][1];
          imgparms[n][3] = xstore[i][j1][2] + zoffs;
          imgparms[n][4] = xstore[i][j2][0];
          imgparms[n][5] = xstore[i][j2][1];
          imgparms[n][6] = xstore[i][j2][2] + zoffs;
          imgparms[n][7] = 0.0;
          ++n;
        }
      }
    }
  }
  numobjs = n;
}

/* ----------------------------------------------------------------------
   current length of lines
------------------------------------------------------------------------- */

double FixGraphicsLines::compute_scalar()
{
  return nvalues;
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image
------------------------------------------------------------------------- */

int FixGraphicsLines::image(int *&objs, double **&parms)
{
  objs = imgobjs;
  parms = imgparms;
  return numobjs;
}
