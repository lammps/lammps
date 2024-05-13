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
   Electronic stopping power
   Contributing authors: K. Avchaciov and T. Metspalu
   Information: k.avchachov@gmail.com
------------------------------------------------------------------------- */

#include "fix_electron_stopping.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"
#include "region.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <exception>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixElectronStopping::FixElectronStopping(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), elstop_ranges(nullptr), idregion(nullptr), region(nullptr), list(nullptr)
{
  scalar_flag = 1;    // Has compute_scalar
  global_freq = 1;    // SeLoss computed every step
  extscalar = 0;      // SeLoss compute_scalar is intensive
  nevery = 1;         // Run fix every step

  // args: 0 = fix ID, 1 = group ID,  2 = "electron/stopping"
  //       3 = Ecut,   4 = file path
  // optional rest: "region" <region name>
  //                "minneigh" <min number of neighbors>

  if (narg < 5) error->all(FLERR, "Illegal fix electron/stopping command: too few arguments");

  Ecut = utils::numeric(FLERR, arg[3], false, lmp);
  if (Ecut <= 0.0) error->all(FLERR, "Illegal fix electron/stopping command: Ecut <= 0");

  int iarg = 5;
  minneigh = 1;
  bool minneighflag = false;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "region") == 0) {
      if (region) error->all(FLERR, "Illegal fix electron/stopping command: region given twice");
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal fix electron/stopping command: region name missing");
      region = domain->get_region_by_id(arg[iarg + 1]);
      if (!region)
        error->all(FLERR, "Region {} for fix electron/stopping does not exist", arg[iarg + 1]);
      idregion = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "minneigh") == 0) {
      if (minneighflag)
        error->all(FLERR, "Illegal fix electron/stopping command: minneigh given twice");
      minneighflag = true;
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal fix electron/stopping command: minneigh number missing");
      minneigh = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (minneigh < 0) error->all(FLERR, "Illegal fix electron/stopping command: minneigh < 0");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix electron/stopping command: unknown argument");
  }

  // Read the input file for energy ranges and stopping powers.
  // First proc 0 reads the file, then bcast to others.
  const int ncol = atom->ntypes + 1;
  if (comm->me == 0) {
    maxlines = 300;
    memory->create(elstop_ranges, ncol, maxlines, "electron/stopping:table");
    read_table(arg[4]);
  }

  MPI_Bcast(&maxlines, 1, MPI_INT, 0, world);
  MPI_Bcast(&table_entries, 1, MPI_INT, 0, world);

  if (comm->me != 0) memory->create(elstop_ranges, ncol, maxlines, "electron/stopping:table");

  MPI_Bcast(&elstop_ranges[0][0], ncol * maxlines, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

FixElectronStopping::~FixElectronStopping()
{
  delete[] idregion;
  memory->destroy(elstop_ranges);
}

/* ---------------------------------------------------------------------- */

int FixElectronStopping::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixElectronStopping::init()
{
  SeLoss_sync_flag = 0;
  SeLoss = 0.0;
  if (idregion) {
    region = domain->get_region_by_id(idregion);
    if (!region) error->all(FLERR, "Region {} for fix electron/stopping does not exist", idregion);
  }

  // need an occasional full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);
}

/* ---------------------------------------------------------------------- */

void FixElectronStopping::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixElectronStopping::post_force(int /*vflag*/)
{
  SeLoss_sync_flag = 0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double dt = update->dt;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  neighbor->build_one(list);
  int *numneigh = list->numneigh;

  for (int i = 0; i < nlocal; ++i) {

    // Do fast checks first, only then the region check
    if (!(mask[i] & groupbit)) continue;

    // Avoid atoms outside bulk material
    if (numneigh[i] < minneigh) continue;

    int itype = type[i];
    double massone = (atom->rmass) ? atom->rmass[i] : atom->mass[itype];
    double v2 = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
    double energy = 0.5 * force->mvv2e * massone * v2;

    if (energy < Ecut) continue;
    if (energy < elstop_ranges[0][0]) continue;
    if (energy > elstop_ranges[0][table_entries - 1])
      error->one(FLERR, "Fix electron/stopping: kinetic energy too high for atom {}: {} vs {}",
                 atom->tag[i], energy, elstop_ranges[0][table_entries - 1]);

    if (region) {
      // Only apply in the given region
      if (region->match(x[i][0], x[i][1], x[i][2]) != 1) continue;
    }

    // Binary search to find correct energy range
    int iup = table_entries - 1;
    int idown = 0;
    while (true) {
      int ihalf = idown + (iup - idown) / 2;
      if (ihalf == idown) break;
      if (elstop_ranges[0][ihalf] < energy)
        idown = ihalf;
      else
        iup = ihalf;
    }

    double Se_lo = elstop_ranges[itype][idown];
    double Se_hi = elstop_ranges[itype][iup];
    double E_lo = elstop_ranges[0][idown];
    double E_hi = elstop_ranges[0][iup];

    // Get electronic stopping with a simple linear interpolation
    double Se = (Se_hi - Se_lo) / (E_hi - E_lo) * (energy - E_lo) + Se_lo;

    double vabs = sqrt(v2);
    double factor = -Se / vabs;

    f[i][0] += v[i][0] * factor;
    f[i][1] += v[i][1] * factor;
    f[i][2] += v[i][2] * factor;

    SeLoss += Se * vabs * dt;    // very rough approx
  }
}

/* ---------------------------------------------------------------------- */

double FixElectronStopping::compute_scalar()
{
  // only sum across procs when changed since last call

  if (SeLoss_sync_flag == 0) {
    MPI_Allreduce(&SeLoss, &SeLoss_all, 1, MPI_DOUBLE, MPI_SUM, world);
    SeLoss_sync_flag = 1;
  }
  return SeLoss_all;
}

/* ----------------------------------------------------------------------
   read electron stopping parameters. only called from MPI rank 0.
   format: energy then one column per atom type
   read as many lines as available.
   energies must be sorted in ascending order.
   ---------------------------------------------------------------------- */

void FixElectronStopping::read_table(const char *file)
{
  const int ncol = atom->ntypes + 1;
  int nlines = 0;
  PotentialFileReader reader(lmp, file, "electron stopping data table");

  try {
    char *line;
    double oldvalue = 0.0;

    while ((line = reader.next_line())) {
      if (nlines >= maxlines) grow_table();
      ValueTokenizer values(line);
      elstop_ranges[0][nlines] = values.next_double();
      if (elstop_ranges[0][nlines] <= oldvalue)
        throw TokenizerException("energy values must be positive and in ascending order", line);

      oldvalue = elstop_ranges[0][nlines];
      for (int i = 1; i < ncol; ++i) elstop_ranges[i][nlines] = values.next_double();

      ++nlines;
    }
  } catch (std::exception &e) {
    error->one(FLERR, "Problem parsing electron stopping data: {}", e.what());
  }
  if (nlines == 0) error->one(FLERR, "Did not find any data in electron/stopping table file");

  table_entries = nlines;
}

/* ---------------------------------------------------------------------- */

void FixElectronStopping::grow_table()
{
  const int ncol = atom->ntypes + 1;
  int new_maxlines = 2 * maxlines;

  double **new_array;
  memory->create(new_array, ncol, new_maxlines, "electron/stopping:table");

  for (int i = 0; i < ncol; i++) memcpy(new_array[i], elstop_ranges[i], maxlines * sizeof(double));

  memory->destroy(elstop_ranges);
  elstop_ranges = new_array;
  maxlines = new_maxlines;
}
