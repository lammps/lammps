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
   Electronic stopping power
   Contributing authors: K. Avchaciov and T. Metspalu
   Information: k.avchachov@gmail.com
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstring>
#include "fix_elstop.h"
#include "mpi.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "force.h"
#include "pair.h"
#include "fix.h"
#include "compute.h"
#include "modify.h"
#include "memory.h"
#include "comm.h"
#include "error.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

FixElstop::FixElstop(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  scalar_flag = 1;  // Has compute_scalar
  global_freq = 1;  // SeLoss computed every step
  extscalar = 0;    // SeLoss compute_scalar is intensive
  nevery = 1;       // Run fix every step


  // Make sure the id for the kinetic energy compute is unique
  // by prepending the ID of this fix.
  int n = strlen(id) + strlen("_ke_atom") + 1;
  id_ke_atom = new char[n];
  strcpy(id_ke_atom, id);
  strcat(id_ke_atom, "_ke_atom");

  char *newarg[3];
  newarg[0] = id_ke_atom;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "ke/atom";
  modify->add_compute(3, newarg);


  // args: 0 = fix ID, 1 = group ID,  2 = "elstop"
  //       3 = Ecut,   4 = file path
  // optional rest: "region" <region name>
  //                "minneigh" <min number of neighbors>

  if (narg < 5)
    error->all(FLERR, "Illegal fix elstop command: too few arguments");

  Ecut = force->numeric(FLERR, arg[3]);
  if (Ecut <= 0.0) error->all(FLERR, "Illegal fix elstop command: Ecut <= 0");

  int iarg = 5;
  iregion = -1;
  minneigh = 1;
  bool minneighflag = false;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "region") == 0) {
      if (iregion >= 0)
         error->all(FLERR, "Illegal fix elstop command: region given twice");
      if (iarg+2 > narg)
        error->all(FLERR, "Illegal fix elstop command: region name missing");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion < 0)
        error->all(FLERR, "Region ID for fix elstop does not exist");
      iarg += 2;
    }
    else if (strcmp(arg[iarg], "minneigh") == 0) {
      if (minneighflag)
         error->all(FLERR, "Illegal fix elstop command: minneigh given twice");
      minneighflag = true;
      if (iarg+2 > narg)
        error->all(FLERR, "Illegal fix elstop command: minneigh number missing");
      minneigh = force->inumeric(FLERR, arg[iarg+1]);
      if (minneigh < 0)
        error->all(FLERR, "Illegal fix elstop command: minneigh < 0");
      iarg += 2;
    }
    else error->all(FLERR, "Illegal fix elstop command: unknown argument");
  }


  // Read the input file for energy ranges and stopping powers.
  // First proc 0 reads the file, then bcast to others.
  const int ncol = atom->ntypes + 1;
  if (comm->me == 0) {
    maxlines = 300;
    memory->create(elstop_ranges, ncol, maxlines, "elstop:tabs");
    read_table(arg[4]);
  }

  MPI_Bcast(&maxlines, 1 , MPI_INT, 0, world);
  MPI_Bcast(&table_entries, 1 , MPI_INT, 0, world);

  if (comm->me != 0)
    memory->create(elstop_ranges, ncol, maxlines, "elstop:tabs");

  MPI_Bcast(&elstop_ranges[0][0], ncol*maxlines, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

FixElstop::~FixElstop()
{
  memory->destroy(elstop_ranges);
  modify->delete_compute(id_ke_atom);
  delete id_ke_atom;
}

/* ---------------------------------------------------------------------- */

int FixElstop::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixElstop::init()
{
  SeLoss_sync_flag = 0;
  SeLoss = 0.0;

  int ikeatom = modify->find_compute(id_ke_atom);
  if (ikeatom < 0)
    error->all(FLERR, "KE compute ID for fix elstop does not exist");
  c_ke = modify->compute[ikeatom];


  // need an occasional full neighbor list
  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void FixElstop::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixElstop::post_force(int /*vflag*/)
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

  c_ke->compute_peratom();
  double *ke = c_ke->vector_atom;

  for (int i = 0; i < nlocal; ++i) {

    // Do fast checks first, only then the region check
    if (!(mask[i] & groupbit)) continue;

    // Avoid atoms outside bulk material
    if (numneigh[i] < minneigh) continue;

    double energy = ke[i];
    if (energy < Ecut) continue;
    if (energy < elstop_ranges[0][0]) continue;
    if (energy > elstop_ranges[0][table_entries - 1])
      error->one(FLERR, "Atom kinetic energy too high for fix elstop");

    if (iregion >= 0) {
      // Only apply in the given region
      if (domain->regions[iregion]->match(x[i][0], x[i][1], x[i][2]) != 1)
        continue;
    }

    // Binary search to find correct energy range
    int iup = table_entries - 1;
    int idown = 0;
    while (true) {
      int ihalf = idown + (iup - idown) / 2;
      if (ihalf == idown) break;
      if (elstop_ranges[0][ihalf] < energy) idown = ihalf;
      else iup = ihalf;
    }

    int itype = type[i];
    double Se_lo = elstop_ranges[itype][idown];
    double Se_hi = elstop_ranges[itype][iup];
    double E_lo = elstop_ranges[0][idown];
    double E_hi = elstop_ranges[0][iup];

    // Get elstop with a simple linear interpolation
    double Se = (Se_hi - Se_lo) / (E_hi - E_lo) * (energy - E_lo) + Se_lo;

    double v2 = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
    double vabs = sqrt(v2);
    double factor = -Se / vabs;

    f[i][0] += v[i][0] * factor;
    f[i][1] += v[i][1] * factor;
    f[i][2] += v[i][2] * factor;

    SeLoss += Se * vabs * dt; // very rough approx
  }
}

/* ---------------------------------------------------------------------- */

double FixElstop::compute_scalar()
{
  // only sum across procs when changed since last call

  if (SeLoss_sync_flag == 0) {
    MPI_Allreduce(&SeLoss, &SeLoss_all, 1, MPI_DOUBLE, MPI_SUM, world);
    SeLoss_sync_flag = 1;
  }
  return SeLoss_all;
}

/* ---------------------------------------------------------------------- */

void FixElstop::read_table(const char *file)
{
  char line[MAXLINE];

  FILE *fp = force->open_potential(file);
  if (fp == NULL) {
    char str[128];
    snprintf(str, 128, "Cannot open stopping range table %s", file);
    error->one(FLERR, str);
  }

  const int ncol = atom->ntypes + 1;

  int l = 0;
  while (true) {
    if (fgets(line, MAXLINE, fp) == NULL) break; // end of file
    if (line[0] == '#') continue; // comment

    char *pch = strtok(line, " \t\n\r");
    if (pch == NULL) continue; // blank line

    if (l >= maxlines) grow_table();

    int i = 0;
    for ( ; i < ncol && pch != NULL; i++) {
      elstop_ranges[i][l] = force->numeric(FLERR, pch);
      pch = strtok(NULL, " \t\n\r");
    }

    if (i != ncol || pch != NULL) // too short or too long
      error->one(FLERR, "fix elstop: Invalid table line");

    if (l >= 1 && elstop_ranges[0][l] <= elstop_ranges[0][l-1])
      error->one(FLERR, "fix elstop: Energies must be in ascending order");

    l++;
  }
  table_entries = l;

  if (table_entries == 0)
    error->one(FLERR, "Did not find any data in elstop table file");

  fclose(fp);
}

/* ---------------------------------------------------------------------- */

void FixElstop::grow_table()
{
  const int ncol = atom->ntypes + 1;
  int new_maxlines = 2 * maxlines;

  double **new_array;
  memory->create(new_array, ncol, new_maxlines, "elstop:tabscopy");

  for (int i = 0; i < ncol; i++)
    memcpy(new_array[i], elstop_ranges[i], maxlines*sizeof(double));

  memory->destroy(elstop_ranges);
  elstop_ranges = new_array;
  maxlines = new_maxlines;
}
