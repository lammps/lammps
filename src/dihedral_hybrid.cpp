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

#include "dihedral_hybrid.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"

#include <cstring>

using namespace LAMMPS_NS;

#define EXTRA 1000

/* ---------------------------------------------------------------------- */

DihedralHybrid::DihedralHybrid(LAMMPS *lmp) : Dihedral(lmp)
{
  writedata = 0;
  nstyles = 0;
  ndihedrallist = nullptr;
  dihedrallist = nullptr;
  maxdihedral = nullptr;
}

/* ---------------------------------------------------------------------- */

DihedralHybrid::~DihedralHybrid()
{
  if (nstyles) {
    for (int i = 0; i < nstyles; i++) delete styles[i];
    delete[] styles;
    for (int i = 0; i < nstyles; i++) delete[] keywords[i];
    delete[] keywords;
  }

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(map);
    delete[] ndihedrallist;
    delete[] maxdihedral;
    for (int i = 0; i < nstyles; i++) memory->destroy(dihedrallist[i]);
    delete[] dihedrallist;
  }
}

/* ---------------------------------------------------------------------- */

void DihedralHybrid::compute(int eflag, int vflag)
{
  int i, j, m, n;

  // save ptrs to original dihedrallist

  int ndihedrallist_orig = neighbor->ndihedrallist;
  int **dihedrallist_orig = neighbor->dihedrallist;

  // if this is re-neighbor step, create sub-style dihedrallists
  // ndihedrallist[] = length of each sub-style list
  // realloc sub-style dihedrallist if necessary
  // load sub-style dihedrallist with 5 values from original dihedrallist

  if (neighbor->ago == 0) {
    for (m = 0; m < nstyles; m++) ndihedrallist[m] = 0;
    for (i = 0; i < ndihedrallist_orig; i++) {
      m = map[dihedrallist_orig[i][4]];
      if (m >= 0) ndihedrallist[m]++;
    }
    for (m = 0; m < nstyles; m++) {
      if (ndihedrallist[m] > maxdihedral[m]) {
        memory->destroy(dihedrallist[m]);
        maxdihedral[m] = ndihedrallist[m] + EXTRA;
        memory->create(dihedrallist[m], maxdihedral[m], 5, "dihedral_hybrid:dihedrallist");
      }
      ndihedrallist[m] = 0;
    }
    for (i = 0; i < ndihedrallist_orig; i++) {
      m = map[dihedrallist_orig[i][4]];
      if (m < 0) continue;
      n = ndihedrallist[m];
      dihedrallist[m][n][0] = dihedrallist_orig[i][0];
      dihedrallist[m][n][1] = dihedrallist_orig[i][1];
      dihedrallist[m][n][2] = dihedrallist_orig[i][2];
      dihedrallist[m][n][3] = dihedrallist_orig[i][3];
      dihedrallist[m][n][4] = dihedrallist_orig[i][4];
      ndihedrallist[m]++;
    }
  }

  // call each sub-style's compute function
  // set neighbor->dihedrallist to sub-style dihedrallist before call
  // accumulate sub-style global/peratom energy/virial in hybrid

  ev_init(eflag, vflag);

  // need to clear per-thread storage here, when using multiple threads
  // with thread-enabled substyles to avoid uninitlialized data access.

  const int nthreads = comm->nthreads;
  if (comm->nthreads > 1) {
    const bigint nall = atom->nlocal + atom->nghost;
    if (eflag_atom) memset(&eatom[0], 0, nall * nthreads * sizeof(double));
    if (vflag_atom) memset(&vatom[0][0], 0, 6 * nall * nthreads * sizeof(double));
  }

  for (m = 0; m < nstyles; m++) {
    neighbor->ndihedrallist = ndihedrallist[m];
    neighbor->dihedrallist = dihedrallist[m];

    styles[m]->compute(eflag, vflag);

    if (eflag_global) energy += styles[m]->energy;
    if (vflag_global)
      for (n = 0; n < 6; n++) virial[n] += styles[m]->virial[n];
    if (eflag_atom) {
      n = atom->nlocal;
      if (force->newton_bond) n += atom->nghost;
      double *eatom_substyle = styles[m]->eatom;
      for (i = 0; i < n; i++) eatom[i] += eatom_substyle[i];
    }
    if (vflag_atom) {
      n = atom->nlocal;
      if (force->newton_bond) n += atom->nghost;
      double **vatom_substyle = styles[m]->vatom;
      for (i = 0; i < n; i++)
        for (j = 0; j < 6; j++) vatom[i][j] += vatom_substyle[i][j];
    }
    if (cvflag_atom) {
      n = atom->nlocal;
      if (force->newton_bond) n += atom->nghost;
      double **cvatom_substyle = styles[m]->cvatom;
      for (i = 0; i < n; i++)
        for (j = 0; j < 9; j++) cvatom[i][j] += cvatom_substyle[i][j];
    }
  }

  // restore ptrs to original dihedrallist

  neighbor->ndihedrallist = ndihedrallist_orig;
  neighbor->dihedrallist = dihedrallist_orig;
}

/* ---------------------------------------------------------------------- */

void DihedralHybrid::allocate()
{
  allocated = 1;
  int np1 = atom->ndihedraltypes + 1;

  memory->create(map, np1, "dihedral:map");
  memory->create(setflag, np1, "dihedral:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;

  ndihedrallist = new int[nstyles];
  maxdihedral = new int[nstyles];
  dihedrallist = new int **[nstyles];
  for (int m = 0; m < nstyles; m++) maxdihedral[m] = 0;
  for (int m = 0; m < nstyles; m++) dihedrallist[m] = nullptr;
}

/* ----------------------------------------------------------------------
   create one dihedral style for each arg in list
------------------------------------------------------------------------- */

void DihedralHybrid::settings(int narg, char **arg)
{
  int i, m;

  if (narg < 1) utils::missing_cmd_args(FLERR, "dihedral_style hybrid", error);

  // delete old lists, since cannot just change settings

  if (nstyles) {
    for (i = 0; i < nstyles; i++) delete styles[i];
    delete[] styles;
    for (i = 0; i < nstyles; i++) delete[] keywords[i];
    delete[] keywords;
  }

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(map);
    delete[] ndihedrallist;
    delete[] maxdihedral;
    for (i = 0; i < nstyles; i++) memory->destroy(dihedrallist[i]);
    delete[] dihedrallist;
  }
  allocated = 0;

  // allocate list of sub-styles

  styles = new Dihedral *[narg];
  keywords = new char *[narg];

  // allocate each sub-style and call its settings() with subset of args
  // allocate uses suffix, but don't store suffix version in keywords,
  //   else syntax in coeff() will not match

  int dummy;
  nstyles = 0;
  i = 0;
  while (i < narg) {
    if (strcmp(arg[i], "hybrid") == 0)
      error->all(FLERR, "Dihedral style hybrid cannot have hybrid as an argument");

    if (strcmp(arg[i], "none") == 0)
      error->all(FLERR, "Dihedral style hybrid cannot have none as an argument");

    for (m = 0; m < nstyles; m++)
      if (strcmp(arg[i], keywords[m]) == 0)
        error->all(FLERR, "Dihedral style hybrid cannot use same dihedral style twice");

    styles[nstyles] = force->new_dihedral(arg[i], 1, dummy);
    keywords[nstyles] = force->store_style(arg[i], 0);

    // determine list of arguments for dihedral style settings
    // by looking for the next known dihedral style name.

    int jarg = i + 1;
    while ((jarg < narg) && !force->dihedral_map->count(arg[jarg]) &&
           !lmp->match_style("dihedral", arg[jarg]))
      jarg++;

    styles[nstyles]->settings(jarg - i - 1, &arg[i + 1]);
    i = jarg;
    nstyles++;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one type
---------------------------------------------------------------------- */

void DihedralHybrid::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->ndihedraltypes, ilo, ihi, error);

  // 2nd arg = dihedral sub-style name
  // allow for "none" or "skip" as valid sub-style name

  int m;
  for (m = 0; m < nstyles; m++)
    if (strcmp(arg[1], keywords[m]) == 0) break;

  int none = 0;
  int skip = 0;
  if (m == nstyles) {
    if (strcmp(arg[1], "none") == 0)
      none = 1;
    else if (strcmp(arg[1], "skip") == 0)
      none = skip = 1;
    else if (strcmp(arg[1], "mbt") == 0)
      error->all(FLERR, "MiddleBondTorsion coeff for hybrid dihedral has invalid format");
    else if (strcmp(arg[1], "ebt") == 0)
      error->all(FLERR, "EndBondTorsion coeff for hybrid dihedral has invalid format");
    else if (strcmp(arg[1], "at") == 0)
      error->all(FLERR, "AngleTorsion coeff for hybrid dihedral has invalid format");
    else if (strcmp(arg[1], "aat") == 0)
      error->all(FLERR, "AngleAngleTorsion coeff for hybrid dihedral has invalid format");
    else if (strcmp(arg[1], "bb13") == 0)
      error->all(FLERR, "BondBond13 coeff for hybrid dihedral has invalid format");
    else
      error->all(FLERR, "Dihedral coeff for hybrid has invalid style");
  }

  // move 1st arg to 2nd arg
  // just copy ptrs, since arg[] points into original input line

  arg[1] = arg[0];

  // invoke sub-style coeff() starting with 1st arg

  if (!none) styles[m]->coeff(narg - 1, &arg[1]);

  // set setflag and which type maps to which sub-style
  // if sub-style is skip: auxiliary class2 setting in data file so ignore
  // if sub-style is none and not skip: set hybrid setflag, wipe out map

  for (int i = ilo; i <= ihi; i++) {
    if (skip)
      continue;
    else if (none) {
      setflag[i] = 1;
      map[i] = -1;
    } else {
      setflag[i] = styles[m]->setflag[i];
      map[i] = m;
    }
  }
}

/* ----------------------------------------------------------------------
   run dihedral style specific initialization
------------------------------------------------------------------------- */

void DihedralHybrid::init_style()
{
  // error if sub-style is not used

  int used;
  for (int istyle = 0; istyle < nstyles; ++istyle) {
    used = 0;
    for (int itype = 1; itype <= atom->ndihedraltypes; ++itype)
      if (map[itype] == istyle) used = 1;
    if (used == 0) error->all(FLERR, "Dihedral hybrid sub-style {} is not used", keywords[istyle]);
  }

  for (int m = 0; m < nstyles; m++)
    if (styles[m]) styles[m]->init_style();
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void DihedralHybrid::write_restart(FILE *fp)
{
  fwrite(&nstyles, sizeof(int), 1, fp);

  int n;
  for (int m = 0; m < nstyles; m++) {
    n = strlen(keywords[m]) + 1;
    fwrite(&n, sizeof(int), 1, fp);
    fwrite(keywords[m], sizeof(char), n, fp);
    styles[m]->write_restart_settings(fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void DihedralHybrid::read_restart(FILE *fp)
{
  int me = comm->me;
  if (me == 0) utils::sfread(FLERR, &nstyles, sizeof(int), 1, fp, nullptr, error);
  MPI_Bcast(&nstyles, 1, MPI_INT, 0, world);
  styles = new Dihedral *[nstyles];
  keywords = new char *[nstyles];

  allocate();

  int n, dummy;
  for (int m = 0; m < nstyles; m++) {
    if (me == 0) utils::sfread(FLERR, &n, sizeof(int), 1, fp, nullptr, error);
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    keywords[m] = new char[n];
    if (me == 0) utils::sfread(FLERR, keywords[m], sizeof(char), n, fp, nullptr, error);
    MPI_Bcast(keywords[m], n, MPI_CHAR, 0, world);
    styles[m] = force->new_dihedral(keywords[m], 0, dummy);
    styles[m]->read_restart_settings(fp);
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double DihedralHybrid::memory_usage()
{
  double bytes = (double) maxeatom * sizeof(double);
  bytes += (double) maxvatom * 6 * sizeof(double);
  bytes += (double) maxcvatom * 9 * sizeof(double);
  for (int m = 0; m < nstyles; m++) bytes += (double) maxdihedral[m] * 5 * sizeof(int);
  for (int m = 0; m < nstyles; m++)
    if (styles[m]) bytes += styles[m]->memory_usage();
  return bytes;
}
