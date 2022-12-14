/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_sna_grid_local.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "sna.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

ComputeSNAGridLocal::ComputeSNAGridLocal(LAMMPS *lmp, int narg, char **arg) :
    ComputeGridLocal(lmp, narg, arg), cutsq(nullptr), radelem(nullptr), wjelem(nullptr)
{
  // skip over arguments used by base class
  // so that argument positions are identical to
  // regular per-atom compute

  arg += nargbase;
  narg -= nargbase;

  // begin code common to all SNAP computes

  double rfac0, rmin0;
  int twojmax, switchflag, bzeroflag, bnormflag, wselfallflag;

  int ntypes = atom->ntypes;
  int nargmin = 6 + 2 * ntypes;

  if (narg < nargmin) error->all(FLERR, "Illegal compute {} command", style);

  // default values

  rmin0 = 0.0;
  switchflag = 1;
  bzeroflag = 1;
  quadraticflag = 0;
  chemflag = 0;
  bnormflag = 0;
  wselfallflag = 0;
  switchinnerflag = 0;
  nelements = 1;

  // process required arguments

  memory->create(radelem, ntypes + 1, "sna/atom:radelem");    // offset by 1 to match up with types
  memory->create(wjelem, ntypes + 1, "sna/atom:wjelem");

  rcutfac = utils::numeric(FLERR, arg[3], false, lmp);
  rfac0 = utils::numeric(FLERR, arg[4], false, lmp);
  twojmax = utils::inumeric(FLERR, arg[5], false, lmp);

  for (int i = 0; i < ntypes; i++) radelem[i + 1] = utils::numeric(FLERR, arg[6 + i], false, lmp);
  for (int i = 0; i < ntypes; i++)
    wjelem[i + 1] = utils::numeric(FLERR, arg[6 + ntypes + i], false, lmp);

  // construct cutsq

  double cut;
  cutmax = 0.0;
  memory->create(cutsq, ntypes + 1, ntypes + 1, "sna/atom:cutsq");
  for (int i = 1; i <= ntypes; i++) {
    cut = 2.0 * radelem[i] * rcutfac;
    if (cut > cutmax) cutmax = cut;
    cutsq[i][i] = cut * cut;
    for (int j = i + 1; j <= ntypes; j++) {
      cut = (radelem[i] + radelem[j]) * rcutfac;
      cutsq[i][j] = cutsq[j][i] = cut * cut;
    }
  }

  // set local input checks

  int sinnerflag = 0;
  int dinnerflag = 0;

  // process optional args

  int iarg = nargmin;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "rmin0") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      rmin0 = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "switchflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      switchflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "bzeroflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      bzeroflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "quadraticflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      quadraticflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "chem") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      chemflag = 1;
      memory->create(map, ntypes + 1, "compute_sna_grid:map");
      nelements = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      for (int i = 0; i < ntypes; i++) {
        int jelem = utils::inumeric(FLERR, arg[iarg + 2 + i], false, lmp);
        if (jelem < 0 || jelem >= nelements) error->all(FLERR, "Illegal compute {} command", style);
        map[i + 1] = jelem;
      }
      iarg += 2 + ntypes;
    } else if (strcmp(arg[iarg], "bnormflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      bnormflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "wselfallflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      wselfallflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "switchinnerflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      switchinnerflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "sinner") == 0) {
      iarg++;
      if (iarg + ntypes > narg) error->all(FLERR, "Illegal compute {} command", style);
      memory->create(sinnerelem, ntypes + 1, "snap:sinnerelem");
      for (int i = 0; i < ntypes; i++)
        sinnerelem[i + 1] = utils::numeric(FLERR, arg[iarg + i], false, lmp);
      sinnerflag = 1;
      iarg += ntypes;
    } else if (strcmp(arg[iarg], "dinner") == 0) {
      iarg++;
      if (iarg + ntypes > narg) error->all(FLERR, "Illegal compute {} command", style);
      memory->create(dinnerelem, ntypes + 1, "snap:dinnerelem");
      for (int i = 0; i < ntypes; i++)
        dinnerelem[i + 1] = utils::numeric(FLERR, arg[iarg + i], false, lmp);
      dinnerflag = 1;
      iarg += ntypes;
    } else
      error->all(FLERR, "Illegal compute {} command", style);
  }

  if (switchinnerflag && !(sinnerflag && dinnerflag))
    error->all(FLERR,
               "Illegal compute {} command: switchinnerflag = 1, missing sinner/dinner keyword",
               style);

  if (!switchinnerflag && (sinnerflag || dinnerflag))
    error->all(FLERR,
               "Illegal compute {} command: switchinnerflag = 0, unexpected sinner/dinner keyword",
               style);

  snaptr = new SNA(lmp, rfac0, twojmax, rmin0, switchflag, bzeroflag, chemflag, bnormflag,
                   wselfallflag, nelements, switchinnerflag);

  ncoeff = snaptr->ncoeff;
  nvalues = ncoeff;
  if (quadraticflag) nvalues += (ncoeff * (ncoeff + 1)) / 2;

  // end code common to all SNAP computes

  size_local_cols = size_local_cols_base + nvalues;
}

/* ---------------------------------------------------------------------- */

ComputeSNAGridLocal::~ComputeSNAGridLocal()
{
  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(cutsq);
  delete snaptr;

  if (chemflag) memory->destroy(map);
}

/* ---------------------------------------------------------------------- */

void ComputeSNAGridLocal::init()
{
  if ((modify->get_compute_by_style("^sna/grid/local$").size() > 1) && (comm->me == 0))
    error->warning(FLERR, "More than one instance of compute sna/grid/local");
  snaptr->init();
}

/* ---------------------------------------------------------------------- */

void ComputeSNAGridLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // compute sna for each gridpoint

  double **const x = atom->x;
  const int *const mask = atom->mask;
  int *const type = atom->type;
  const int ntotal = atom->nlocal + atom->nghost;

  // insure rij, inside, and typej are of size jnum

  snaptr->grow_rij(ntotal);

  int igrid = 0;
  for (int iz = nzlo; iz <= nzhi; iz++)
    for (int iy = nylo; iy <= nyhi; iy++)
      for (int ix = nxlo; ix <= nxhi; ix++) {
        double xgrid[3];
        grid2x(ix, iy, iz, xgrid);
        const double xtmp = xgrid[0];
        const double ytmp = xgrid[1];
        const double ztmp = xgrid[2];

        // currently, all grid points are type 1
        // not clear what a better choice would be

        const int itype = 1;
        int ielem = 0;
        if (chemflag) ielem = map[itype];

        // rij[][3] = displacements between atom I and those neighbors
        // inside = indices of neighbors of I within cutoff
        // typej = types of neighbors of I within cutoff

        int ninside = 0;
        for (int j = 0; j < ntotal; j++) {

          // check that j is in compute group

          if (!(mask[j] & groupbit)) continue;

          const double delx = xtmp - x[j][0];
          const double dely = ytmp - x[j][1];
          const double delz = ztmp - x[j][2];
          const double rsq = delx * delx + dely * dely + delz * delz;
          int jtype = type[j];
          int jelem = 0;
          if (chemflag) jelem = map[jtype];
          if (rsq < cutsq[jtype][jtype] && rsq > 1e-20) {
            snaptr->rij[ninside][0] = delx;
            snaptr->rij[ninside][1] = dely;
            snaptr->rij[ninside][2] = delz;
            snaptr->inside[ninside] = j;
            snaptr->wj[ninside] = wjelem[jtype];
            snaptr->rcutij[ninside] = 2.0 * radelem[jtype] * rcutfac;
            if (switchinnerflag) {
              snaptr->sinnerij[ninside] = sinnerelem[jelem];
              snaptr->dinnerij[ninside] = dinnerelem[jelem];
            }
            if (chemflag)
              snaptr->element[ninside] = jelem;    // element index for multi-element snap
            ninside++;
          }
        }

        snaptr->compute_ui(ninside, ielem);
        snaptr->compute_zi();
        snaptr->compute_bi(ielem);

        // linear contributions

        for (int icoeff = 0; icoeff < ncoeff; icoeff++)
          alocal[igrid][size_local_cols_base + icoeff] = snaptr->blist[icoeff];

        // quadratic contributions

        if (quadraticflag) {
          int ncount = ncoeff;
          for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
            double bveci = snaptr->blist[icoeff];
            alocal[igrid][size_local_cols_base + ncount++] = 0.5 * bveci * bveci;
            for (int jcoeff = icoeff + 1; jcoeff < ncoeff; jcoeff++)
              alocal[igrid][size_local_cols_base + ncount++] = bveci * snaptr->blist[jcoeff];
          }
        }
        igrid++;
      }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double ComputeSNAGridLocal::memory_usage()
{
  double nbytes = snaptr->memory_usage();    // SNA object
  int n = atom->ntypes + 1;
  nbytes += (double) n * sizeof(int);    // map

  return nbytes;
}
