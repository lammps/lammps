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

#include "compute_gaussian_grid_local.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathConst::MY_2PI;
using MathSpecial::powint;

ComputeGaussianGridLocal::ComputeGaussianGridLocal(LAMMPS *lmp, int narg, char **arg) :
    ComputeGridLocal(lmp, narg, arg), cutsq(nullptr), radelem(nullptr),
    sigmaelem(nullptr), prefacelem(nullptr), argfacelem(nullptr)
{
  // skip over arguments used by base class
  // so that argument positions are identical to
  // regular per-atom compute

  arg += nargbase;
  narg -= nargbase;

  //double rfac0, rmin0;
  //int twojmax, switchflag, bzeroflag, bnormflag, wselfallflag;

  int ntypes = atom->ntypes;
  int nargmin = 4 + 2 * ntypes;

  if (narg < nargmin) error->all(FLERR, "Illegal compute {} command", style);

  // process required arguments

  memory->create(radelem, ntypes + 1, "gaussian/atom:radelem");    // offset by 1 to match up with types
  memory->create(sigmaelem, ntypes + 1, "gaussian/atom:sigmaelem");
  memory->create(prefacelem, ntypes + 1, "gaussian/atom:prefacelem");
  memory->create(argfacelem, ntypes + 1, "gaussian/atom:argfacelem");

  rcutfac = utils::numeric(FLERR, arg[3], false, lmp);

  for (int i = 0; i < ntypes; i++) radelem[i + 1] = utils::numeric(FLERR, arg[4 + i], false, lmp);
  for (int i = 0; i < ntypes; i++)
    sigmaelem[i + 1] = utils::numeric(FLERR, arg[ntypes + 4 + i], false, lmp);

  // construct cutsq
  double cut;
  cutmax = 0.0;
  memory->create(cutsq, ntypes + 1, ntypes + 1, "gaussian/atom:cutsq");
  for (int i = 1; i <= ntypes; i++) {
    cut = 2.0 * radelem[i] * rcutfac;
    if (cut > cutmax) cutmax = cut;
    cutsq[i][i] = cut * cut;
    for (int j = i + 1; j <= ntypes; j++) {
      cut = (radelem[i] + radelem[j]) * rcutfac;
      cutsq[i][j] = cutsq[j][i] = cut * cut;
    }
  }

  size_local_cols = size_local_cols_base + ntypes;

  // pre-compute coefficients
  for (int i = 0; i < ntypes; i++) {
    prefacelem[i + 1] = 1.0/powint(sigmaelem[i + 1] * sqrt(MY_2PI), 3);
    argfacelem[i + 1] = 1.0/(2.0 * sigmaelem[i + 1] * sigmaelem[i + 1]);
  }
}

/* ---------------------------------------------------------------------- */

ComputeGaussianGridLocal::~ComputeGaussianGridLocal()
{
  if (copymode) return;
  memory->destroy(radelem);
  memory->destroy(sigmaelem);
  memory->destroy(prefacelem);
  memory->destroy(argfacelem);
  memory->destroy(cutsq);
}

/* ---------------------------------------------------------------------- */

void ComputeGaussianGridLocal::init()
{
  if ((modify->get_compute_by_style("^gaussian/grid/local$").size() > 1) && (comm->me == 0))
    error->warning(FLERR, "More than one instance of compute gaussian/grid/local");
}

/* ---------------------------------------------------------------------- */

void ComputeGaussianGridLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // compute gaussian for each gridpoint

  double **const x = atom->x;
  const int *const mask = atom->mask;
  int *const type = atom->type;
  const int ntotal = atom->nlocal + atom->nghost;

  int igrid = 0;
  for (int iz = nzlo; iz <= nzhi; iz++)
    for (int iy = nylo; iy <= nyhi; iy++)
      for (int ix = nxlo; ix <= nxhi; ix++) {
        double xgrid[3];
        grid2x(ix, iy, iz, xgrid);
        const double xtmp = xgrid[0];
        const double ytmp = xgrid[1];
        const double ztmp = xgrid[2];

        // Zeroing out the components, which are filled as a sum.
        for (int icol = size_local_cols_base; icol < size_local_cols; icol++){
          alocal[igrid][icol] = 0.0;
        }

        for (int j = 0; j < ntotal; j++) {

          // check that j is in compute group

          if (!(mask[j] & groupbit)) continue;

          const double delx = xtmp - x[j][0];
          const double dely = ytmp - x[j][1];
          const double delz = ztmp - x[j][2];
          const double rsq = delx * delx + dely * dely + delz * delz;
          int jtype = type[j];
          if (rsq < cutsq[jtype][jtype]) {
            int icol = size_local_cols_base + jtype - 1;
            alocal[igrid][icol] += prefacelem[jtype] * exp(-rsq * argfacelem[jtype]);
          }
        }
        igrid++;
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double ComputeGaussianGridLocal::memory_usage()
{
  int n = atom->ntypes + 1;
  int nbytes = (double) n * sizeof(int);    // map

  return nbytes;
}
