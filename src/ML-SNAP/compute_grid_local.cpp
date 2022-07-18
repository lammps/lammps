/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_grid_local.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>

// For the subdomain test below; grid-points and subdomain boundaries
// sometimes differ by minimal amounts (in the order of 2e-17).
static constexpr double EPSILON = 1.0e-10;

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeGridLocal::ComputeGridLocal(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), alocal(nullptr)
{
  if (narg < 6) error->all(FLERR, "Illegal compute grid/local command");

  local_flag = 1;
  size_local_cols = 0;
  size_local_rows = 0;
  extarray = 0;

  int iarg0 = 3;
  int iarg = iarg0;
  if (strcmp(arg[iarg], "grid") == 0) {
    if (iarg + 4 > narg) error->all(FLERR, "Illegal compute grid/local command");
    nx = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
    ny = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
    nz = utils::inumeric(FLERR, arg[iarg + 3], false, lmp);
    if (nx <= 0 || ny <= 0 || nz <= 0)
      error->all(FLERR, "All grid/local dimensions must be positive");
    iarg += 4;
  } else
    error->all(FLERR, "Illegal compute grid/local command");

  nargbase = iarg - iarg0;

  size_local_cols_base = 6;
  gridlocal_allocated = 0;
}

/* ---------------------------------------------------------------------- */

ComputeGridLocal::~ComputeGridLocal()
{
  deallocate();
}

/* ---------------------------------------------------------------------- */

void ComputeGridLocal::setup()
{
  deallocate();
  set_grid_global();
  set_grid_local();
  allocate();
  assign_coords();
}

/* ----------------------------------------------------------------------
   convert global array indexes to box coords
------------------------------------------------------------------------- */

void ComputeGridLocal::grid2x(int ix, int iy, int iz, double *x)
{
  x[0] = ix * delx;
  x[1] = iy * dely;
  x[2] = iz * delz;

  if (triclinic) domain->lamda2x(x, x);
}

/* ----------------------------------------------------------------------
   convert global array indexes to lamda coords; for orthorombic
   cells defaults to grid2x.
------------------------------------------------------------------------- */

void ComputeGridLocal::grid2lamda(int ix, int iy, int iz, double *x)
{
  x[0] = ix * delx;
  x[1] = iy * dely;
  x[2] = iz * delz;
}

/* ----------------------------------------------------------------------
   create arrays
------------------------------------------------------------------------- */

void ComputeGridLocal::allocate()
{
  if (nxlo <= nxhi && nylo <= nyhi && nzlo <= nzhi) {
    gridlocal_allocated = 1;
    memory->create(alocal, size_local_rows, size_local_cols, "compute/grid/local:alocal");
    array_local = alocal;
  }
}

/* ----------------------------------------------------------------------
   free arrays
------------------------------------------------------------------------- */

void ComputeGridLocal::deallocate()
{
  if (gridlocal_allocated) {
    gridlocal_allocated = 0;
    memory->destroy(alocal);
  }
  array_local = nullptr;
}

/* ----------------------------------------------------------------------
   set global grid
------------------------------------------------------------------------- */

void ComputeGridLocal::set_grid_global()
{
  // calculate grid layout

  triclinic = domain->triclinic;

  if (triclinic == 0) {
    prd = domain->prd;
    boxlo = domain->boxlo;
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    prd = domain->prd_lamda;
    boxlo = domain->boxlo_lamda;
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];

  delxinv = nx / xprd;
  delyinv = ny / yprd;
  delzinv = nz / zprd;

  delx = 1.0 / delxinv;
  dely = 1.0 / delyinv;
  delz = 1.0 / delzinv;
}

/* ----------------------------------------------------------------------
   set local subset of grid that I own
   n xyz lo/hi = 3d brick that I own (inclusive)
------------------------------------------------------------------------- */

void ComputeGridLocal::set_grid_local()
{
  // nx,ny,nz = extent of global grid
  // indices into the global grid range from 0 to N-1 in each dim
  // if grid point is inside my sub-domain I own it,
  //   this includes sub-domain lo boundary but excludes hi boundary
  // ixyz lo/hi = inclusive lo/hi bounds of global grid sub-brick I own
  // if proc owns no grid cells in a dim, then ilo > ihi
  // if 2 procs share a boundary a grid point is exactly on,
  //   the 2 equality if tests insure a consistent decision
  //   as to which proc owns it

  double xfraclo, xfrachi, yfraclo, yfrachi, zfraclo, zfrachi;

  if (comm->layout != Comm::LAYOUT_TILED) {
    xfraclo = comm->xsplit[comm->myloc[0]];
    xfrachi = comm->xsplit[comm->myloc[0] + 1];
    yfraclo = comm->ysplit[comm->myloc[1]];
    yfrachi = comm->ysplit[comm->myloc[1] + 1];
    zfraclo = comm->zsplit[comm->myloc[2]];
    zfrachi = comm->zsplit[comm->myloc[2] + 1];
  } else {
    xfraclo = comm->mysplit[0][0];
    xfrachi = comm->mysplit[0][1];
    yfraclo = comm->mysplit[1][0];
    yfrachi = comm->mysplit[1][1];
    zfraclo = comm->mysplit[2][0];
    zfrachi = comm->mysplit[2][1];
  }

  nxlo = static_cast<int>(xfraclo * nx);
  if (1.0 * nxlo != xfraclo * nx) nxlo++;
  nxhi = static_cast<int>(xfrachi * nx);
  if (1.0 * nxhi == xfrachi * nx) nxhi--;

  nylo = static_cast<int>(yfraclo * ny);
  if (1.0 * nylo != yfraclo * ny) nylo++;
  nyhi = static_cast<int>(yfrachi * ny);
  if (1.0 * nyhi == yfrachi * ny) nyhi--;

  nzlo = static_cast<int>(zfraclo * nz);
  if (1.0 * nzlo != zfraclo * nz) nzlo++;
  nzhi = static_cast<int>(zfrachi * nz);
  if (1.0 * nzhi == zfrachi * nz) nzhi--;

  size_local_rows = (nxhi - nxlo + 1) * (nyhi - nylo + 1) * (nzhi - nzlo + 1);
}

/* ----------------------------------------------------------------------
   copy coords to local array
------------------------------------------------------------------------- */

void ComputeGridLocal::assign_coords()
{
  int igrid = 0;
  for (int iz = nzlo; iz <= nzhi; iz++)
    for (int iy = nylo; iy <= nyhi; iy++)
      for (int ix = nxlo; ix <= nxhi; ix++) {
        alocal[igrid][0] = ix;
        alocal[igrid][1] = iy;
        alocal[igrid][2] = iz;
        double xgrid[3];

        // for triclinic: create gridpoint in lamda coordinates and transform after check.
        // for orthorombic: create gridpoint in box coordinates.

        if (triclinic)
          grid2lamda(ix, iy, iz, xgrid);
        else
          grid2x(ix, iy, iz, xgrid);

        // ensure gridpoint is not strictly outside subdomain

        if ((sublo[0] - xgrid[0]) > EPSILON || (xgrid[0] - subhi[0]) > EPSILON ||
            (sublo[1] - xgrid[1]) > EPSILON || (xgrid[1] - subhi[1]) > EPSILON ||
            (sublo[2] - xgrid[2]) > EPSILON || (xgrid[2] - subhi[2]) > EPSILON)
          error->one(FLERR, "Invalid gridpoint position in compute grid/local");

        // convert lamda to x, y, z, after sudomain check

        if (triclinic) domain->lamda2x(xgrid, xgrid);

        alocal[igrid][3] = xgrid[0];
        alocal[igrid][4] = xgrid[1];
        alocal[igrid][5] = xgrid[2];
        igrid++;
      }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeGridLocal::memory_usage()
{
  double nbytes = (double) size_local_rows * size_local_cols * sizeof(double);    // gridlocal
  return nbytes;
}
