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

#include "compute_grid.h"

#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeGrid::ComputeGrid(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), grid(nullptr), gridall(nullptr), gridlocal(nullptr)
{
  if (narg < 6) error->all(FLERR, "Illegal compute grid command");

  array_flag = 1;
  size_array_cols = 0;
  size_array_rows = 0;
  extarray = 0;

  int iarg0 = 3;
  int iarg = iarg0;
  if (strcmp(arg[iarg], "grid") == 0) {
    if (iarg + 4 > narg) error->all(FLERR, "Illegal compute grid command");
    nx = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
    ny = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
    nz = utils::inumeric(FLERR, arg[iarg + 3], false, lmp);
    if (nx <= 0 || ny <= 0 || nz <= 0) error->all(FLERR, "All grid dimensions must be positive");
    iarg += 4;
  } else
    error->all(FLERR, "Illegal compute grid command");

  nargbase = iarg - iarg0;

  size_array_rows = nx * ny * nz;
  size_array_cols_base = 3;
  gridlocal_allocated = 0;
}

/* ---------------------------------------------------------------------- */

ComputeGrid::~ComputeGrid()
{
  deallocate();
}

/* ---------------------------------------------------------------------- */

void ComputeGrid::setup()
{
  deallocate();
  set_grid_global();
  set_grid_local();
  allocate();
}

/* ----------------------------------------------------------------------
   convert global array index to box coords
------------------------------------------------------------------------- */

void ComputeGrid::grid2x(int igrid, double *x)
{
  int iz = igrid / (nx * ny);
  igrid -= iz * (nx * ny);
  int iy = igrid / nx;
  igrid -= iy * nx;
  int ix = igrid;

  x[0] = ix * delx;
  x[1] = iy * dely;
  x[2] = iz * delz;

  if (triclinic) domain->lamda2x(x, x);
}

/* ----------------------------------------------------------------------
   copy coords to global array
------------------------------------------------------------------------- */

void ComputeGrid::assign_coords_all()
{
  double x[3];
  for (int igrid = 0; igrid < size_array_rows; igrid++) {
    grid2x(igrid, x);
    gridall[igrid][0] = x[0];
    gridall[igrid][1] = x[1];
    gridall[igrid][2] = x[2];
  }
}

/* ----------------------------------------------------------------------
   create arrays
------------------------------------------------------------------------- */

void ComputeGrid::allocate()
{
  // allocate arrays

  memory->create(grid, size_array_rows, size_array_cols, "grid:grid");
  memory->create(gridall, size_array_rows, size_array_cols, "grid:gridall");
  if (nxlo <= nxhi && nylo <= nyhi && nzlo <= nzhi) {
    gridlocal_allocated = 1;
    memory->create4d_offset(gridlocal, size_array_cols, nzlo, nzhi, nylo, nyhi, nxlo, nxhi,
                            "grid:gridlocal");
  }
  array = gridall;
}

/* ----------------------------------------------------------------------
   free arrays
------------------------------------------------------------------------- */

void ComputeGrid::deallocate()
{
  memory->destroy(grid);
  memory->destroy(gridall);
  if (gridlocal_allocated) {
    gridlocal_allocated = 0;
    memory->destroy4d_offset(gridlocal, nzlo, nylo, nxlo);
  }
  array = nullptr;
}

/* ----------------------------------------------------------------------
   set global grid
------------------------------------------------------------------------- */

void ComputeGrid::set_grid_global()
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

void ComputeGrid::set_grid_local()
{
  // nx,ny,nz = extent of global grid
  // indices into the global grid range from 0 to N-1 in each dim
  // if grid point is inside my sub-domain I own it,
  //   this includes sub-domain lo boundary but excludes hi boundary
  // ixyz lo/hi = inclusive lo/hi bounds of global grid sub-brick I own
  // if proc owns no grid cells in a dim, then ilo > ihi
  // if 2 procs share a boundary a grid point is exactly on,
  //   the 2 equality if tests ensure a consistent decision
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

  ngridlocal = (nxhi - nxlo + 1) * (nyhi - nylo + 1) * (nzhi - nzlo + 1);
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeGrid::memory_usage()
{
  double nbytes = (double) size_array_rows * size_array_cols * sizeof(double);    // grid
  nbytes += (double) size_array_rows * size_array_cols * sizeof(double);          // gridall
  nbytes += (double) size_array_cols * ngridlocal * sizeof(double);               // gridlocal
  return nbytes;
}
