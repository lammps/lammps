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

#include "compute_grid.h"
#include <mpi.h>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeGrid::ComputeGrid(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal compute grid command");

  array_flag = 1;
  size_array_cols = 0;
  size_array_rows = 0;
  extarray = 0;

  int iarg0 = 3;
  int iarg = iarg0;
  if (strcmp(arg[iarg],"grid") == 0) {
    if (iarg+4 > narg) error->all(FLERR,"Illegal compute grid command");
    nx = force->inumeric(FLERR,arg[iarg+1]);
    ny = force->inumeric(FLERR,arg[iarg+2]);
    nz = force->inumeric(FLERR,arg[iarg+3]);
    if (nx <= 0 || ny <= 0 || nz <= 0)
      error->all(FLERR,"All grid dimensions must be positive");
    iarg += 4;
  } else error->all(FLERR,"Illegal compute grid command");

  nargbase = iarg - iarg0;

  size_array_rows = ngrid = nx*ny*nz;
  size_array_cols_base = 3;
}

/* ---------------------------------------------------------------------- */

ComputeGrid::~ComputeGrid()
{
  memory->destroy(grid);
  memory->destroy(grid_local);
}

/* ---------------------------------------------------------------------- */

void ComputeGrid::init()
{
}

/* ---------------------------------------------------------------------- */

void ComputeGrid::setup()
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

  delxinv = nx/xprd;
  delyinv = ny/yprd;
  delzinv = nz/zprd;
  
  delx = 1.0/delxinv;
  dely = 1.0/delyinv;
  delz = 1.0/delzinv;

  allocate();
  assign_grid_coords();
  assign_grid_local();
}

/* ----------------------------------------------------------------------
   convert global array index to box coords
------------------------------------------------------------------------- */

void ComputeGrid::grid2x(int igrid, double *x)
{
  int iz = igrid / (nx*ny);
  igrid -= iz * (nx*ny);
  int iy = igrid / nx;
  igrid -= iy * nx;
  int ix = igrid;

  x[0] = ix*delx;
  x[1] = iy*dely;
  x[2] = iz*delz;

  if (triclinic) domain->lamda2x(x, x);

}

/* ----------------------------------------------------------------------
   check if grid point is local
------------------------------------------------------------------------- */

int ComputeGrid::check_grid_local(int igrid)
{
  double x[3];

  int iz = igrid / (nx*ny);
  igrid -= iz * (nx*ny);
  int iy = igrid / nx;
  igrid -= iy * nx;
  int ix = igrid;

  x[0] = ix*delx;
  x[1] = iy*dely;
  x[2] = iz*delz;

  int islocal = 
    x[0] >= sublo[0] && x[0] < subhi[0] &&
    x[1] >= sublo[1] && x[1] < subhi[1] &&
    x[2] >= sublo[2] && x[2] < subhi[2];

  return islocal;
}

/* ----------------------------------------------------------------------
   copy coords to global array
------------------------------------------------------------------------- */

void ComputeGrid::assign_grid_coords()
{
  double x[3];
  for (int igrid = 0; igrid < ngrid; igrid++) {
    grid2x(igrid,x);
    grid[igrid][0] = x[0];
    grid[igrid][1] = x[1];
    grid[igrid][2] = x[2];
  }
}

/* ----------------------------------------------------------------------
   copy coords to global array
------------------------------------------------------------------------- */

void ComputeGrid::assign_grid_local()
{
  double x[3];
  for (int igrid = 0; igrid < ngrid; igrid++) {
    if (check_grid_local(igrid))
      grid_local[igrid] = 1;
    else {
      grid_local[igrid] = 0;
      memset(grid[igrid],0,size_array_cols);
    }
  }
}

/* ----------------------------------------------------------------------
   free and reallocate arrays
------------------------------------------------------------------------- */

void ComputeGrid::allocate()
{
  // grow global array if necessary

  memory->destroy(grid);
  memory->destroy(grid_local);
  memory->create(grid,size_array_rows,size_array_cols,"grid:grid");
  memory->create(gridall,size_array_rows,size_array_cols,"grid:gridall");
  memory->create(grid_local,size_array_rows,"grid:grid_local");
  array = gridall;
}
/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeGrid::memory_usage()
{
  double nbytes = size_array_rows*size_array_cols * 
    sizeof(double);                             // grid
  nbytes += size_array_rows*sizeof(int); // grid_local
  return nbytes;
}
