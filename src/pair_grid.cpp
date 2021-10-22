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

#include "pair_grid.h"
#include <mpi.h>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGrid::PairGrid(LAMMPS *lmp) :
  Pair(lmp), gridlocal(nullptr), alocal(nullptr)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;


  ndesc = 0;
  ngridlocal = 0;

  ndesc_base = 6;
  gridlocal_allocated = 0;
  beta_max = 0;
  beta = nullptr;
}

/* ---------------------------------------------------------------------- */

PairGrid::~PairGrid()
{
  if (copymode) return;

  memory->destroy(beta);

  deallocate_grid();
}

/* ---------------------------------------------------------------------- */

void PairGrid::init()
{
}

/* ---------------------------------------------------------------------- */

void PairGrid::setup()
{
  printf("Hello, world! C\n");  
  deallocate_grid();
  set_grid_global();
  set_grid_local();
  allocate_grid();
  assign_coords();
}

/* ----------------------------------------------------------------------
   convert global array indexes to box coords
------------------------------------------------------------------------- */

void PairGrid::grid2x(int ix, int iy, int iz, double *x)
{
  x[0] = ix*delx;
  x[1] = iy*dely;
  x[2] = iz*delz;

  if (triclinic) domain->lamda2x(x, x);
}

/* ----------------------------------------------------------------------
   create arrays
------------------------------------------------------------------------- */

void PairGrid::allocate_grid()
{
  if (nxlo <= nxhi && nylo <= nyhi && nzlo <= nzhi) {
    gridlocal_allocated = 1;
    memory->create4d_offset(gridlocal,ndesc,nzlo,nzhi,nylo,nyhi,
			    nxlo,nxhi,"pair/grid:gridlocal");
    memory->create(alocal, ngridlocal, ndesc, "pair/grid:alocal");
  }
}

/* ----------------------------------------------------------------------
   free arrays
------------------------------------------------------------------------- */

void PairGrid::deallocate_grid()
{
  if (gridlocal_allocated) {
    gridlocal_allocated = 0;
    memory->destroy4d_offset(gridlocal,nzlo,nylo,nxlo);
    memory->destroy(alocal);
  }
}

/* ----------------------------------------------------------------------
   set global grid
------------------------------------------------------------------------- */

void PairGrid::set_grid_global()
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
}

/* ----------------------------------------------------------------------
   set local subset of grid that I own
   n xyz lo/hi = 3d brick that I own (inclusive)
------------------------------------------------------------------------- */

void PairGrid::set_grid_local()
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

  double xfraclo,xfrachi,yfraclo,yfrachi,zfraclo,zfrachi;

  if (comm->layout != Comm::LAYOUT_TILED) {
    xfraclo = comm->xsplit[comm->myloc[0]];
    xfrachi = comm->xsplit[comm->myloc[0]+1];
    yfraclo = comm->ysplit[comm->myloc[1]];
    yfrachi = comm->ysplit[comm->myloc[1]+1];
    zfraclo = comm->zsplit[comm->myloc[2]];
    zfrachi = comm->zsplit[comm->myloc[2]+1];
  } else {
    xfraclo = comm->mysplit[0][0];
    xfrachi = comm->mysplit[0][1];
    yfraclo = comm->mysplit[1][0];
    yfrachi = comm->mysplit[1][1];
    zfraclo = comm->mysplit[2][0];
    zfrachi = comm->mysplit[2][1];
  }

  nxlo = static_cast<int> (xfraclo * nx);
  if (1.0*nxlo != xfraclo*nx) nxlo++;
  nxhi = static_cast<int> (xfrachi * nx);
  if (1.0*nxhi == xfrachi*nx) nxhi--;

  nylo = static_cast<int> (yfraclo * ny);
  if (1.0*nylo != yfraclo*ny) nylo++;
  nyhi = static_cast<int> (yfrachi * ny);
  if (1.0*nyhi == yfrachi*ny) nyhi--;

  nzlo = static_cast<int> (zfraclo * nz);
  if (1.0*nzlo != zfraclo*nz) nzlo++;
  nzhi = static_cast<int> (zfrachi * nz);
  if (1.0*nzhi == zfrachi*nz) nzhi--;

  ngridlocal = (nxhi - nxlo + 1) * (nyhi - nylo + 1) * (nzhi - nzlo + 1);
}

/* ----------------------------------------------------------------------
   copy coords to local array
------------------------------------------------------------------------- */

void PairGrid::assign_coords()
{
  int igrid = 0;
  for (int iz = nzlo; iz <= nzhi; iz++)
    for (int iy = nylo; iy <= nyhi; iy++)
      for (int ix = nxlo; ix <= nxhi; ix++) {
	alocal[igrid][0] = ix;
	alocal[igrid][1] = iy;
	alocal[igrid][2] = iz;
	double xgrid[3];
	grid2x(ix, iy, iz, xgrid);
	alocal[igrid][3] = xgrid[0];
	alocal[igrid][4] = xgrid[1];
	alocal[igrid][5] = xgrid[2];
	igrid++;
      }
}

/* ----------------------------------------------------------------------
   copy the 4d gridlocal array values to the 2d local array
------------------------------------------------------------------------- */

void PairGrid::copy_gridlocal_to_local_array()
{
  int igrid = 0;
  for (int iz = nzlo; iz <= nzhi; iz++)
    for (int iy = nylo; iy <= nyhi; iy++)
      for (int ix = nxlo; ix <= nxhi; ix++) {
	for (int icol = ndesc_base; icol < ndesc; icol++)
	  alocal[igrid][icol] = gridlocal[icol][iz][iy][ix];
	igrid++;
      }
}

/* ----------------------------------------------------------------------
   get beta from someplace
------------------------------------------------------------------------- */

void PairGrid::compute_beta()
{
  int igrid = 0;
  for (int iz = nzlo; iz <= nzhi; iz++)
    for (int iy = nylo; iy <= nyhi; iy++)
      for (int ix = nxlo; ix <= nxhi; ix++) {
	for (int icol = ndesc_base; icol < ndesc; icol++)
	  beta[igrid][icol] = 1.0;
	igrid++;
      }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairGrid::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGrid::settings(int narg, char ** arg)
{
  if (narg < 4) error->all(FLERR,"Illegal pair style command");
  int iarg0 = 0;
  int iarg = iarg0;
  if (strcmp(arg[iarg],"grid") == 0) {
    if (iarg+4 > narg) error->all(FLERR,"Illegal pair grid command");
    nx = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
    ny = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
    nz = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
    if (nx <= 0 || ny <= 0 || nz <= 0)
      error->all(FLERR,"All grid/local dimensions must be positive");
    iarg += 4;
  } else error->all(FLERR,"Illegal pair grid command");
  nargbase = iarg - iarg0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGrid::coeff(int narg, char **arg)
{
  if (!allocated) allocate();
  if (narg != 2) error->all(FLERR,"Incorrect args for pair coefficients");
 //  if (narg != 2 + atom->ntypes) error->all(FLERR,"Incorrect args for pair coefficients");

  //  map_element2type(narg-4,arg+4);
  //  map_element2type(0,nullptr);

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGrid::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style grid requires newton pair on");

  // no neighbor list

  // int irequest = neighbor->request(this,instance_me);
  // neighbor->requests[irequest]->half = 0;
  // neighbor->requests[irequest]->full = 1;

}

/* ----------------------------------------------------------------------
   return maximum force cut off distance
------------------------------------------------------------------------- */

double PairGrid::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  return cutmax;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double PairGrid::memory_usage()
{
  int nbytes = ndesc*ngridlocal*sizeof(double); // gridlocal
  return nbytes;
}
