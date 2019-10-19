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
  extarray = 1;

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

  size_array_rows = nx*ny*nz;
}

/* ---------------------------------------------------------------------- */

ComputeGrid::~ComputeGrid()
{
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
  } else {
    prd = domain->prd_lamda;
    boxlo = domain->boxlo_lamda;
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

  // sufficient conditions for stencil bounding rcut

  // require |delz*mz|^2 <= rcut^2
  // require |dely*my|^2 <= rcut^2 + |delyz*mz_max|^2
  // require |delx*mx|^2 <= rcut^2 + |delxz*mz_max|^2 + |delxy*my_max|^2

  double delxy = domain->xy/ny;
  double delxz = domain->xz/nz;
  double delyz = domain->yz/nz;
  
  if (!triclinic) {
    mz = cutmax*delzinv + 1;
    my = sqrt(cutmax*cutmax + pow(delyz*mz,2))*delyinv + 1;
    mx = sqrt(cutmax*cutmax + pow(delxz*mz,2) 
                            + pow(delxy*my,2))*delxinv + 1;
  } else {
    double delxinvtmp = nx/domain->xprd;
    double delyinvtmp = ny/domain->yprd;
    double delzinvtmp = nz/domain->zprd;
    mz = cutmax*delzinvtmp + 1;
    my = sqrt(cutmax*cutmax + pow(delyz*mz,2))*delyinvtmp + 1;
    mx = sqrt(cutmax*cutmax + pow(delxz*mz,2) 
                            + pow(delxy*my,2))*delxinvtmp + 1;
  }

  printf("mx = %d\n",mx);
  printf("my = %d\n",my);
  printf("mz = %d\n",mz);

  // size global grid to accomodate periodic interactions
  
  nxfull = nx + 2*mx;
  nyfull = ny + 2*my;
  nzfull = nz + 2*mz;
  nxyfull = nxfull * nyfull;

  printf("nxfull = %d\n",nxfull);
  printf("nyfull = %d\n",nyfull);
  printf("nzfull = %d\n",nzfull);

  x0full = boxlo[0] - mx*delx;
  y0full = boxlo[1] - my*dely;
  z0full = boxlo[2] - mz*delz;

  allocate();
}

/* ----------------------------------------------------------------------
   convert grid index to box coords
------------------------------------------------------------------------- */

void ComputeGrid::igridfull2x(int igrid, double *x)
{
  int iz = igrid / nxyfull;
  igrid -= iz*nxyfull;
  int iy = igrid / nxfull;
  igrid -= igrid*nxfull;
  int ix = igrid;

  x[0] = x0full+ix*delx;
  x[1] = y0full+iy*dely;
  x[2] = z0full+iz*delz;

  if (triclinic) domain->lamda2x(x, x);

}

/* ----------------------------------------------------------------------
   gather global array from full grid
------------------------------------------------------------------------- */

void ComputeGrid::gather_global_array()
{
  int iarray;
  memset(&array[0][0],0,size_array_rows*size_array_cols*sizeof(double));

  for (int igrid = 0; igrid < ngridfull; igrid++) {

    // inefficient, should exploit shared ix structure

    iarray = igridfull2iarray(igrid);
    for (int icol = 0; icol < size_array_cols; icol++)
      array[iarray][icol] += gridfull[igrid][icol];
  }
}

/* ----------------------------------------------------------------------
   convert full grid index to compute array index
   inefficient, should exploit shared ix structure
------------------------------------------------------------------------- */

int ComputeGrid::igridfull2iarray(int igrid)
{
  int iz = igrid / nxyfull;
  igrid -= iz*nxyfull;
  int iy = igrid / nxfull;
  igrid -= igrid*nxfull;
  int ix = igrid;

  ix -= mx;
  iy -= my;
  iz -= mz;

  while (ix < 0) ix += nx;
  while (iy < 0) iy += ny;
  while (iz < 0) iz += nz;

  while (ix >= nx) ix -= nx;
  while (iy >= ny) iy -= ny;
  while (iz >= nz) iz -= nz;

  int iarray = (iz * ny + iy) * nx + ix;

  return iarray;
}

/* ----------------------------------------------------------------------
   free and reallocate arrays
------------------------------------------------------------------------- */

void ComputeGrid::allocate()
{
  ngridfull = nxfull*nyfull*nzfull;

  // grow global array if necessary

  memory->destroy(array);
  memory->create(array,size_array_rows,size_array_cols,"sna/grid:array");
  memory->create(gridfull,ngridfull,size_array_cols,"sna/grid:gridfull");
}
/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeGrid::memory_usage()
{
  int nbytes = 0;
  return nbytes;
}
