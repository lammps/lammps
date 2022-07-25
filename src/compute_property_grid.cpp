/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_property_grid.h"

#include "domain.h"
#include "error.h"
#include "grid2d.h"
#include "grid3d.h"
#include "memory.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

enum { ID, X, Y, Z, XS, YS, ZS, XC, YC, ZC, XSC, YSC, ZSC };

#define DELTA 10000

/* ---------------------------------------------------------------------- */

ComputePropertyGrid::ComputePropertyGrid(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), pack_choice(nullptr), 
  grid2d(nullptr), grid3d(nullptr), 
  vec2d(nullptr), array2d(nullptr), vec3d(nullptr), array3d(nullptr)
{
  if (narg < 7) error->all(FLERR, "Illegal compute property/grid command");

  pergrid_flag = 1;

  dimension = domain->dimension;

  nx = utils::inumeric(FLERR,arg[3],false,lmp);
  ny = utils::inumeric(FLERR,arg[4],false,lmp);
  nz = utils::inumeric(FLERR,arg[5],false,lmp);

  if (dimension == 2 && nz != 1) 
    error->all(FLERR,"Compute property/grid for 2d requires nz = 1");

  if (nx <= 0 || ny <= 0 || nz <= 0) 
    error->all(FLERR, "Illegal compute property/grid command");
 
  nvalues = narg - 6;
  pack_choice = new FnPtrPack[nvalues];

  for (int iarg = 6; iarg < narg; iarg++) {
    int jarg = iarg - 6;

    if (strcmp(arg[iarg], "id") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_id;

    } else if (strcmp(arg[iarg], "ix") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_ix;
    } else if (strcmp(arg[iarg], "iy") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_iy;
    } else if (strcmp(arg[iarg], "iz") == 0) {
      if (dimension == 2) 
        error->all(FLERR,"Compute property/grid for 2d cannot use z coord");
      pack_choice[jarg] = &ComputePropertyGrid::pack_iz;

    } else if (strcmp(arg[iarg], "x") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_x;
    } else if (strcmp(arg[iarg], "y") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_y;
    } else if (strcmp(arg[iarg], "z") == 0) {
      if (dimension == 2) 
        error->all(FLERR,"Compute property/grid for 2d cannot use z coord");
      pack_choice[jarg] = &ComputePropertyGrid::pack_z;

    } else if (strcmp(arg[iarg], "xs") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_xs;
    } else if (strcmp(arg[iarg], "ys") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_ys;
    } else if (strcmp(arg[iarg], "zs") == 0) {
      if (dimension == 2) 
        error->all(FLERR,"Compute property/grid for 2d cannot use z coord");
      pack_choice[jarg] = &ComputePropertyGrid::pack_zs;

    } else if (strcmp(arg[iarg], "xc") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_xc;
    } else if (strcmp(arg[iarg], "yc") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_yc;
    } else if (strcmp(arg[iarg], "zc") == 0) {
      if (dimension == 2) 
        error->all(FLERR,"Compute property/grid for 2d cannot use z coord");
      pack_choice[jarg] = &ComputePropertyGrid::pack_zc;

    } else if (strcmp(arg[iarg], "xsc") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_xsc;
    } else if (strcmp(arg[iarg], "ysc") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_ysc;
    } else if (strcmp(arg[iarg], "zsc") == 0) {
      if (dimension == 2) 
        error->all(FLERR,"Compute property/grid for 2d cannot use z coord");
      pack_choice[jarg] = &ComputePropertyGrid::pack_zsc;

    } else error->all(FLERR, "Illegal compute property/grid command");
  }

  // instantiate the Grid class and allocate per-grid memory
  // NOTE: need new memory create methods for 2d

  if (dimension == 2) {
    
  } else {
    grid3d = new Grid3d(lmp, world, nx, ny, nz, 0, 0.0, 0.0,
                        nxlo_in, nxhi_in, nylo_in, nyhi_in, nzlo_in, nzhi_in, 
                        nxlo_out, nxhi_out, nylo_out, nyhi_out, 
                        nzlo_out, nzhi_out);
    if (nvalues == 1)
      memory->create3d_offset(vec3d, nzlo_out, nzhi_out, nylo_out, 
                              nyhi_out, nxlo_out,
                              nxhi_out, "property/grid:vec3d");
    else
      memory->create4d_offset_last(array3d, nzlo_out, nzhi_out, nylo_out, 
                                   nyhi_out, nxlo_out,
                                   nxhi_out, nvalues, "property/grid:array3d");
  }
}

/* ---------------------------------------------------------------------- */

ComputePropertyGrid::~ComputePropertyGrid()
{
  delete[] pack_choice;

  delete grid2d;
  delete grid3d;
  //memory->destroy2d_offset(vec2d);
  //memory->destroy2d_offset(array2d);
  memory->destroy3d_offset(vec3d,nzlo_out,nylo_out,nxlo_out);
  memory->destroy4d_offset_last(array3d,nzlo_out,nylo_out,nxlo_out);
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::compute_pergrid()
{
  invoked_pergrid = update->ntimestep;

  // set current size for portion of grid on each proc
  // may change between compute invocations due to load balancing
  
  if (dimension == 2)
    grid2d->query_bounds(nxlo_in,nxhi_in,nylo_in,nyhi_in);
  else
    grid3d->query_bounds(nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in);

  // reallocate data vector or array if changed




  // fill data vector or array with values for my grid pts

  if (nvalues == 1) {
    (this->*pack_choice[0])(0);
  } else {
    for (int n = 0; n < nvalues; n++) (this->*pack_choice[n])(n);
  }
}

/* ----------------------------------------------------------------------
   return index of grid associated with name
   this class can store M named grids, indexed 0 to M-1
   also set dim for 2d vs 3d grid
------------------------------------------------------------------------- */

int ComputePropertyGrid::get_grid_by_name(char *name, int &dim)
{
  if (strcmp(name,"grid") == 0) {
    dim = dimension;
    return 0;
  }

  return -1;
}

/* ----------------------------------------------------------------------
   return ptr to Grid data struct for grid with index
   this class can store M named grids, indexed 0 to M-1
------------------------------------------------------------------------- */

void *ComputePropertyGrid::get_grid_by_index(int index)
{
  if (index == 0) {
    if (dimension == 2) return grid2d;
    else return grid3d;
  }
  return nullptr;
}

/* ----------------------------------------------------------------------
   return index of data associated with name in grid with index igrid
   this class can store M named grids, indexed 0 to M-1
   each grid can store G named data sets, indexed 0 to G-1
   a data set name can be associated with multiple grids
   also set ncol for data set, 0 = vector, 1-N for array with N columns
   vector = single value per grid pt, array = N values per grid pt
------------------------------------------------------------------------- */

int ComputePropertyGrid::get_griddata_by_name(int igrid, char *name, int &ncol)
{
  if (igrid == 0 && strcmp(name,"data") == 0) {
    if (nvalues == 1) ncol = 0;
    else ncol = nvalues;
    return 0;
  }

  return -1;
}

/* ----------------------------------------------------------------------
   return ptr to multidim data array associated with index
   this class can store G named data sets, indexed 0 to M-1
------------------------------------------------------------------------- */

void *ComputePropertyGrid::get_griddata_by_index(int index)
{
  if (index == 0) {
    if (dimension == 2) {
      if (nvalues == 1) return vec2d;
      else return array2d;
    } else {
      if (nvalues == 1) return vec3d;
      else return array3d;
    }
  }
  return nullptr;
}

/* ----------------------------------------------------------------------
   memory usage of grid data
------------------------------------------------------------------------- */

double ComputePropertyGrid::memory_usage()
{
  double bytes = 0.0;
  //double bytes = (double) nmax * nvalues * sizeof(double);
  //bytes += (double) nmax * 2 * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   one method for every keyword compute property/grid can output
------------------------------------------------------------------------- */

void ComputePropertyGrid::pack_id(int n)
{
  if (dimension == 2) {
    if (nvalues == 0) {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          vec2d[iy][ix] = iy*nx + ix + 1;
    } else {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          array2d[iy][ix][n] = iy*nx + ix + 1;
    }
  } else if (dimension == 3) {
    if (nvalues == 0) {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d[iz][iy][ix] = iz*ny*nx + iy*nx + ix + 1;
    } else {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            array3d[iz][iy][ix][n] = iz*ny*nx + iy*nx + ix + 1;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_ix(int n)
{
  if (dimension == 2) {
    if (nvalues == 0) {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          vec2d[iy][ix] = ix + 1;
    } else {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          array2d[iy][ix][n] = ix + 1;
    }
  } else if (dimension == 3) {
    if (nvalues == 0) {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d[iz][iy][ix] = ix + 1;
    } else {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            array3d[iz][iy][ix][n] = ix + 1;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_iy(int n)
{
  if (dimension == 2) {
    if (nvalues == 0) {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          vec2d[iy][ix] = iy + 1;
    } else {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          array2d[iy][ix][n] = iy + 1;
    }
  } else if (dimension == 3) {
    if (nvalues == 0) {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d[iz][iy][ix] = iy + 1;
    } else {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            array3d[iz][iy][ix][n] = iy + 1;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_iz(int n)
{
  if (nvalues == 0) {
    for (int iz = nzlo_in; iz <= nzhi_in; iz++)
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          vec3d[iz][iy][ix] = iz + 1;
  } else {
    for (int iz = nzlo_in; iz <= nzhi_in; iz++)
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          array3d[iz][iy][ix][n] = iz + 1;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_x(int n)
{
  double boxlo,dx;

  if (dimension == 2) {
    grid2d->query_box(0,boxlo,dx);
    if (nvalues == 0) {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          vec2d[iy][ix] = boxlo + ix*dx;
    } else {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          array2d[iy][ix][n] = boxlo + ix*dx;
    }
  } else if (dimension == 3) {
    grid3d->query_box(0,boxlo,dx);
    if (nvalues == 0) {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d[iz][iy][ix] = boxlo + ix*dx;
    } else {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            array3d[iz][iy][ix][n] = boxlo + ix*dx;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_y(int n)
{
  double boxlo,dy;

  if (dimension == 2) {
    grid2d->query_box(1,boxlo,dy);
    if (nvalues == 0) {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          vec2d[iy][ix] = boxlo + iy*dy;
    } else {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          array2d[iy][ix][n] = boxlo + iy*dy;
    }
  } else if (dimension == 3) {
    grid3d->query_box(1,boxlo,dy);
    if (nvalues == 0) {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d[iz][iy][ix] = boxlo + iy*dy;
    } else {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            array3d[iz][iy][ix][n] = boxlo + iy*dy;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_z(int n)
{
  double boxlo,dz;
  grid3d->query_box(2,boxlo,dz);

  if (nvalues == 0) {
    for (int iz = nzlo_in; iz <= nzhi_in; iz++)
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          vec3d[iz][iy][ix] = boxlo + iz*dz;
  } else {
    for (int iz = nzlo_in; iz <= nzhi_in; iz++)
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          array3d[iz][iy][ix][n] = boxlo + iz*dz;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_xs(int n)
{
  double dx = 1.0/nx;

  if (dimension == 2) {
    if (nvalues == 0) {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          vec2d[iy][ix] = ix*dx;
    } else {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          array2d[iy][ix][n] = ix*dx;
    }
  } else if (dimension == 3) {
    if (nvalues == 0) {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d[iz][iy][ix] = ix*dx;
    } else {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            array3d[iz][iy][ix][n] = ix*dx;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_ys(int n)
{
  double dy = 1.0/ny;

  if (dimension == 2) {
    if (nvalues == 0) {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          vec2d[iy][ix] = iy*dy;
    } else {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          array2d[iy][ix][n] = iy*dy;
    }
  } else if (dimension == 3) {
    if (nvalues == 0) {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d[iz][iy][ix] = iy*dy;
    } else {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            array3d[iz][iy][ix][n] = iy*dy;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_zs(int n)
{
  double dz = 1.0/nz;

  if (nvalues == 0) {
    for (int iz = nzlo_in; iz <= nzhi_in; iz++)
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          vec3d[iz][iy][ix] = iz*dz;
  } else {
    for (int iz = nzlo_in; iz <= nzhi_in; iz++)
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          array3d[iz][iy][ix][n] = iz*dz;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_xc(int n)
{
  double boxlo,dx;

  if (dimension == 2) {
    grid2d->query_box(0,boxlo,dx);
    if (nvalues == 0) {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          vec2d[iy][ix] = boxlo + (ix+0.5)*dx;
    } else {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          array2d[iy][ix][n] = boxlo + (ix+0.5)*dx;
    }
  } else if (dimension == 3) {
    grid3d->query_box(0,boxlo,dx);
    if (nvalues == 0) {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d[iz][iy][ix] = boxlo + (ix+0.5)*dx;
    } else {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            array3d[iz][iy][ix][n] = boxlo + (ix+0.5)*dx;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_yc(int n)
{
  double boxlo,dy;

  if (dimension == 2) {
    grid2d->query_box(1,boxlo,dy);
    if (nvalues == 0) {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          vec2d[iy][ix] = boxlo + (iy+0.5)*dy;
    } else {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          array2d[iy][ix][n] = boxlo + (iy+0.5)*dy;
    }
  } else if (dimension == 3) {
    grid3d->query_box(1,boxlo,dy);
    if (nvalues == 0) {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d[iz][iy][ix] = boxlo + (iy+0.5)*dy;
    } else {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            array3d[iz][iy][ix][n] = boxlo + (iy+0.5)*dy;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_zc(int n)
{
  double boxlo,dz;
  grid3d->query_box(2,boxlo,dz);

  if (nvalues == 0) {
    for (int iz = nzlo_in; iz <= nzhi_in; iz++)
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          vec3d[iz][iy][ix] = boxlo + (iz+0.5)*dz;
  } else {
    for (int iz = nzlo_in; iz <= nzhi_in; iz++)
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          array3d[iz][iy][ix][n] = boxlo + (iz+0.5)*dz;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_xsc(int n)
{
  double dx = 1.0/nx;

  if (dimension == 2) {
    if (nvalues == 0) {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          vec2d[iy][ix] = (ix+0.5)*dx;
    } else {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          array2d[iy][ix][n] = (ix+0.5)*dx;
    }
  } else if (dimension == 3) {
    if (nvalues == 0) {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d[iz][iy][ix] = (ix+0.5)*dx;
    } else {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            array3d[iz][iy][ix][n] = (ix+0.5)*dx;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_ysc(int n)
{
  double dy = 1.0/ny;

  if (dimension == 2) {
    if (nvalues == 0) {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          vec2d[iy][ix] = (iy+0.5)*dy;
    } else {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          array2d[iy][ix][n] = (iy+0.5)*dy;
    }
  } else if (dimension == 3) {
    if (nvalues == 0) {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d[iz][iy][ix] = (iy+0.5)*dy;
    } else {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            array3d[iz][iy][ix][n] = (iy+0.5)*dy;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_zsc(int n)
{
  double dz = 1.0/nz;

  if (nvalues == 0) {
    for (int iz = nzlo_in; iz <= nzhi_in; iz++)
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          vec3d[iz][iy][ix] = (iz+0.5)*dz;
  } else {
    for (int iz = nzlo_in; iz <= nzhi_in; iz++)
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++)
          array3d[iz][iy][ix][n] = (iz+0.5)*dz;
  }
}
