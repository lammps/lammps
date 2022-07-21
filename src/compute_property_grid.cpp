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
#include "memory.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

enum { ID, X, Y, Z, XS, YS, ZS, XC, YC, ZC, XSC, YSC, ZSC };

#define DELTA 10000

/* ---------------------------------------------------------------------- */

ComputePropertyGrid::ComputePropertyGrid(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), pack_choice(nullptr)
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
    if (strcmp(arg[iarg], "id") == 0) {
      pack_choice[iarg] = &ComputePropertyGrid::pack_id;

    } else if (strcmp(arg[iarg], "x") == 0) {
      pack_choice[iarg] = &ComputePropertyGrid::pack_x;
    } else if (strcmp(arg[iarg], "y") == 0) {
      pack_choice[iarg] = &ComputePropertyGrid::pack_y;
    } else if (strcmp(arg[iarg], "z") == 0) {
      if (dimension == 2) 
        error->all(FLERR,"Compute property/grid for 2d cannot use z coord");
      pack_choice[iarg] = &ComputePropertyGrid::pack_z;

    } else if (strcmp(arg[iarg], "xs") == 0) {
      pack_choice[iarg] = &ComputePropertyGrid::pack_xs;
    } else if (strcmp(arg[iarg], "ys") == 0) {
      pack_choice[iarg] = &ComputePropertyGrid::pack_ys;
    } else if (strcmp(arg[iarg], "zs") == 0) {
      if (dimension == 2) 
        error->all(FLERR,"Compute property/grid for 2d cannot use z coord");
      pack_choice[iarg] = &ComputePropertyGrid::pack_zs;

    } else if (strcmp(arg[iarg], "xc") == 0) {
      pack_choice[iarg] = &ComputePropertyGrid::pack_xc;
    } else if (strcmp(arg[iarg], "yc") == 0) {
      pack_choice[iarg] = &ComputePropertyGrid::pack_yc;
    } else if (strcmp(arg[iarg], "zc") == 0) {
      if (dimension == 2) 
        error->all(FLERR,"Compute property/grid for 2d cannot use z coord");
      pack_choice[iarg] = &ComputePropertyGrid::pack_zc;

    } else if (strcmp(arg[iarg], "xsc") == 0) {
      pack_choice[iarg] = &ComputePropertyGrid::pack_xsc;
    } else if (strcmp(arg[iarg], "ysc") == 0) {
      pack_choice[iarg] = &ComputePropertyGrid::pack_ysc;
    } else if (strcmp(arg[iarg], "zsc") == 0) {
      if (dimension == 2) 
        error->all(FLERR,"Compute property/grid for 2d cannot use z coord");
      pack_choice[iarg] = &ComputePropertyGrid::pack_zsc;

    } else error->all(FLERR, "Illegal compute property/grid command");
  }

  
}

/* ---------------------------------------------------------------------- */

ComputePropertyGrid::~ComputePropertyGrid()
{
  delete[] pack_choice;
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::compute_pergrid()
{
  invoked_pergrid = update->ntimestep;

  // fill vector or array with local values

  /*
  if (nvalues == 1) {
    buf = vlocal;
    (this->*pack_choice[0])(0);
  } else {
    if (alocal) buf = &alocal[0][0];
    for (int n = 0; n < nvalues; n++) (this->*pack_choice[n])(n);
  }
  */
}

/* ----------------------------------------------------------------------
   return index of grid associated with name
   this class can store M named grids, indexed 0 to M-1
   also set dim for 2d vs 3d grid
------------------------------------------------------------------------- */

int ComputePropertyGrid::get_grid_by_name(char *name, int &dim)
{
  if (strcmp(name,"grid") == 0) {
    dim = 3;
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
  if (index == 0) return gc;
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
    ncol = nvalues;
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
  if (index == 0) return data;
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
   packed into buf starting at n with stride nvalues
------------------------------------------------------------------------- */

void ComputePropertyGrid::pack_id(int n)
{
  int i;

  /*
  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    buf[n] = tag[i];
    n += nvalues;
  }
  */
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_x(int n)
{
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_y(int n)
{
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_z(int n)
{
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_xs(int n)
{
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_ys(int n)
{
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_zs(int n)
{
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_xc(int n)
{
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_yc(int n)
{
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_zc(int n)
{
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_xsc(int n)
{
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_ysc(int n)
{
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_zsc(int n)
{
}
