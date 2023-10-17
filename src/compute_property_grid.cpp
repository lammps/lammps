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

#include "compute_property_grid.h"

#include "comm.h"
#include "domain.h"
#include "error.h"
#include "grid2d.h"
#include "grid3d.h"
#include "memory.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

enum { LOW, CTR };
enum { UNSCALED, SCALED };

#define DELTA 10000

/* ---------------------------------------------------------------------- */

ComputePropertyGrid::ComputePropertyGrid(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), grid2d(nullptr), grid3d(nullptr), vec2d(nullptr), vec3d(nullptr),
    array2d(nullptr), array3d(nullptr), pack_choice(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "compute property/grid", error);

  pergrid_flag = 1;

  dimension = domain->dimension;

  nxgrid = utils::inumeric(FLERR, arg[3], false, lmp);
  nygrid = utils::inumeric(FLERR, arg[4], false, lmp);
  nzgrid = utils::inumeric(FLERR, arg[5], false, lmp);

  if (dimension == 2 && nzgrid != 1)
    error->all(FLERR, "Compute property/grid for 2d requires nz = 1");

  if (nxgrid <= 0 || nygrid <= 0 || nzgrid <= 0)
    error->all(FLERR, "All compute property/grid grid counts must be > 0");

  nvalues = narg - 6;
  pack_choice = new FnPtrPack[nvalues];

  for (int iarg = 6; iarg < narg; iarg++) {
    int jarg = iarg - 6;

    if (strcmp(arg[iarg], "id") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_id;
    } else if (strcmp(arg[iarg], "proc") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_proc;

    } else if (strcmp(arg[iarg], "ix") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_indices<0>;
    } else if (strcmp(arg[iarg], "iy") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_indices<1>;
    } else if (strcmp(arg[iarg], "iz") == 0) {
      if (dimension == 2) error->all(FLERR, "Compute property/grid for 2d cannot use z coord");
      pack_choice[jarg] = &ComputePropertyGrid::pack_indices<2>;

    } else if (strcmp(arg[iarg], "x") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_coords<LOW, UNSCALED, 0>;
    } else if (strcmp(arg[iarg], "y") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_coords<LOW, UNSCALED, 1>;
    } else if (strcmp(arg[iarg], "z") == 0) {
      if (dimension == 2) error->all(FLERR, "Compute property/grid for 2d cannot use z coord");
      pack_choice[jarg] = &ComputePropertyGrid::pack_coords<LOW, UNSCALED, 2>;

    } else if (strcmp(arg[iarg], "xs") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_coords<LOW, SCALED, 0>;
    } else if (strcmp(arg[iarg], "ys") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_coords<LOW, SCALED, 1>;
    } else if (strcmp(arg[iarg], "zs") == 0) {
      if (dimension == 2) error->all(FLERR, "Compute property/grid for 2d cannot use z coord");
      pack_choice[jarg] = &ComputePropertyGrid::pack_coords<LOW, SCALED, 2>;

    } else if (strcmp(arg[iarg], "xc") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_coords<CTR, UNSCALED, 0>;
    } else if (strcmp(arg[iarg], "yc") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_coords<CTR, UNSCALED, 1>;
    } else if (strcmp(arg[iarg], "zc") == 0) {
      if (dimension == 2) error->all(FLERR, "Compute property/grid for 2d cannot use z coord");
      pack_choice[jarg] = &ComputePropertyGrid::pack_coords<CTR, UNSCALED, 2>;

    } else if (strcmp(arg[iarg], "xsc") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_coords<CTR, SCALED, 0>;
    } else if (strcmp(arg[iarg], "ysc") == 0) {
      pack_choice[jarg] = &ComputePropertyGrid::pack_coords<CTR, SCALED, 1>;
    } else if (strcmp(arg[iarg], "zsc") == 0) {
      if (dimension == 2) error->all(FLERR, "Compute property/grid for 2d cannot use z coord");
      pack_choice[jarg] = &ComputePropertyGrid::pack_coords<CTR, SCALED, 2>;

    } else
      error->all(FLERR, "Unknown compute property/grid keyword: {}", arg[iarg]);
  }

  // initial setup of distributed grid

  allocate_grid();
}

/* ---------------------------------------------------------------------- */

ComputePropertyGrid::~ComputePropertyGrid()
{
  delete[] pack_choice;

  deallocate_grid();
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::init()
{
  triclinic = domain->triclinic;
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::compute_pergrid()
{
  invoked_pergrid = update->ntimestep;

  // fill data vector or array with values for my grid pts

  if (nvalues == 1) {
    (this->*pack_choice[0])(0);
  } else {
    for (int n = 0; n < nvalues; n++) (this->*pack_choice[n])(n);
  }
}

/* ----------------------------------------------------------------------
   subset of grid assigned to each proc may have changed
   deallocate and reallocate Grid class and local data structs
   called by load balancer when proc subdomains are adjusted
---------------------------------------------------------------------- */

void ComputePropertyGrid::reset_grid()
{
  deallocate_grid();
  allocate_grid();
}

/* ----------------------------------------------------------------------
   return index of grid associated with name
   this class can store M named grids, indexed 0 to M-1
   also set dim for 2d vs 3d grid
   return -1 if grid name not found
------------------------------------------------------------------------- */

int ComputePropertyGrid::get_grid_by_name(const std::string &name, int &dim)
{
  if (name == "grid") {
    dim = dimension;
    return 0;
  }

  return -1;
}

/* ----------------------------------------------------------------------
   return ptr to Grid data struct for grid with index
   this class can store M named grids, indexed 0 to M-1
   return nullptr if index is invalid
------------------------------------------------------------------------- */

void *ComputePropertyGrid::get_grid_by_index(int index)
{
  if (index == 0) {
    if (dimension == 2)
      return grid2d;
    else
      return grid3d;
  }

  return nullptr;
}

/* ----------------------------------------------------------------------
   return index of data associated with name in grid with index igrid
   this class can store M named grids, indexed 0 to M-1
   each grid can store G named data sets, indexed 0 to G-1
     a data set name can be associated with multiple grids
   set ncol for data set, 0 = vector, 1-N for array with N columns
     vector = single value per grid pt, array = N values per grid pt
   return -1 if data name not found
------------------------------------------------------------------------- */

int ComputePropertyGrid::get_griddata_by_name(int igrid, const std::string &name, int &ncol)
{
  if ((igrid == 0) && (name == "data")) {
    if (nvalues == 1)
      ncol = 0;
    else
      ncol = nvalues;
    return 0;
  }

  return -1;
}

/* ----------------------------------------------------------------------
   return ptr to multidim data array associated with index
   this class can store G named data sets, indexed 0 to M-1
   return nullptr if index is invalid
------------------------------------------------------------------------- */

void *ComputePropertyGrid::get_griddata_by_index(int index)
{
  if (index == 0) {
    if (dimension == 2) {
      if (nvalues == 1)
        return vec2d;
      else
        return array2d;
    } else {
      if (nvalues == 1)
        return vec3d;
      else
        return array3d;
    }
  }

  return nullptr;
}

/* ----------------------------------------------------------------------
   instantiate the Grid class and allocate local per-grid memory
---------------------------------------------------------------------- */

void ComputePropertyGrid::allocate_grid()
{
  if (dimension == 2) {
    grid2d = new Grid2d(lmp, world, nxgrid, nygrid);
    grid2d->setup_grid(nxlo_in, nxhi_in, nylo_in, nyhi_in, nxlo_out, nxhi_out, nylo_out, nyhi_out);

    if (nvalues == 1)
      memory->create2d_offset(vec2d, nylo_out, nyhi_out, nxlo_out, nxhi_out, "property/grid:vec2d");
    else
      memory->create3d_offset_last(array2d, nylo_out, nyhi_out, nxlo_out, nxhi_out, nvalues,
                                   "property/grid:array2d");
    ngridout = (nxhi_out - nxlo_out + 1) * (nyhi_out - nylo_out + 1);

  } else {
    grid3d = new Grid3d(lmp, world, nxgrid, nygrid, nzgrid);
    grid3d->setup_grid(nxlo_in, nxhi_in, nylo_in, nyhi_in, nzlo_in, nzhi_in,
                       nxlo_out, nxhi_out, nylo_out, nyhi_out, nzlo_out, nzhi_out);
    if (nvalues == 1)
      memory->create3d_offset(vec3d, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                              "property/grid:vec3d");
    else
      memory->create4d_offset_last(array3d, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out,
                                   nxhi_out, nvalues, "property/grid:array3d");
    ngridout = (nxhi_out - nxlo_out + 1) * (nyhi_out - nylo_out + 1) * (nzhi_out - nzlo_out + 1);
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::deallocate_grid()
{
  delete grid2d;
  delete grid3d;
  memory->destroy2d_offset(vec2d, nylo_out, nxlo_out);
  memory->destroy3d_offset_last(array2d, nylo_out, nxlo_out);
  memory->destroy3d_offset(vec3d, nzlo_out, nylo_out, nxlo_out);
  memory->destroy4d_offset_last(array3d, nzlo_out, nylo_out, nxlo_out);
}

/* ----------------------------------------------------------------------
   memory usage of grid data
------------------------------------------------------------------------- */

double ComputePropertyGrid::memory_usage()
{
  double bytes = (double) ngridout * nvalues * sizeof(double);
  return bytes;
}

// ----------------------------------------------------------------------
// pack methods for all values
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   grid point IDs
------------------------------------------------------------------------- */

void ComputePropertyGrid::pack_id(int n)
{
  if (dimension == 2) {
    if (nvalues == 1) {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++) vec2d[iy][ix] = iy * nxgrid + ix + 1;
    } else {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++) array2d[iy][ix][n] = iy * nxgrid + ix + 1;
    }
  } else if (dimension == 3) {
    if (nvalues == 1) {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d[iz][iy][ix] = iz * nygrid * nxgrid + iy * nxgrid + ix + 1;
    } else {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            array3d[iz][iy][ix][n] = iz * nygrid * nxgrid + iy * nxgrid + ix + 1;
    }
  }
}

/* ----------------------------------------------------------------------
   grid point owning processor
------------------------------------------------------------------------- */

void ComputePropertyGrid::pack_proc(int n)
{
  int me = comm->me;

  if (dimension == 2) {
    if (nvalues == 1) {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++) vec2d[iy][ix] = me;
    } else {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++) array2d[iy][ix][n] = me;
    }
  } else if (dimension == 3) {
    if (nvalues == 1) {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            vec3d[iz][iy][ix] = me;
    } else {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            array3d[iz][iy][ix][n] = me;
    }
  }
}

/* ----------------------------------------------------------------------
   grid indices via templating
------------------------------------------------------------------------- */

template <int IDIM> void ComputePropertyGrid::pack_indices(int n)
{
  if (dimension == 2) {
    if (nvalues == 1) {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
          if (IDIM == 0) vec2d[iy][ix] = ix + 1;
          if (IDIM == 1) vec2d[iy][ix] = iy + 1;
        }
    } else {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
          if (IDIM == 0) array2d[iy][ix][n] = ix + 1;
          if (IDIM == 1) array2d[iy][ix][n] = iy + 1;
        }
    }

  } else if (dimension == 3) {
    if (nvalues == 1) {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
            if (IDIM == 0) vec3d[iz][iy][ix] = ix + 1;
            if (IDIM == 1) vec3d[iz][iy][ix] = iy + 1;
            if (IDIM == 2) vec3d[iz][iy][ix] = iz + 1;
          }
    } else {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
            if (IDIM == 0) array3d[iz][iy][ix][n] = ix + 1;
            if (IDIM == 1) array3d[iz][iy][ix][n] = iy + 1;
            if (IDIM == 2) array3d[iz][iy][ix][n] = iz + 1;
          }
    }
  }
}

/* ----------------------------------------------------------------------
   grid point coords
   LOW/CTR, SCALED/UNSCALED, DIM = 0/1 via templating
   dimension and orthogonal/tricilic in code logic
------------------------------------------------------------------------- */

template <int POS, int MODE, int IDIM> void ComputePropertyGrid::pack_coords(int n)
{
  double boxlo, delta;
  double lamda[3], xone[3];

  // 2d grid

  if (dimension == 2) {

    // for coords which are orthogonal OR scaled

    if (!triclinic || MODE == SCALED) {

      if (MODE == UNSCALED) {
        boxlo = domain->boxlo[IDIM];
        if (IDIM == 0) delta = domain->prd[IDIM] / nxgrid;
        if (IDIM == 1) delta = domain->prd[IDIM] / nygrid;
      }
      if (MODE == SCALED) {
        boxlo = 0.0;
        if (IDIM == 0) delta = 1.0 / nxgrid;
        if (IDIM == 1) delta = 1.0 / nygrid;
      }

      if (nvalues == 1) {
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
            if (POS == LOW) {
              if (IDIM == 0) vec2d[iy][ix] = boxlo + ix * delta;
              if (IDIM == 1) vec2d[iy][ix] = boxlo + iy * delta;
            }
            if (POS == CTR) {
              if (IDIM == 0) vec2d[iy][ix] = boxlo + (ix + 0.5) * delta;
              if (IDIM == 1) vec2d[iy][ix] = boxlo + (iy + 0.5) * delta;
            }
          }

      } else {
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
            if (POS == LOW) {
              if (IDIM == 0) array2d[iy][ix][n] = boxlo + ix * delta;
              if (IDIM == 1) array2d[iy][ix][n] = boxlo + iy * delta;
            }
            if (POS == CTR) {
              if (IDIM == 0) array2d[iy][ix][n] = boxlo + (ix + 0.5) * delta;
              if (IDIM == 1) array2d[iy][ix][n] = boxlo + (iy + 0.5) * delta;
            }
          }
      }

    // only for coords which are triclinic AND unscaled

    } else {

      double dx = 1.0 / nxgrid;
      double dy = 1.0 / nygrid;
      lamda[2] = 0.0;

      if (nvalues == 1) {
        for (int iy = nylo_in; iy <= nyhi_in; iy++) {
          if (POS == LOW) lamda[1] = iy * dy;
          else lamda[1] = (iy + 0.5) * dy;
          for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
            if (POS == LOW) lamda[0] = ix * dx;
            else lamda[0] = (ix + 0.5) * dx;
            domain->lamda2x(lamda, xone);
            if (IDIM == 0) vec2d[iy][ix] = xone[0];
            if (IDIM == 1) vec2d[iy][ix] = xone[1];
          }
        }

      } else {
        for (int iy = nylo_in; iy <= nyhi_in; iy++) {
          if (POS == LOW) lamda[1] = iy * dy;
          else lamda[1] = (iy + 0.5) * dy;
          for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
            if (POS == LOW) lamda[0] = ix * dx;
            else lamda[0] = (ix + 0.5) * dx;
            domain->lamda2x(lamda, xone);
            if (IDIM == 0) array2d[iy][ix][n] = xone[0];
            if (IDIM == 1) array2d[iy][ix][n] = xone[1];
          }
        }
      }
    }

  // 3d grid

  } else if (dimension == 3) {

    // for coords which are orthogonal OR scaled

    if (!triclinic || MODE == SCALED) {

      if (MODE == UNSCALED) {
        boxlo = domain->boxlo[IDIM];
        if (IDIM == 0) delta = domain->prd[IDIM] / nxgrid;
        if (IDIM == 1) delta = domain->prd[IDIM] / nygrid;
        if (IDIM == 2) delta = domain->prd[IDIM] / nzgrid;
      }
      if (MODE == SCALED) {
        boxlo = 0.0;
        if (IDIM == 0) delta = 1.0 / nxgrid;
        if (IDIM == 1) delta = 1.0 / nygrid;
        if (IDIM == 2) delta = 1.0 / nzgrid;
      }

      if (nvalues == 1) {
        for (int iz = nzlo_in; iz <= nzhi_in; iz++)
          for (int iy = nylo_in; iy <= nyhi_in; iy++)
            for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
              if (POS == LOW) {
                if (IDIM == 0) vec3d[iz][iy][ix] = boxlo + ix * delta;
                if (IDIM == 1) vec3d[iz][iy][ix] = boxlo + iy * delta;
                if (IDIM == 2) vec3d[iz][iy][ix] = boxlo + iz * delta;
              }
              if (POS == CTR) {
                if (IDIM == 0) vec3d[iz][iy][ix] = boxlo + (ix + 0.5) * delta;
                if (IDIM == 1) vec3d[iz][iy][ix] = boxlo + (iy + 0.5) * delta;
                if (IDIM == 2) vec3d[iz][iy][ix] = boxlo + (iz + 0.5) * delta;
              }
            }

      } else {
        for (int iz = nzlo_in; iz <= nzhi_in; iz++)
          for (int iy = nylo_in; iy <= nyhi_in; iy++)
            for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
              if (POS == LOW) {
                if (IDIM == 0) array3d[iz][iy][ix][n] = boxlo + ix * delta;
                if (IDIM == 1) array3d[iz][iy][ix][n] = boxlo + iy * delta;
                if (IDIM == 2) array3d[iz][iy][ix][n] = boxlo + iz * delta;
              }
              if (POS == CTR) {
                if (IDIM == 0) array3d[iz][iy][ix][n] = boxlo + (ix + 0.5) * delta;
                if (IDIM == 1) array3d[iz][iy][ix][n] = boxlo + (iy + 0.5) * delta;
                if (IDIM == 2) array3d[iz][iy][ix][n] = boxlo + (iz + 0.5) * delta;
              }
            }
      }

    // only for coords which are triclinic AND unscaled

    } else {

      double dx = 1.0 / nxgrid;
      double dy = 1.0 / nygrid;
      double dz = 1.0 / nzgrid;

      if (nvalues == 1) {
        for (int iz = nzlo_in; iz <= nzhi_in; iz++) {
          if (POS == LOW) lamda[2] = iz * dz;
          else lamda[2] = (iz + 0.5) * dz;
          for (int iy = nylo_in; iy <= nyhi_in; iy++) {
            if (POS == LOW) lamda[1] = iy * dy;
            else lamda[1] = (iy + 0.5) * dy;
            for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
              if (POS == LOW) lamda[0] = ix * dx;
              else lamda[0] = (ix + 0.5) * dx;
              domain->lamda2x(lamda, xone);
              if (IDIM == 0) vec3d[iz][iy][ix] = xone[0];
              if (IDIM == 1) vec3d[iz][iy][ix] = xone[1];
              if (IDIM == 2) vec3d[iz][iy][ix] = xone[2];
            }
          }
        }

      } else {
        for (int iz = nzlo_in; iz <= nzhi_in; iz++) {
          if (POS == LOW) lamda[2] = iz * dz;
          else lamda[2] = (iz + 0.5) * dz;
          for (int iy = nylo_in; iy <= nyhi_in; iy++) {
            if (POS == LOW) lamda[1] = iy * dy;
            else lamda[1] = (iy + 0.5) * dy;
            for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
              if (POS == LOW) lamda[0] = ix * dx;
              else lamda[0] = (ix + 0.5) * dx;
              domain->lamda2x(lamda, xone);
              if (IDIM == 0) array3d[iz][iy][ix][n] = xone[0];
              if (IDIM == 1) array3d[iz][iy][ix][n] = xone[1];
              if (IDIM == 2) array3d[iz][iy][ix][n] = xone[2];
            }
          }
        }
      }
    }
  }
}
