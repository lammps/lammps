/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_COMPUTE_GRID_H
#define LMP_COMPUTE_GRID_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeGrid : public Compute {
 public:
  ComputeGrid(class LAMMPS *, int, char **);
  ~ComputeGrid() override;
  void setup() override;
  void compute_array() override = 0;

  double memory_usage() override;

 protected:
  int nx, ny, nz;                            // global grid dimensions
  int nxlo, nxhi, nylo, nyhi, nzlo, nzhi;    // local grid bounds, inclusive
  int ngridlocal;                            // number of local grid points
  int nvalues;                               // number of values per grid point
  double **grid;                             // global grid
  double **gridall;                          // global grid summed over procs
  double ****gridlocal;                      // local grid
  int triclinic;                             // triclinic flag
  double *boxlo, *prd;                       // box info (units real/ortho or reduced/tri)
  double *sublo, *subhi;                     // subdomain info (units real/ortho or reduced/tri)
  double delxinv, delyinv, delzinv;          // inverse grid spacing
  double delx, dely, delz;                   // grid spacing
  int nargbase;                              // number of base class args
  double cutmax;                             // largest cutoff distance
  int size_array_cols_base;                  // number of columns used for coords, etc.
  int gridlocal_allocated;                   // shows if gridlocal allocated

  void allocate();               // create arrays
  void deallocate();             // free arrays
  void grid2x(int, double *);    // convert grid point to coord
  void assign_coords_all();      // assign coords for global grid
  void set_grid_global();        // set global grid
  void set_grid_local();         // set bounds for local grid
};

}    // namespace LAMMPS_NS

#endif
