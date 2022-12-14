/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_COMPUTE_GRID_LOCAL_H
#define LMP_COMPUTE_GRID_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeGridLocal : public Compute {
 public:
  ComputeGridLocal(class LAMMPS *, int, char **);
  ~ComputeGridLocal() override;
  void setup() override;
  void compute_local() override = 0;

  double memory_usage() override;

 protected:
  int nx, ny, nz;                            // global grid dimensions
  int nxlo, nxhi, nylo, nyhi, nzlo, nzhi;    // local grid bounds, inclusive
  int nvalues;                               // number of values per grid point
  double **alocal;                           // pointer to Compute::array_local
  int triclinic;                             // triclinic flag
  double *boxlo, *prd;                       // box info (units real/ortho or reduced/tri)
  double *sublo, *subhi;                     // subdomain info (units real/ortho or reduced/tri)
  double delxinv, delyinv, delzinv;          // inverse grid spacing
  double delx, dely, delz;                   // grid spacing
  int nargbase;                              // number of base class args
  double cutmax;                             // largest cutoff distance
  int size_local_cols_base;                  // number of columns used for coords, etc.
  int gridlocal_allocated;                   // shows if gridlocal allocated

  void allocate();                             // create arrays
  void deallocate();                           // free arrays
  void grid2x(int, int, int, double *);        // convert global indices to coordinates
  void grid2lamda(int, int, int, double *);    // convert global indices to lamda coordinates
  void set_grid_global();                      // set global grid
  void set_grid_local();                       // set bounds for local grid
  void assign_coords();                        // assign coords for grid
};

}    // namespace LAMMPS_NS

#endif
