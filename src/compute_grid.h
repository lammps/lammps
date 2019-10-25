/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
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
  virtual ~ComputeGrid();
  void init();
  void setup();
  virtual void compute_array() = 0;

  double memory_usage();

 protected:
  int nx, ny, nz;                      // grid dimensions
  int ngrid;                           // number of grid points
  int nvalues;                         // number of values per grid point
  double **grid;                       // global grid
  double **gridall;                    // global grid summed over procs
  int triclinic;                       // triclinic flag
  double *boxlo, *prd;                 // box info (units real/ortho or reduced/tri)
  double *sublo, *subhi;               // subdomain info (units real/ortho or reduced/tri)
  double delxinv,delyinv,delzinv;      // inverse grid spacing
  double delx,dely,delz;               // grid spacing
  int nargbase;                        // number of base class args 
  double cutmax;                       // largest cutoff distance
  int size_array_cols_base;            // number of columns used for coords, etc.
  int *grid_local;                     // local flag for each grid point
  void allocate();
  void grid2x(int, double*);           // convert grid point to coord
  void assign_grid_coords();           // assign coords for grid
  void assign_grid_local();            // set local flag for each grid point
  int check_grid_local(int);           // check if grid point is local
 private:
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
