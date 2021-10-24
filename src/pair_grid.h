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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(grid,PairGrid);
// clang-format on
#else

#ifndef LMP_PAIR_GRID_H
#define LMP_PAIR_GRID_H

#include "pair.h"

namespace LAMMPS_NS {

class PairGrid : public Pair {
 public:
  PairGrid(class LAMMPS *);
  virtual ~PairGrid();
  virtual void init_style(){};
  void setup();
  virtual void compute(int, int) {
    printf("DANGER! This function should always be overridden by child\n");
  };

  void settings(int, char **);
  virtual void coeff(int, char **);
  virtual double init_one(int, int);
  double memory_usage();
  
 protected:
  int nx, ny, nz;                      // global grid dimensions
  int nxlo, nxhi, nylo, nyhi, nzlo, nzhi; // local grid bounds, inclusive
  int ngridlocal;                      // number of local grid points
  int nvalues;                         // number of values per grid point
  double ****gridlocal;                // local grid
  double **alocal;                     // pointer to Compute::array_local
  int triclinic;                       // triclinic flag
  double *boxlo, *prd;                 // box info (units real/ortho or reduced/tri)
  double *sublo, *subhi;               // subdomain info (units real/ortho or reduced/tri)
  double delxinv,delyinv,delzinv;      // inverse grid spacing
  double delx,dely,delz;               // grid spacing
  int nargbase;                        // number of base class args 
  double cutmax;                       // largest cutoff distance
  int ndesc;                           // number of descriptors
  int ndesc_base;                      // number of columns used for coords, etc.
  int gridlocal_allocated;             // shows if gridlocal allocated
  double **beta;                       // betas for all local grid points in list
  int beta_max;                        // length of beta

  void allocate();                     // allocate pairstyle arrays
  void allocate_grid();                // create grid arrays
  void deallocate_grid();              // free grid arrays
  void grid2x(int, int, int, double*); // convert global indices to coordinates
  void set_grid_global();              // set global grid
  void set_grid_local();               // set bounds for local grid
  void assign_coords();                // assign coords for grid
  void copy_gridlocal_to_local_array();// copy 4d gridlocal array to 2d local array
  void compute_beta();                 // get betas from someplace
 private:
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
