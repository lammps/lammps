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
  int nxfull, nyfull, nzfull;          // grid dimensions with ghost points
  int nxyfull;                         // nx_full*ny_full 
  int ngridfull;                       // number of full grid points
  double **gridfull;                   // full grid points
  int mx, my, mz;                      // cutmax stencil dimensions
  int triclinic;                       // triclinic flag
  double *boxlo, *prd;                 // box info (units real/ortho or reduced/tri)
  double delxinv,delyinv,delzinv;      // inverse grid spacing
  double delx,dely,delz;               // grid spacing
  double x0full, y0full, z0full;       // origin of full grid
  int nargbase;                        // number of base class args 
  double cutmax;                       // largest cutoff distance
  virtual void allocate();
  void igridfull2x(int, double*);      // convert full grid point to coord
  void gather_global_array();          // gather global array from full grid
  void copy_local_grid();              // copy local grid to global array
  int igridfull2iarray(int);           // convert full grid index to compute array index
  int iarray2igridfull(int);           // convert compute array index to full grid index

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
