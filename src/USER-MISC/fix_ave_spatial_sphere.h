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

/* Adapted from fix ave/spatial by Niall Jackson <niall.jackson@gmail.com>*/

#ifdef FIX_CLASS

FixStyle(ave/spatial/sphere,FixAveSpatialSphere)

#else

#ifndef LMP_FIX_AVE_SPATIAL_SPHERE_H
#define LMP_FIX_AVE_SPATIAL_SPHERE_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixAveSpatialSphere : public Fix {
 public:
  FixAveSpatialSphere(class LAMMPS *, int, char **);
  ~FixAveSpatialSphere();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  double memory_usage();
  double compute_array(int,int);
  void reset_timestep(bigint);
  
 private:
  int me; //rank of the MPI Process
  int nfreq; //fix nfreq parameter
  int nrepeat; //fix nrepeat parameter
  bigint nvalid; //the next timestep on which I'll be invoked
  
  //controls for various optional behaviours
  int normflag; //what type of normalisation to do
  int norm;
  int scaleflag; //default units (lattice, box, etc.)
  double scale; //coordinate scaling factor
  int regionflag; //restricted to a particular region?
  char *idregion; //name of the region to use
  class Region *region; //pointer to the region
  FILE *fp; //pointer for the output file
  long int filepos; //file position pointer
  int ave; //averaging mode
  int nwindow; //number of averaging windows
  int iwindow; //current window
  int window_limit; 
  int overwrite; //continuously overwrite the output (requires ave running)
  
  //used to keep track of which per-atom values we deal with
  int nvalues; //number of variables to average
  int *argindex; //if the variable is an array, this is the offset in that array
  int *value2index;
  int *which; //what sort of variable is arg i? Variable, compute, etc?
  char **ids; //names of the variables
  int maxvar; //current size of the varatom array
  double* varatom; //contains the peratom values of a variable
  
  //details of the sphere and the bins
  int maxbin; //current number of bins in memory (relevant if box changes)
  int nbins; //number of spherical bins
  int maxatom; //current size of the bin array
  int *bin; //stores the bin of each atom
  double *coord; //values of r at the mid points of the bins
  double *binvol; //volumes of the bins
  double *count_one, *count_many, *count_sum, *count_total; //bin populations
  double **values_one, **values_many, **values_sum, **values_total; //accumulated bin values
  double origin[3]; //origin coordinates of the sphere
  int origin_type[3]; //are origin coordinates constant or variable?
  char *origin_ids[3]; //store the names of variables used to access the origin
  int origin_index[3]; //indices for the origin variables
  int origin_val2idx[3]; //compute/variable indices
  double r_min, r_minsq; //minimum radius, and its square
  double r_max, r_maxsq; //maximum radius, and its square
  double deltar, inv_deltar; //radial width of a bin, and its inverse (and their squares)
  
  void setup_bins(); //create the bin arrays
  void set_bin_volumes(); //calculate the volume of each bin
  void bin_atoms(); //put the atom into bins based on their coordinates
  bigint nextvalid(); //return the next timestep on which this is invoked
  
  //NEED TO BE CATEGORISED
  int irepeat;
  double **count_list;
  double ***values_list;
  
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix ave/spatial/spherical does not exist

Self-explanatory.

E: Cannot open fix ave/spatial/spherical file %s

The specified file cannot be opened.  Check that the path and name are
correct.

E: Compute ID for fix ave/spatial/spherical does not exist

Self-explanatory.

E: Fix ave/spatial/spherical compute does not calculate per-atom values

A compute used by fix ave/spatial/spherical must generate per-atom values.

E: Fix ave/spatial/spherical compute does not calculate a per-atom vector

A compute used by fix ave/spatial/spherical must generate per-atom values.

E: Fix ave/spatial/spherical compute does not calculate a per-atom array

Self-explanatory.

E: Fix ave/spatial/spherical compute vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix ID for fix ave/spatial/spherical does not exist

Self-explanatory.

E: Fix ave/spatial/spherical fix does not calculate per-atom values

A fix used by fix ave/spatial/spherical must generate per-atom values.

E: Fix ave/spatial/spherical fix does not calculate a per-atom vector

A fix used by fix ave/spatial/spherical must generate per-atom values.

E: Fix ave/spatial/spherical fix does not calculate a per-atom array

Self-explanatory.

E: Fix ave/spatial/spherical fix vector is accessed out-of-range

The index for the vector is out of bounds.

E: Variable name for fix ave/spatial/spherical does not exist

Self-explanatory.

E: Fix ave/spatial/spherical variable is not atom-style variable

A variable used by fix ave/spatial/spherical must generate per-atom values.

E: Fix ave/spatial/spherical for triclinic boxes requires units reduced

Self-explanatory.

E: Fix ave/spatial/spherical requires reduced units if the box changes size.

If the box size changes, only the units reduced option can be
used.

E: Fix for fix ave/spatial/spherical not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/spatial/spherical is
requesting a value on a non-allowed timestep.

E: Fix ave/spatial/spherical missed timestep

You cannot reset the timestep to a value beyond where the fix
expects to next perform averaging.

*/
