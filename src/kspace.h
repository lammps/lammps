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

#ifndef LMP_KSPACE_H
#define LMP_KSPACE_H

#include "pointers.h"

namespace LAMMPS_NS {

class KSpace : protected Pointers {
  friend class ThrOMP;
 public:
  double energy;
  double virial[6];
  double g_ewald;
  int nx_pppm,ny_pppm,nz_pppm;
 
  KSpace(class LAMMPS *, int, char **);
  virtual ~KSpace() {}
  void modify_params(int, char **);
  void *extract(const char *);

  virtual void init() = 0;
  virtual void setup() = 0;
  virtual void compute(int, int) = 0;
  virtual void timing(int, double &, double &) {}
  virtual double memory_usage() {return 0.0;}

 protected:
  double slab_volfactor;
  int gridflag,gewaldflag;
  int order;
  int slabflag;
  double scale;
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Bad kspace_modify slab parameter

Kspace_modify value for the slab/volume keyword must be >= 2.0.

W: Kspace_modify slab param < 2.0 may cause unphysical behavior

The kspace_modify slab parameter should be larger to insure periodic
grids padded with empty space do not overlap.

*/
