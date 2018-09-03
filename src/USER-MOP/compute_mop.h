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

/*------------------------------------------------------------------------
  Contributing Authors : Romain Vermorel (LFCR), Laurent Joly (ULyon)
  --------------------------------------------------------------------------*/

#ifdef COMPUTE_CLASS

ComputeStyle(mop,ComputeMop)

#else

#ifndef LMP_COMPUTE_MOP_H
#define LMP_COMPUTE_MOP_H

#include "compute.h"

namespace LAMMPS_NS {

  class ComputeMop : public Compute {
  public:
    ComputeMop(class LAMMPS *, int, char **);
    virtual ~ComputeMop();
    void init();
    void init_list(int, class NeighList *);
    void compute_vector();

  private:

    void compute_pairs();

    int me,nvalues,dir;
    int *which;

    double *values_local,*values_global;
    double pos,pos1,dt,nktv2p,ftm2v;
    double area;
    class NeighList *list;

  };

}

#endif
#endif

/* ERROR/WARNING messages:

   E: Illegal ... command

   Self-explanatory.  Check the input script syntax and compare to the
   documentation for the command.  You can use -echo screen as a
   command-line option when running LAMMPS to see the offending line.

   E: Pair style does not support compute mop

   The pair style does not have a single() function, so it can
   not be invoked by compute mop/spatial.



*/

