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

ComputeStyle(mop/profile,ComputeMopProfile)

#else

#ifndef LMP_COMPUTE_MOP_PROFILE_H
#define LMP_COMPUTE_MOP_PROFILE_H

#include "compute.h"

namespace LAMMPS_NS {

  class ComputeMopProfile : public Compute {
  public:
    ComputeMopProfile(class LAMMPS *, int, char **);
    virtual ~ComputeMopProfile();
    void init();
    void init_list(int, class NeighList *);
    void compute_array();

  private:

    void compute_pairs();
    void setup_bins();

    int me,nvalues,dir;
    int *which;

    int originflag;
    double origin,delta,offset,invdelta;
    int nbins;
    double **coord,**coordp;
    double **values_local,**values_global;

    int ndim;
    double dt,nktv2p,ftm2v;
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

   E: Pair style does not support compute mop/profile

   The pair style does not have a single() function, so it can
   not be invoked by compute mop/profile.



*/

