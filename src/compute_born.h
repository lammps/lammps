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
  Contributing Authors : Germain Clavier (UCA)
  --------------------------------------------------------------------------*/

#ifdef COMPUTE_CLASS

ComputeStyle(born,ComputeBorn)

#else

#ifndef LMP_COMPUTE_BORN_H
#define LMP_COMPUTE_BORN_H

#include "compute.h"

namespace LAMMPS_NS {

  class ComputeBorn : public Compute {
  public:
    ComputeBorn(class LAMMPS *, int, char **);
    virtual ~ComputeBorn();
    void init();
    void init_list(int, class NeighList *);
    void compute_vector();

  private:

    void compute_pairs();
    void compute_bonds();
    void compute_angles();
    void compute_dihedrals();
    void compute_impropers();

    int me,nvalues;
    int *which;

    int pairflag, bondflag, angleflag;
    int dihedflag, impflag;

    int const albe[21][2] = {
      {0,0},  // C11
      {1,1},  // C22
      {2,2},  // C33
      {1,2},  // C44
      {0,2},  // C55
      {0,1},  // C66
      {0,1},  // C12
      {0,2},  // C13
      {0,3},  // C14
      {0,4},  // C15
      {0,5},  // C16
      {1,2},  // C23
      {1,3},  // C24
      {1,4},  // C25
      {1,5},  // C26
      {2,3},  // C34
      {2,4},  // C35
      {2,5},  // C36
      {3,4},  // C45
      {3,5},  // C46
      {4,5}   // C56
    };

    int const albemunu[21][4] = {
      {0,0,0,0},  // C11
      {1,1,1,1},  // C22
      {2,2,2,2},  // C33
      {1,2,1,2},  // C44
      {0,2,0,2},  // C55
      {0,1,0,1},  // C66
      {0,0,1,1},  // C12
      {0,0,2,2},  // C13
      {0,0,1,2},  // C14
      {0,0,0,2},  // C15
      {0,0,0,1},  // C16
      {1,1,2,2},  // C23
      {1,1,1,2},  // C24
      {1,1,0,2},  // C25
      {1,1,0,1},  // C26
      {2,2,1,2},  // C34
      {2,2,0,2},  // C35
      {2,2,0,1},  // C36
      {1,2,0,2},  // C45
      {1,2,0,1},  // C46
      {0,1,0,2}   // C56
    };

    double *values_local,*values_global;
    double pos,pos1,dt,nktv2p,ftm2v;
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

 E: ... style does not support compute born

 Some component of the force field (pair, bond, angle...) does not provide
 a function to return the Born term contribution.
 */

