/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */

/*------------------------------------------------------------------------
  Contributing Authors : Germain Clavier (TUe)
  --------------------------------------------------------------------------*/

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(born/matrix,ComputeBornMatrix);
// clang-format on
#else

#ifndef LMP_COMPUTE_BORN_MATRIX_H
#define LMP_COMPUTE_BORN_MATRIX_H

#include "compute.h"

namespace LAMMPS_NS {

  class ComputeBornMatrix : public Compute {
  public:
    ComputeBornMatrix(class LAMMPS *, int, char **);
    virtual ~ComputeBornMatrix();
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
    int dihedflag, impflag, kspaceflag;

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

 E: ... style does not support compute born/matrix

 Some component of the force field (pair, bond, angle...) does not provide
 a function to return the Born term contribution.
 */

