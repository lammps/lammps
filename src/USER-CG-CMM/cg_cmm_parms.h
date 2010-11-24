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

/* ----------------------------------------------------------------------
   Common parameters for the CMM coarse grained MD potentials.
   Contributing author: Axel Kohlmeyer (Temple U) <akohlmey@gmail.com>
------------------------------------------------------------------------- */

#ifndef LMP_CG_CMM_PARMS_H
#define LMP_CG_CMM_PARMS_H

namespace LAMMPS_NS {

  class CGCMMParms {
    public:

    // CG type flags. list of supported LJ exponent combinations
    enum {CG_NOT_SET=0, CG_LJ9_6, CG_LJ12_4, CG_LJ12_6, NUM_CG_TYPES,
          CG_COUL_NONE, CG_COUL_CUT, CG_COUL_DEBYE, CG_COUL_LONG};

    int find_cg_type(const char *);

    protected:
    // coarse grain flags
    static const char * const cg_type_list[NUM_CG_TYPES];
    static const double cg_prefact[NUM_CG_TYPES];
    static const double cg_pow1[NUM_CG_TYPES] ;
    static const double cg_pow2[NUM_CG_TYPES] ;

  };

}
#endif
