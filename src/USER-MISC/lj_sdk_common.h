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
#if 0

/* ----------------------------------------------------------------------
   Common data for the Shinoda, DeVane, Klein (SDK) coars grain model
   Contributing author: Axel Kohlmeyer (Temple U) <akohlmey@gmail.com>
------------------------------------------------------------------------- */

#ifndef LMP_LJ_SDK_COMMON_H
#define LMP_LJ_SDK_COMMON_H

#include "string.h"

namespace LAMMPS_NS {

  namespace LJSDKParms {

    // CG type flags. list of supported LJ exponent combinations
    enum {CG_NOT_SET=0, CG_LJ9_6, CG_LJ12_4, CG_LJ12_6, NUM_CG_TYPES};

    static int find_cg_type(const char *label) {
      for (int i=0; i < NUM_CG_TYPES; ++i) {
	if (strcmp(label,cg_type_list[i]) == 0) {
	  return i;
	}
      }
      return CG_NOT_SET;
    }

    static const char * const cg_type_list[] = {"none", "lj9_6", "lj12_4", "lj12_6"};
    static const double cg_prefact[] = {0.0, 6.75,  2.59807621135332, 4.0};
    static const double cg_pow1[]    = {0.0, 9.00, 12.0,             12.0};
    static const double cg_pow2[]    = {0.0, 6.00,  4.0,              6.0};

  }

  class PairPairLJSDKCommon : public Pair {
   public:

    PairLJSDKCommon(class LAMMPS *);
    virtual ~PairPairLJSDKCommon();

    virtual void settings(int, char **);
    virtual void coeff(int, char **);
    virtual void init_style();
    virtual void init_list(int, class NeighList *);
    virtual double init_one(int, int);

    virtual void write_restart(FILE *);
    virtual void read_restart(FILE *);
    virtual void write_restart_settings(FILE *);
    virtual void read_restart_settings(FILE *);

    virtual double memory_usage();

   protected:

    // coarse grain flags
    int **cg_type;
    
    // lennard jones parameters
    double cut_lj_global, **cutsq, **cut_lj;
    double **epsilon, **sigma;
    double **lj1, **lj2, **lj3, **lj4, **offset;

    // coulomb parameters
    int allocated_coul; // 0/1 = whether coulomb arrays are allocated
    double cut_coul_global, cut_coulsq_global, g_ewald;

    // tables
    double tabinnersq;
    double *rtable,*drtable,*ftable,*dftable,*ctable,*dctable;
    double *etable,*detable,*ptable,*dptable,*vtable,*dvtable;
    int ncoulshiftbits,ncoulmask;
    
    // methods
    virtual void allocate();

   private:

    // disable default constructor
    PairLJSDKCommon();

}
#endif


#endif
