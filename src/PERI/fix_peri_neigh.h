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

#ifdef FIX_CLASS

FixStyle(PERI_NEIGH,FixPeriNeigh)

#else

#ifndef LMP_FIX_PERI_NEIGH_H
#define LMP_FIX_PERI_NEIGH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPeriNeigh : public Fix {
  friend class PairPeriPMB;
  friend class PairPeriPMBOMP;
  friend class PairPeriLPS;
  friend class PairPeriVES;   //NEW
  friend class PairPeriLPSOMP;
  friend class ComputeDamageAtom;

 public:
  FixPeriNeigh(class LAMMPS *,int, char **);
  virtual ~FixPeriNeigh();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void min_setup(int);

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  void write_restart(FILE *);
  void restart(char *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);


 protected:
  int first;                 // flag for first time initialization
  int maxpartner;            // max # of peridynamic neighs for any atom
  int *npartner;             // # of neighbors for each atom
  int **partner;             // neighs for each atom, stored as global IDs
  double **deviatorextention; // Deviatoric extention     
  double **deviatorBackextention; // Deviatoric back extention 
  double **r0;               // initial distance to partners
  double **r1;               // Instanteneous distance to partners *** NEW ***
  double *thetaOld;          // Dilatation Old one
  double *vinter;            // sum of vfrac for bonded neighbors
  double *wvolume;           // weighted volume of particle
  int isPMB;
  int isLPS;
  int isVES;

  class NeighList *list;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Duplicate particle in PeriDynamic bond - simulation box is too small

This is likely because your box length is shorter than 2 times
the bond length.

*/
