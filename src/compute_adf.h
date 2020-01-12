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

#ifdef COMPUTE_CLASS

ComputeStyle(adf,ComputeADF)

#else

#ifndef LMP_COMPUTE_ADF_H
#define LMP_COMPUTE_ADF_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeADF : public Compute {
 public:
  ComputeADF(class LAMMPS *, int, char **);
  ~ComputeADF();
  void init();
  void init_list(int, class NeighList *);
  void compute_array();

 private:
  int nbin;                            // # of adf bins
  int ntriples;                        // # of adf triples
  double deltax,deltaxinv;             // bin width and inverse-width
  int *ilo,*ihi,*jlo,*jhi,*klo,*khi;
  double **hist;         // histogram bins
  double **histall;      // summed histogram bins across all procs

  double *rcutinnerj, *rcutinnerk; // list of inner cutoffs
  double *rcutouterj, *rcutouterk; // list of outer cutoffs

  class NeighList *list;   // full neighbor list

  int *iatomcount;               // local number of central atoms
  int *iatomcountall;            // total number of central atoms
  int **iatomflag;               // 1 if type is central atom in ADF
  int *maxjatom, *maxkatom;      // allocated size jatom, katom neighlist
  int *numjatom, *numkatom;      // actual size of jatom, katom neighlist
  int **neighjatom, **neighkatom;// list of short neighbor lists
  int **jatomflag, **katomflag;  // 1 if type is neighbor atom in ADF

  int *maxjkatom;                // allocated size short neighlist
  int *numjkatom;                // actual size of short neighlist
  int **neighjkatom;             // list of short neighbor lists
  int **bothjkatom;              // 1 if atom is in both jatom and katom lists
  double ***delrjkatom;          // list of 4-vectors: delx, dely, delx, and 1/r

  double rad2deg;        // conversion factor from radians to degrees
  int ordinate_style;    // DEGREE, RADIAN, or COSINE
  int cutflag;          // 1 if at least one outer cutoff specified
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute adf requires a pair style be defined or outer cutoff specified

Self-explanatory.

E: Compute adf cutoff exceeds ghost atom range - use comm_modify cutoff command

Self-explanatary.

*/
