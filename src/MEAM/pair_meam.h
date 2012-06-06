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

#ifdef PAIR_CLASS

PairStyle(meam,PairMEAM)

#else

#ifndef LMP_PAIR_MEAM_H
#define LMP_PAIR_MEAM_H

extern "C" {
  void meam_setup_global_(int *, int *, double *, int *, double *, double *,
                         double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *,
                         double *, double *, int *);
  void meam_setup_param_(int *, double *, int *, int *, int *);
  void meam_setup_done_(double *);

  void meam_dens_init_(int *, int *, int *, int *, int *,
                       double *, int *, int *, int *, int *,
                       double *, double *, double *, double *,
                       double *, double *,
                       double *, double *, double *, double *, double *,
                       int *);

  void meam_dens_final_(int *, int *, int *, int *, int *, double *, double *,
                        int *, int *, int *,
                        double *, double *, double *, double *,
                        double *, double *, double *,
                        double *, double *, double *, double *,
                        double *, double *,
                        double *, double *, double *, double *, int *);

  void meam_force_(int *, int *, int *, int *, int *, int *,
                   double *, double *, int *, int *, int *,
                   double *, int *, int *, int *, int *, double *, double *,
                   double *, double *, double *, double *, double *, double *,
                   double *, double *, double *, double *, double *, double *,
                   double *, double *, double *, double *, double *, double *, int *);

  void meam_cleanup_();
}


#include "pair.h"

namespace LAMMPS_NS {

class PairMEAM : public Pair {
 public:
  PairMEAM(class LAMMPS *);
  ~PairMEAM();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);

  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  double *mass;                 // mass of each element

  int *map;                     // mapping from atom types to elements
  int *fmap;                    // Fortran version of map array for MEAM lib

  int maxneigh;
  double *scrfcn,*dscrfcn,*fcpair;

  int nmax;
  double *rho,*rho0,*rho1,*rho2,*rho3,*frhop;
  double *gamma,*dgamma1,*dgamma2,*dgamma3,*arho2b;
  double **arho1,**arho2,**arho3,**arho3b,**t_ave,**tsq_ave;

  void allocate();
  void read_files(char *, char *);
  void neigh_strip(int, int *, int *, int **);
  void neigh_f2c(int, int *, int *, int **);
  void neigh_c2f(int, int *, int *, int **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: MEAM library error %d

A call to the MEAM Fortran library returned an error.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style MEAM requires newton pair on

See the newton command.  This is a restriction to use the MEAM
potential.

E: Cannot open MEAM potential file %s

The specified MEAM potential file cannot be opened.  Check that the
path and name are correct.

E: Incorrect format in MEAM potential file

Incorrect number of words per line in the potential file.

E: Unrecognized lattice type in MEAM file 1

The lattice type in an entry of the MEAM library file is not
valid.

E: Did not find all elements in MEAM library file

The requested elements were not found in the MEAM file.

E: Keyword %s in MEAM parameter file not recognized

Self-explanatory.

E: Unrecognized lattice type in MEAM file 2

The lattice type in an entry of the MEAM parameter file is not
valid.

*/
