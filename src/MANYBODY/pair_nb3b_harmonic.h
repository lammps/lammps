/* --*- c++ -*- ---------------------------------------------------------
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

PairStyle(nb3b/harmonic,PairNb3bHarmonic)

#else

#ifndef LMP_PAIR_NB3B_HARMONIC_H
#define LMP_PAIR_NB3B_HARMONIC_H

#include "pair.h"

namespace LAMMPS_NS {

class PairNb3bHarmonic : public Pair {
 public:
  PairNb3bHarmonic(class LAMMPS *);
  virtual ~PairNb3bHarmonic();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();

 protected:
  struct Param {
    double k_theta, theta0, cutoff;
    double cut,cutsq;
    int ielement,jelement,kelement;
  };
  
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  int ***elem2param;            // mapping from element triplets to parameters
  int *map;                     // mapping from atom types to elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  Param *params;                // parameter set for an I-J-K interaction

  void allocate();
  void read_file(char *);
  void setup();
  void twobody(Param *, double, double &, int, double &);
  void threebody(Param *, Param *, Param *, double, double, double *, double *,
		 double *, double *, int, double &);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style COMB requires atom IDs

This is a requirement to use the AIREBO potential.

E: Pair style COMB requires newton pair on

See the newton command.  This is a restriction to use the COMB
potential.

E: Pair style COMB requires atom attribute q

Self-explanatory.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open COMB potential file %s

The specified COMB potential file cannot be opened.  Check that the
path and name are correct.

E: Incorrect format in COMB potential file

Incorrect number of words per line in the potential file.

E: Illegal COMB parameter

One or more of the coefficients defined in the potential file is
invalid.

E: Potential file has duplicate entry

The potential file for a SW or Tersoff potential has more than
one entry for the same 3 ordered elements.

E: Potential file is missing an entry

The potential file for a SW or Tersoff potential does not have a
needed entry.

W: Pair COMB charge %.10f with force %.10f hit min barrier

Something is possibly wrong with your model.

W: Pair COMB charge %.10f with force %.10f hit max barrier

Something is possibly wrong with your model.

E: Neighbor list overflow, boost neigh_modify one

There are too many neighbors of a single atom.  Use the neigh_modify
command to increase the max number of neighbors allowed for one atom.
You may also want to boost the page size.

*/
