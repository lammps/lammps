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

#ifdef PAIR_CLASS

PairStyle(list,PairList)

#else

#ifndef LMP_PAIR_LIST_H
#define LMP_PAIR_LIST_H

#include "pair.h"

namespace LAMMPS_NS {

class PairList : public Pair {
 public:
  PairList(class LAMMPS *);
  virtual ~PairList();

  virtual void compute(int, int);
  virtual void settings(int, char **);
  virtual void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);
  virtual double memory_usage();

 protected:
  void allocate();

  enum { NONE=0, HARM, MORSE, LJ126 };

  // potential specific parameters
  struct harm_p  { double k, r0;          };
  struct morse_p { double d0, alpha, r0;  };
  struct lj126_p { double epsilon, sigma; };

  union parm_u { 
    struct harm_p harm;
    struct morse_p morse;
    struct lj126_p lj126;
  };

  typedef struct {
    tagint id1,id2;        // global atom ids
    double cutsq;       // cutoff**2 for this pair
    double offset;      // energy offset
    union parm_u parm;  // parameters for style
  } list_parm_t;    

 protected:
  double cut_global;    // global cutoff distance
  int *style;           // list of styles for pair interactions
  list_parm_t *params;  // lisf of pair interaction parameters
  int npairs;           // # of atom pairs in global list
  int check_flag;       // 1 if checking for missing pairs
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Not all pairs processed in pair_style list

Not all interacting pairs for which coefficients were found. This can be intentional
and then you need to set the 'nocheck' option. If not, it usually means that the 
communication cutoff is too small. This can be ameliorated by either increasing
the cutoff in the pair_style command or the communication cutoff.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot open pair list file

Self-explanatory.  The file with the list of pairs cannot be open for reading.
Check the path and permissions.

E: Incorrectly formatted ...

Self-explanatory.  The content of the pair list file does not match the documented
format. Please re-read the documentation and carefully compare it to your file.

E: Unknown pair list potential style

Self-explanatory.  You requested a potential type that is not yet implemented or have a typo.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style list requires atom IDs

Self-explanatory.  The pairs in the list are identified via atom IDs, so they need to be present.

E: Pair style list requires an atom map

Self-explanatory.  Atoms are looked up via an atom map. Create one using the atom_style map command.

*/
