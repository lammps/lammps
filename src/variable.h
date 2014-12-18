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

#ifndef LMP_VARIABLE_H
#define LMP_VARIABLE_H

#include "stdlib.h"
#include "pointers.h"

namespace LAMMPS_NS {

class Variable : protected Pointers {
 public:
  Variable(class LAMMPS *);
  ~Variable();
  void set(int, char **);
  void set(char *, int, char **);
  int next(int, char **);
  int find(char *);
  int equalstyle(int);
  int atomstyle(int);
  char *retrieve(char *);
  double compute_equal(int);
  double compute_equal(char *);
  void compute_atom(int, int, double *, int, int);
  tagint int_between_brackets(char *&, int);
  double evaluate_boolean(char *);

  void equal_save(int, char *&);
  void equal_restore(int, char *);
  void equal_override(int, double);

  unsigned int data_mask(int ivar);
  unsigned int data_mask(char *str);

 private:
  int me;
  int nvar;                // # of defined variables
  int maxvar;              // max # of variables following lists can hold
  char **names;            // name of each variable
  int *style;              // style of each variable
  int *num;                // # of values for each variable
  int *which;              // next available value for each variable
  int *pad;                // 1 = pad loop/uloop variables with 0s, 0 = no pad
  class VarReader **reader;   // variable that reads from file
  char ***data;            // str value of each variable's values

  int *eval_in_progress;   // flag if evaluation of variable is in progress

  class RanMars *randomequal;   // random number generator for equal-style vars
  class RanMars *randomatom;    // random number generator for atom-style vars

  int precedence[17];      // precedence level of math operators
                           // set length to include up to OR in enum

  struct Tree {            // parse tree for atom-style variables
    double value;          // single scalar  
    double *array;         // per-atom or per-type list of doubles
    int *iarray;           // per-atom list of ints
    bigint *barray;        // per-atom list of bigints
    int type;              // operation, see enum{} in variable.cpp
    int nstride;           // stride between atoms if array is a 2d array
    int selfalloc;         // 1 if array is allocated here, else 0
    int ivalue1,ivalue2;   // extra values for needed for gmask,rmask,grmask
    int nextra;            // # of additional args beyond first 2
    Tree *first,*second;   // ptrs further down tree for first 2 args
    Tree **extra;          // ptrs further down tree for nextra args
  };

  void remove(int);
  void grow();
  void copy(int, char **, char **);
  double evaluate(char *, Tree **);
  double collapse_tree(Tree *);
  double eval_tree(Tree *, int);
  void free_tree(Tree *);
  int find_matching_paren(char *, int, char *&);
  int math_function(char *, char *, Tree **, Tree **, int &, double *, int &);
  int group_function(char *, char *, Tree **, Tree **, int &, double *, int &);
  int region_function(char *);
  int special_function(char *, char *, Tree **, Tree **,
                       int &, double *, int &);
  void peratom2global(int, char *, double *, int, tagint,
                      Tree **, Tree **, int &, double *, int &);
  int is_atom_vector(char *);
  void atom_vector(char *, Tree **, Tree **, int &);
  int is_constant(char *);
  double constant(char *);
  int parse_args(char *, char **);
  char *find_next_comma(char *);
  void print_tree(Tree *, int);
};

class VarReader : protected Pointers {
 public:
  class FixStore *fix;
  char *id_fix;

  VarReader(class LAMMPS *, char *, char *, int);
  ~VarReader();
  int read_scalar(char *);
  int read_peratom();

 private:
  int me,style;
  FILE *fp;
  char *buffer;
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: World variable count doesn't match # of partitions

A world-style variable must specify a number of values equal to the
number of processor partitions.

E: Universe/uloop variable count < # of partitions

A universe or uloop style variable must specify a number of values >= to the
number of processor partitions.

E: All universe/uloop variables must have same # of values

Self-explanatory.

E: Cannot redefine variable as a different style

An equal-style variable can be re-defined but only if it was
originally an equal-style variable.

E: File variable could not read value

Check the file assigned to the variable.

E: Atomfile variable could not read values

Check the file assigned to the variable.

E: Variable name must be alphanumeric or underscore characters

Self-explanatory.

E: Invalid variable in next command

Self-explanatory.

E: All variables in next command must be same style

Self-explanatory.

E: Invalid variable style with next command

Variable styles {equal} and {world} cannot be used in a next
command.

E: Next command must list all universe and uloop variables

This is to insure they stay in sync.

E: Invalid syntax in variable formula

Self-explanatory.

E: Variable evaluation before simulation box is defined

Cannot evaluate a compute or fix or atom-based value in a variable
before the simulation has been setup.

E: Invalid compute ID in variable formula

The compute is not recognized.

E: Compute used in variable between runs is not current

Computes cannot be invoked by a variable in between runs.  Thus they
must have been evaluated on the last timestep of the previous run in
order for their value(s) to be accessed.  See the doc page for the
variable command for more info.

E: Variable formula compute vector is accessed out-of-range

Self-explanatory.

E: Variable formula compute array is accessed out-of-range

Self-explanatory.

E: Per-atom compute in equal-style variable formula

Equal-style variables cannot use per-atom quantities.

E: Mismatched compute in variable formula

A compute is referenced incorrectly or a compute that produces per-atom
values is used in an equal-style variable formula.

E: Invalid fix ID in variable formula

The fix is not recognized.

E: Fix in variable not computed at compatible time

Fixes generate their values on specific timesteps.  The variable is
requesting the values on a non-allowed timestep.

E: Variable formula fix vector is accessed out-of-range

Self-explanatory.

E: Variable formula fix array is accessed out-of-range

Self-explanatory.

E: Per-atom fix in equal-style variable formula

Equal-style variables cannot use per-atom quantities.

E: Mismatched fix in variable formula

A fix is referenced incorrectly or a fix that produces per-atom
values is used in an equal-style variable formula.

E: Invalid variable name in variable formula

Variable name is not recognized.

E: Variable has circular dependency

A circular dependency is when variable "a" in used by variable "b" and
variable "b" is also used by varaible "a".  Circular dependencies with
longer chains of dependence are also not allowed.

E: Invalid variable evaluation in variable formula

A variable used in a formula could not be evaluated.

E: Atom-style variable in equal-style variable formula

Atom-style variables generate one value per atom which is not allowed
in an equal-style variable.

E: Atomfile-style variable in equal-style variable formula

Self-explanatory.

E: Mismatched variable in variable formula

A variable is referenced incorrectly or an atom-style variable that
produces per-atom values is used in an equal-style variable
formula.

E: Invalid math/group/special function in variable formula

Self-explanatory.

E: Invalid thermo keyword in variable formula

The keyword is not recognized.

E: Divide by 0 in variable formula

Self-explanatory.

E: Modulo 0 in variable formula

Self-explanatory.

E: Power by 0 in variable formula

Self-explanatory.

E: Sqrt of negative value in variable formula

Self-explanatory.

E: Log of zero/negative value in variable formula

Self-explanatory.

E: Arcsin of invalid value in variable formula

Argument of arcsin() must be between -1 and 1.

E: Arccos of invalid value in variable formula

Argument of arccos() must be between -1 and 1.

E: Invalid math function in variable formula

Self-explanatory.

E: Non digit character between brackets in variable

Self-explantory.

E: Mismatched brackets in variable

Self-explanatory.

E: Empty brackets in variable

There is no variable syntax that uses empty brackets.  Check
the variable doc page.

E: Index between variable brackets must be positive

Self-explanatory.

E: Cannot use ramp in variable formula between runs

This is because the ramp() function is time dependent.

E: Cannot use vdisplace in variable formula between runs

This is a function of elapsed time.

E: Cannot use swiggle in variable formula between runs

This is a function of elapsed time.

E: Cannot use cwiggle in variable formula between runs

This is a function of elapsed time.

E: Group ID in variable formula does not exist

Self-explanatory.

E: Invalid group function in variable formula

Group function is not recognized.

E: Region ID in variable formula does not exist

Self-explanatory.

E: Invalid special function in variable formula

Self-explanatory.

E: Gmask function in equal-style variable formula

Gmask is per-atom operation.

E: Rmask function in equal-style variable formula

Rmask is per-atom operation.

E: Grmask function in equal-style variable formula

Grmask is per-atom operation.

E: Variable ID in variable formula does not exist

Self-explanatory.

E: Atomfile variable in equal-style variable formula

Self-explanatory.

E: Invalid variable style in special function next

Only file-style or atomfile-style variables can be used with next().

E: Indexed per-atom vector in variable formula without atom map

Accessing a value from an atom vector requires the ability to lookup
an atom index, which is provided by an atom map.  An atom map does not
exist (by default) for non-molecular problems.  Using the atom_modify
map command will force an atom map to be created.

E: Variable uses atom property that isn't allocated

Self-explanatory.

E: Invalid atom vector in variable formula

The atom vector is not recognized.

E: Atom vector in equal-style variable formula

Atom vectors generate one value per atom which is not allowed
in an equal-style variable.

E: Expected floating point parameter in variable definition

The quantity being read is a non-numeric value.

E: Expected integer parameter in variable definition

The quantity being read is a floating point or non-numeric value.

E: Invalid Boolean syntax in if command

Self-explanatory.

E: Cannot open file variable file %s

The specified file cannot be opened.  Check that the path and name are
correct.

E: Cannot use atomfile-style variable unless atom map exists

Self-explanatory.  See the atom_modify command to create a map.

E: Invalid atom ID in variable file

Self-explanatory.

*/
