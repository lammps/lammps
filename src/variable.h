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

#include <cstdio>
#include "pointers.h"

namespace LAMMPS_NS {

class Variable : protected Pointers {
 friend class Info;
 public:
  Variable(class LAMMPS *);
  ~Variable();
  void set(int, char **);
  void set(char *, int, char **);
  int set_string(char *, char *);
  int next(int, char **);

  int find(char *);
  void set_arrays(int);
  void python_command(int, char **);

  int equalstyle(int);
  int atomstyle(int);
  int vectorstyle(int);
  char *pythonstyle(char *, char *);
  int internalstyle(int);

  char *retrieve(char *);
  double compute_equal(int);
  double compute_equal(char *);
  void compute_atom(int, int, double *, int, int);
  int compute_vector(int, double **);
  void internal_set(int, double);

  tagint int_between_brackets(char *&, int);
  double evaluate_boolean(char *);

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
  double *dvalue;          // single numeric value for internal variables

  struct VecVar {
    int n,nmax;
    bigint currentstep;
    double *values;
  };
  VecVar *vecs;

  int *eval_in_progress;       // flag if evaluation of variable is in progress
  int treetype;                // ATOM or VECTOR flag for formula evaluation

  class RanMars *randomequal;   // random number generator for equal-style vars
  class RanMars *randomatom;    // random number generator for atom-style vars

  int precedence[18];      // precedence level of math operators
                           // set length to include up to XOR in enum

  struct Tree {            // parse tree for atom-style or vector-style vars
    double value;          // single scalar
    double *array;         // per-atom or per-type list of doubles
    int *iarray;           // per-atom list of ints
    bigint *barray;        // per-atom list of bigints
    int type;              // operation, see enum{} in variable.cpp
    int nvector;           // length of array for vector-style variable
    int nstride;           // stride between atoms if array is a 2d array
    int selfalloc;         // 1 if array is allocated here, else 0
    int ivalue1,ivalue2;   // extra values needed for gmask,rmask,grmask
    int nextra;            // # of additional args beyond first 2
    Tree *first,*second;   // ptrs further down tree for first 2 args
    Tree **extra;          // ptrs further down tree for nextra args

    Tree() :
      array(NULL), iarray(NULL), barray(NULL),
      selfalloc(0), ivalue1(0), ivalue2(0), nextra(0),
      first(NULL), second(NULL), extra(NULL) {}
  };

  int compute_python(int);
  void remove(int);
  void grow();
  void copy(int, char **, char **);
  double evaluate(char *, Tree **, int);
  double collapse_tree(Tree *);
  double eval_tree(Tree *, int);
  int size_tree_vector(Tree *);
  int compare_tree_vector(int, int);
  void free_tree(Tree *);
  int find_matching_paren(char *, int, char *&, int);
  int math_function(char *, char *, Tree **, Tree **,
                    int &, double *, int &, int);
  int group_function(char *, char *, Tree **, Tree **,
                     int &, double *, int &, int);
  int region_function(char *, int);
  int special_function(char *, char *, Tree **, Tree **,
                       int &, double *, int &, int);
  void peratom2global(int, char *, double *, int, tagint,
                      Tree **, Tree **, int &, double *, int &);
  int is_atom_vector(char *);
  void atom_vector(char *, Tree **, Tree **, int &);
  int is_constant(char *);
  double constant(char *);
  int parse_args(char *, char **);
  char *find_next_comma(char *);
  void print_var_error(const char *, int, const char *, int, int global=1);
  void print_tree(Tree *, int);
};

class VarReader : protected Pointers {
 public:
  class FixStore *fixstore;
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

E: Cannot open temporary file for world counter.

Self-explanatory.

E: All universe/uloop variables must have same # of values

Self-explanatory.

E: Cannot redefine variable as a different style

An equal-style variable can be re-defined but only if it was
originally an equal-style variable.

E: File variable could not read value

Check the file assigned to the variable.

E: Atomfile variable could not read values

Check the file assigned to the variable.

E: LAMMPS is not built with Python embedded

This is done by including the PYTHON package before LAMMPS is built.
This is required to use python-style variables.

E: Variable name '%s' must have only alphanumeric characters or underscore

UNDOCUMENTED

E: Invalid variable '%s' in next command

UNDOCUMENTED

E: All variables in next command must have same style

UNDOCUMENTED

E: Invalid variable style with next command

Variable styles {equal} and {world} cannot be used in a next
command.

E: Incorrect conversion in format string

A format style variable was not using either a %f, a %g, or a %e conversion.

E: Next command must list all universe and uloop variables

This is to insure they stay in sync.

E: Python variable '%s' does not match Python function

UNDOCUMENTED

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

E: Variable name between brackets must be alphanumeric or underscore characters

Self-explanatory.

E: Non digit character between brackets in variable

Self-explanatory.

E: Mismatched brackets in variable

Self-explanatory.

E: Empty brackets in variable

There is no variable syntax that uses empty brackets.  Check
the variable doc page.

E: Invalid variable name in variable formula

Variable name is not recognized.

E: Invalid variable evaluation in variable formula

A variable used in a formula could not be evaluated.

E: Index between variable brackets must be positive

Self-explanatory.

E: Indexed per-atom vector in variable formula without atom map

Accessing a value from an atom vector requires the ability to lookup
an atom index, which is provided by an atom map.  An atom map does not
exist (by default) for non-molecular problems.  Using the atom_modify
map command will force an atom map to be created.

E: Variable atom ID is too large

Specified ID is larger than the maximum allowed atom ID.

E: Variable uses atom property that isn't allocated

Self-explanatory.

E: Invalid atom vector in variable formula

The atom vector is not recognized.

E: Atom vector in equal-style variable formula

Atom vectors generate one value per atom which is not allowed
in an equal-style variable.

E: Too many args in variable function

More args are used than any variable function allows.

E: Invalid Boolean syntax in if command

Self-explanatory.

E: Cannot open file variable file %s

The specified file cannot be opened.  Check that the path and name are
correct.

E: Cannot use atomfile-style variable unless atom map exists

Self-explanatory.  See the atom_modify command to create a map.

E: Invalid atom ID in variable file

Self-explanatory.

U: Variable name must be alphanumeric or underscore characters

Self-explanatory.

U: Invalid variable in next command

Self-explanatory.

U: All variables in next command must be same style

Self-explanatory.

U: Variable has circular dependency

A circular dependency is when variable "a" in used by variable "b" and
variable "b" is also used by variable "a".  Circular dependencies with
longer chains of dependence are also not allowed.

U: Python variable does not match Python function

This matching is defined by the python-style variable and the python
command.

U: Python variable has no function

No python command was used to define the function associated with the
python-style variable.

U: Invalid syntax in variable formula

Self-explanatory.

U: Variable evaluation before simulation box is defined

Cannot evaluate a compute or fix or atom-based value in a variable
before the simulation has been setup.

U: Invalid compute ID in variable formula

The compute is not recognized.

U: Compute used in variable between runs is not current

Computes cannot be invoked by a variable in between runs.  Thus they
must have been evaluated on the last timestep of the previous run in
order for their value(s) to be accessed.  See the doc page for the
variable command for more info.

U: Variable formula compute vector is accessed out-of-range

Self-explanatory.

U: Variable formula compute array is accessed out-of-range

Self-explanatory.

U: Per-atom compute in equal-style variable formula

Equal-style variables cannot use per-atom quantities.

U: Mismatched compute in variable formula

A compute is referenced incorrectly or a compute that produces per-atom
values is used in an equal-style variable formula.

U: Invalid fix ID in variable formula

The fix is not recognized.

U: Fix in variable not computed at compatible time

Fixes generate their values on specific timesteps.  The variable is
requesting the values on a non-allowed timestep.

U: Variable formula fix vector is accessed out-of-range

Self-explanatory.

U: Variable formula fix array is accessed out-of-range

Self-explanatory.

U: Per-atom fix in equal-style variable formula

Equal-style variables cannot use per-atom quantities.

U: Mismatched fix in variable formula

A fix is referenced incorrectly or a fix that produces per-atom
values is used in an equal-style variable formula.

U: Atom-style variable in equal-style variable formula

Atom-style variables generate one value per atom which is not allowed
in an equal-style variable.

U: Atomfile-style variable in equal-style variable formula

Self-explanatory.

U: Mismatched variable in variable formula

A variable is referenced incorrectly or an atom-style variable that
produces per-atom values is used in an equal-style variable
formula.

U: Invalid math/group/special function in variable formula

Self-explanatory.

U: Invalid thermo keyword in variable formula

The keyword is not recognized.

U: Cannot use ramp in variable formula between runs

This is because the ramp() function is time dependent.

U: Cannot use vdisplace in variable formula between runs

This is a function of elapsed time.

U: Cannot use swiggle in variable formula between runs

This is a function of elapsed time.

U: Cannot use cwiggle in variable formula between runs

This is a function of elapsed time.

U: Group ID in variable formula does not exist

Self-explanatory.

U: Invalid group function in variable formula

Group function is not recognized.

U: Region ID in variable formula does not exist

Self-explanatory.

U: Invalid special function in variable formula

Self-explanatory.

U: Gmask function in equal-style variable formula

Gmask is per-atom operation.

U: Rmask function in equal-style variable formula

Rmask is per-atom operation.

U: Grmask function in equal-style variable formula

Grmask is per-atom operation.

U: Variable ID in variable formula does not exist

Self-explanatory.

U: Atomfile variable in equal-style variable formula

Self-explanatory.

U: Invalid variable style in special function next

Only file-style or atomfile-style variables can be used with next().

U: Invalid is_active() function in variable formula

Self-explanatory.

U: Invalid is_available() function in variable formula

Self-explanatory.

U: Invalid is_defined() function in variable formula

Self-explanatory.

*/
