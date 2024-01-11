/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_VARIABLE_H
#define LMP_VARIABLE_H

#include "pointers.h"

namespace LAMMPS_NS {
class Region;

class Variable : protected Pointers {
  friend class Info;

 public:
  Variable(class LAMMPS *);
  ~Variable() override;
  void set(int, char **);
  void set(const std::string &);
  void set(char *, int, char **);
  int set_string(const char *, const char *);
  int next(int, char **);

  int find(const char *);
  void set_arrays(int);
  void python_command(int, char **);
  void purge_atomfile();

  int equalstyle(int);
  int atomstyle(int);
  int vectorstyle(int);
  char *pythonstyle(char *, char *);
  int internalstyle(int);

  char *retrieve(const char *);
  double compute_equal(int);
  double compute_equal(const std::string &);
  void compute_atom(int, int, double *, int, int);
  int compute_vector(int, double **);
  void internal_set(int, double);

  tagint int_between_brackets(char *&, int);
  double evaluate_boolean(char *);

 public:
  int nvar;        // # of defined variables
  char **names;    // name of each variable

  // must match "varstyles" array in info.cpp
  enum {
    INDEX,
    LOOP,
    WORLD,
    UNIVERSE,
    ULOOP,
    STRING,
    GETENV,
    SCALARFILE,
    ATOMFILE,
    FORMAT,
    EQUAL,
    ATOM,
    VECTOR,
    PYTHON,
    TIMER,
    INTERNAL
  };
  static constexpr int VALUELENGTH = 64;

 private:
  int me;
  int maxvar;                  // max # of variables following lists can hold
  int *style;                  // style of each variable
  int *num;                    // # of values for each variable
  int *which;                  // next available value for each variable
  int *pad;                    // 1 = pad loop/uloop variables with 0s, 0 = no pad
  class VarReader **reader;    // variable that reads from file
  char ***data;                // str value of each variable's values
  double *dvalue;              // single numeric value for internal variables

  struct VecVar {
    int n, nmax;
    int dynamic;
    bigint currentstep;
    double *values;
  };
  VecVar *vecs;

  int *eval_in_progress;    // flag if evaluation of variable is in progress
  int treetype;             // ATOM or VECTOR flag for formula evaluation

  class RanMars *randomequal;    // random number generator for equal-style vars
  class RanMars *randomatom;     // random number generator for atom-style vars

  int precedence[18];    // precedence level of math operators
                         // set length to include up to XOR in enum

  struct Tree {              // parse tree for atom-style or vector-style vars
    double value;            // single scalar
    double *array;           // per-atom or per-type list of doubles
    int *iarray;             // per-atom list of ints
    bigint *barray;          // per-atom list of bigints
    int type;                // operation, see enum{} in variable.cpp
    int nvector;             // length of array for vector-style variable
    int nstride;             // stride between atoms if array is a 2d array
    int selfalloc;           // 1 if array is allocated here, else 0
    int ivalue;              // extra value needed for gmask, grmask
    int nextra;              // # of additional args beyond first 2
    Region *region;          // region pointer for rmask, grmask
    Tree *first, *second;    // ptrs further down tree for first 2 args
    Tree **extra;            // ptrs further down tree for nextra args

    Tree() :
        array(nullptr), iarray(nullptr), barray(nullptr), selfalloc(0), ivalue(0), nextra(0),
        region(nullptr), first(nullptr), second(nullptr), extra(nullptr)
    {
    }
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
  int math_function(char *, char *, Tree **, Tree **, int &, double *, int &, int);
  int group_function(char *, char *, Tree **, Tree **, int &, double *, int &, int);
  Region *region_function(char *, int);
  int special_function(char *, char *, Tree **, Tree **, int &, double *, int &, int);
  int feature_function(char *, char *, Tree **, Tree **, int &, double *, int &, int);
  void peratom2global(int, char *, double *, int, tagint, Tree **, Tree **, int &, double *, int &);
  void custom2global(int *, double *, int, tagint, Tree **, Tree **, int &, double *, int &);
  int is_atom_vector(char *);
  void atom_vector(char *, Tree **, Tree **, int &);
  int parse_args(char *, char **);
  void parse_vector(int, char *);
  char *find_next_comma(char *);
  void print_var_error(const std::string &, int, const std::string &, int, int global = 1);
  void print_tree(Tree *, int);
};

class VarReader : protected Pointers {
 public:
  class FixStoreAtom *fixstore;
  char *id_fix;

  VarReader(class LAMMPS *, char *, char *, int);
  ~VarReader() override;
  int read_scalar(char *);
  int read_peratom();

 private:
  int me, style;
  FILE *fp;
  char *buffer;
};

}    // namespace LAMMPS_NS

#endif
