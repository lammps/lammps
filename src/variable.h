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

#ifndef LMP_VARIABLE_H
#define LMP_VARIABLE_H

#include "pointers.h"

namespace LAMMPS_NS {

class Variable : protected Pointers {
 public:
  Variable(class LAMMPS *);
  ~Variable();
  void set(int, char **);
  void set(char *, char *);
  int next(int, char **);
  int find(char *);
  int equalstyle(int);
  int atomstyle(int);
  char *retrieve(char *);
  double compute_equal(int);
  void compute_atom(int, int, double *, int, int);
  int int_between_brackets(char *&);
  double evaluate_boolean(char *);

 private:
  int me;
  int nvar;                // # of defined variables
  int maxvar;              // max # of variables arrays can hold
  char **names;            // name of each variable
  int *style;              // style of each variable
  int *num;                // # of values for each variable
  int *which;              // next available value for each variable
  int *pad;                // 1 = pad loop/uloop variables with 0s, 0 = no pad
  char ***data;            // str value of each variable's values
  double PI;

  class RanMars *randomequal;   // random number generator for equal-style vars
  class RanMars *randomatom;    // random number generator for atom-style vars

  int precedence[16];      // precedence level of math operators
                           // set length to include OR in enum

  struct Tree {            // parse tree for atom-style variables
    double value;
    double *array;
    int *iarray;
    int type;
    int nstride;
    int ivalue1,ivalue2;
    Tree *left,*middle,*right;
  };

  void remove(int);
  void extend();
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
  void peratom2global(int, char *, double *, int, int,
		      Tree **, Tree **, int &, double *, int &);
  int is_atom_vector(char *);
  void atom_vector(char *, Tree **, Tree **, int &);
  int is_constant(char *);
  double constant(char *);
  double numeric(char *);
  int inumeric(char *);
  void print_tree(Tree *, int);
};

}

#endif
