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

#ifndef VARIABLE_H
#define VARIABLE_H

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
  int peratom(int);
  char *retrieve(char *);

  void build_parse_tree(int);
  void evaluate_parse_tree(int, double *);
  void free_parse_tree();

 private:
  int me;
  int nvar;                // # of defined variables
  int maxvar;              // max # of variables arrays can hold
  char **names;            // name of each variable
  int *style;              // style of each variable
  int *num;                // # of values for each variable
  int *index;              // next available value for each variable
  char ***data;            // str value of each variable's values

  struct Tree {
    double value;
    double *array;
    int nstride;
    int type;
    Tree *left,*right;
  };

  Tree *ptree;             // parse tree for an ATOM variable

  void copy(int, char **, char **);
  double evaluate(char *, Tree *);
  void remove(int);
  double eval_tree(Tree *, int);
  void free_tree(Tree *);
};

}

#endif
