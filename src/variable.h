/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef VARIABLE_H
#define VARIABLE_H

#include "lammps.h"

class Variable : public LAMMPS {
 public:
  Variable();
  ~Variable();
  void set(int, char **);
  void set(char *, char *);
  int next(int, char **);
  char *retrieve(char *);

 private:
  int me;
  int nvar;                // # of defined variables
  int maxvar;              // max # of variables arrays can hold
  char **names;            // name of each variable
  int *style;              // style of each variable
  int *num;                // # of values for each variable
  int *index;              // next available value for each variable
  char ***data;            // str value of each variable's values

  int find(char *);
  void copy(int, char **, char **);
  char *evaluate(char *);
  void remove(int);
};

#endif
