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

#ifndef LMP_PYTHON_H
#define LMP_PYTHON_H

#include "pointers.h"

namespace LAMMPS_NS {

class Python : protected Pointers {
 public:
  int python_exists;

  Python(class LAMMPS *);
  ~Python();
  void command(int, char **);
  void invoke_function(int, char *);
  int find(char *);
  int variable_match(char *, char *, int);

 private:
  int ninput,noutput;
  char **istr;
  char *ostr,*format;
  void *pyMain;
  
  struct PyFunc {
    char *name;
    int ninput,noutput;
    int *itype,*ivarflag;
    int *ivalue;
    double *dvalue;
    char **svalue;
    int otype;
    char *ovarname;
    void *pFunc;
  };

  PyFunc *pfuncs;
  int nfunc;

  int create_entry(char *);
  void deallocate(int);
};

}

#endif
