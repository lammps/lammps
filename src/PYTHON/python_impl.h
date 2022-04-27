/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_PYTHON_IMPL_H
#define LMP_PYTHON_IMPL_H

#include "lmppython.h"
#include "pointers.h"

namespace LAMMPS_NS {

class PythonImpl : protected Pointers, public PythonInterface {

 public:
  PythonImpl(class LAMMPS *);
  ~PythonImpl() override;
  void command(int, char **) override;
  void invoke_function(int, char *) override;
  int find(const char *) override;
  int variable_match(const char *, const char *, int) override;
  char *long_string(int) override;
  int execute_string(char *) override;
  int execute_file(char *) override;
  bool has_minimum_version(int major, int minor) override;
  static void finalize();

 private:
  void *pyMain;

  struct PyFunc {
    char *name;
    int ninput, noutput;
    int *itype, *ivarflag;
    int *ivalue;
    double *dvalue;
    char **svalue;
    int otype;
    char *ovarname;
    char *longstr;
    int length_longstr;
    void *pFunc;
  };

  PyFunc *pfuncs;
  int nfunc;

  int create_entry(char *, int, int, int, char **, char *, char *);
  void deallocate(int);
};

}    // namespace LAMMPS_NS

#endif
