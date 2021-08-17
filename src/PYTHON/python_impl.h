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
  bool external_interpreter;

  PythonImpl(class LAMMPS *);
  ~PythonImpl();
  void command(int, char **);
  void invoke_function(int, char *);
  int find(const char *);
  int variable_match(const char *, const char *, int);
  char *long_string(int);
  int execute_string(char *);
  int execute_file(char *);
  bool has_minimum_version(int major, int minor);

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

/* ERROR/WARNING messages:

E: Could not initialize embedded Python

The main module in Python was not accessible.

E: Invalid python command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Python invoke of undefined function

Cannot invoke a function that has not been previously defined.

E: Python variable does not match Python function

This matching is defined by the python-style variable and the python
command.

E: Could not process Python source command

UNDOCUMENTED

E: Could not open Python file

The specified file of Python code cannot be opened.  Check that the
path and name are correct.

E: Could not process Python file

The Python code in the specified file was not run successfully by
Python, probably due to errors in the Python code.

E: Could not process Python string

The Python code in the here string was not run successfully by Python,
probably due to errors in the Python code.

E: Could not find Python function

The provided Python code was run successfully, but it not
define a callable function with the required name.

E: Python function is not callable

The provided Python code was run successfully, but it not
define a callable function with the required name.

E: Could not create Python function arguments

This is an internal Python error, possibly because the number
of inputs to the function is too large.

E: Could not evaluate Python function input variable

Self-explanatory.

E: Unsupported variable type

UNDOCUMENTED

E: Python function evaluation failed

The Python function did not run successfully and/or did not return a
value (if it is supposed to return a value).  This is probably due to
some error condition in the function.

E: Python command length keyword cannot be used unless output is a string

UNDOCUMENTED

U: Cannot embed Python when also extending Python with LAMMPS

When running LAMMPS via Python through the LAMMPS library interface
you cannot also user the input script python command.

*/
