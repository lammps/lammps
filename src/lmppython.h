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

#ifndef LMP_PYTHON_H
#define LMP_PYTHON_H

#include "pointers.h"

namespace LAMMPS_NS {

class PythonInterface {
 public:
  virtual ~PythonInterface() = default;
  virtual void command(int, char **) = 0;
  virtual void invoke_function(int, char *) = 0;
  virtual int find(const char *) = 0;
  virtual int variable_match(const char *, const char *, int) = 0;
  virtual char *long_string(int ifunc) = 0;
  virtual int execute_string(char *) = 0;
  virtual int execute_file(char *) = 0;
  virtual bool has_minimum_version(int major, int minor) = 0;
};

class Python : protected Pointers {
 public:
  Python(class LAMMPS *);
  ~Python() override;

  void command(int, char **);
  void invoke_function(int, char *);
  int find(const char *);
  int variable_match(const char *, const char *, int);
  char *long_string(int ifunc);
  int execute_string(char *);
  int execute_file(char *);
  bool has_minimum_version(int major, int minor);

  bool is_enabled() const;
  void init();
  static void finalize();

 private:
  PythonInterface *impl;
};

}    // namespace LAMMPS_NS
#endif
