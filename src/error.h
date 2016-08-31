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

#ifndef LMP_ERROR_H
#define LMP_ERROR_H

#include "pointers.h"
#include <string>
#include <exception>

namespace LAMMPS_NS {

class LAMMPSException : public std::exception
{
public:
  std::string message;

  LAMMPSException(std::string msg) : message(msg) {
  }

  ~LAMMPSException() throw() {
  }

  virtual const char * what() const throw() {
    return message.c_str();
  }
};

class LAMMPSAbortException : public LAMMPSException {
public:
  MPI_Comm universe;

  LAMMPSAbortException(std::string msg, MPI_Comm universe) :
    LAMMPSException(msg),
    universe(universe)
  {
  }
};

class Error : protected Pointers {
  char * last_error_message;

 public:
  Error(class LAMMPS *);

  void universe_all(const char *, int, const char *);
  void universe_one(const char *, int, const char *);
  void universe_warn(const char *, int, const char *);

  void all(const char *, int, const char *);
  void one(const char *, int, const char *);
  void warning(const char *, int, const char *, int = 1);
  void message(const char *, int, const char *, int = 1);
  void done(int = 0); // 1 would be fully backwards compatible

  char * get_last_error() const;
  void   set_last_error(const char * msg);
};

}

#endif

/* ERROR/WARNING messages:

*/
