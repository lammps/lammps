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

#ifndef LMP_EXCEPTIONS_H
#define LMP_EXCEPTIONS_H

#include <mpi.h>
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

enum ErrorType {
   ERROR_NONE   = 0,
   ERROR_NORMAL = 1,
   ERROR_ABORT  = 2
};

}

#endif
