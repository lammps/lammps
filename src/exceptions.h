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

#ifndef LMP_EXCEPTIONS_H
#define LMP_EXCEPTIONS_H

#include <exception>
#include <mpi.h>
#include <string>

namespace LAMMPS_NS {

class LAMMPSException : public std::exception {
 public:
  std::string message;

  LAMMPSException(const std::string &msg) : message(msg) {}

  const char *what() const noexcept override { return message.c_str(); }
};

class LAMMPSAbortException : public LAMMPSException {
 public:
  MPI_Comm universe;

  LAMMPSAbortException(const std::string &msg, MPI_Comm _universe) :
      LAMMPSException(msg), universe(_universe)
  {
  }
};

enum ErrorType { ERROR_NONE = 0, ERROR_NORMAL = 1, ERROR_ABORT = 2 };

}    // namespace LAMMPS_NS

#endif
