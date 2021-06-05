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

/* ----------------------------------------------------------------------
   Contributing author: Richard Berger (Temple U)
------------------------------------------------------------------------- */

#ifndef LMP_FILE_WRITER_H
#define LMP_FILE_WRITER_H

#include <string>

namespace LAMMPS_NS {

class FileWriter {
 public:
  FileWriter() = default;
  virtual ~FileWriter() = default;
  virtual void open(const std::string &path, bool append = false) = 0;
  virtual void close() = 0;
  virtual void flush() = 0;
  virtual size_t write(const void *buffer, size_t length) = 0;
  virtual bool isopen() const = 0;
};

class FileWriterException : public std::exception {
  std::string message;

 public:
  FileWriterException(const std::string &msg) : message(msg) {}

  ~FileWriterException() throw() {}

  virtual const char *what() const throw() { return message.c_str(); }
};

}    // namespace LAMMPS_NS

#endif
