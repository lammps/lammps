/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributed by Norbert Podhorszki (Oak Ridge National Laboratory)
------------------------------------------------------------------------- */

#ifdef READER_CLASS
// clang-format off
ReaderStyle(adios, ReaderADIOS);
// clang-format on
#else

#ifndef LMP_READER_ADIOS_H
#define LMP_READER_ADIOS_H

#include "reader.h"

#include <map>
#include <string>
#include <vector>

namespace LAMMPS_NS {
class ReadADIOSInternal;

class ReaderADIOS : public Reader {
 public:
  ReaderADIOS(class LAMMPS *);
  virtual ~ReaderADIOS();

  virtual void settings(int, char **);

  virtual int read_time(bigint &);
  virtual void skip();
  virtual bigint read_header(double[3][3], int &, int &, int, int, int *, char **, int, int, int &,
                             int &, int &, int &);
  virtual void read_atoms(int, int, double **);

  virtual void open_file(const char *);
  virtual void close_file();

 private:
  int *fieldindex;         // mapping of input fields to dump
  uint64_t nAtomsTotal;    // current number of atoms in entire dump step
  uint64_t nAtoms;         // current number of atoms for this process
                           // (Sum(nAtoms)=nAtomsTotal)
  uint64_t atomOffset;     // starting atom position for this process to read

  bigint nstep;    // current (time) step number
  bigint nid;      // current atom id.

  int me;
  ReadADIOSInternal *internal;

  int find_label(const std::string &label, const std::map<std::string, int> &labels);
};

}    // namespace LAMMPS_NS

#endif
#endif
