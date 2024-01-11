/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributed by Timothy Sirk
------------------------------------------------------------------------- */

#ifndef LMP_READER_H
#define LMP_READER_H

#include "pointers.h"

namespace LAMMPS_NS {

class Reader : protected Pointers {
 public:
  enum { ID, TYPE, X, Y, Z, VX, VY, VZ, Q, IX, IY, IZ, FX, FY, FZ };
  enum { UNSET, NOSCALE_NOWRAP, NOSCALE_WRAP, SCALE_NOWRAP, SCALE_WRAP };

  Reader(class LAMMPS *);
  ~Reader() override;

  virtual void settings(int, char **);

  virtual int read_time(bigint &) = 0;
  virtual void skip() = 0;
  virtual bigint read_header(double[3][3], int &, int &, int, int, int *, char **, int, int, int &,
                             int &, int &, int &) = 0;
  virtual void read_atoms(int, int, double **) = 0;

  virtual void open_file(const std::string &);
  virtual void close_file();

 protected:
  FILE *fp;           // pointer to opened file or pipe
  bool compressed;    // flag for dump file compression
  bool binary;        // flag for (native) binary files
};

}    // namespace LAMMPS_NS

#endif
