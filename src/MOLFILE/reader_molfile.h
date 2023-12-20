/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributed by Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef READER_CLASS
// clang-format off
ReaderStyle(molfile,ReaderMolfile);
// clang-format on
#else

#ifndef LMP_READER_MOLFILE_H
#define LMP_READER_MOLFILE_H

#include "reader.h"

namespace LAMMPS_NS {

class ReaderMolfile : public Reader {
 public:
  ReaderMolfile(class LAMMPS *);
  ~ReaderMolfile() override;

  void settings(int, char **) override;

  int read_time(bigint &) override;
  void skip() override;
  bigint read_header(double[3][3], int &, int &, int, int, int *, char **, int, int, int &, int &,
                     int &, int &) override;
  void read_atoms(int, int, double **) override;

  void open_file(const std::string &) override;
  void close_file() override;

 private:
  int *fieldindex;    // mapping of input fields to dump

  class MolfileInterface *mf;
  float *coords;    // pointer to temporary coordinate storage
  float *vels;      // pointer to temporary velocity storage
  int *types;       // pointer to temporary type info storage
  float cell[6];    // box info (stored as, a, b, c, alpha, beta, gamma)
  int natoms;       // current number of atoms
  int needvels;     // 1 if velocities are required, otherwise 0

  bigint nstep;    // current (time) step number
  bigint nid;      // current atom id.

  int me;
};

}    // namespace LAMMPS_NS

#endif
#endif
