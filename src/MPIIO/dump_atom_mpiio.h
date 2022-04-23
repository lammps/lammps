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

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(atom/mpiio,DumpAtomMPIIO);
// clang-format on
#else

#ifndef LMP_DUMP_ATOM_MPIIO_H
#define LMP_DUMP_ATOM_MPIIO_H

#include "dump_atom.h"

namespace LAMMPS_NS {

class DumpAtomMPIIO : public DumpAtom {

 public:
  DumpAtomMPIIO(class LAMMPS *, int, char **);
  ~DumpAtomMPIIO() override;

 protected:
  bigint
      sumFileSize;    // size in bytes of the file up through this rank offset from the end of the header data
  char *headerBuffer;    // buffer for holding header data

  MPI_File mpifh;
  MPI_Offset mpifo, offsetFromHeader, headerSize, currentFileSize;
  int performEstimate;    // switch for write_data and write_header methods to use for gathering data and detemining filesize for preallocation vs actually writing the data
  char *filecurrent;      // name of file for this round (with % and * replaced)

  void openfile() override;
  void write_header(bigint) override;
  void write() override;
  void write_data(int, double *) override;

  void init_style() override;
  typedef void (DumpAtomMPIIO::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;    // ptr to write header functions
  void header_binary(bigint);
  void header_binary_triclinic(bigint);
  void header_item(bigint);
  void header_item_triclinic(bigint);

#if defined(_OPENMP)
  int convert_string_omp(int, double *);     // multithreaded version of convert_string
  int convert_image_omp(int, double *);      // multithreaded version of convert_image
  int convert_noimage_omp(int, double *);    // multithreaded version of convert_noimage
#endif

  int convert_string(int, double *) override;
  typedef void (DumpAtomMPIIO::*FnPtrData)(int, double *);
  FnPtrData write_choice;    // ptr to write data functions
  void write_binary(int, double *);
  void write_string(int, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
