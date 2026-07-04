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

#ifndef LMP_STL_READER_H
#define LMP_STL_READER_H

#include "pointers.h"

#include <exception>
#include <string>
#include <utility>
#include <vector>

namespace LAMMPS_NS {

// exception thrown for malformed or unreadable STL files

class STLReaderException : public std::exception {
  std::string message;

 public:
  explicit STLReaderException(std::string msg) : message(std::move(msg)) {}
  [[nodiscard]] const char *what() const noexcept override { return message.c_str(); }
};

class STLReader : protected Pointers {
 public:
  // a single STL triangle: outward facet normal + 3 vertices (x,y,z each)

  struct Triangle {
    double normal[3];
    double vert[3][3];
  };

  STLReader(class LAMMPS *);
  ~STLReader() override;

  // read triangle vertices from an STL file on rank 0 and broadcast to all ranks
  // caller_tris is filled with 9 doubles per triangle (3 vertices, normal omitted)
  // returns the number of triangles.  Used by fix surface/global.

  int read_file(const char *filename, double **&caller_tris);

  // stand-alone parser (no MPI, no other LAMMPS state): auto-detects ASCII vs
  // binary STL format and returns the list of triangles.  Throws an
  // STLReaderException on any malformed or unreadable input.  This routine is
  // covered by the unit tests and is shared by all STL consumers in LAMMPS.

  static std::vector<Triangle> parse(const std::string &filename, std::string *title = nullptr);

 private:
  static std::vector<Triangle> parse_text(const std::string &filename, std::string &title);
  static std::vector<Triangle> parse_binary(FILE *fp, std::string &title);

  int ntris, maxtris;
  double **tris;
};

}    // namespace LAMMPS_NS

#endif
