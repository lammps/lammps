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
#include <utility>

namespace LAMMPS_NS {

/** Exception thrown by the STLReader class for malformed or unreadable files */

class STLReaderException : public std::exception {
  std::string message;

 public:
  explicit STLReaderException(std::string msg) : message(std::move(msg)) {}
  [[nodiscard]] const char *what() const noexcept override { return message.c_str(); }
};

/** Class for reading triangle meshes from files in STL format

\verbatim embed:rst

A file in `STL (stereolithography) format
<https://en.wikipedia.org/wiki/STL_(file_format)>`_ stores a triangle
mesh as a flat list of triangles.  Each triangle consists of an outward
facing normal vector and the coordinates of its three vertices.  There
are two variants of the format: a plain text ("ASCII") version and a
binary version.  This class transparently supports both: the format is
detected from the file contents and *not* from the file name.

Most functionality is in the static :cpp:func:`STLReader::parse()
<LAMMPS_NS::STLReader::parse>` function, which is independent from any
LAMMPS instance and thus can also be called from code that has no access
to one.  It reads and parses the entire file and returns the triangles
as a vector of :cpp:struct:`STLReader::Triangle
<LAMMPS_NS::STLReader::Triangle>` data structures.  All error conditions
are reported by throwing an :cpp:class:`STLReaderException
<LAMMPS_NS::STLReaderException>`, so the calling code decides how errors
are handled.  Since the file is read in its entirety by every MPI
process calling it, this is best suited for small to medium sized
meshes.

.. code-block:: c++
   :caption: Use of the STLReader class in the create_atoms command

   #include "stl_reader.h"

   std::vector<STLReader::Triangle> triangles;
   try {
     std::string title;
     triangles = STLReader::parse(filename, &title);
     if (comm->me == 0)
       utils::logmesg(lmp, "Reading STL object {} with {} triangles from file {}\n",
                      title.empty() ? "(unnamed)" : title, triangles.size(), filename);
   } catch (std::exception &e) {
     error->all(FLERR, Error::NOLASTLINE, "{}", e.what());
   }

   for (const auto &tri : triangles) {
     // process tri.normal[3] and tri.vert[3][3]
   }

The :cpp:func:`STLReader::read_file() <LAMMPS_NS::STLReader::read_file>`
member function is an alternative for the case where the mesh should be
read on MPI rank 0 only and then communicated to all MPI processes.  It
requires a class instance and thus a LAMMPS instance, aborts LAMMPS with
an error message when the file cannot be read, and stores the vertices in
an array that is managed by the class instance.

\endverbatim
*/

class STLReader : protected Pointers {
 public:
  /** Storage for a single triangle of an STL mesh */

  struct Triangle {
    double normal[3];     //!< outward facing normal vector of the facet
    double vert[3][3];    //!< x-, y-, and z-coordinates of the 3 vertices
  };

  STLReader(class LAMMPS *);
  ~STLReader() override;

  int read_file(const char *filename, double **&caller_tris);

  static std::vector<Triangle> parse(const std::string &filename, std::string *title_out = nullptr);

 private:
  static std::vector<Triangle> parse_text(const std::string &filename, std::string &title);
  static std::vector<Triangle> parse_binary(FILE *fp, std::string &title);

  int ntris, maxtris;
  double **tris;
};

}    // namespace LAMMPS_NS

#endif
