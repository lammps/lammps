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

#ifndef LMP_VTK_WRITER_H
#define LMP_VTK_WRITER_H

#include <cstdio>
#include <exception>
#include <string>
#include <utility>
#include <vector>

namespace LAMMPS_NS {

// exception thrown for inconsistent data or unwritable files

class VTKWriterException : public std::exception {
  std::string message;

 public:
  explicit VTKWriterException(std::string msg) : message(std::move(msg)) {}
  [[nodiscard]] const char *what() const noexcept override { return message.c_str(); }
};

/** Write data files in the formats used by the VTK visualization toolkit

\verbatim embed:rst

This class writes the subset of the VTK file formats that LAMMPS needs,
without requiring the VTK library itself.  Both the "legacy" format and the
XML format are supported, each with ASCII or binary encoding.  Reading VTK
files is not supported.

Usage is always the same: create a writer for the desired format, select
exactly one dataset type, attach any number of data arrays, and write.  All
inconsistencies (mismatched array lengths, missing dataset, unwritable file)
are reported by throwing a ``VTKWriterException``, so callers are expected to
catch it and turn it into a LAMMPS error.

.. code-block:: c++

   VTKWriter writer(VTKWriter::XML, false);
   writer.set_polydata(coords);
   writer.add_point_array("id", 1, ids);
   try {
     writer.write("dump.vtp");
   } catch (VTKWriterException &e) {
     error->one(FLERR, e.what());
   }

\endverbatim */

class VTKWriter {
 public:
  enum Flavor { LEGACY, XML };

  VTKWriter(Flavor flavor, bool binary);

  /** set the title written into the header of legacy files (ignored for XML) */

  void set_title(const std::string &title);

  // dataset selection.  exactly one of these must be called before write().
  // coordinates are passed as 3 doubles per point.

  /** points with one vertex cell each, written as POLYDATA / PolyData */

  void set_polydata(const std::vector<double> &xyz);

  /** points with one vertex cell each, written as UNSTRUCTURED_GRID */

  void set_unstructured_grid(const std::vector<double> &xyz);

  /** a single hexahedron cell from 8 corners, written as UNSTRUCTURED_GRID.
      corners follow the VTK ordering, i.e. the bottom face counterclockwise
      followed by the top face counterclockwise */

  void set_hexahedron(const double corners[8][3]);

  /** grid with non-uniform spacing given by the coordinates along each axis */

  void set_rectilinear_grid(const std::vector<double> &xc, const std::vector<double> &yc,
                            const std::vector<double> &zc);

  /** grid with uniform spacing, written as STRUCTURED_POINTS / ImageData */

  void set_image_data(const int dim[3], const double origin[3], const double spacing[3]);

  // data arrays.  "data" holds ncomp consecutive values per point or cell.

  void add_point_array(const std::string &name, int ncomp, const std::vector<int> &data);
  void add_point_array(const std::string &name, int ncomp, const std::vector<double> &data);
  void add_point_array(const std::string &name, const std::vector<std::string> &data);
  void add_cell_array(const std::string &name, int ncomp, const std::vector<int> &data);
  void add_cell_array(const std::string &name, int ncomp, const std::vector<double> &data);
  void add_cell_array(const std::string &name, const std::vector<std::string> &data);

  /** mark a previously added array as the default one for coloring.  it is
      written as SCALARS in legacy files and referenced by the Scalars
      attribute in XML files */

  void set_active_scalars(const std::string &name);

  /** number of points implied by the selected dataset */

  int number_of_points() const { return npoints; }

  /** number of cells implied by the selected dataset */

  int number_of_cells() const { return ncells; }

  void write(const std::string &filename);
  void write(FILE *fp);

 private:
  enum Dataset { NONE, POLYDATA, UNSTRUCTURED, RECTILINEAR, IMAGE };
  enum DataType { INT, DOUBLE, STRING };

  struct DataArray {
    std::string name;
    DataType type;
    int ncomp;
    std::vector<int> ivalues;
    std::vector<double> dvalues;
    std::vector<std::string> svalues;
  };

  Flavor flavor;
  bool binary;
  Dataset dataset;
  std::string title;
  std::string scalars;

  int npoints, ncells;
  std::vector<double> points;
  std::vector<double> xcoord, ycoord, zcoord;
  int dims[3];
  double origin[3], spacing[3];
  int celltype;

  std::vector<DataArray> point_arrays, cell_arrays;

  void set_vertex_cells(const std::vector<double> &xyz, Dataset type);
  void add_array(std::vector<DataArray> &arrays, int nitems, const char *kind, DataArray &&array);

  // legacy format

  void write_legacy(FILE *fp);
  void write_legacy_arrays(FILE *fp, const std::vector<DataArray> &arrays, const char *keyword,
                           int nitems);
  void write_legacy_array_data(FILE *fp, const DataArray &array);
  void write_legacy_cells(FILE *fp, const char *keyword);

  // XML format

  void write_xml(FILE *fp);
  void write_xml_arrays(FILE *fp, const std::vector<DataArray> &arrays, const char *tag,
                        int indent);
  void write_xml_array(FILE *fp, const DataArray &array, int indent);
  void write_xml_data_array(FILE *fp, const char *type, const std::string &name, int ncomp,
                            const std::string &payload, int indent);
  void write_xml_cells(FILE *fp, const char *tag, int indent);
};

}    // namespace LAMMPS_NS

#endif
