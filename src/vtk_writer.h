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

// exception thrown for inconsistent data or files that cannot be written

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
inconsistencies (mismatched array lengths, missing dataset, a file that
cannot be opened) are reported by throwing a ``VTKWriterException``, so callers are expected to
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
  enum Precision { SINGLE, DOUBLE };

  /** Coordinates larger than this many length units keep fewer than 3 digits
      after the decimal point when stored in single precision, which has a
      relative resolution of about 1.2e-7.  single_precision_resolution()
      scales this with force->angstrom so that the same physical resolution
      applies in every unit system. */

  static constexpr double SINGLE_PRECISION_LIMIT = 1.0e4;

  /** Resolution of single precision storage at the given coordinate magnitude.
   *
   * Callers pass the value reported by max_single_precision_value() and warn
   * with the returned resolution when it is not zero.  Keeping this test next
   * to SINGLE_PRECISION_LIMIT ensures all callers warn consistently.
   *
   * \param  maxcoord  largest coordinate magnitude written in single precision
   * \param  angstrom  length units per Angstrom, that is force->angstrom
   * \return absolute resolution at that magnitude when it warrants a warning, 0.0 otherwise */

  static double single_precision_resolution(double maxcoord, double angstrom);

  /** XML type name for floating point data at the given precision.
   *
   * Exposed so that the parallel summary files written by dump vtk always
   * declare exactly the types that the piece files contain.
   *
   * \param  precision  SINGLE or DOUBLE
   * \return "Float32" or "Float64" */

  static const char *xml_real_type(Precision precision);

  /** Byte order attribute written into XML files on this machine.
   *
   * \return "LittleEndian" or "BigEndian" */

  static const char *xml_byte_order();

  /** Create a writer for one file.
   *
   * All floating point data, that is point coordinates and data arrays alike,
   * is stored in single precision by default, which is what visualization
   * programs expect.  See max_single_precision_value() for how to tell when
   * that is not enough for the coordinates.
   *
   * \param  flavor     LEGACY for the simple legacy format, XML for the XML format
   * \param  binary     true to store the data as binary numbers instead of as text
   * \param  precision  precision used for all floating point data */

  VTKWriter(Flavor flavor, bool binary, Precision precision = SINGLE);

  /** Set the title written into the header of legacy files, ignored for XML.
   *
   * \param  title  descriptive text, truncated to the 256 characters the format allows */

  void set_title(const std::string &title);

  // dataset selection.  exactly one of these must be called before write().
  // coordinates are passed as 3 doubles per point.  coordinate and data
  // vectors are taken by value: pass with std::move() to avoid copying data
  // that is no longer needed after the call.

  /** Select points with one vertex cell each, written as POLYDATA / PolyData.
   *
   * \param  xyz  coordinates, 3 values per point */

  void set_polydata(std::vector<double> xyz);

  /** Select points with one vertex cell each, written as UNSTRUCTURED_GRID.
   *
   * \param  xyz  coordinates, 3 values per point */

  void set_unstructured_grid(std::vector<double> xyz);

  /** Select a single hexahedron cell, written as UNSTRUCTURED_GRID.
   *
   * \param  corners  the 8 corners in VTK ordering, that is the bottom face
   *                  counterclockwise followed by the top face counterclockwise */

  void set_hexahedron(const double corners[8][3]);

  /** Select a grid with non-uniform spacing.
   *
   * \param  xc  coordinates of the cell boundaries along x
   * \param  yc  coordinates of the cell boundaries along y
   * \param  zc  coordinates of the cell boundaries along z */

  void set_rectilinear_grid(std::vector<double> xc, std::vector<double> yc, std::vector<double> zc);

  /** Select a grid with uniform spacing, written as STRUCTURED_POINTS / ImageData.
   *
   * \param  dim      number of grid points in each dimension
   * \param  origin   coordinates of the first grid point
   * \param  spacing  distance between grid points in each dimension */

  void set_image_data(const int dim[3], const double origin[3], const double spacing[3]);

  // data arrays.  "data" holds ncomp consecutive values per point or cell
  // and is taken by value, so callers can std::move() it in.

  void add_point_array(const std::string &name, int ncomp, std::vector<int> data);
  void add_point_array(const std::string &name, int ncomp, std::vector<double> data);
  void add_point_array(const std::string &name, std::vector<std::string> data);
  void add_cell_array(const std::string &name, int ncomp, std::vector<int> data);
  void add_cell_array(const std::string &name, int ncomp, std::vector<double> data);
  void add_cell_array(const std::string &name, std::vector<std::string> data);

  /** Mark a previously added array as the default one for coloring.
   *
   * The array is written as SCALARS in legacy files and referenced by the
   * Scalars attribute in XML files.
   *
   * \param  name  name of an array that was added before */

  void set_active_scalars(const std::string &name);

  /** largest absolute value among the coordinates that are written in single
      precision, or 0.0 if there are none.  callers use this to warn when the
      resolution of single precision is no longer sufficient, which happens
      for very large coordinates because they need absolute resolution.  data
      arrays are not tracked: their values only need the relative resolution
      that single precision always provides.
   *
   * \return largest coordinate magnitude written in single precision, 0.0 if there is none */

  [[nodiscard]] double max_single_precision_value() const { return maxsingle; }

  /** Number of points implied by the selected dataset.
   *
   * \return number of points */

  [[nodiscard]] int number_of_points() const { return npoints; }

  /** Number of cells implied by the selected dataset.
   *
   * \return number of cells */

  [[nodiscard]] int number_of_cells() const { return ncells; }

  void write(const std::string &filename);
  void write(FILE *fp);

 private:
  enum Dataset { NONE, POLYDATA, UNSTRUCTURED, RECTILINEAR, IMAGE };
  enum DataType { TYPE_INT, TYPE_DOUBLE, TYPE_STRING };

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
  Precision prec;
  double maxsingle;
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

  void set_vertex_cells(std::vector<double> &&xyz, Dataset type);
  void add_array(std::vector<DataArray> &arrays, int nitems, const char *kind, DataArray &&array);

  // coordinates and floating point data arrays honor the selected precision.
  // only coordinates are tracked for the single precision warning, because
  // only they need absolute resolution.

  void track_single(const std::vector<double> &values);
  void write_legacy_reals(FILE *fp, const std::vector<double> &values);
  [[nodiscard]] std::string xml_reals(const std::vector<double> &values, int indent) const;
  [[nodiscard]] const char *legacy_real_type() const;
  [[nodiscard]] const char *xml_real_type() const;

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
