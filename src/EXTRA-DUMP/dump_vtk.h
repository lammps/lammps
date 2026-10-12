/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   This file initially came from LIGGGHTS (www.liggghts.com)
   Copyright (2014) DCS Computing GmbH, Linz
   Copyright (2015) Johannes Kepler University Linz

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(vtk,DumpVTK);
// clang-format on
#else

#ifndef LMP_DUMP_VTK_H
#define LMP_DUMP_VTK_H

#include "dump_custom.h"
#include "vtk_writer.h"

#include <map>
#include <set>

namespace LAMMPS_NS {

/**
 * @brief DumpVTK class
 *        write atom data to vtk files.
 *
 * Similar to the DumpCustom class but uses the built-in VTKWriter class to write data
 * to vtk simple legacy or xml format depending on the filename extension specified.
 * (Since this conflicts with the way binary output is specified, dump_modify allows to
 * set the binary flag for this dump command explicitly).
 * In contrast to DumpCustom class the attributes to be packed are stored in a std::map
 * to avoid duplicate entries and enforce correct ordering of vector components (except
 * for computes and fixes - these have to be given in the right order in the input script).
 * (Note: std::map elements are sorted by their keys.)
 * This dump command does not support compressed files, buffering or custom format strings,
 * multiproc is only supported by the xml formats.
 */

class DumpVTK : public DumpCustom {
 public:
  DumpVTK(class LAMMPS *, int, char **);
  ~DumpVTK() override;

  void write() override;

 protected:
  char *label;    // string for dump file header

  int vtk_file_format;    // which vtk file format to write (vtk, vtp, vtu ...)

  std::map<int, int> field2index;    // which compute,fix,variable calcs this field
  std::map<int, int> argindex;       // index into compute,fix scalar_atom,vector_atom
                                     // 0 for scalar_atom, 1-N for vector_atom values

  // private methods

  void init_style() override;
  void write_header(bigint) override;
  int count() override;
  void pack(tagint *) override;
  void write_data(int, double *) override;

  int parse_vtk_fields(int, char **);
  void identify_vectors();
  int modify_param(int, char **) override;

  using FnPtrWrite = void (DumpVTK::*)(int, double *);
  FnPtrWrite write_choice;    // ptr to write data functions
  void write_vtk(int, double *);
  void write_vtp(int, double *);
  void write_vtu(int, double *);
  void write_xml_snapshot(int, double *, bool unstructured);    // shared body of vtp/vtu
  void write_pvtk();                       // write parallel .pvtp/.pvtu summary file
  std::string pvtk_piece_filename(int);    // per-proc piece file name as referenced in summary

  void write_points(VTKWriter::Flavor, bool unstructured);    // write the atom data file
  void check_coordinate_precision(double);    // warn if single precision is too coarse
  void write_domain(VTKWriter::Flavor);       // write the box data file

  using FnPtrPack = void (DumpVTK::*)(int);
  std::map<int, FnPtrPack> pack_choice;    // ptrs to pack functions
  std::map<int, int> vtype;                // data type
  std::map<int, std::string> name;         // attribute labels
  std::set<int> vector_set;                // set of vector attributes
  int current_pack_choice_key;

  // data collected for the current snapshot.  one entry per attribute that
  // is not a point coordinate, in the order the values appear in buf.

  struct VTKArray {
    std::string name;
    int type;                            // Dump::INT, Dump::DOUBLE or Dump::STRING
    int ncomp;                           // 1 or 3
    std::vector<double> values;          // for INT and DOUBLE
    std::vector<std::string> strings;    // for STRING
  };

  std::vector<double> points;    // coordinates of the dumped atoms
  std::vector<VTKArray> myarrays;

  int n_calls_;
  int precision_warned;              // 1 after the single precision warning was printed
  VTKWriter::Precision writeprec;    // precision of floating point output, dump_modify double
  char *filecurrent;
  char *domainfilecurrent;
  char *parallelfilecurrent;
  char *multiname_ex;

  void setFileCurrent();
  void buf2arrays(int, double *);    // transfer data from buf array to vtk arrays
  void reset_vtk_data_containers();

  // customize by adding a method prototype
  void pack_vtk_compute(int);
  void pack_vtk_fix(int);
  void pack_vtk_variable(int);
  void pack_vtk_custom(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
