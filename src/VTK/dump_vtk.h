/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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
#include <map>
#include <set>

#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>

class vtkAbstractArray;
class vtkRectilinearGrid;
class vtkUnstructuredGrid;

namespace LAMMPS_NS {

/**
 * @brief DumpVTK class
 *        write atom data to vtk files.
 *
 * Similar to the DumpCustom class but uses the vtk library to write data to vtk simple
 * legacy or xml format depending on the filename extension specified. (Since this
 * conflicts with the way binary output is specified, dump_modify allows to set the
 * binary flag for this dump command explicitly).
 * In contrast to DumpCustom class the attributes to be packed are stored in a std::map
 * to avoid duplicate entries and enforce correct ordering of vector components (except
 * for computes and fixes - these have to be given in the right order in the input script).
 * (Note: std::map elements are sorted by their keys.)
 * This dump command does not support compressed files, buffering or custom format strings,
 * multiproc is only supported by the xml formats, multifile option has to be used.
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
  double memory_usage() override;

  int parse_fields(int, char **);
  void identify_vectors();
  int add_compute(const char *);
  int add_fix(const char *);
  int add_variable(const char *);
  int add_custom(const char *, int);
  int modify_param(int, char **) override;

  typedef void (DumpVTK::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;    // ptr to write header functions
  void header_vtk(bigint);

  typedef void (DumpVTK::*FnPtrWrite)(int, double *);
  FnPtrWrite write_choice;    // ptr to write data functions
  void write_vtk(int, double *);
  void write_vtp(int, double *);
  void write_vtu(int, double *);

  void prepare_domain_data(vtkRectilinearGrid *);
  void prepare_domain_data_triclinic(vtkUnstructuredGrid *);
  void write_domain_vtk();
  void write_domain_vtk_triclinic();
  void write_domain_vtr();
  void write_domain_vtu_triclinic();

  typedef void (DumpVTK::*FnPtrPack)(int);
  std::map<int, FnPtrPack> pack_choice;    // ptrs to pack functions
  std::map<int, int> vtype;                // data type
  std::map<int, std::string> name;         // attribute labels
  std::set<int> vector_set;                // set of vector attributes
  int current_pack_choice_key;

  // vtk data containers
  vtkSmartPointer<vtkPoints> points;
  vtkSmartPointer<vtkCellArray> pointsCells;
  std::map<int, vtkSmartPointer<vtkAbstractArray>> myarrays;

  int n_calls_;
  double (*boxcorners)[3];    // corners of triclinic domain box
  char *filecurrent;
  char *domainfilecurrent;
  char *parallelfilecurrent;
  char *multiname_ex;

  void setFileCurrent();
  void buf2arrays(int, double *);    // transfer data from buf array to vtk arrays
  void reset_vtk_data_containers();

  // customize by adding a method prototype
  void pack_compute(int);
  void pack_fix(int);
  void pack_variable(int);
  void pack_custom(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
