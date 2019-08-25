/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
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

DumpStyle(vtk,DumpVTK)

#else

#ifndef LMP_DUMP_VTK_H
#define LMP_DUMP_VTK_H

#include "dump_custom.h"
#include <map>
#include <set>
#include <string>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>

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
  virtual ~DumpVTK();

  virtual void write();
 protected:
  char *label;               // string for dump file header

  int vtk_file_format;       // which vtk file format to write (vtk, vtp, vtu ...)

  std::map<int, int> field2index; // which compute,fix,variable calcs this field
  std::map<int, int> argindex;    // index into compute,fix scalar_atom,vector_atom
                                  // 0 for scalar_atom, 1-N for vector_atom values

  // private methods

  virtual void init_style();
  virtual void write_header(bigint);
  int count();
  void pack(tagint *);
  virtual void write_data(int, double *);
  bigint memory_usage();

  int parse_fields(int, char **);
  void identify_vectors();
  int add_compute(char *);
  int add_fix(char *);
  int add_variable(char *);
  int add_custom(char *, int);
  virtual int modify_param(int, char **);

  typedef void (DumpVTK::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_vtk(bigint);

  typedef void (DumpVTK::*FnPtrWrite)(int, double *);
  FnPtrWrite write_choice;             // ptr to write data functions
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
  std::map<int, FnPtrPack> pack_choice;  // ptrs to pack functions
  std::map<int, int> vtype;              // data type
  std::map<int, std::string> name;       // attribute labels
  std::set<int> vector_set;              // set of vector attributes
  int current_pack_choice_key;

  // vtk data containers
  vtkSmartPointer<vtkPoints> points;
  vtkSmartPointer<vtkCellArray> pointsCells;
  std::map<int, vtkSmartPointer<vtkAbstractArray> > myarrays;

  int n_calls_;
  double (*boxcorners)[3]; // corners of triclinic domain box
  char *filecurrent;
  char *domainfilecurrent;
  char *parallelfilecurrent;
  char *multiname_ex;

  void setFileCurrent();
  void buf2arrays(int, double *); // transfer data from buf array to vtk arrays
  void reset_vtk_data_containers();

  // customize by adding a method prototype
  void pack_compute(int);
  void pack_fix(int);
  void pack_variable(int);
  void pack_custom(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: No dump custom arguments specified

The dump custom command requires that atom quantities be specified to
output to dump file.

E: Invalid attribute in dump custom command

Self-explanatory.

E: Dump_modify format string is too short

There are more fields to be dumped in a line of output than your
format string specifies.

E: Could not find dump custom compute ID

Self-explanatory.

E: Could not find dump custom fix ID

Self-explanatory.

E: Dump custom and fix not computed at compatible times

The fix must produce per-atom quantities on timesteps that dump custom
needs them.

E: Could not find dump custom variable name

Self-explanatory.

E: Could not find custom per-atom property ID

Self-explanatory.

E: Region ID for dump custom does not exist

Self-explanatory.

E: Compute used in dump between runs is not current

The compute was not invoked on the current timestep, therefore it
cannot be used in a dump between runs.

E: Threshhold for an atom property that isn't allocated

A dump threshold has been requested on a quantity that is
not defined by the atom style used in this simulation.

E: Dumping an atom property that isn't allocated

The chosen atom style does not define the per-atom quantity being
dumped.

E: Dump custom compute does not compute per-atom info

Self-explanatory.

E: Dump custom compute does not calculate per-atom vector

Self-explanatory.

E: Dump custom compute does not calculate per-atom array

Self-explanatory.

E: Dump custom compute vector is accessed out-of-range

Self-explanatory.

E: Dump custom fix does not compute per-atom info

Self-explanatory.

E: Dump custom fix does not compute per-atom vector

Self-explanatory.

E: Dump custom fix does not compute per-atom array

Self-explanatory.

E: Dump custom fix vector is accessed out-of-range

Self-explanatory.

E: Dump custom variable is not atom-style variable

Only atom-style variables generate per-atom quantities, needed for
dump output.

E: Custom per-atom property ID is not floating point

Self-explanatory.

E: Custom per-atom property ID is not integer

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Dump_modify region ID does not exist

Self-explanatory.

E: Dump modify element names do not match atom types

Number of element names must equal number of atom types.

E: Invalid attribute in dump modify command

Self-explanatory.

E: Could not find dump modify compute ID

Self-explanatory.

E: Dump modify compute ID does not compute per-atom info

Self-explanatory.

E: Dump modify compute ID does not compute per-atom vector

Self-explanatory.

E: Dump modify compute ID does not compute per-atom array

Self-explanatory.

E: Dump modify compute ID vector is not large enough

Self-explanatory.

E: Could not find dump modify fix ID

Self-explanatory.

E: Dump modify fix ID does not compute per-atom info

Self-explanatory.

E: Dump modify fix ID does not compute per-atom vector

Self-explanatory.

E: Dump modify fix ID does not compute per-atom array

Self-explanatory.

E: Dump modify fix ID vector is not large enough

Self-explanatory.

E: Could not find dump modify variable name

Self-explanatory.

E: Dump modify variable is not atom-style variable

Self-explanatory.

E: Could not find dump modify custom atom floating point property ID

Self-explanatory.

E: Could not find dump modify custom atom integer property ID

Self-explanatory.

E: Invalid dump_modify threshold operator

Operator keyword used for threshold specification in not recognized.

*/
