/* ----------------------------------------------------------------------
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

/* ----------------------------------------------------------------------
   Contributing authors:
   Daniel Queteschiner (DCS, JKU)
   Christoph Kloss (DCS)
   Richard Berger (JKU)
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(custom/vtk,DumpCustomVTK)

#else

#ifndef LMP_DUMP_CUSTOM_VTK_H
#define LMP_DUMP_CUSTOM_VTK_H

#include "dump.h"
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
 * @brief DumpCustomVTK class
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
class DumpCustomVTK : public Dump {
 public:
  DumpCustomVTK(class LAMMPS *, int, char **);
  virtual ~DumpCustomVTK();

  virtual void write();
 protected:
  int nevery;                // dump frequency for output
  char *label;               // string for dump file header
  int iregion;               // -1 if no region, else which region
  char *idregion;            // region ID
  int nthresh;               // # of defined thresholds
  int *thresh_array;         // array to threshold on for each nthresh
  int *thresh_op;            // threshold operation for each nthresh
  double *thresh_value;      // threshold value for each nthresh

  int vtk_file_format;       // which vtk file format to write (vtk, vtp, vtu ...)

  int nchoose;               // # of selected atoms
  int maxlocal;              // size of atom selection and variable arrays
  int *choose;               // local indices of selected atoms
  double *dchoose;           // value for each atom to threshhold against
  int *clist;                // compressed list of indices of selected atoms

  int nfield;                // # of keywords listed by user
  int ioptional;             // index of start of optional args

  std::map<int, int> field2index; // which compute,fix,variable calcs this field
  std::map<int, int> argindex;    // index into compute,fix scalar_atom,vector_atom
                                  // 0 for scalar_atom, 1-N for vector_atom values

  int ncompute;              // # of Compute objects used by dump
  char **id_compute;         // their IDs
  class Compute **compute;   // list of ptrs to the Compute objects

  int nfix;                  // # of Fix objects used by dump
  char **id_fix;             // their IDs
  class Fix **fix;           // list of ptrs to the Fix objects

  int nvariable;             // # of Variables used by dump
  char **id_variable;        // their names
  int *variable;             // list of indices for the Variables
  double **vbuf;             // local storage for variable evaluation

  int ncustom;               // # of custom atom properties
  char **id_custom;          // their names
  int *flag_custom;          // their data type

  int ntypes;                // # of atom types
  char **typenames;          // array of element names for each type

  // private methods

  virtual void init_style();
  virtual void write_header(bigint);
  int count();
  void pack(int *);
  virtual void write_data(int, double *);
  bigint memory_usage();

  int parse_fields(int, char **);
  void identify_vectors();
  int add_compute(char *);
  int add_fix(char *);
  int add_variable(char *);
  int add_custom(char *, int);
  virtual int modify_param(int, char **);

  typedef void (DumpCustomVTK::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_vtk(bigint);

  typedef void (DumpCustomVTK::*FnPtrWrite)(int, double *);
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

  typedef void (DumpCustomVTK::*FnPtrPack)(int);
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

  void pack_id(int);
  void pack_molecule(int);
  void pack_proc(int);
  void pack_procp1(int);
  void pack_type(int);
  void pack_mass(int);

  void pack_x(int);
  void pack_y(int);
  void pack_z(int);
  void pack_xs(int);
  void pack_ys(int);
  void pack_zs(int);
  void pack_xs_triclinic(int);
  void pack_ys_triclinic(int);
  void pack_zs_triclinic(int);
  void pack_xu(int);
  void pack_yu(int);
  void pack_zu(int);
  void pack_xu_triclinic(int);
  void pack_yu_triclinic(int);
  void pack_zu_triclinic(int);
  void pack_xsu(int);
  void pack_ysu(int);
  void pack_zsu(int);
  void pack_xsu_triclinic(int);
  void pack_ysu_triclinic(int);
  void pack_zsu_triclinic(int);
  void pack_ix(int);
  void pack_iy(int);
  void pack_iz(int);

  void pack_vx(int);
  void pack_vy(int);
  void pack_vz(int);
  void pack_fx(int);
  void pack_fy(int);
  void pack_fz(int);
  void pack_q(int);
  void pack_mux(int);
  void pack_muy(int);
  void pack_muz(int);
  void pack_mu(int);
  void pack_radius(int);
  void pack_diameter(int);

  void pack_omegax(int);
  void pack_omegay(int);
  void pack_omegaz(int);
  void pack_angmomx(int);
  void pack_angmomy(int);
  void pack_angmomz(int);
  void pack_tqx(int);
  void pack_tqy(int);
  void pack_tqz(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: No dump custom arguments specified

The dump custom command requires that atom quantities be specified to
output to dump file.

E: Invalid attribute in dump custom command

Self-explantory.

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

A dump threshhold has been requested on a quantity that is
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

Self-explantory.

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

E: Invalid dump_modify threshhold operator

Operator keyword used for threshold specification in not recognized.

*/
