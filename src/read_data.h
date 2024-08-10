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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(read_data,ReadData);
// clang-format on
#else

#ifndef LMP_READ_DATA_H
#define LMP_READ_DATA_H

#include "command.h"

namespace LAMMPS_NS {
class Fix;
class ReadData : public Command {
 public:
  ReadData(class LAMMPS *);
  ~ReadData() override;
  void command(int, char **) override;
  static bool is_data_section(const std::string &);

 private:
  int me, compressed;
  char *line, *keyword, *buffer, *style;
  FILE *fp;
  char **coeffarg;
  int ncoeffarg, maxcoeffarg;
  std::string argoffset1, argoffset2;

  bigint id_offset, mol_offset;

  int nlocal_previous;
  bigint natoms;
  bigint nbonds, nangles, ndihedrals, nimpropers;
  int ntypes, nbondtypes, nangletypes, ndihedraltypes, nimpropertypes;

  bigint nellipsoids;
  class AtomVecEllipsoid *avec_ellipsoid;
  bigint nlines;
  class AtomVecLine *avec_line;
  bigint ntris;
  class AtomVecTri *avec_tri;
  bigint nbodies;
  class AtomVecBody *avec_body;

  // type labels

  class LabelMap *lmap;

  // box info read from file

  int triclinic, triclinic_general;
  int xloxhi_flag, yloyhi_flag, zlozhi_flag, tilt_flag;
  int avec_flag, bvec_flag, cvec_flag, abc_origin_flag;

  double boxlo[3], boxhi[3];
  double xy, xz, yz;
  double avec[3], bvec[3], cvec[3];
  double abc_origin[3];

  // optional args

  int addflag, offsetflag, shiftflag, coeffflag, settypeflag;
  int tlabelflag, blabelflag, alabelflag, dlabelflag, ilabelflag;
  tagint addvalue;
  int toffset, boffset, aoffset, doffset, ioffset;
  double shift[3];
  int extra_atom_types, extra_bond_types, extra_angle_types;
  int extra_dihedral_types, extra_improper_types;
  int groupbit;

  int nfix;
  Fix **fix_index;
  char **fix_header;
  char **fix_section;

  // methods

  void open(const std::string &);
  void scan(int &, int &, int &, int &);
  int reallocate(int **, int, int);
  void header(int);
  void parse_keyword(int);
  void skip_lines(bigint);
  void parse_coeffs(char *, const char *, int, int, int, int *);
  int style_match(const char *, const char *);

  void atoms();
  void velocities();

  void bonds(int);
  void bond_scan(int, char *, int *);
  void angles(int);
  void dihedrals(int);
  void impropers(int);

  void bonus(bigint, class AtomVec *, const char *);
  void bodies(int, class AtomVec *);

  void mass();
  void paircoeffs();
  void pairIJcoeffs();
  void bondcoeffs();
  void anglecoeffs(int);
  void dihedralcoeffs(int);
  void impropercoeffs(int);
  void typelabels(int);

  void fix(Fix *, char *);
};

}    // namespace LAMMPS_NS

#endif
#endif
