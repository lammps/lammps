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

/* ----------------------------------------------------------------------
   Contributing Author: Jacob Gissinger (jgissing@stevens.edu)
------------------------------------------------------------------------- */

#ifndef LMP_REACTION_PARSER_H
#define LMP_REACTION_PARSER_H

#include "pointers.h"    // IWYU pragma: export
#include "reaction.h"

#include <array>
#include <string>
#include <vector>

namespace LAMMPS_NS {

class ReactionParser : protected Pointers {
public:
  ReactionParser(LAMMPS *lmp) : Pointers(lmp) {}

  void parse_reaction(char **, int, int, Reaction &);
  int get_chirality(double[12]);                           // get handedness given an ordered set of coordinates

private:
  static constexpr int MAXLINE = 1024;                     // max length of line read from files
  static constexpr int MAXCONARGS = 14;                    // max # of arguments for any type of constraint + rxnID

  FILE *fp;

  void validate_variable_keyword(const char *, int);
  void read_map_file(Reaction &);
  void EdgeIDs(char *, Reaction &, int);
  void Equivalences(char *, Reaction &, int);
  void DeleteAtoms(char *, Reaction &, int);
  void CreateAtoms(char *, Reaction &, int);
  void CustomCharges(int, Reaction &);
  void ChiralCenters(char *, Reaction &, int);
  void ReadWildcards(char *, Reaction &, int);
  void ReadConstraints(char *, Reaction &);
  void readID(char *, Reaction::Constraint &, Reaction &, int);
  void readline(char *);
  void parse_keyword(int, char *, char *);
  int firstint(char *, const char *);
};

}    // namespace LAMMPS_NS

#endif
