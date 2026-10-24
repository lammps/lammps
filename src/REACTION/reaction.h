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

#ifndef LMP_REACTION_H
#define LMP_REACTION_H

#include "pointers.h"    // IWYU pragma: export

#include <array>
#include <string>
#include <vector>

namespace LAMMPS_NS {
static constexpr int MAXNAME = 256;                      // max character length of react-ID

enum class Molecule_Keys { OFF, INTER, INTRA };          // values for molecule_keyword

struct ReactionTemplate {
  std::string name;                                      // name of reaction template
  class Molecule *reactant;                              // pre-reacted molecule template
  class Molecule *product;                               // post-reacted molecule template
  int iatomtype, jatomtype;                              // initiator atom types
  int ibonding, jbonding;                                // initiator atom IDs
  int closeneigh;                                        // indicates if initiator atoms are 1-2, 1-3, or 1-4 neighbors
  std::string constraintstr;                             // boolean logic for multiple constraints
  std::string mapfilename;                               // name of the map file
  int custom_charges_fragid;                             // molecule fragment defining which atomic charges to update
  int rescale_charges_flag;                              // if nonzero, indicates number of atoms whose charges are updated
  double mol_total_charge;                               // sum of charges of post-reaction atoms whose charges are updated
  int nnewmolids;                                        // number of unique new molids needed for each reaction
  int create_atoms_flag, modify_create_fragid;           // molecule fragment defining which atoms are created
  double overlapsq;                                      // squared distance for overlap detection during atom creation

  struct ReactionAtomFlags {
    int edge;                                            // true if atom in molecule template has incorrect valency
    int landlocked;                                      // true if atom is at least three bonds away from edge atoms
    bool wildcard;                                       // true if atom type contains a wildcard
    int recharged;                                       // true if atom whose charge should be updated
    int deleted;                                         // true if atom in pre-reacted template to delete
    int created;                                         // true if atom in post-reacted template to create
    int newmolid;                                        // for molmap option: mol IDs in post, but not in pre, re-indexed from 1
    std::array<int, 6> chiral;                           // pre-react chiral atoms. 1) flag 2) orientation 3-4) ordered atom types
    std::array<int, 2> amap;                             // atom map: clmn 1 = product atom IDs, clmn 2: reactant atom IDs
    std::array<int, 2> ramap;                            // reverse amap
  };
  std::vector<ReactionAtomFlags> atoms;

  struct Constraint {
    int ID;
    enum class Type { DISTANCE, ANGLE, DIHEDRAL, ARRHENIUS, RMSD, CUSTOM } type;
    struct Distance { double rminsq, rmaxsq; } distance;
    struct Angle { double amin, amax; } angle;
    struct Dihedral { double amin, amax, amin2, amax2; } dihedral;
    struct RMSD { double rmsdmax; } rmsd;
    struct Arrhenius { double A, n, E_a, seed; class RanMars *rrhandom; } arrhenius;
    struct Custom { std::string str; } custom;
    enum class IDType { ATOM, FRAG };
    static constexpr int MAXCONIDS = 4;
    std::array<int, MAXCONIDS> ids;
    std::array<IDType, MAXCONIDS> idtypes{};
    bool satisfied;
  };
  std::vector<Constraint> constraints;
};

struct Reaction : public ReactionTemplate {
  int ID;                                                // indexed from 0
  int nevery, groupbits;
  double rminsq, rmaxsq;
  double fraction;
  int reaction_count, reaction_count_total;
  int local_rxn_count, ghostly_rxn_count;
  int nlocalkeep, nghostlykeep;
  int seed, limit_duration;
  int stabilize_steps_flag;
  Molecule_Keys molecule_keyword;
  int v_nevery, v_rmin, v_rmax, v_prob;                  // ID of variable, -1 if static
  std::vector<std::array<tagint, 2>> attempts;           // stores sim atom IDs of initiator atoms
};

}    // namespace LAMMPS_NS

#endif
