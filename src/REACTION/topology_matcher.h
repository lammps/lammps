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

#ifndef LMP_TOPOLOGY_MATCHER_H
#define LMP_TOPOLOGY_MATCHER_H

#include "pointers.h"    // IWYU pragma: export
#include "reaction.h"
#include "reaction_constraints.h"

#include <vector>

namespace LAMMPS_NS {

class TopologyMatcher : protected Pointers {
public:
  TopologyMatcher(class LAMMPS *);
  ~TopologyMatcher();

  ReactionConstraints *rxn_constraints;

  bool match_topology(std::vector<tagint> &, Reaction &, std::array<tagint, 2>);

private:
  static constexpr int MAXGUESS = 20;                      // max # of guesses allowed by superimpose algorithm

  enum class Status { ACCEPT, REJECT, PROCEED,
                      CONTINUE, GUESSFAIL, RESTORE };      // values for superimpose algorithm status
  Status status;

  struct Superimpose {
    int avail_guesses;                                     // num of restore points available
    std::vector<int> guess_branch;                         // used when there is more than two choices when guessing
    struct StatePoint {
      int pion, neigh, trace, glove_counter;
      std::vector<tagint> glove, pioneer_count, pioneers;
    } sp;
  };
  std::vector<Superimpose::StatePoint> restore_pts;

  void make_a_guess(Superimpose &, Reaction &rxn);
  void neighbor_loop(Superimpose &, Reaction &rxn);
  void check_a_neighbor(Superimpose &, Reaction &rxn);
  void crosscheck_the_neighbor(Superimpose &, Reaction &rxn);
  void inner_crosscheck_loop(Superimpose &, Reaction &rxn);
  int ring_check(Reaction &rxn, const std::vector<tagint> &);
  bool compare_atomtype(int, Reaction &rxn, int);
};

}    // namespace LAMMPS_NS

#endif
