/* ----------------------------------------------------------------------
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

#include "topology_matcher.h"

#include "error.h"
#include "molecule.h"

#include <algorithm>

using namespace LAMMPS_NS;

TopologyMatcher::TopologyMatcher(LAMMPS *lmp, bool constraints_flag, bool limit_tags_flag) :
  Pointers(lmp), constraints_flag(constraints_flag), limit_tags_flag(limit_tags_flag),
  default_topo(lmp), topo(&default_topo)
{

  // checking reaction constraints inside topology_matcher

  if (constraints_flag) rxn_constraints = new ReactionConstraints(lmp);

  status = Status::PROCEED;
}

TopologyMatcher::~TopologyMatcher()
{
  if (constraints_flag) delete rxn_constraints;
}

/* ----------------------------------------------------------------------
  match_topology arguments:
    super: structure to hold global indices needed for subroutines
    rxn: reaction instance to be matched (or not)
    rxn_attempt: atom IDs of two simulation particles corresponding to initiator atoms
  Output: magic string "SUCCESS" for successful topology match,
          else empty string for no topology match, else error message
  Overall summary:
  Set up global variables. Loop through all pairs; loop through Pioneers
  until Superimpose Algorithm is completed for each pair.
  various statuses of superimpose algorithm:
    ACCEPT: site successfully matched to pre-reacted template
    REJECT: site does not match pre-reacted template
    PROCEED: normal execution (non-guessing mode)
    CONTINUE: a neighbor has been assigned, skip to next neighbor
    GUESSFAIL: a guess has failed (if no more restore points, status = 'REJECT')
    RESTORE: restore mode, load most recent restore point
------------------------------------------------------------------------- */

std::string TopologyMatcher::match_topology(std::vector<tagint> &outglove, Reaction &rxn, std::array<tagint, 2> rxn_attempt)
{
  status = Status::PROCEED;

  // quick description of important global indices you'll see floating about:
  // pion: the pioneer loop index
  // neigh: in the first neighbor index
  // trace: retraces the first neighbors. once you choose a first neighbor, you then check for other neighbors of same type
  // pioneers: during Superimpose Algorithm, atoms which have been assigned, but whose first neighbors haven't
  // glove: global IDs. index indicates is pre-reaction ID-1, value is mapped sim atom ID
  // glove_counter: used to determine when to terminate Superimpose Algorithm

  Superimpose super;
  Superimpose::StatePoint &sp = super.sp;

  sp.glove.resize(rxn.reactant->natoms);
  sp.pioneers.resize(rxn.reactant->natoms);
  sp.pioneer_count.resize(rxn.reactant->natoms);

  restore_pts.resize(MAXGUESS);
  for (auto &restore_pt : restore_pts) {
    restore_pt.glove.resize(rxn.reactant->natoms);
    restore_pt.pioneers.resize(rxn.reactant->natoms);
    restore_pt.pioneer_count.resize(rxn.reactant->natoms);
  }

  sp.pion = sp.neigh = sp.trace = sp.glove_counter = 0;
  std::fill(sp.glove.begin(), sp.glove.end(), 0);

  std::vector<int> &guess_branch = super.guess_branch;
  guess_branch.resize(MAXGUESS, 0);
  std::fill(guess_branch.begin(), guess_branch.end(), 0);

  sp.glove[rxn.ibonding-1] = rxn_attempt[0];
  sp.glove_counter++;
  sp.glove[rxn.jbonding-1] = rxn_attempt[1];
  sp.glove_counter++;

  // special case, only two atoms in reaction templates
  // then: bonding reactant nspecials guaranteed to be equal, and either 0 or 1
  if (sp.glove_counter == rxn.reactant->natoms) {
    tagint atom1 = sp.glove[rxn.ibonding-1];
    tagint atom2 = sp.glove[rxn.jbonding-1];
    if ( (topo->nxspecial(atom1) == rxn.reactant->nspecial[rxn.ibonding-1][0] &&
          topo->nxspecial(atom2) == topo->nxspecial(atom1)) &&
         (topo->nxspecial(atom1) == 0 ||
          topo->xspecial(atom1,0) == atom->tag[atom2]) &&
         constraints_check(rxn, sp.glove) ) status = Status::ACCEPT;
    else status = Status::REJECT;
  }

  super.avail_guesses = 0;

  std::fill(sp.pioneer_count.begin(), sp.pioneer_count.end(), 0);

  for (int i = 0; i < rxn.reactant->nspecial[rxn.ibonding-1][0]; i++)
    sp.pioneer_count[rxn.reactant->special[rxn.ibonding-1][i]-1]++;

  for (int i = 0; i < rxn.reactant->nspecial[rxn.jbonding-1][0]; i++)
    sp.pioneer_count[rxn.reactant->special[rxn.jbonding-1][i]-1]++;


  int hang_catch = 0;
  while (status != Status::ACCEPT && status != Status::REJECT) {

    std::fill(sp.pioneers.begin(), sp.pioneers.end(), 0);

    int npioneers = 0;
    for (int i = 0; i < rxn.reactant->natoms; i++) {
      if (sp.glove[i] != 0 && sp.pioneer_count[i] < rxn.reactant->nspecial[i][0] && rxn.atoms[i].edge == 0) {

        // make sure all neighbors aren't already assigned
        // an issue discovered for coarse-grained example
        int assigned_count = 0;
        int nfirst_neighs = rxn.reactant->nspecial[i][0];
        for (int j = 0; j < nfirst_neighs; j++) {
          for (int k = 0; k < rxn.reactant->natoms; k++) {
            if (topo->xspecial(sp.glove[i],j) == sp.glove[k]) {
              assigned_count++;
              break;
            }
          }
        }
        if (assigned_count == nfirst_neighs) continue;

        sp.pioneers[i] = 1;
        npioneers++;
      }
    }
    if (npioneers == 0) status = Status::GUESSFAIL;

    // run through the pioneers
    // due to use of restore points, 'pion' index can change in loop
    for (sp.pion = 0; sp.pion < rxn.reactant->natoms; sp.pion++) {
      if (sp.pioneers[sp.pion] || status == Status::GUESSFAIL) {
        auto err_msg = make_a_guess(super, rxn);
        if (!err_msg.empty()) return err_msg;
        if (status == Status::ACCEPT || status == Status::REJECT) break;
      }
    }

    hang_catch++;
    // let's go ahead and catch the simplest of hangs
    //if (hang_catch > rxn.reactant->natoms*4)
    if (hang_catch > atom->nlocal*30) {
      return "Excessive iteration of superimpose algorithm. "
          "Please check that all pre-reaction template atoms are linked to an initiator atom, "
          "via at least one path that does not involve edge atoms.";
    }
  }
  if (status == Status::ACCEPT) {
    outglove = sp.glove;
    return "SUCCESS";
  } return {};
}

/* ----------------------------------------------------------------------
  Screen for obvious algorithm fails. This is the return point when a guess
  has failed: check for available restore points.
------------------------------------------------------------------------- */

std::string TopologyMatcher::make_a_guess(Superimpose &super, Reaction &rxn)
{
  Superimpose::StatePoint &sp = super.sp;
  int &avail_guesses = super.avail_guesses;

  int nfirst_neighs = rxn.reactant->nspecial[sp.pion][0];

  // per-atom property indicating if in bond/react master group
  int *i_limit_tags;
  if (limit_tags_flag) {
    int flag,cols;
    int index1 = atom->find_custom("limit_tags",flag,cols); // TODO: this needs to be an input
    i_limit_tags = atom->ivector[index1];
  }

  if (status == Status::GUESSFAIL && avail_guesses == 0) {
    status = Status::REJECT;
    return {};
  }

  if (status == Status::GUESSFAIL && avail_guesses > 0) {
    // load restore point
    for (int i = 0; i < rxn.reactant->natoms; i++) {
      sp.glove[i] = restore_pts[avail_guesses-1].glove[i];
      sp.pioneer_count[i] = restore_pts[avail_guesses-1].pioneer_count[i];
      sp.pioneers[i] = restore_pts[avail_guesses-1].pioneers[i];
    }
    sp.pion = restore_pts[avail_guesses-1].pion;
    sp.neigh = restore_pts[avail_guesses-1].neigh;
    sp.trace = restore_pts[avail_guesses-1].trace;
    sp.glove_counter = restore_pts[avail_guesses-1].glove_counter;
    status = Status::RESTORE;
    auto err_msg = neighbor_loop(super, rxn);
    if (!err_msg.empty()) return err_msg;
    if (status != Status::PROCEED) return {};
  }

  nfirst_neighs = rxn.reactant->nspecial[sp.pion][0];

  //  check if any of first neighbors are in bond_react_MASTER_group
  //  if so, this constitutes a fail
  //  because still undergoing a previous reaction!
  //  could technically fail unnecessarily during a wrong guess if near edge atoms
  //  we accept this temporary and infrequent decrease in reaction occurrences

  for (int i = 0; i < topo->nxspecial(sp.glove[sp.pion]); i++) {
    if (atom->map(topo->xspecial(sp.glove[sp.pion],i)) < 0) {
      return "Superimpose algorithm needs ghost atoms from further away"; // parallel issues.
    }
    if (limit_tags_flag && i_limit_tags[(int)atom->map(topo->xspecial(sp.glove[sp.pion],i))] != 0) {
      status = Status::GUESSFAIL;
      return {};
    }
  }

  // check for same number of neighbors between unreacted mol and simulation
  if (nfirst_neighs != topo->nxspecial(sp.glove[sp.pion])) {
    status = Status::GUESSFAIL;
    return {};
  }

  // check if all neigh atom types are the same between simulation and unreacted mol
  std::multiset<int> mol_types, lcl_types;
  for (int i = 0; i < nfirst_neighs; i++) {
    int imolatom = (int) rxn.reactant->special[sp.pion][i]-1;
    if (!rxn.atoms[imolatom].wildcard) mol_types.insert(rxn.reactant->type[imolatom]);
    lcl_types.insert(topo->type(topo->xspecial(sp.glove[sp.pion],i)));
  }

  if (!std::includes(lcl_types.begin(), lcl_types.end(), mol_types.begin(), mol_types.end())) {
    status = Status::GUESSFAIL;
    return {};
  }

  // okay everything seems to be in order. let's assign some ID pairs!!!
  auto err_msg = neighbor_loop(super, rxn);
  if (!err_msg.empty()) return err_msg;
  return {};
}

/* ----------------------------------------------------------------------
  Loop through all First Bonded Neighbors of the current Pioneer.
  Prepare appropriately if we are in Restore Mode.
------------------------------------------------------------------------- */

std::string TopologyMatcher::neighbor_loop(Superimpose &super, Reaction &rxn)
{
  Superimpose::StatePoint &sp = super.sp;

  int nfirst_neighs = rxn.reactant->nspecial[sp.pion][0];

  if (status == Status::RESTORE) {
    auto err_msg = check_a_neighbor(super, rxn);
    if (!err_msg.empty()) return err_msg;
    return {};
  }

  for (sp.neigh = 0; sp.neigh < nfirst_neighs; sp.neigh++) {
    if (sp.glove[(int)rxn.reactant->special[sp.pion][sp.neigh]-1] == 0) {
      auto err_msg = check_a_neighbor(super, rxn);
      if (!err_msg.empty()) return err_msg;
    }
  }
  return {};
  // status should still = PROCEED
}

/* ----------------------------------------------------------------------
  Check if we can assign this First Neighbor to pre-reacted template
  without guessing. If so, do it! If not, call crosscheck_the_nieghbor().
------------------------------------------------------------------------- */

std::string TopologyMatcher::check_a_neighbor(Superimpose &super, Reaction &rxn)
{
  Superimpose::StatePoint &sp = super.sp;

  int nfirst_neighs = rxn.reactant->nspecial[sp.pion][0];

  if (status != Status::RESTORE) {
    // special consideration for hydrogen atoms (and all first neighbors bonded to no other atoms) (and aren't edge atoms)
    if (rxn.reactant->nspecial[(int)rxn.reactant->special[sp.pion][sp.neigh]-1][0] == 1 && rxn.atoms[(int)rxn.reactant->special[sp.pion][sp.neigh]-1].edge == 0) {

      for (int i = 0; i < nfirst_neighs; i++) {

        int checktype = topo->type(topo->xspecial(sp.glove[sp.pion],i));
        int reactant_atom = (int) rxn.reactant->special[sp.pion][sp.neigh];
        if (compare_atomtype(checktype, rxn, reactant_atom) &&
            topo->nxspecial(topo->xspecial(sp.glove[sp.pion],i)) == 1) {

          int already_assigned = 0;
          for (int j = 0; j < rxn.reactant->natoms; j++) {
            if (sp.glove[j] == topo->xspecial(sp.glove[sp.pion],i)) {
              already_assigned = 1;
              break;
            }
          }

          if (already_assigned == 0) {
            sp.glove[(int)rxn.reactant->special[sp.pion][sp.neigh]-1] = topo->xspecial(sp.glove[sp.pion],i);

            //another check for ghost atoms. perhaps remove the one in make_a_guess
            if (atom->map(sp.glove[(int)rxn.reactant->special[sp.pion][sp.neigh]-1]) < 0) {
              return "Superimpose algorithm needs ghost atoms from further away";
            }

            for (int j = 0; j < rxn.reactant->nspecial[rxn.reactant->special[sp.pion][sp.neigh]-1][0]; j++) {
              sp.pioneer_count[rxn.reactant->special[rxn.reactant->special[sp.pion][sp.neigh]-1][j]-1]++;
            }

            sp.glove_counter++;
            if (sp.glove_counter == rxn.reactant->natoms) {
              if (ring_check(rxn, sp.glove) && constraints_check(rxn, sp.glove)) status = Status::ACCEPT;
              else status = Status::GUESSFAIL;
              return {};
            }
            // status should still == PROCEED
            return {};
          }
        }
      }
      // we are here if no matching atom found
      status = Status::GUESSFAIL;
      return {};
    }
  }

  crosscheck_the_neighbor(super, rxn);
  if (status != Status::PROCEED) {
    if (status == Status::CONTINUE)
      status = Status::PROCEED;
    return {};
  }

  // finally ready to match non-duplicate, non-edge atom IDs!!

  for (int i = 0; i < nfirst_neighs; i++) {

    int checktype = topo->type(topo->xspecial(sp.glove[sp.pion],i));
    int reactant_atom = (int) rxn.reactant->special[sp.pion][sp.neigh];
    if (compare_atomtype(checktype, rxn, reactant_atom)) {
      int already_assigned = 0;

      //check if a first neighbor of the pioneer is already assigned to pre-reacted template
      for (int j = 0; j < rxn.reactant->natoms; j++) {
        if (sp.glove[j] == topo->xspecial(sp.glove[sp.pion],i)) {
          already_assigned = 1;
          break;
        }
      }

      if (already_assigned == 0) {
        sp.glove[(int)rxn.reactant->special[sp.pion][sp.neigh]-1] = topo->xspecial(sp.glove[sp.pion],i);

        //another check for ghost atoms. perhaps remove the one in make_a_guess
        if (atom->map(sp.glove[(int)rxn.reactant->special[sp.pion][sp.neigh]-1]) < 0) {
          return "Superimpose algorithm needs ghost atoms from further away";
        }

        for (int ii = 0; ii < rxn.reactant->nspecial[rxn.reactant->special[sp.pion][sp.neigh]-1][0]; ii++) {
          sp.pioneer_count[rxn.reactant->special[rxn.reactant->special[sp.pion][sp.neigh]-1][ii]-1]++;
        }

        sp.glove_counter++;
        if (sp.glove_counter == rxn.reactant->natoms) {
          if (ring_check(rxn, sp.glove) && constraints_check(rxn, sp.glove)) status = Status::ACCEPT;
          else status = Status::GUESSFAIL;
          return {};
          // will never complete here when there are edge atoms
          // ...actually that could be wrong if people get creative...shouldn't affect anything
        }
        // status should still = PROCEED
        return {};
      }
    }
  }
  return {};
  // status is still 'PROCEED' if we are here!
}

/* ----------------------------------------------------------------------
  Check if there a viable guess to be made. If so, prepare to make a
  guess by recording a restore point.
------------------------------------------------------------------------- */

std::string TopologyMatcher::crosscheck_the_neighbor(Superimpose &super, Reaction &rxn)
{
  Superimpose::StatePoint &sp = super.sp;
  int &avail_guesses = super.avail_guesses;

  int nfirst_neighs = rxn.reactant->nspecial[sp.pion][0];

  if (status == Status::RESTORE) {
    inner_crosscheck_loop(super, rxn);
    return {};
  }

  for (sp.trace = 0; sp.trace < nfirst_neighs; sp.trace++) {
    int reactant_atom1 = (int) rxn.reactant->special[sp.pion][sp.trace];
    int atom1type = rxn.reactant->type[reactant_atom1-1];
    int reactant_atom2 = (int) rxn.reactant->special[sp.pion][sp.neigh];
    int atom2type = rxn.reactant->type[reactant_atom2-1];
    if (sp.neigh != sp.trace && (compare_atomtype(atom2type, rxn, reactant_atom1) || compare_atomtype(atom1type, rxn, reactant_atom2)) &&
        sp.glove[rxn.reactant->special[sp.pion][sp.trace]-1] == 0) {

      if (avail_guesses == MAXGUESS) {
        error->warning(FLERR,"Superimpose algorithm failed because MAXGUESS set too small. ask developer for info");
        status = Status::GUESSFAIL;
        return {};
      }
      avail_guesses++;
      for (int i = 0; i < rxn.reactant->natoms; i++) {
        restore_pts[avail_guesses-1].glove[i] = sp.glove[i];
        restore_pts[avail_guesses-1].pioneer_count[i] = sp.pioneer_count[i];
        restore_pts[avail_guesses-1].pioneers[i] = sp.pioneers[i];
      }
      restore_pts[avail_guesses-1].pion = sp.pion;
      restore_pts[avail_guesses-1].neigh = sp.neigh;
      restore_pts[avail_guesses-1].trace = sp.trace;
      restore_pts[avail_guesses-1].glove_counter = sp.glove_counter;

      inner_crosscheck_loop(super, rxn);
      return {};
    }
  }
  return {};
  // status is still 'PROCEED' if we are here!
}

/* ----------------------------------------------------------------------
  We are ready to make a guess. If there are multiple possible choices
  for this guess, keep track of these.
------------------------------------------------------------------------- */

std::string TopologyMatcher::inner_crosscheck_loop(Superimpose &super, Reaction &rxn)
{
  Superimpose::StatePoint &sp = super.sp;
  int &avail_guesses = super.avail_guesses;
  std::vector<int> &guess_branch = super.guess_branch;

  // arbitrarily limited to 5 identical first neighbors
  tagint tag_choices[5];
  int nfirst_neighs = rxn.reactant->nspecial[sp.pion][0];

  int num_choices = 0;
  for (int i = 0; i < nfirst_neighs; i++) {
    int checktype = topo->type(topo->xspecial(sp.glove[sp.pion],i));
    int reactant_atom = (int) rxn.reactant->special[sp.pion][sp.neigh];
    if (compare_atomtype(checktype, rxn, reactant_atom)) {
      if (num_choices == 5) { // here failed because too many identical first neighbors. but really no limit if situation arises
        status = Status::GUESSFAIL;
        return {};
      }
      tag_choices[num_choices++] = topo->xspecial(sp.glove[sp.pion],i);
    }
  }

  // guess branch is for when multiple identical neighbors. then we guess each one in turn
  // guess_branch must work even when avail_guesses = 0. so index accordingly!
  // ...actually, avail_guesses should never be zero here anyway
  if (guess_branch[avail_guesses-1] == 0) guess_branch[avail_guesses-1] = num_choices;

  for (int i=1; i < num_choices; ++i) {
    tagint hold = tag_choices[i];
    int j = i - 1;
    while ((j >= 0) && (tag_choices[j] > hold)) {
      tag_choices[j+1] = tag_choices[j];
      --j;
    }
    tag_choices[j+1] = hold;
  }

  for (int i = guess_branch[avail_guesses-1]-1; i >= 0; i--) {
    int already_assigned = 0;
    for (int j = 0; j < rxn.reactant->natoms; j++) {
      if (sp.glove[j] == tag_choices[i]) {
        already_assigned = 1;
        break;
      }
    }
    if (already_assigned == 1) {
      guess_branch[avail_guesses-1]--;
      if (guess_branch[avail_guesses-1] == 0) {
        status = Status::REJECT;
        return {};
      }
    } else {
      sp.glove[rxn.reactant->special[sp.pion][sp.neigh]-1] = tag_choices[i];
      guess_branch[avail_guesses-1]--;
      break;
    }
  }

  //another check for ghost atoms. perhaps remove the one in make_a_guess
  if (atom->map(sp.glove[(int)rxn.reactant->special[sp.pion][sp.neigh]-1]) < 0) {
    return "Superimpose algorithm needs ghost atoms from further away";
  }

  if (guess_branch[avail_guesses-1] == 0) avail_guesses--;

  for (int i = 0; i < rxn.reactant->nspecial[rxn.reactant->special[sp.pion][sp.neigh]-1][0]; i++) {
    sp.pioneer_count[rxn.reactant->special[rxn.reactant->special[sp.pion][sp.neigh]-1][i]-1]++;
  }
  sp.glove_counter++;
  if (sp.glove_counter == rxn.reactant->natoms) {
    if (ring_check(rxn, sp.glove) && constraints_check(rxn, sp.glove)) status = Status::ACCEPT;
    else status = Status::GUESSFAIL;
    return {};
  }
  status = Status::CONTINUE;
  return {};
}

/* ----------------------------------------------------------------------
  Check constraints if check_constraint flag is true
------------------------------------------------------------------------- */

bool TopologyMatcher::constraints_check(Reaction &rxn, const std::vector<tagint> &glove)
{
  if (!constraints_flag) return true;
  else return rxn_constraints->check(rxn, glove);
}

/* ----------------------------------------------------------------------
  Check that newly assigned atoms have correct bonds
  Necessary for certain ringed structures
------------------------------------------------------------------------- */

bool TopologyMatcher::ring_check(Reaction &rxn, const std::vector<tagint> &glove)
{
  // ring_check can be made more efficient by re-introducing 'frozen' atoms
  // 'frozen' atoms have been assigned and also are no longer pioneers

  // double check the number of neighbors match for all non-edge atoms
  // otherwise, atoms at 'end' of symmetric ring can behave like edge atoms
  for (int i = 0; i < rxn.reactant->natoms; i++)
    if (rxn.atoms[i].edge == 0 &&
        rxn.reactant->nspecial[i][0] != topo->nxspecial(glove[i]))
      return false;

  for (int i = 0; i < rxn.reactant->natoms; i++) {
    for (int j = 0; j < rxn.reactant->nspecial[i][0]; j++) {
      int ring_fail = 1;
      int ispecial = rxn.reactant->special[i][j];
      for (int k = 0; k < topo->nxspecial(glove[i]); k++) {
        if (topo->xspecial(glove[i],k) == glove[ispecial-1]) {
          ring_fail = false;
          break;
        }
      }
      if (ring_fail == 1) return false;
    }
  }
  return true;
}

/* ----------------------------------------------------------------------
check if an atom type matches that of a reactant atom
------------------------------------------------------------------------- */

bool TopologyMatcher::compare_atomtype(int checktype, Reaction &rxn, int reactant_atom)
{
  int iatom = reactant_atom - 1; // index of reactant atom
  if (checktype == rxn.reactant->type[iatom] || rxn.atoms[iatom].wildcard) return true;
  else return false;
}
