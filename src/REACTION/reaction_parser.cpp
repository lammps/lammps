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

#include "reaction_parser.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "group.h"
#include "input.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "molecule.h"
#include "tokenizer.h"
#include "variable.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

void ReactionParser::parse_reaction(char **arg, int iarg, int rxn_narg, Reaction &rxn) {

  // set defaults

  rxn.fraction = 1.0;
  rxn.seed = 12345;
  rxn.stabilize_steps_flag = 0;
  rxn.custom_charges_fragid = -1;
  rxn.rescale_charges_flag = 0;
  rxn.create_atoms_flag = 0;
  rxn.modify_create_fragid = -1;
  rxn.overlapsq = 0.0;
  rxn.mol_total_charge = 0.0;
  rxn.molecule_keyword = Molecule_Keys::OFF;
  rxn.limit_duration = 60;
  rxn.reaction_count = 0;
  rxn.local_rxn_count = 0;
  rxn.ghostly_rxn_count = 0;
  rxn.reaction_count_total = 0;
  rxn.v_rmin = rxn.v_rmax = -1;
  rxn.v_nevery = rxn.v_prob = -1;

  rxn.name = arg[iarg++];
  if (rxn.name.size()+1 > MAXNAME) error->all(FLERR,"Reaction name (react-ID) is too long (limit: 255 characters)");

  int groupid = group->find(arg[iarg++]);
  if (groupid == -1) error->all(FLERR,"Could not find fix group ID");
  rxn.groupbits = group->bitmask[groupid];

  if (strncmp(arg[iarg],"v_",2) == 0) {
    rxn.v_nevery = input->variable->find(&arg[iarg][2]);
    validate_variable_keyword(&arg[iarg][2], rxn.v_nevery);
  } else {
    rxn.nevery = utils::inumeric(FLERR,arg[iarg],false,lmp);
    if (rxn.nevery <= 0) error->all(FLERR,"Illegal fix bond/react command: "
                                     "'Nevery' must be a positive integer");
  }
  iarg++;

  double cutoff;
  if (strncmp(arg[iarg],"v_",2) == 0) {
    rxn.v_rmin = input->variable->find(&arg[iarg][2]);
    validate_variable_keyword(&arg[iarg][2], rxn.v_rmin);
    cutoff = input->variable->compute_equal(rxn.v_rmin);
  } else cutoff = utils::numeric(FLERR,arg[iarg],false,lmp);
    if (cutoff < 0.0) error->all(FLERR,"Illegal fix bond/react command: "
                                 "'Rmin' cannot be negative");
    rxn.rminsq = cutoff*cutoff;
  iarg++;

  if (strncmp(arg[iarg],"v_",2) == 0) {
    rxn.v_rmax = input->variable->find(&arg[iarg][2]);
    validate_variable_keyword(&arg[iarg][2], rxn.v_rmax);
    cutoff = input->variable->compute_equal(rxn.v_rmax);
  } else cutoff = utils::numeric(FLERR,arg[iarg],false,lmp);
    if (cutoff < 0.0) error->all(FLERR,"Illegal fix bond/react command:"
                                 "'Rmax' cannot be negative");
    rxn.rmaxsq = cutoff*cutoff;
  iarg++;

  int mol_idx = atom->find_molecule(arg[iarg++]);
  if (mol_idx == -1) error->all(FLERR,"Pre-reaction molecule template ID for "
                                           "fix bond/react does not exist");
  rxn.reactant = atom->molecules[mol_idx];
  mol_idx = atom->find_molecule(arg[iarg++]);
  if (mol_idx == -1) error->all(FLERR,"Post-reaction molecule template ID for "
                                         "fix bond/react does not exist");
  rxn.product = atom->molecules[mol_idx];

  // store map file name
  rxn.mapfilename = arg[iarg];
  iarg++;

  while (iarg < rxn_narg) {
    if (strcmp(arg[iarg],"prob") == 0) {
      if (iarg+3 > rxn_narg) error->all(FLERR,"Illegal fix bond/react command: "
                                    "'prob' keyword has too few arguments");
      // check if probability is a variable
      if (strncmp(arg[iarg+1],"v_",2) == 0) {
        rxn.v_prob = input->variable->find(&arg[iarg+1][2]);
        validate_variable_keyword(&arg[iarg+1][2], rxn.v_prob);
        rxn.fraction = input->variable->compute_equal(rxn.v_prob);
      } else {
        // otherwise probability should be a number
        rxn.fraction = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      }
      rxn.seed = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      if (rxn.fraction < 0.0 || rxn.fraction > 1.0)
        error->all(FLERR,"Illegal fix bond/react command: "
                   "probability fraction must between 0 and 1, inclusive");
      if (rxn.seed <= 0) error->all(FLERR,"Illegal fix bond/react command: "
                                     "probability seed must be positive");
      iarg += 3;
    } else if (strcmp(arg[iarg],"stabilize_steps") == 0) {
      if (iarg+2 > rxn_narg) error->all(FLERR,"Illegal fix bond/react command: "
                                    "'stabilize_steps' has too few arguments");
      rxn.limit_duration = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      rxn.stabilize_steps_flag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"custom_charges") == 0) {
      if (iarg+2 > rxn_narg) error->all(FLERR,"Illegal fix bond/react command: "
                                    "'custom_charges' has too few arguments");
      if (strcmp(arg[iarg+1],"no") == 0) rxn.custom_charges_fragid = -1; //default
      else {
        rxn.custom_charges_fragid = rxn.reactant->findfragment(arg[iarg+1]);
        if (rxn.custom_charges_fragid < 0) error->one(FLERR,"Fix bond/react: Molecule fragment for "
                                                       "'custom_charges' keyword does not exist");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"rescale_charges") == 0) {
      if (iarg+2 > rxn_narg) error->all(FLERR,"Illegal fix bond/react command: "
                                    "'rescale_charges' has too few arguments");
      if (strcmp(arg[iarg+1],"no") == 0) rxn.rescale_charges_flag = 0; //default
      else if (strcmp(arg[iarg+1],"yes") == 0) {
        if (!atom->q_flag) error->all(FLERR,"Illegal fix bond/react command: cannot use "
                                    "'rescale_charges' without atomic charges enabled");
        if (!rxn.product->qflag) error->all(FLERR,"Illegal fix bond/react command: cannot use "
                                    "'rescale_charges' without Charges section in post-reaction template");
        rxn.rescale_charges_flag = 1; // overloaded below to also indicate number of atoms to update
      } else error->one(FLERR,"Bond/react: Illegal option for 'rescale_charges' keyword");
      iarg += 2;
    } else if (strcmp(arg[iarg],"molecule") == 0) {
      if (iarg+2 > rxn_narg) error->all(FLERR,"Illegal fix bond/react command: "
                                    "'molecule' has too few arguments");
      if (strcmp(arg[iarg+1],"off") == 0) rxn.molecule_keyword = Molecule_Keys::OFF; //default
      else if (strcmp(arg[iarg+1],"inter") == 0) rxn.molecule_keyword = Molecule_Keys::INTER;
      else if (strcmp(arg[iarg+1],"intra") == 0) rxn.molecule_keyword = Molecule_Keys::INTRA;
      else error->one(FLERR,"Fix bond/react: Illegal option for 'molecule' keyword");
      iarg += 2;
    } else if (strcmp(arg[iarg],"modify_create") == 0) {
      if (iarg++ > rxn_narg) error->all(FLERR,"Illegal fix bond/react command: "
                                    "'modify_create' has too few arguments");
      while (iarg < rxn_narg && strcmp(arg[iarg],"react") != 0) {
        if (strcmp(arg[iarg],"fit") == 0) {
          if (iarg+2 > rxn_narg) error->all(FLERR,"Illegal fix bond/react command: "
                                        "'modify_create' has too few arguments");
          if (strcmp(arg[iarg+1],"all") == 0) rxn.modify_create_fragid = -1; //default
          else {
            rxn.modify_create_fragid = rxn.product->findfragment(arg[iarg+1]);
            if (rxn.modify_create_fragid < 0) error->one(FLERR,"Fix bond/react: Molecule fragment for "
                                                           "'modify_create' keyword does not exist");
          }
          iarg += 2;
        } else if (strcmp(arg[iarg],"overlap") == 0) {
          if (iarg+2 > rxn_narg) error->all(FLERR,"Illegal fix bond/react command: "
                                        "'modify_create' has too few arguments");
          rxn.overlapsq = utils::numeric(FLERR,arg[iarg+1],false,lmp);
          rxn.overlapsq *= rxn.overlapsq;
          iarg += 2;
        } else break;
      }
    } else if (strcmp(arg[iarg],"rate_limit") == 0) {
      error->all(FLERR,"Fix bond/react: 'rate_limit' as an 'individual keyword' has been deprecated. "
                       "Please use the 'rate_limit' common keyword instead, which can be applied to one or more reactions.");
    } else if (strcmp(arg[iarg],"max_rxn") == 0) {
      error->all(FLERR,"Fix bond/react: 'max_rxn' as an 'individual keyword' has been deprecated. "
                       "Please use the 'max_rxn' common keyword instead, which can be applied to one or more reactions.");
    } else error->all(FLERR,"Illegal fix bond/react command: unknown keyword");
  }

  // set more defaults
  rxn.nnewmolids = 0;
  rxn.atoms.resize(rxn.product->natoms);
  int idx = 1;
  for (auto &atm : rxn.atoms) {
    atm.edge = 0;
    atm.wildcard = false;
    atm.recharged = 1; // update all partial charges by default
    atm.deleted = 0;
    atm.created = 0;
    atm.newmolid = 0;
    atm.chiral.fill(0);
    // default amap to their own molecule template atom ID
    // all but created atoms will be updated
    atm.amap.fill(idx++);
  }

  // read map file

  fp = fopen(rxn.mapfilename.c_str(),"r");
  if (fp == nullptr) error->one(FLERR, "Fix bond/react: Cannot open map file {}", rxn.mapfilename);
  rxn.reactant->check_attributes();
  rxn.product->check_attributes();
  read_map_file(rxn);
  fclose(fp);
  rxn.iatomtype = rxn.reactant->type[rxn.ibonding-1];
  rxn.jatomtype = rxn.reactant->type[rxn.jbonding-1];
  //find_landlocked_atoms(rxn);
  if (rxn.custom_charges_fragid >= 0) CustomCharges(rxn.custom_charges_fragid,rxn);
};

/* ----------------------------------------------------------------------
return handedness (1 or -1) of a chiral center, given ordered set of coordinates
------------------------------------------------------------------------- */

int ReactionParser::get_chirality(double four_coords[12])
{
  // define oriented plane with first three coordinates
  double vec1[3],vec2[3],vec3[3],vec4[3],mean3[3],dot;

  for (int i = 0; i < 3; i++) {
    vec1[i] = four_coords[i]-four_coords[i+3];
    vec2[i] = four_coords[i+3]-four_coords[i+6];
  }

  MathExtra::cross3(vec1,vec2,vec3);

  for (int i = 0; i < 3; i++) {
    mean3[i] = (four_coords[i] + four_coords[i+3] +
                four_coords[i+6])/3;
    vec4[i] = four_coords[i+9] - mean3[i];
  }

  dot = MathExtra::dot3(vec3,vec4);
  dot = dot/fabs(dot);
  return (int) dot;
}

/* ----------------------------------------------------------------------
add equal-style variable to keyword argument list
------------------------------------------------------------------------- */

void ReactionParser::validate_variable_keyword(const char *myarg, int var_id)
{
  if (var_id < 0)
    error->all(FLERR,"Fix bond/react: Variable name {} does not exist",myarg);
  if (!input->variable->equalstyle(var_id))
    error->all(FLERR,"Fix bond/react: Variable {} is not equal-style",myarg);
}

/* ----------------------------------------------------------------------
read map file
------------------------------------------------------------------------- */

void ReactionParser::read_map_file(Reaction &rxn)
{
  int nedge, nequivalent, nchiral, nwild, ndelete, ncreate;
  nedge = nchiral = nwild = ndelete = ncreate = 0;
  char line[MAXLINE] = {'\0'};
  char keyword[MAXLINE] = {'\0'};
  char *eof,*ptr;

  // skip 1st line of file
  eof = fgets(line,MAXLINE,fp);
  if (eof == nullptr) error->one(FLERR,"Fix bond/react: Unexpected end of superimpose file");

  // read header lines
  // skip blank lines or lines that start with "#"
  // stop when read an unrecognized line

  while (true) {

    readline(line);

    // trim anything from '#' onward
    // if line is blank, continue

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    if (strspn(line," \t\n\r") == strlen(line)) continue;

    if (strstr(line,"edgeIDs")) {
      nedge = firstint(line, "Map file header is incorrectly formatted");
      if ((nedge < 0) || (nedge > rxn.reactant->natoms))
        error->one(FLERR,"Fix bond/react: Invalid number of edgeIDs in map file");
    }
    else if (strstr(line,"equivalences")) {
      nequivalent = firstint(line, "Map file header is incorrectly formatted");
      if (nequivalent != rxn.reactant->natoms)
        error->one(FLERR,"Fix bond/react: Number of equivalences in map file must "
                   "equal number of atoms in reaction templates");
    }
    else if (strstr(line,"deleteIDs")) {
      ndelete = firstint(line, "Map file header is incorrectly formatted");
      if ((ndelete < 0) || (ndelete > rxn.reactant->natoms))
        error->one(FLERR,"Fix bond/react: Invalid number of deleteIDs in map file");
    } else if (strstr(line,"createIDs")) {
      ncreate = firstint(line, "Map file header is incorrectly formatted");
      if ((ncreate < 0) || (ncreate > rxn.product->natoms))
        error->one(FLERR,"Fix bond/react: Invalid number of createIDs in map file");
    } else if (strstr(line,"chiralIDs")) {
      nchiral = firstint(line, "Map file header is incorrectly formatted");
      if ((nchiral < 0) || (nchiral > rxn.reactant->natoms))
        error->one(FLERR,"Fix bond/react: Invalid number of chiralIDs in map file");
    } else if (strstr(line,"wildcards")) {
      nwild = firstint(line, "Map file header is incorrectly formatted");
      if ((nwild < 0) || (nwild > rxn.reactant->natoms))
        error->one(FLERR,"Fix bond/react: Invalid number of wildcards in map file");
    } else if (strstr(line,"constraints")) {
      int nconstraints = firstint(line, "Map file header is incorrectly formatted");
      if ((nconstraints < 0) || (nconstraints > 4096))
        error->one(FLERR,"Fix bond/react: Invalid number of constraints in map file");
      rxn.constraints.resize(nconstraints);
      for (int i = 0; i < nconstraints; i++) rxn.constraints[i].ID = i;
    } else break;
  }

  if (ncreate == 0 && rxn.reactant->natoms != rxn.product->natoms)
    error->all(FLERR,"Fix bond/react: Reaction templates must contain the same number of atoms");
  else if (ncreate > 0 && rxn.reactant->natoms + ncreate != rxn.product->natoms)
    error->all(FLERR,"Fix bond/react: Incorrect number of created atoms");

  // grab keyword and skip next line

  parse_keyword(0,line,keyword);
  readline(line);

  // loop over sections of superimpose file

  int equivflag = 0, bondflag = 0;
  while (strlen(keyword)) {
    if (strcmp(keyword,"InitiatorIDs") == 0 || strcmp(keyword,"BondingIDs") == 0) {
      if (strcmp(keyword,"BondingIDs") == 0)
        if (comm->me == 0) error->warning(FLERR,"Fix bond/react: The BondingIDs section title has been deprecated. Please use InitiatorIDs instead.");
      bondflag = 1;
      readline(line);
      rxn.ibonding = firstint(line, "InitiatorIDs section is incorrectly formatted");
      if (rxn.ibonding > rxn.reactant->natoms)
        error->one(FLERR,"Fix bond/react: Invalid template atom ID in map file");
      readline(line);
      rxn.jbonding = firstint(line, "InitiatorIDs section is incorrectly formatted");
      if (rxn.jbonding > rxn.reactant->natoms)
        error->one(FLERR,"Fix bond/react: Invalid template atom ID in map file");
    } else if (strcmp(keyword,"EdgeIDs") == 0) {
      EdgeIDs(line, rxn, nedge);
    } else if (strcmp(keyword,"Equivalences") == 0) {
      equivflag = 1;
      Equivalences(line, rxn, nequivalent);
    } else if (strcmp(keyword,"DeleteIDs") == 0) {
      DeleteAtoms(line, rxn, ndelete);
    } else if (strcmp(keyword,"CreateIDs") == 0) {
      CreateAtoms(line, rxn, ncreate);
    } else if (strcmp(keyword,"ChiralIDs") == 0) {
      ChiralCenters(line, rxn, nchiral);
    } else if (strcmp(keyword,"Wildcards") == 0) {
      ReadWildcards(line, rxn, nwild);
    } else if (strcmp(keyword,"Constraints") == 0) {
      ReadConstraints(line, rxn);
    } else error->one(FLERR,"Fix bond/react: Unknown section in map file");

    parse_keyword(1,line,keyword);

  }

  // error check
  for (int i = 0; i < rxn.reactant->natoms; i++) {
    int my_equiv = rxn.atoms[i].ramap[1];
    if (rxn.atoms[my_equiv-1].created == 1)
      error->all(FLERR,"Fix bond/react: Created atoms cannot also be listed in Equivalences section\n");
  }

  // error check
  if (bondflag == 0 || equivflag == 0)
    error->all(FLERR,"Fix bond/react: Map file missing InitiatorIDs or Equivalences section\n");
}

void ReactionParser::EdgeIDs(char *line, Reaction &rxn, int nedge)
{
  // puts a 1 at edge(edgeID)

  int tmp;
  for (int i = 0; i < nedge; i++) {
    readline(line);
    tmp = firstint(line, "EdgeIDs section is incorrectly formatted");
    if ((tmp < 1) || (tmp > rxn.reactant->natoms))
      error->one(FLERR,"Fix bond/react: Invalid template atom ID in map file");
    rxn.atoms[tmp-1].edge = 1;
  }
}

void ReactionParser::Equivalences(char *line, Reaction &rxn, int nequivalent)
{
  int tmp1 = 0, tmp2 = 0;
  for (int i = 0; i < nequivalent; i++) {
    readline(line);
    try {
      ValueTokenizer values(line);
      tmp1 = values.next_int();
      tmp2 = values.next_int();
    } catch (TokenizerException &) {
      error->one(FLERR, "Equivalences section is incorrectly formatted");
    }
    if ((tmp1 < 1) || (tmp1 > rxn.reactant->natoms) || (tmp2 < 1) || (tmp2 > rxn.product->natoms))
      error->one(FLERR,"Fix bond/react: Invalid template atom ID in map file");
    //equivalences is-> clmn 1: post-reacted, clmn 2: pre-reacted
    rxn.atoms[tmp2-1].amap[0] = tmp2;
    rxn.atoms[tmp2-1].amap[1] = tmp1;
    //reverse_equiv is-> clmn 1: pre-reacted, clmn 2: post-reacted
    rxn.atoms[tmp1-1].ramap[0] = tmp1;
    rxn.atoms[tmp1-1].ramap[1] = tmp2;
  }
  // sanity check for one-to-one mapping for equivalences
  for (int i = 0; i < rxn.product->natoms; i++) {
    if (rxn.atoms[i].created == 1) continue;
    for (int j = i+1; j < rxn.product->natoms; j++) {
      if (rxn.atoms[j].created == 1) continue;
      if (rxn.atoms[i].amap[0] == rxn.atoms[j].amap[0] ||
          rxn.atoms[i].amap[1] == rxn.atoms[j].amap[1]) {
        error->one(FLERR,"Fix bond/react: Repeated atoms IDs in Equivalences section");
      }
    }
  }
}

void ReactionParser::DeleteAtoms(char *line, Reaction &rxn, int ndelete)
{
  int tmp;
  for (int i = 0; i < ndelete; i++) {
    readline(line);
    tmp = firstint(line, "DeleteIDs section is incorrectly formatted");
    if ((tmp < 1) || (tmp > rxn.reactant->natoms))
      error->one(FLERR,"Fix bond/react: Invalid template atom ID in map file");
    rxn.atoms[tmp-1].deleted = 1;
  }
}

void ReactionParser::CreateAtoms(char *line, Reaction &rxn, int ncreate)
{
  rxn.create_atoms_flag = 1;
  int tmp;
  for (int i = 0; i < ncreate; i++) {
    readline(line);
    tmp = firstint(line, "CreateIDs section is incorrectly formatted");
    if ((tmp < 1) || (tmp > rxn.product->natoms))
      error->one(FLERR, Error::NOLASTLINE, "Fix bond/react: Invalid atom ID in CreateIDs section of map file");
    rxn.atoms[tmp-1].created = 1;
  }
  if (rxn.product->xflag == 0)
    error->one(FLERR, Error::NOLASTLINE,
               "Fix bond/react: 'Coords' section required in post-reaction template when creating new atoms");
  if (atom->rmass_flag && !rxn.product->rmassflag)
    error->one(FLERR, Error::NOLASTLINE,
               "Fix bond/react: 'Masses' section required in post-reaction template when creating new atoms "
               "and per-atom masses are defined.");
}

void ReactionParser::CustomCharges(int ifragment, Reaction &rxn)
{
  for (int i = 0; i < rxn.reactant->natoms; i++)
    if (rxn.reactant->fragmentmask[ifragment][i])
      rxn.atoms[i].recharged = 1;
    else
      rxn.atoms[i].recharged = 0;
}

void ReactionParser::ChiralCenters(char *line, Reaction &rxn, int nchiral)
{
  int tmp;
  for (int i = 0; i < nchiral; i++) {
    readline(line);
    tmp = firstint(line, "ChiralIDs section is incorrectly formatted");
    if ((tmp < 1) || (tmp > rxn.reactant->natoms))
      error->one(FLERR,"Fix bond/react: Invalid template atom ID in map file");
    rxn.atoms[tmp-1].chiral[0] = 1;
    if (rxn.reactant->xflag == 0)
      error->one(FLERR,"Fix bond/react: Molecule template 'Coords' section required for chiralIDs keyword");
    if ((int) rxn.reactant->nspecial[tmp-1][0] != 4)
      error->one(FLERR,"Fix bond/react: Chiral atoms must have exactly four first neighbors");
    for (int j = 0; j < 4; j++) {
      for (int k = j+1; k < 4; k++) {
        if (rxn.reactant->type[rxn.reactant->special[tmp-1][j]-1] ==
            rxn.reactant->type[rxn.reactant->special[tmp-1][k]-1])
          error->one(FLERR,"Fix bond/react: First neighbors of chiral atoms must be of mutually different types");
      }
    }
    // record order of atom types, and coords
    double my4coords[12] = {0.0};
    for (int j = 0; j < 4; j++) {
      rxn.atoms[tmp-1].chiral[j+2] = rxn.reactant->type[rxn.reactant->special[tmp-1][j]-1];
      for (int k = 0; k < 3; k++) {
        my4coords[3*j+k] = rxn.reactant->x[rxn.reactant->special[tmp-1][j]-1][k];
      }
    }
    // get orientation
    rxn.atoms[tmp-1].chiral[1] = get_chirality(my4coords);
  }
}

void ReactionParser::ReadWildcards(char *line, Reaction &rxn, int nwild)
{
  int tmp;
  for (int i = 0; i < nwild; i++) {
    readline(line);
    tmp = firstint(line, "Wildcards section is incorrectly formatted");
    if ((tmp < 1) || (tmp > rxn.reactant->natoms))
      error->one(FLERR,"Bond/react: Invalid template atom ID in map file");
    rxn.atoms[tmp-1].wildcard = true;
  }
}

void ReactionParser::ReadConstraints(char *line, Reaction &rxn)
{
  int rv;
  double tmp[MAXCONARGS];
  char **strargs,*ptr,*lptr;
  memory->create(strargs,MAXCONARGS,MAXLINE,"bond/react:strargs");
  auto *constraint_type = new char[MAXLINE];
  rxn.constraintstr = "("; // string for boolean constraint logic
  for (auto &constraint : rxn.constraints) {
    readline(line);
    // find left parentheses, add to constraintstr, and update line
    for (int j = 0; j < (int)strlen(line); j++) {
      if (line[j] == '(') rxn.constraintstr += "(";
      if (isalpha(line[j])) {
        line = line + j;
        break;
      }
    }
    // 'C' indicates where to sub in next constraint
    rxn.constraintstr += "C";
    // special consideration for 'custom' constraint
    // find final double quote, or skip two words
    lptr = line;
    if ((ptr = strrchr(lptr,'\"'))) lptr = ptr+1;
    else {
      while (lptr[0] != ' ') lptr++; // skip first 'word'
      while (lptr[0] == ' ' || lptr[0] == '\t') lptr++; // skip blanks
      while (lptr[0] != ' ') lptr++; // skip second 'word'
    }
    // find right parentheses
    for (int j = 0; j < (int)strlen(lptr); j++)
      if (lptr[j] == ')') rxn.constraintstr += ")";
    // find logic symbols, and trim line via ptr
    if ((ptr = strstr(lptr,"&&"))) {
      rxn.constraintstr += "&&";
      *ptr = '\0';
    } else if ((ptr = strstr(lptr,"||"))) {
      rxn.constraintstr += "||";
      *ptr = '\0';
    } else if (constraint.ID+1 < (int)rxn.constraints.size()) {
      rxn.constraintstr += "&&";
    }
    if ((ptr = strchr(lptr,')')))
      *ptr = '\0';
    rv = sscanf(line,"%s",constraint_type);
    if (rv != 1) error->one(FLERR, "Constraints section is incorrectly formatted");
    if (strcmp(constraint_type,"distance") == 0) {
      constraint.type = Reaction::Constraint::Type::DISTANCE;
      rv = sscanf(line,"%*s %s %s %lg %lg",strargs[0],strargs[1],&tmp[0],&tmp[1]);
      if (rv != 4) error->one(FLERR, "Distance constraint is incorrectly formatted");
      readID(strargs[0], constraint, rxn, 0);
      readID(strargs[1], constraint, rxn, 1);
      // cutoffs
      constraint.distance.rminsq = tmp[0]*tmp[0]; // using square of distance
      constraint.distance.rmaxsq = tmp[1]*tmp[1];
    } else if (strcmp(constraint_type,"angle") == 0) {
      constraint.type = Reaction::Constraint::Type::ANGLE;
      rv = sscanf(line,"%*s %s %s %s %lg %lg",strargs[0],strargs[1],strargs[2],&tmp[0],&tmp[1]);
      if (rv != 5) error->one(FLERR, "Angle constraint is incorrectly formatted");
      readID(strargs[0], constraint, rxn, 0);
      readID(strargs[1], constraint, rxn, 1);
      readID(strargs[2], constraint, rxn, 2);
      constraint.angle.amin = tmp[0]/180.0 * MY_PI;
      constraint.angle.amax = tmp[1]/180.0 * MY_PI;
    } else if (strcmp(constraint_type,"dihedral") == 0) {
      constraint.type = Reaction::Constraint::Type::DIHEDRAL;
      tmp[2] = 181.0; // impossible range
      tmp[3] = 182.0;
      rv = sscanf(line,"%*s %s %s %s %s %lg %lg %lg %lg",strargs[0],strargs[1],
             strargs[2],strargs[3],&tmp[0],&tmp[1],&tmp[2],&tmp[3]);
      if (rv != 6 && rv != 8) error->one(FLERR, "Dihedral constraint is incorrectly formatted");
      readID(strargs[0], constraint, rxn, 0);
      readID(strargs[1], constraint, rxn, 1);
      readID(strargs[2], constraint, rxn, 2);
      readID(strargs[3], constraint, rxn, 3);
      constraint.dihedral.amin = tmp[0]/180.0 * MY_PI;
      constraint.dihedral.amax = tmp[1]/180.0 * MY_PI;
      constraint.dihedral.amin2 = tmp[2]/180.0 * MY_PI;
      constraint.dihedral.amax2 = tmp[3]/180.0 * MY_PI;
    } else if (strcmp(constraint_type,"arrhenius") == 0) {
      constraint.type = Reaction::Constraint::Type::ARRHENIUS;
      rv = sscanf(line,"%*s %lg %lg %lg %lg",&tmp[0],&tmp[1],&tmp[2],&tmp[3]);
      if (rv != 4) error->one(FLERR, "Arrhenius constraint is incorrectly formatted");
      constraint.arrhenius.A = tmp[0];
      constraint.arrhenius.n = tmp[1];
      constraint.arrhenius.E_a = tmp[2];
      constraint.arrhenius.seed = tmp[3];
    } else if (strcmp(constraint_type,"rmsd") == 0) {
      constraint.type = Reaction::Constraint::Type::RMSD;
      strcpy(strargs[0],"0");
      rv = sscanf(line,"%*s %lg %s",&tmp[0],strargs[0]);
      if (rv != 1 && rv != 2) error->one(FLERR, "RMSD constraint is incorrectly formatted");
      constraint.rmsd.rmsdmax = tmp[0];
      constraint.ids[0] = -1; // optional molecule fragment
      if (isalpha(strargs[0][0])) {
        int ifragment = rxn.reactant->findfragment(strargs[0]);
        if (ifragment < 0) error->one(FLERR,"Fix bond/react: Molecule fragment does not exist");
        else constraint.ids[0] = ifragment;
      }
    } else if (strcmp(constraint_type,"custom") == 0) {
      constraint.type = Reaction::Constraint::Type::CUSTOM;
      std::vector<std::string> args = utils::split_words(line);
      constraint.custom.str = args[1];
    } else error->one(FLERR,"Fix bond/react: Illegal constraint type in 'Constraints' section of map file");
  }
  rxn.constraintstr += ")"; // close boolean constraint logic string
  delete[] constraint_type;
  memory->destroy(strargs);
}

/* ----------------------------------------------------------------------
if ID starts with character, assume it is a pre-reaction molecule fragment ID
otherwise, it is a pre-reaction atom ID
---------------------------------------------------------------------- */

void ReactionParser::readID(char *strarg, Reaction::Constraint &constraint, Reaction &rxn, int i)
{
  if (isalpha(strarg[0])) {
    constraint.idtypes[i] = Reaction::Constraint::IDType::FRAG; // fragment vs. atom ID flag
    int ifragment = rxn.reactant->findfragment(strarg);
    if (ifragment < 0)
      error->one(FLERR,"Fix bond/react: Molecule fragment {} does not exist", strarg);
    constraint.ids[i] = ifragment;
  } else {
    constraint.idtypes[i] = Reaction::Constraint::IDType::ATOM; // fragment vs. atom ID flag
    int iatom = utils::inumeric(FLERR, strarg, true, lmp);
    if (iatom > rxn.reactant->natoms)
      error->one(FLERR,"Fix bond/react: Invalid template atom ID {} in map file", strarg);
    constraint.ids[i] = iatom;
  }
}

void ReactionParser::readline(char *line)
{
  int n;
  if (comm->me == 0) {
    if (fgets(line,MAXLINE,fp) == nullptr) n = 0;
    else n = strlen(line) + 1;
  }
  MPI_Bcast(&n,1,MPI_INT,0,world);
  if (n == 0) error->all(FLERR,"Fix bond/react: Unexpected end of map file");
  MPI_Bcast(line,n,MPI_CHAR,0,world);
}

void ReactionParser::parse_keyword(int flag, char *line, char *keyword)
{
  if (flag) {

    // read upto non-blank line plus 1 following line
    // eof is set to 1 if any read hits end-of-file

    int eof = 0;
    if (comm->me == 0) {
      if (fgets(line,MAXLINE,fp) == nullptr) eof = 1;
      while (eof == 0 && strspn(line," \t\n\r") == strlen(line)) {
        if (fgets(line,MAXLINE,fp) == nullptr) eof = 1;
      }
      if (fgets(keyword,MAXLINE,fp) == nullptr) eof = 1;
    }

    // if eof, set keyword empty and return

    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) {
      keyword[0] = '\0';
      return;
    }

    // bcast keyword line to all procs

    int n;
    if (comm->me == 0) n = strlen(line) + 1;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);
  }

  // copy non-whitespace portion of line into keyword

  int start = strspn(line," \t\n\r");
  int stop = strlen(line) - 1;
  while (line[stop] == ' ' || line[stop] == '\t'
         || line[stop] == '\n' || line[stop] == '\r') stop--;
  line[stop+1] = '\0';
  strcpy(keyword,&line[start]);
}

// parse the first integer on a line, error out with errmsg if not a number

int ReactionParser::firstint(char *line, const char *errmsg)
{
  try {
    return ValueTokenizer(line).next_int();
  } catch (TokenizerException &) {
    error->one(FLERR, errmsg);
  }
  return -1;
}
